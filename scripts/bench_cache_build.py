#!/usr/bin/env python3
"""Benchmark VEPYR cache building (build_cache) parallelism.

Etap 0 of the cache-parallelization plan. Measures wall-clock + peak RSS of
``vepyr.build_cache`` for a single-chromosome (or arbitrary subset) raw Ensembl
VEP cache, across a swept knob (``partitions`` today; ``build_concurrency`` once
the scheduler lands). Each configuration runs in a fresh subprocess with a clean
output directory, mirroring vepyr/e2e-testing/scripts/benchmark_parallel_lookup.py.

Outputs (in --results-dir):
  runs.jsonl           one JSON object per run
  runs.csv             flat table
  summary.json         environment + per-config median/min/max + speedup
  RESULTS.md           human-readable table + scaling curve

The current vepyr.build_cache signature has NO chromosome filter, so we build a
chr-subset raw cache by copying the requested per-chromosome directories (plus
info.txt / chr_synonyms.txt) into a temporary local_cache.
"""

from __future__ import annotations

import argparse
import csv
import inspect
import json
import os
import platform
import resource
import shutil
import statistics
import subprocess
import sys
import time
from pathlib import Path

# Knobs that build_cache may accept for scaling. "partitions" exists today;
# "build_concurrency" is added by the unified scheduler (Etap 3).
SUPPORTED_KNOBS = ("partitions", "build_concurrency")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument(
        "--vep-cache-root",
        type=Path,
        default=Path.home() / ".vep" / "homo_sapiens" / "115_GRCh38",
        help="Raw Ensembl VEP cache dir containing info.txt + per-chrom dirs.",
    )
    parser.add_argument(
        "--chromosomes",
        default="1",
        help="Comma-separated chromosomes to build (subset cache). e.g. '1' or '22'.",
    )
    parser.add_argument(
        "--prepared-cache",
        type=Path,
        help="Use an already-prepared subset local_cache instead of copying.",
    )
    parser.add_argument("--cache-type", default="ensembl")
    parser.add_argument("--release", type=int, default=115)
    parser.add_argument("--assembly", default="GRCh38")
    parser.add_argument(
        "--knob",
        choices=SUPPORTED_KNOBS,
        default="partitions",
        help="Which build_cache parameter to sweep for the scaling curve.",
    )
    parser.add_argument(
        "--values",
        default="1,2,4,8,16",
        help="Comma-separated knob values (e.g. concurrency/partition levels).",
    )
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument(
        "--fixed-partitions",
        type=int,
        help="Hold build_cache 'partitions' fixed while sweeping a different knob.",
    )
    parser.add_argument(
        "--keep-output",
        action="store_true",
        help="Keep built cache dirs (needed to use one as a correctness baseline).",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        help="Where to write built caches (default: <results-dir>/caches).",
    )
    # --single-run internals (one build, prints JSON)
    parser.add_argument("--single-run", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--local-cache", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--out-cache-dir", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--knob-value", type=int, help=argparse.SUPPRESS)
    return parser.parse_args()


def parse_chromosomes(value: str) -> list[str]:
    chroms = [c.strip().removeprefix("chr") for c in value.split(",") if c.strip()]
    if not chroms:
        raise ValueError("--chromosomes must contain at least one chromosome")
    return chroms


def parse_values(value: str) -> list[int]:
    values = [int(v) for v in value.split(",")]
    if not values or any(v <= 0 for v in values):
        raise ValueError("--values must contain positive integers")
    return values


def max_rss_bytes() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform == "darwin":
        return int(value)  # darwin reports bytes
    return int(value) * 1024  # linux reports KiB


def prepare_subset_cache(vep_root: Path, chroms: list[str], dest: Path) -> Path:
    """Copy info.txt + requested per-chrom dirs into a fresh subset local_cache."""
    if dest.exists():
        shutil.rmtree(dest)
    dest.mkdir(parents=True)
    for meta in ("info.txt", "chr_synonyms.txt"):
        src = vep_root / meta
        if src.exists():
            shutil.copy2(src, dest / meta)
    wanted = set(chroms)
    for chrom in chroms:
        src_dir = vep_root / chrom
        if not src_dir.is_dir():
            raise FileNotFoundError(f"chromosome dir not found: {src_dir}")
        shutil.copytree(src_dir, dest / chrom)
    present = sorted(p.name for p in dest.iterdir() if p.is_dir())
    if set(present) != wanted:
        raise RuntimeError(f"prepared cache chroms {present} != requested {sorted(wanted)}")
    return dest


def git_value(*args: str) -> str:
    try:
        return subprocess.run(
            ["git", *args], check=True, capture_output=True, text=True
        ).stdout.strip()
    except Exception:
        return "unknown"


def run_once(args: argparse.Namespace) -> int:
    """Build the cache once, print a JSON measurement line."""
    import vepyr

    knob_kwargs: dict[str, int] = {}
    sig = inspect.signature(vepyr.build_cache)
    if args.knob not in sig.parameters:
        raise SystemExit(
            f"build_cache does not accept knob '{args.knob}'. "
            f"Available: {sorted(sig.parameters)}"
        )
    knob_kwargs[args.knob] = args.knob_value
    if args.fixed_partitions is not None and args.knob != "partitions":
        knob_kwargs["partitions"] = args.fixed_partitions

    started = time.perf_counter()
    files = vepyr.build_cache(
        args.release,
        str(args.out_cache_dir),
        cache_type=args.cache_type,
        assembly=args.assembly,
        local_cache=str(args.local_cache),
        show_progress=False,
        overwrite=True,
        **knob_kwargs,
    )
    elapsed = time.perf_counter() - started

    total_rows = sum(rows for _, rows in files)
    total_bytes = 0
    for path, _ in files:
        try:
            total_bytes += Path(path).stat().st_size
        except OSError:
            pass
    print(
        json.dumps(
            {
                "elapsed_seconds": round(elapsed, 6),
                "max_rss_bytes": max_rss_bytes(),
                "num_files": len(files),
                "total_rows": total_rows,
                "total_parquet_bytes": total_bytes,
            },
            sort_keys=True,
        )
    )
    return 0


def write_csv(path: Path, rows: list[dict]) -> None:
    fields = [
        "knob",
        "knob_value",
        "repeat",
        "elapsed_seconds",
        "max_rss_bytes",
        "num_files",
        "total_rows",
        "total_parquet_bytes",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows({f: row.get(f, "") for f in fields} for row in rows)


def summarize(rows: list[dict]) -> list[dict]:
    summaries = []
    for value in sorted({row["knob_value"] for row in rows}):
        selected = [r for r in rows if r["knob_value"] == value]
        elapsed = [r["elapsed_seconds"] for r in selected]
        rss = [r["max_rss_bytes"] for r in selected]
        rowcounts = {r["total_rows"] for r in selected}
        summaries.append(
            {
                "knob_value": value,
                "runs": len(selected),
                "median_seconds": round(statistics.median(elapsed), 6),
                "min_seconds": round(min(elapsed), 6),
                "max_seconds": round(max(elapsed), 6),
                "median_max_rss_bytes": int(statistics.median(rss)),
                "total_rows": sorted(rowcounts),
                "rows_consistent": len(rowcounts) == 1,
            }
        )
    baseline = next((s for s in summaries if s["knob_value"] == 1), None)
    base_t = baseline["median_seconds"] if baseline else None
    for s in summaries:
        s["speedup_vs_1"] = (
            round(base_t / s["median_seconds"], 3) if base_t else None
        )
    return summaries


def write_results_md(path: Path, summaries: list[dict], meta: dict) -> None:
    lines = [
        f"# Cache build benchmark — {meta['chromosomes']}",
        "",
        f"- created: {meta['created_at']}",
        f"- knob swept: `{meta['knob']}`",
        f"- platform: {meta['platform']}",
        f"- git: {meta['git_branch']}@{meta['git_commit'][:10]}",
        f"- raw cache: {meta['local_cache']}",
        "",
        f"| {meta['knob']} | median s | min s | max s | speedup | RSS GiB | rows |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for s in summaries:
        rows_repr = s["total_rows"][0] if s["rows_consistent"] else f"!{s['total_rows']}"
        lines.append(
            f"| {s['knob_value']} | {s['median_seconds']:.2f} | {s['min_seconds']:.2f} "
            f"| {s['max_seconds']:.2f} | {s['speedup_vs_1'] or '-'}x "
            f"| {s['median_max_rss_bytes'] / 1024**3:.2f} | {rows_repr} |"
        )
    if not all(s["rows_consistent"] for s in summaries):
        lines += ["", "> WARNING: row counts differ across configs — output not stable!"]
    path.write_text("\n".join(lines) + "\n")


def orchestrate(args: argparse.Namespace) -> int:
    chroms = parse_chromosomes(args.chromosomes)
    values = parse_values(args.values)
    if args.repeats <= 0:
        raise ValueError("--repeats must be positive")

    args.results_dir.mkdir(parents=True, exist_ok=True)
    output_root = args.output_root or (args.results_dir / "caches")
    output_root.mkdir(parents=True, exist_ok=True)

    if args.prepared_cache is not None:
        local_cache = args.prepared_cache
        if not local_cache.is_dir():
            raise FileNotFoundError(local_cache)
    else:
        local_cache = prepare_subset_cache(
            args.vep_cache_root, chroms, args.results_dir / "raw_subset"
        )
    print(f"raw subset cache: {local_cache} (chroms={chroms})", flush=True)

    rows: list[dict] = []
    raw_path = args.results_dir / "runs.jsonl"
    csv_path = args.results_dir / "runs.csv"

    scheduled = [(v, r) for v in values for r in range(1, args.repeats + 1)]
    with raw_path.open("w") as raw_handle:
        for knob_value, repeat in scheduled:
            out_dir = output_root / f"{args.knob}{knob_value}_r{repeat}"
            if out_dir.exists():
                shutil.rmtree(out_dir)
            command = [
                sys.executable,
                str(Path(__file__).resolve()),
                "--single-run",
                "--results-dir", str(args.results_dir),
                "--local-cache", str(local_cache),
                "--out-cache-dir", str(out_dir),
                "--cache-type", args.cache_type,
                "--release", str(args.release),
                "--assembly", args.assembly,
                "--knob", args.knob,
                "--knob-value", str(knob_value),
            ]
            if args.fixed_partitions is not None:
                command += ["--fixed-partitions", str(args.fixed_partitions)]
            completed = subprocess.run(command, check=True, capture_output=True, text=True)
            measurement = json.loads(completed.stdout.strip().splitlines()[-1])
            row = {"knob": args.knob, "knob_value": knob_value, "repeat": repeat, **measurement}
            rows.append(row)
            raw_handle.write(json.dumps(row, sort_keys=True) + "\n")
            raw_handle.flush()
            write_csv(csv_path, rows)
            print(
                f"{args.knob}={knob_value} repeat={repeat}: "
                f"{measurement['elapsed_seconds']:.2f}s, "
                f"rss={measurement['max_rss_bytes'] / 1024**3:.2f}GiB, "
                f"rows={measurement['total_rows']}",
                flush=True,
            )
            if not args.keep_output:
                shutil.rmtree(out_dir, ignore_errors=True)

    summaries = summarize(rows)
    meta = {
        "created_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "platform": platform.platform(),
        "python": sys.version,
        "git_commit": git_value("rev-parse", "HEAD"),
        "git_branch": git_value("branch", "--show-current"),
        "local_cache": str(local_cache),
        "chromosomes": ",".join(chroms),
        "knob": args.knob,
        "values": values,
        "repeats": args.repeats,
        "fixed_partitions": args.fixed_partitions,
    }
    (args.results_dir / "summary.json").write_text(
        json.dumps({**meta, "results": summaries}, indent=2) + "\n"
    )
    write_results_md(args.results_dir / "RESULTS.md", summaries, meta)
    write_csv(csv_path, rows)
    print(f"\nwrote results to {args.results_dir}", flush=True)
    return 0


def main() -> int:
    args = parse_args()
    return run_once(args) if args.single_run else orchestrate(args)


if __name__ == "__main__":
    raise SystemExit(main())
