#!/usr/bin/env python3
"""Re-measure the VEPYR plugin series on chr1/chr22: five repeats, all with measure_once.py.

Added 2026-10-25. In the September series (runs/5repeats_m5) repeat 1 came from
`vepyr-diffly benchmark-plugins` (time including interpreter start-up, sampled RSS of the
process tree) and repeats 2-5 from measure_once.py. This run measures all five repeats with
measure_once.py under the same settings as repeats 2-5 of run_5repeats_benchmark.py:
same Python environment and vepyr build, workers=1, everything + hgvs + reference FASTA,
cache/116_GRCh38_merged, the same per-plugin plugin_cache directories (symlinks into
plugin_cache/), repeat as the outer loop and configuration as the inner loop, chr22 before chr1.
Inputs are the chr1/chr22 slices recreated by step S0 of the completion run with the same
bcftools command as on 2026-09-24 (the September files were deleted).

Usage: python3 -I tools/rerun_plugin_series_5x.py [output dir]
"""

from __future__ import annotations

import json
import statistics
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path("/Users/lukaszjezapkowicz/workspace/data_vepyr")
PY = ROOT / ".venv-native/bin/python3"
MEASURE_SCRIPT = ROOT / "tools/measure_once.py"
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "runs/plugin_series_m5_20261025"
TMP = ROOT / "output/plugin_series_tmp"
CACHE_DIR = ROOT / "cache/116_GRCh38_merged"
FASTA = ROOT / "input/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
INPUTS = {c: ROOT / f"runs/m5_completion_20261024/inputs/HG002_normalized_{c}.vcf.gz" for c in ("chr22", "chr1")}
PLUGIN_ROOT = ROOT / "runs/chapter5-m5-20260924T061203Z/correctness/{chrom}_individual_plugins/{plugin}/plugin_cache"
CONFIGS = ["baseline", "clinvar", "spliceai", "cadd", "alphamissense", "dbnsfp"]
REPEATS = range(1, 6)
PAUSE_S = {"chr22": 2, "chr1": 3}


def run_single(chrom: str, config: str, out_vcf: Path) -> dict:
    cmd = [
        str(PY), str(MEASURE_SCRIPT),
        "--input-vcf", str(INPUTS[chrom]),
        "--cache-dir", str(CACHE_DIR),
        "--reference-fasta", str(FASTA),
        "--workers", "1",
        "--output-vcf", str(out_vcf),
    ]
    if config != "baseline":
        cmd += ["--plugin-cache-root", str(PLUGIN_ROOT).format(chrom=chrom, plugin=config)]
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
    data = json.loads(proc.stdout.strip().splitlines()[-1])
    data["wall_seconds"] = time.perf_counter() - t0
    return data


def summarise(runs: list[dict]) -> dict:
    out: dict = {}
    for chrom in INPUTS:
        base = None
        for config in CONFIGS:
            rs = [r for r in runs if r["chrom"] == chrom and r["config"] == config]
            t = [r["annotation_seconds"] for r in rs]
            gib = [r["max_rss_bytes"] / 1024**3 for r in rs]
            mean_t = statistics.mean(t)
            if config == "baseline":
                base = mean_t
            out[f"{chrom}/{config}"] = {
                "n": len(rs),
                "mean_s": round(mean_t, 4),
                "min_s": round(min(t), 4),
                "max_s": round(max(t), 4),
                "rsd_pct": round(100 * statistics.stdev(t) / mean_t, 2) if len(t) > 1 else None,
                "overhead_pct": None if config == "baseline" else round(100 * (mean_t - base) / base, 2),
                "mean_peak_rss_gib": round(statistics.mean(gib), 4),
            }
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    TMP.mkdir(parents=True, exist_ok=True)
    for d in (OUT, TMP):
        (d / ".metadata_never_index").touch()
    raw = OUT / "raw_runs.jsonl"
    if raw.exists():
        sys.exit(f"{raw} exists; refusing to append to an earlier series")
    for chrom in INPUTS:
        for config in CONFIGS[1:]:
            target = Path(str(PLUGIN_ROOT).format(chrom=chrom, plugin=config)) / "plugin" / config
            if not target.resolve().is_dir():
                sys.exit(f"missing plugin cache {target}")

    runs = []
    for chrom in INPUTS:
        for rep in REPEATS:
            for config in CONFIGS:
                print(f"{chrom} {config} repeat {rep}/5 ...", end="", flush=True)
                res = run_single(chrom, config, TMP / f"{chrom}_{config}_r{rep}.vcf")
                res.update(chrom=chrom, config=config, repeat=rep,
                           timestamp=datetime.now(timezone.utc).isoformat())
                runs.append(res)
                with raw.open("a", encoding="utf-8") as f:
                    f.write(json.dumps(res, sort_keys=True) + "\n")
                print(f" {res['annotation_seconds']:.2f} s, {res['max_rss_bytes'] / 1024**3:.2f} GiB", flush=True)
                time.sleep(PAUSE_S[chrom])

    summary = summarise(runs)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
