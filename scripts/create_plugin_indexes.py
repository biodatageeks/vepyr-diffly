from __future__ import annotations

import argparse
import gzip
import shutil
import struct
import subprocess
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class IndexSpec:
    label: str
    pattern: str
    kind: str  # "vcf" | "tsv"


PLUGIN_INDEX_SPECS: dict[str, list[IndexSpec]] = {
    "clinvar": [IndexSpec("clinvar", "clinvar.vcf.gz", "vcf")],
    "spliceai": [
        IndexSpec(
            "spliceai",
            "spliceai_scores.masked.snv.ensembl_mane.grch38.110.vcf.gz",
            "vcf",
        )
    ],
    "dbnsfp": [IndexSpec("dbnsfp", "dbNSFP*.gz", "tsv")],
    "alphamissense": [IndexSpec("alphamissense", "AlphaMissense*.tsv.gz", "tsv")],
    "cadd": [
        IndexSpec("cadd_snv", "whole_genome_SNVs.tsv.gz", "tsv"),
        IndexSpec("cadd_indel", "gnomad.genomes.r4.0.indel.tsv.gz", "tsv"),
    ],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create .tbi indexes for plugin source files. BGZF inputs are indexed directly; "
            "plain gzip inputs can optionally be recompressed to sibling .bgz files first."
        )
    )
    parser.add_argument(
        "--plugins-dir",
        type=Path,
        default=Path("plugins"),
        help="Directory containing plugin source files (default: ./plugins).",
    )
    parser.add_argument(
        "--plugins",
        default="clinvar,spliceai,cadd,alphamissense,dbnsfp",
        help="Comma-separated plugin names to process.",
    )
    parser.add_argument(
        "--recompress-plain-gzip",
        action="store_true",
        help=(
            "For plain gzip inputs, create a sibling BGZF file (*.bgz) and index that file. "
            "The original .gz is left untouched."
        ),
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force tabix reindexing even when a .tbi file already exists.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print planned actions without writing files.",
    )
    return parser.parse_args()


def is_bgzf(path: Path) -> bool:
    with path.open("rb") as handle:
        header = handle.read(18)
    if len(header) < 18:
        return False
    if header[0:2] != b"\x1f\x8b":
        return False
    flg = header[3]
    if not (flg & 0x04):
        return False
    xlen = struct.unpack("<H", header[10:12])[0]
    extra = header[12 : 12 + xlen]
    i = 0
    while i + 4 <= len(extra):
        si1 = extra[i : i + 1]
        si2 = extra[i + 1 : i + 2]
        slen = struct.unpack("<H", extra[i + 2 : i + 4])[0]
        data_start = i + 4
        data_end = data_start + slen
        if data_end > len(extra):
            break
        if si1 == b"B" and si2 == b"C":
            return True
        i = data_end
    return False


def resolve_plugin_names(raw: str) -> list[str]:
    names = [name.strip().lower() for name in raw.split(",") if name.strip()]
    if not names:
        raise SystemExit("no plugins requested")
    unknown = sorted(set(names) - set(PLUGIN_INDEX_SPECS))
    if unknown:
        raise SystemExit(
            f"unsupported plugin name(s): {', '.join(unknown)}; supported: {', '.join(sorted(PLUGIN_INDEX_SPECS))}"
        )
    return names


def resolve_inputs(plugins_dir: Path, plugin_names: list[str]) -> list[tuple[str, IndexSpec, Path]]:
    resolved: list[tuple[str, IndexSpec, Path]] = []
    for plugin_name in plugin_names:
        for spec in PLUGIN_INDEX_SPECS[plugin_name]:
            matches = sorted(plugins_dir.glob(spec.pattern))
            if not matches:
                raise SystemExit(
                    f"missing input for {plugin_name}: expected a file matching '{spec.pattern}' in {plugins_dir}"
                )
            if len(matches) > 1:
                raise SystemExit(
                    f"ambiguous input for {plugin_name}: pattern '{spec.pattern}' matched {len(matches)} files"
                )
            resolved.append((plugin_name, spec, matches[0]))
    return resolved


def ensure_tool(name: str) -> str:
    path = shutil.which(name)
    if path is None:
        raise SystemExit(f"required tool not found on PATH: {name}")
    return path


def run_with_progress(
    cmd: list[str],
    *,
    label: str,
    action: str,
    stdout=None,
    progress_interval_seconds: float = 5.0,
) -> float:
    started = time.monotonic()
    proc = subprocess.Popen(cmd, stdout=stdout, stderr=subprocess.PIPE, text=True)
    last_reported = started
    while True:
        returncode = proc.poll()
        now = time.monotonic()
        if returncode is not None:
            break
        if now - last_reported >= progress_interval_seconds:
            elapsed = now - started
            print(f"[{label}] {action}: still running ({elapsed:.1f}s elapsed)")
            last_reported = now
        time.sleep(0.2)

    stderr = proc.stderr.read() if proc.stderr is not None else ""
    if returncode != 0:
        raise subprocess.CalledProcessError(returncode, cmd, stderr=stderr)

    elapsed = time.monotonic() - started
    print(f"[{label}] {action}: done in {elapsed:.1f}s")
    return elapsed


def bgzf_sibling_path(path: Path) -> Path:
    if path.name.endswith(".gz"):
        return path.with_name(path.name[:-3] + ".bgz")
    return path.with_suffix(path.suffix + ".bgz")


def recompress_to_bgzf(path: Path, dry_run: bool) -> Path:
    bgzip_bin = ensure_tool("bgzip")
    output_path = bgzf_sibling_path(path)
    print(f"[{path.name}] recompress plain gzip -> {output_path.name}")
    if dry_run:
        return output_path

    with tempfile.NamedTemporaryFile(
        prefix=output_path.name + ".", suffix=".tmp", dir=output_path.parent, delete=False
    ) as temp_handle:
        temp_path = Path(temp_handle.name)
    try:
        started = time.monotonic()
        sink = temp_path.open("wb")
        proc = subprocess.Popen(
            [bgzip_bin, "-c"],
            stdin=subprocess.PIPE,
            stdout=sink,
            stderr=subprocess.PIPE,
        )
        assert proc.stdin is not None
        last_reported = started
        with gzip.open(path, "rb") as source, sink:
            while chunk := source.read(1024 * 1024):
                proc.stdin.write(chunk)
                now = time.monotonic()
                if now - last_reported >= 5.0:
                    elapsed = now - started
                    print(
                        f"[{path.name}] bgzip -> {output_path.name}: still running "
                        f"({elapsed:.1f}s elapsed)"
                    )
                    last_reported = now
        proc.stdin.close()
        stderr = proc.stderr.read().decode() if proc.stderr is not None else ""
        returncode = proc.wait()
        if returncode != 0:
            raise subprocess.CalledProcessError(
                returncode, [bgzip_bin, "-c"], stderr=stderr
            )
        elapsed = time.monotonic() - started
        print(f"[{path.name}] bgzip -> {output_path.name}: done in {elapsed:.1f}s")
        temp_path.replace(output_path)
    finally:
        try:
            temp_path.unlink()
        except FileNotFoundError:
            pass
    return output_path


def tabix_command(path: Path, kind: str, force: bool) -> list[str]:
    cmd = [ensure_tool("tabix")]
    if force:
        cmd.append("-f")
    if kind == "vcf":
        cmd.extend(["-p", "vcf", str(path)])
        return cmd
    if kind == "tsv":
        cmd.extend(["-s", "1", "-b", "2", "-e", "2", "-c", "#", str(path)])
        return cmd
    raise ValueError(f"unsupported index kind: {kind}")


def needs_reindex(path: Path) -> bool:
    index_path = path.with_suffix(path.suffix + ".tbi")
    if not index_path.exists():
        return True
    return index_path.stat().st_mtime < path.stat().st_mtime


def process_input(spec: IndexSpec, path: Path, *, recompress_plain_gzip: bool, force: bool, dry_run: bool) -> None:
    target_path = path
    bgzf = is_bgzf(path)
    if not bgzf:
        if not recompress_plain_gzip:
            print(
                f"[{spec.label}] skip: {path.name} is plain gzip, not BGZF. "
                "Use --recompress-plain-gzip to create a sibling .bgz file first."
            )
            return
        target_path = recompress_to_bgzf(path, dry_run=dry_run)

    if not force and not dry_run and not needs_reindex(target_path):
        print(f"[{spec.label}] up-to-date: {target_path.name}.tbi already exists")
        return

    cmd = tabix_command(target_path, spec.kind, force=True if dry_run else force)
    print(f"[{spec.label}] index: {' '.join(cmd)}")
    if not dry_run:
        run_with_progress(
            cmd,
            label=spec.label,
            action=f"tabix {target_path.name}",
        )


def main() -> int:
    args = parse_args()
    plugin_names = resolve_plugin_names(args.plugins)
    inputs = resolve_inputs(args.plugins_dir, plugin_names)
    total = len(inputs)
    for index, (_, spec, path) in enumerate(inputs, start=1):
        print(f"[{index}/{total}] {spec.label}: {path.name}")
        process_input(
            spec,
            path,
            recompress_plain_gzip=args.recompress_plain_gzip,
            force=args.force,
            dry_run=args.dry_run,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
