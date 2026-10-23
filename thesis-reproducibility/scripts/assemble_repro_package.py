#!/usr/bin/env python3
"""Assemble the thesis reproducibility package (small artifacts only).

Copies the measurement artifacts referenced by the thesis tables into
<target>/thesis-reproducibility/, writes MANIFEST.sha256 over the copied files
and LARGE_ARTIFACTS.tsv for large outputs that are deliberately not included.
Idempotent: re-running refreshes the copies (used again after the October 2026
completion run).

Usage: python3 assemble_repro_package.py <worktree-of-vepyr-diffly>
"""
from __future__ import annotations

import hashlib
import shutil
import sys
from pathlib import Path

D = Path("/Users/lukaszjezapkowicz/workspace/data_vepyr")
T = Path("/Users/lukaszjezapkowicz/Desktop/magisterka/praca/magisterska_jezapkowicz/runs")
E2E = D / "tools/vepyr/e2e-testing/reports"
SCRATCH_FIG = Path("/private/tmp/claude-501/-Users-lukaszjezapkowicz/f509acfb-c3f0-42c9-81e1-ec0bc03e8e1c/scratchpad/fig_plugin_overhead.py")
MAX_FILE = 5_000_000
SKIP_SUFFIX = (".vcf", ".vcf.gz", ".parquet", ".bgz", ".tbi", ".csi", ".gz")

# (destination inside the package, source, glob patterns or None for the whole tree)
COPIES = [
    ("results/annotation-2026-09/chapter5-m5-20260924T061203Z", D / "runs/chapter5-m5-20260924T061203Z",
     ["environment.txt", "series.log", "status.txt", "benchmarks/**/*", "correctness/commands.log",
      "correctness/correctness.log", "correctness/chr1_individual_plugins/benchmark_plugins_summary.*",
      "correctness/chr22_individual_plugins/benchmark_plugins_summary.*"]),
    ("results/annotation-2026-09/5repeats_m5", D / "runs/5repeats_m5", None),
    ("results/annotation-2026-09/vep_5repeats_m5", D / "runs/vep_5repeats_m5", None),
    ("results/annotation-2026-09/vep_chr1_warmup", D / "runs/vep_chr1_warmup", None),
    ("results/annotation-2026-10/plugin_scaling_m5_5repeats_20261013", T / "plugin_scaling_m5_5repeats_20261013", None),
    ("results/concordance-2026-09/e2e-reports", E2E, ["fast_chr1_merged_116_*", "fast_chr22_merged_116_*",
                                                      "fast_chr1_merged_plugins_116_*", "fast_chr22_merged_plugins_116_*"]),
    ("results/vep-full-genome-2026-10/vep_full_repeats_m5_20261023", D / "runs/vep_full_repeats_m5_20261023", None),
    ("results/cache-build-2026-10/parallel_full_cache/queue_20_10_repeats", T / "parallel_full_cache/queue_20_10_repeats", None),
    ("results/cache-build-2026-10/full_cache_20", T / "full_cache_20", None),
    ("results/cache-build-2026-10/full_core_5x", T / "full_core_5x", None),
    ("results/completion-2026-10/m5_completion_20261024", D / "runs/m5_completion_20261024", None),
    ("scripts", D / "tools", ["measure_once.py", "run_5repeats_benchmark.py", "run_vep_5repeats_suite.py",
                              "run_vep_chr1_baseline.py", "run_full_vep_bench.py", "run_vep_full_repeats.sh",
                              "m5_completion_batch.py", "m5_progress.py", "assemble_repro_package.py"]),
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def wanted(p: Path) -> bool:
    return p.is_file() and not p.name.endswith(SKIP_SUFFIX) and p.stat().st_size <= MAX_FILE \
        and "__pycache__" not in p.parts and not p.name.startswith(".")


def main() -> int:
    pkg = Path(sys.argv[1]) / "thesis-reproducibility"
    large = []
    for dest, src, patterns in COPIES:
        if not src.exists():
            print(f"skip (missing): {src}")
            continue
        files = []
        for pat in (patterns or ["**/*"]):
            files += [p for p in src.glob(pat)]
        for p in sorted(set(files)):
            if not p.is_file():
                continue
            if not wanted(p):
                if p.stat().st_size > MAX_FILE or p.name.endswith(SKIP_SUFFIX):
                    large.append((str(p), p.stat().st_size))
                continue
            out = pkg / dest / p.relative_to(src)
            out.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(p, out)
    if SCRATCH_FIG.exists():
        (pkg / "scripts").mkdir(parents=True, exist_ok=True)
        shutil.copy2(SCRATCH_FIG, pkg / "scripts/fig_plugin_overhead.py")
    with (pkg / "LARGE_ARTIFACTS.tsv").open("w") as f:
        f.write("path_on_author_machine\tbytes\n")
        for path, size in sorted(large):
            f.write(f"{path}\t{size}\n")
    files = sorted(p for p in pkg.rglob("*") if p.is_file() and p.name != "MANIFEST.sha256")
    with (pkg / "MANIFEST.sha256").open("w") as f:
        for p in files:
            f.write(f"{sha256(p)}  {p.relative_to(pkg)}\n")
    total = sum(p.stat().st_size for p in files)
    print(f"package: {len(files)} files, {total / 1e6:.1f} MB; large artifacts listed: {len(large)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
