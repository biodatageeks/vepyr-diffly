#!/usr/bin/env python3
"""Mean peak RSS (GiB) of the VEPYR plugin series on chr1/chr22, as reported in the thesis.

Repeat 1 of these series comes from `vepyr-diffly benchmark-plugins` (peak_rss_mb, MiB of the
sampled process tree); summary_5repeats.json stored it as MiB/1000 under `max_rss_gb`.
Repeats 2-5 come from measure_once.py (ru_maxrss bytes). Both are converted to GiB here.

Usage: python3 plugin_rss_gib.py [<5repeats_m5 dir> <chapter5-m5 correctness dir>]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

D = Path("/Users/lukaszjezapkowicz/workspace/data_vepyr/runs")
REPEATS = Path(sys.argv[1]) if len(sys.argv) > 1 else D / "5repeats_m5"
REPEAT1 = Path(sys.argv[2]) if len(sys.argv) > 2 else D / "chapter5-m5-20260924T061203Z/correctness"


def main() -> None:
    summary = json.loads((REPEATS / "summary_5repeats.json").read_text())
    for chrom in ("chr1", "chr22"):
        bench = json.loads((REPEAT1 / f"{chrom}_individual_plugins/benchmark_plugins_summary.json").read_text())
        first = {(r["plugin"] or "baseline"): r["peak_rss_mb"] / 1024 for r in bench["results"]}
        for cfg, runs in summary[f"{chrom}_plugins"].items():
            rest = [r["max_rss_bytes"] / 1024**3 for r in runs if "max_rss_bytes" in r]
            values = [first[cfg], *rest]
            print(f"{chrom} {cfg:14s} n={len(values)} mean_peak_rss_gib={sum(values) / len(values):.2f}")


if __name__ == "__main__":
    main()
