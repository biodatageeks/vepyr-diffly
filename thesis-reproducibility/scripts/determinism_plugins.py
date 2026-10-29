#!/usr/bin/env python3
"""Determinism of VEPYR with all five plugins over the workers setting (thesis item B16).

The thesis shows one body checksum for workers 1-10 in the baseline scenario (chr1 and
chr1-22). Here the HG002 chr22 and chr1 slices are annotated with all five plugins
(plugin cache of the plugin series, everything, hgvs, FASTA) for workers 1-10, REPEATS
times in shuffled order; every body SHA-256 must equal the VEP 116 five-plugin reference
(results/references-2026-10-25). The Linux side of B16 comes from vepyr_in_vm.py (body
equal to the native output).

Usage: determinism_plugins.py OUT_DIR [REPEATS, default 2]
Run with the .venv-native interpreter.
"""
from __future__ import annotations

import json
import random
import subprocess
import sys
import time
from pathlib import Path

D = Path("/Users/lukaszjezapkowicz/workspace/data_vepyr")
PKG = Path(__file__).resolve().parents[1]
INPUT = {"chr22": D / "runs/m5_completion_20261024/inputs/HG002_normalized_chr22.vcf.gz",
         "chr1": D / "runs/m5_completion_20261024/inputs/HG002_normalized_chr1.vcf.gz"}
CHILD = r'''
import hashlib, json, sys, time, vepyr
t0 = time.perf_counter()
vepyr.annotate(vcf=sys.argv[1], cache_dir=sys.argv[2], everything=True, hgvs=True, reference_fasta=sys.argv[3],
               plugin_cache_root=sys.argv[4], workers=int(sys.argv[5]), output_vcf=sys.argv[6], show_progress=False)
dt = time.perf_counter() - t0
h = hashlib.sha256(); n = 0
with open(sys.argv[6], "rb") as f:
    for line in f:
        if not line.startswith(b"#"):
            h.update(line); n += 1
print(json.dumps({"seconds": dt, "body_sha256": h.hexdigest(), "records": n}))
'''


def main() -> int:
    out = Path(sys.argv[1]).resolve()
    reps = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    out.mkdir(parents=True, exist_ok=True)
    ref = {c: (PKG / f"results/references-2026-10-25/{c}_reference_body.sha256").read_text().split()[0] for c in INPUT}
    rng = random.Random(20261030)
    cells = [(c, w, r) for r in range(1, reps + 1) for c in INPUT for w in range(1, 11)]
    for r in range(1, reps + 1):
        block = [x for x in cells if x[2] == r]
        rng.shuffle(block)
        for chrom, w, rep in block:
            vcf = out / "out.vcf"
            p = subprocess.run([sys.executable, "-c", CHILD, str(INPUT[chrom]), str(D / "cache/116_GRCh38_merged"),
                                str(D / "input/Homo_sapiens.GRCh38.dna.primary_assembly.fa"),
                                str(D / "cache/plugin_cache_116"), str(w), str(vcf)], capture_output=True, text=True)
            vcf.unlink(missing_ok=True)
            row = {"chrom": chrom, "workers": w, "repeat": rep, "returncode": p.returncode,
                   "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
            if p.returncode == 0:
                row.update(json.loads(p.stdout.strip().splitlines()[-1]))
                row["equal_vep_reference"] = row["body_sha256"] == ref[chrom]
            else:
                row["stderr_tail"] = p.stderr[-1000:]
            with (out / "runs.jsonl").open("a") as f:
                f.write(json.dumps(row) + "\n")
            print(json.dumps({k: v for k, v in row.items() if k != "stderr_tail"}), flush=True)
    rows = [json.loads(l) for l in (out / "runs.jsonl").read_text().splitlines()]
    summ = {c: {"runs": sum(r["chrom"] == c for r in rows),
                "distinct_bodies": sorted({r.get("body_sha256") for r in rows if r["chrom"] == c}),
                "all_equal_vep_reference": all(r.get("equal_vep_reference") for r in rows if r["chrom"] == c)}
            for c in INPUT}
    (out / "summary.json").write_text(json.dumps(summ, indent=1) + "\n")
    print(json.dumps(summ, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
