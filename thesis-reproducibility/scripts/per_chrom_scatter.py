#!/usr/bin/env python3
"""Scatter-gather baseline (thesis, D10): annotate chr1-22 as independent per-chromosome
jobs on one machine, P jobs at a time, longest chromosome first.

Modes:
  vep    one `docker run` of VEP 116 per chromosome (same options as the full-input
         series, no --fork), P concurrent containers in the Colima VM;
  vepyr  one Python process per chromosome running vepyr.annotate(workers=1).
For each repeat: makespan (start of the first job to the end of the last), per-job
wall time, and the SHA-256 of the per-chromosome bodies concatenated in chromosome
order, which must equal the body of the single full-input run (VEP repeat 2:
5fee6224...).

Usage: per_chrom_scatter.py MODE OUT_DIR [REPEATS] [P]
"""
from __future__ import annotations

import concurrent.futures as cf
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

D = Path("/Users/lukaszjezapkowicz/workspace/data_vepyr")
IMAGE = "ensemblorg/ensembl-vep:release_116.0"
VEP_CACHE = D / "homo_sapiens_merged/116_GRCh38"
FASTA_DIR = D / "input"
FASTA = "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
FULL = D / "input/HG002_normalized.vcf.gz"
PYN = D / ".venv-native/bin/python3"
CORE = D / "cache/116_GRCh38_merged"
EXPECTED_BODY = "5fee6224bf135d80075712141cc06194e0730fa253ff8b131a9d25d661647ee9"
CHROMS = [f"chr{i}" for i in range(1, 23)]


def sh(cmd, **kw):
    return subprocess.run(cmd, check=True, capture_output=True, text=True, **kw)


def slice_inputs(dst: Path) -> dict:
    dst.mkdir(parents=True, exist_ok=True)
    counts = {}
    for c in CHROMS:
        out = dst / f"{c}.vcf.gz"
        if not out.exists():
            sh(["bcftools", "view", "-r", c, "-Oz", "-o", str(out), str(FULL)])
            sh(["bcftools", "index", "-f", "-t", str(out)])
        counts[c] = int(subprocess.run(f"bcftools view -H {out} | wc -l", shell=True, capture_output=True, text=True).stdout)
    return counts


def vep_job(c: str, inputs: Path, out: Path) -> dict:
    t = out / f"{c}.time.txt"
    cmd = ["/usr/bin/time", "-p", "-o", str(t), "docker", "run", "--rm", "--user", "501:20", "--env", "HOME=/tmp",
           "-v", f"{VEP_CACHE}:/opt/vep/.vep/homo_sapiens_merged/116_GRCh38:ro",
           "-v", f"{FASTA_DIR}:/ref:ro", "-v", f"{inputs}:/input:ro", "-v", f"{out}:/output",
           IMAGE, "vep", "--dir", "/opt/vep/.vep", "--cache", "--merged", "--offline", "--assembly", "GRCh38",
           "--input_file", f"/input/{c}.vcf.gz", "--output_file", f"/output/{c}.vcf", "--vcf", "--force_overwrite",
           "--no_stats", "--everything", "--hgvs", "--fasta", f"/ref/{FASTA}"]
    t0 = time.perf_counter(); p = subprocess.run(cmd, capture_output=True, text=True); t1 = time.perf_counter()
    return {"chrom": c, "start": t0, "end": t1, "wall_s": t1 - t0, "returncode": p.returncode, "stderr_tail": p.stderr[-500:]}


VEPYR_CHILD = r'''
import sys, vepyr
vepyr.annotate(vcf=sys.argv[1], cache_dir=sys.argv[2], everything=True, hgvs=True, reference_fasta=sys.argv[3],
               workers=1, output_vcf=sys.argv[4], show_progress=False)
'''


def vepyr_job(c: str, inputs: Path, out: Path) -> dict:
    cmd = [str(PYN), "-c", VEPYR_CHILD, str(inputs / f"{c}.vcf.gz"), str(CORE), str(FASTA_DIR / FASTA), str(out / f"{c}.vcf")]
    t0 = time.perf_counter(); p = subprocess.run(cmd, capture_output=True, text=True); t1 = time.perf_counter()
    return {"chrom": c, "start": t0, "end": t1, "wall_s": t1 - t0, "returncode": p.returncode, "stderr_tail": p.stderr[-500:]}


def concat_sha(out: Path) -> tuple[str, int]:
    h, n = hashlib.sha256(), 0
    for c in CHROMS:
        with (out / f"{c}.vcf").open("rb") as f:
            for line in f:
                if not line.startswith(b"#"):
                    h.update(line); n += 1
    return h.hexdigest(), n


def main() -> int:
    mode, root = sys.argv[1], Path(sys.argv[2])
    repeats = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    p = int(sys.argv[4]) if len(sys.argv) > 4 else 10
    root.mkdir(parents=True, exist_ok=True)
    inputs = D / "runs/scatter_inputs_chr1_22"
    counts = slice_inputs(inputs)
    order = sorted(CHROMS, key=lambda c: -counts[c])
    job = vep_job if mode == "vep" else vepyr_job
    for rep in range(1, repeats + 1):
        out = root / f"{mode}_r{rep}"; out.mkdir(exist_ok=True)
        t0 = time.perf_counter()
        with cf.ThreadPoolExecutor(max_workers=p) as ex:
            res = list(ex.map(lambda c: job(c, inputs, out), order))
        makespan = time.perf_counter() - t0
        ok = all(r["returncode"] == 0 for r in res)
        sha, n = concat_sha(out) if ok else (None, 0)
        rec = {"mode": mode, "repeat": rep, "processes": p, "makespan_s": makespan, "ok": ok, "body_sha256": sha,
               "records": n, "matches_full_vep": sha == EXPECTED_BODY,
               "jobs": [{k: (v - t0 if k in ("start", "end") else v) for k, v in r.items()} for r in res]}
        with (root / "runs.jsonl").open("a") as f:
            f.write(json.dumps(rec) + "\n")
        print(json.dumps({k: rec[k] for k in ("mode", "repeat", "makespan_s", "ok", "records", "matches_full_vep")}), flush=True)
        for c in CHROMS:
            (out / f"{c}.vcf").unlink(missing_ok=True)
        if not ok:
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
