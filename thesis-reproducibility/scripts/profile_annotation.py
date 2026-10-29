#!/usr/bin/env python3
"""CPU accounting and stage profile of VEPYR annotation vs workers (thesis, 2026-10-29).
The SHA-256 of every output body is recorded (outside the timed call), which extends
the workers=1 vs workers=8 determinism check to all worker counts and repeats.

For every (input, workers) pair the annotation runs in a fresh Python process (same
call as measure_once.py: everything, hgvs, reference FASTA, plain VCF output) and
records wall time, user and system CPU time of the process (all threads), peak RSS,
block reads/writes and context switches from getrusage. One untimed warm-up run per
input precedes the series and is reported separately. Timing runs use no profiling
hooks; afterwards one run per workers value is repeated with VEP_PROFILE=1, whose
stderr (per-contig pipeline_profile and vcf_sink_profile lines) is kept.

Usage: profile_annotation.py OUT_DIR
Environment: PYN (.venv-native python with vepyr 0.7.0), CACHE, FASTA, CHR1, FULL.
"""
from __future__ import annotations

import json
import os
import random
import subprocess
import sys
import time
from pathlib import Path

D = Path("/Users/lukaszjezapkowicz/workspace/data_vepyr")
PYN = os.environ.get("PYN", str(D / ".venv-native/bin/python3"))
CACHE = os.environ.get("CACHE", str(D / "cache/116_GRCh38_merged"))
FASTA = os.environ.get("FASTA", str(D / "input/Homo_sapiens.GRCh38.dna.primary_assembly.fa"))
INPUTS = {
    "chr1": os.environ.get("CHR1", str(D / "runs/m5_completion_20261024/inputs/HG002_normalized_chr1.vcf.gz")),
    "full": os.environ.get("FULL", str(D / "input/HG002_normalized.vcf.gz")),
}
WORKERS = [1, 2, 4, 8, 10]
REPEATS = {"chr1": 3, "full": 1}

CHILD = r'''
import json, resource, sys, time, vepyr
vcf, out, workers, cache, fasta = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5]
r0 = resource.getrusage(resource.RUSAGE_SELF)
t0 = time.perf_counter()
vepyr.annotate(vcf=vcf, cache_dir=cache, everything=True, hgvs=True, reference_fasta=fasta,
               workers=workers, output_vcf=out, show_progress=False)
wall = time.perf_counter() - t0
r1 = resource.getrusage(resource.RUSAGE_SELF)
d = lambda a: getattr(r1, a) - getattr(r0, a)
print(json.dumps({"wall_s": wall, "user_s": d("ru_utime"), "sys_s": d("ru_stime"),
                  "maxrss_bytes": r1.ru_maxrss, "inblock": d("ru_inblock"), "oublock": d("ru_oublock"),
                  "nvcsw": d("ru_nvcsw"), "nivcsw": d("ru_nivcsw")}))
'''


def body_sha256(path: Path) -> tuple[str, int]:
    """SHA-256 of the VCF body (header lines skipped), computed after timing."""
    import hashlib
    h, n = hashlib.sha256(), 0
    with path.open("rb") as f:
        for line in f:
            if not line.startswith(b"#"):
                h.update(line)
                n += 1
    return h.hexdigest(), n


def therm() -> str:
    return subprocess.run(["pmset", "-g", "therm"], capture_output=True, text=True).stdout.strip()


def run(out_dir: Path, name: str, workers: int, tag: str, profile: bool) -> dict:
    out_vcf = out_dir / "tmp_output.vcf"
    env = dict(os.environ)
    if profile:
        env["VEP_PROFILE"] = "1"
    p = subprocess.run([PYN, "-c", CHILD, INPUTS[name], str(out_vcf), str(workers), CACHE, FASTA],
                       capture_output=True, text=True, env=env)
    rec = {"input": name, "workers": workers, "tag": tag, "profile": profile, "returncode": p.returncode,
           "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    if p.returncode == 0:
        rec.update(json.loads(p.stdout.strip().splitlines()[-1]))
        rec["output_bytes"] = out_vcf.stat().st_size
        rec["body_sha256"], rec["records"] = body_sha256(out_vcf)
        rec["cpu_s"] = rec["user_s"] + rec["sys_s"]
        rec["avg_busy_cores"] = rec["cpu_s"] / rec["wall_s"]
    if profile or p.returncode != 0:
        (out_dir / f"stderr_{name}_w{workers}_{tag}.txt").write_text(p.stderr)
    out_vcf.unlink(missing_ok=True)
    return rec


def main() -> int:
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    jl = out / "runs.jsonl"
    rng = random.Random(20261029)
    (out / "therm_start.txt").write_text(therm())

    def log(rec):
        with jl.open("a") as f:
            f.write(json.dumps(rec) + "\n")
        print(json.dumps({k: rec.get(k) for k in ("input", "workers", "tag", "wall_s", "cpu_s", "avg_busy_cores")}), flush=True)

    for name in ("chr1", "full"):
        log(run(out, name, 4, "warmup", False))
        for rep in range(1, REPEATS[name] + 1):
            order = WORKERS[:]
            rng.shuffle(order)
            for w in order:
                log(run(out, name, w, f"r{rep}", False))
                time.sleep(3)
        if name == "chr1":
            for w in WORKERS:
                log(run(out, name, w, "profile", True))
                time.sleep(3)
    (out / "therm_end.txt").write_text(therm())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
