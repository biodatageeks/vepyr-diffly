#!/usr/bin/env python3
"""VEP counterpart of the independent-process test (thesis item D21, 2026-10-30).

As step S6 of m5_completion_batch.py did for VEPYR (tab:multiproc-read): n independent
annotations of the same HG002 chr22 baseline started together, here n `docker run`
commands of VEP 116 without --fork (options of the timing series), all in the Colima VM
and reading the same cache. Group time = from starting the first command to the end of
the last (perf_counter, container start included), for n = 1, 2, 4, 8 and five repeats in
the order of S6 (n ascending within a repeat). Reported per n: mean group time, T_n/T_1
and the throughput n / (T_n/T_1); per container the /usr/bin/time real time.

Usage: vep_multiproc.py OUT_DIR [REPEATS]   (Colima running)
"""
from __future__ import annotations

import json
import statistics
import subprocess
import sys
import time
from pathlib import Path

D = Path("/Users/lukaszjezapkowicz/workspace/data_vepyr")
IMAGE = "ensemblorg/ensembl-vep:release_116.0"
CHR22 = D / "runs/m5_completion_20261024/inputs/HG002_normalized_chr22.vcf.gz"
FASTA = "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
NS = (1, 2, 4, 8)


def cmd(out: Path, i: int) -> list[str]:
    return ["/usr/bin/time", "-p", "-o", str(out / f"time_{i}.txt"), "docker", "run", "--rm", "--user", "501:20",
            "--env", "HOME=/tmp",
            "-v", f"{D}/homo_sapiens_merged/116_GRCh38:/opt/vep/.vep/homo_sapiens_merged/116_GRCh38:ro",
            "-v", f"{D}/input:/ref:ro", "-v", f"{CHR22.parent}:/input:ro", "-v", f"{out}:/output", IMAGE, "vep",
            "--dir", "/opt/vep/.vep", "--cache", "--merged", "--offline", "--assembly", "GRCh38",
            "--input_file", f"/input/{CHR22.name}", "--output_file", f"/output/vep_{i}.vcf", "--vcf",
            "--force_overwrite", "--no_stats", "--everything", "--hgvs", "--fasta", f"/ref/{FASTA}"]


def main() -> int:
    out = Path(sys.argv[1]).resolve()
    reps = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    out.mkdir(parents=True, exist_ok=True)
    raw = out / "raw_runs.jsonl"
    done = {(r["n"], r["repeat"]) for r in map(json.loads, raw.read_text().splitlines())} if raw.exists() else set()
    for rep in range(1, reps + 1):
        for n in NS:
            if (n, rep) in done:
                continue
            run = out / f"n{n}_r{rep}"
            run.mkdir(exist_ok=True)
            t0 = time.perf_counter()
            procs = [subprocess.Popen(cmd(run, i), stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
                     for i in range(n)]
            errs = [p.communicate()[1] for p in procs]
            group = time.perf_counter() - t0
            if any(p.returncode for p in procs):
                raise SystemExit(f"VEP failed: {[e[-500:] for e in errs]}")
            per = [float(next(l.split()[1] for l in (run / f"time_{i}.txt").read_text().splitlines()
                              if l.startswith("real"))) for i in range(n)]
            records = [sum(1 for l in open(run / f"vep_{i}.vcf") if not l.startswith("#")) for i in range(n)]
            for i in range(n):
                (run / f"vep_{i}.vcf").unlink()
            r = {"n": n, "repeat": rep, "group_wall_s": group, "container_real_s": per, "records": records,
                 "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
            with raw.open("a") as f:
                f.write(json.dumps(r) + "\n")
            print(json.dumps(r), flush=True)
            time.sleep(3)
    runs = list(map(json.loads, raw.read_text().splitlines()))
    t1 = statistics.mean(r["group_wall_s"] for r in runs if r["n"] == 1)
    summ = {}
    for n in NS:
        g = [r["group_wall_s"] for r in runs if r["n"] == n]
        summ[n] = {"repeats": len(g), "mean_group_s": statistics.mean(g), "min_s": min(g), "max_s": max(g),
                   "sd_s": statistics.stdev(g) if len(g) > 1 else None,
                   "relative_group_time": statistics.mean(g) / t1, "throughput": n / (statistics.mean(g) / t1)}
    (out / "summary.json").write_text(json.dumps(summ, indent=1) + "\n")
    print(json.dumps(summ, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
