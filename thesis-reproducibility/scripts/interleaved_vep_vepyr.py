#!/usr/bin/env python3
"""Same-session VEP vs VEPYR timing (thesis, D5): the chr22 (and optionally chr1)
baseline annotated alternately by VEP 116 (one process, no --fork, docker run timed
with /usr/bin/time -p as in the September series) and VEPYR (workers=1, timed with
time.perf_counter around vepyr.annotate as in measure_once.py), REPEATS times, so
both tools see the same machine state and the same core cache build. CPU time (D6): for VEP
the container's cgroup cpu.stat read after vep exits (all VEP processes; the command runs
under `sh -c` for that), for VEPYR getrusage(RUSAGE_SELF) around vepyr.annotate. The body SHA-256
(lines not starting with #) of every output is recorded after the timing, as a sanity check.

Usage: interleaved_vep_vepyr.py OUT_DIR CHROM [REPEATS]
"""
import hashlib, json, subprocess, sys, time
from pathlib import Path

D = Path("/Users/lukaszjezapkowicz/workspace/data_vepyr")
IMAGE = "ensemblorg/ensembl-vep:release_116.0"
INPUT = {"chr22": D / "runs/m5_completion_20261024/inputs/HG002_normalized_chr22.vcf.gz",
         "chr1": D / "runs/m5_completion_20261024/inputs/HG002_normalized_chr1.vcf.gz"}
PYN = D / ".venv-native/bin/python3"
CHILD = r'''
import json, resource, sys, time, vepyr
u0 = resource.getrusage(resource.RUSAGE_SELF)
t0 = time.perf_counter()
vepyr.annotate(vcf=sys.argv[1], cache_dir=sys.argv[2], everything=True, hgvs=True, reference_fasta=sys.argv[3],
               workers=1, output_vcf=sys.argv[4], show_progress=False)
dt = time.perf_counter() - t0
u1 = resource.getrusage(resource.RUSAGE_SELF)
print(json.dumps({"annotation_seconds": dt, "user_s": u1.ru_utime - u0.ru_utime, "sys_s": u1.ru_stime - u0.ru_stime}))
'''
# vep runs unchanged; afterwards the container's own cgroup reports the CPU time of every process in it
WRAP = 'vep "$@"; rc=$?; echo CPUSTAT; cat /sys/fs/cgroup/cpu.stat; exit $rc'


def body_sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for line in f:
            if not line.startswith(b"#"):
                h.update(line)
    return h.hexdigest()


def vep(chrom, out):
    inp = INPUT[chrom]
    t = out / "vep.time.txt"
    cmd = ["/usr/bin/time", "-p", "-o", str(t), "docker", "run", "--rm", "--user", "501:20", "--env", "HOME=/tmp",
           "-v", f"{D}/homo_sapiens_merged/116_GRCh38:/opt/vep/.vep/homo_sapiens_merged/116_GRCh38:ro",
           "-v", f"{D}/input:/ref:ro", "-v", f"{inp.parent}:/input:ro", "-v", f"{out}:/output", IMAGE, "sh", "-c", WRAP, "vep",
           "--dir", "/opt/vep/.vep", "--cache", "--merged", "--offline", "--assembly", "GRCh38",
           "--input_file", f"/input/{inp.name}", "--output_file", "/output/vep.vcf", "--vcf", "--force_overwrite",
           "--no_stats", "--everything", "--hgvs", "--fasta", "/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa"]
    p = subprocess.run(cmd, capture_output=True, text=True)
    real = float(next(l.split()[1] for l in t.read_text().splitlines() if l.startswith("real")))
    body = body_sha256(out / "vep.vcf") if p.returncode == 0 else None
    (out / "vep.vcf").unlink(missing_ok=True)
    stat = dict(l.split() for l in p.stdout.split("CPUSTAT\n", 1)[-1].splitlines() if len(l.split()) == 2)
    cpu = {k: int(stat[f"{k}_usec"]) / 1e6 for k in ("usage", "user", "system") if f"{k}_usec" in stat}
    return {"tool": "vep", "returncode": p.returncode, "seconds": real, "cgroup_cpu_s": cpu, "body_sha256": body}


def vepyr(chrom, out):
    p = subprocess.run([str(PYN), "-c", CHILD, str(INPUT[chrom]), str(D / "cache/116_GRCh38_merged"),
                        str(D / "input/Homo_sapiens.GRCh38.dna.primary_assembly.fa"), str(out / "vepyr.vcf")],
                       capture_output=True, text=True)
    body = body_sha256(out / "vepyr.vcf") if p.returncode == 0 else None
    (out / "vepyr.vcf").unlink(missing_ok=True)
    r = json.loads(p.stdout.strip().splitlines()[-1]) if p.returncode == 0 else {}
    return {"tool": "vepyr", "returncode": p.returncode, "seconds": r.get("annotation_seconds"),
            "user_s": r.get("user_s"), "sys_s": r.get("sys_s"), "body_sha256": body}


def main():
    out, chrom = Path(sys.argv[1]), sys.argv[2]
    reps = int(sys.argv[3]) if len(sys.argv) > 3 else 5
    out.mkdir(parents=True, exist_ok=True)
    for rep in range(1, reps + 1):
        order = [vep, vepyr] if rep % 2 else [vepyr, vep]
        for fn in order:
            r = fn(chrom, out); r.update(chrom=chrom, repeat=rep, utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
            with (out / "runs.jsonl").open("a") as f:
                f.write(json.dumps(r) + "\n")
            print(json.dumps(r), flush=True)
            if r["returncode"] != 0:
                return 1
            time.sleep(5)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
