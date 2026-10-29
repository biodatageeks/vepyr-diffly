#!/usr/bin/env python3
"""Same-session VEP vs VEPYR timing (thesis, D5): the chr22 (and optionally chr1)
baseline annotated alternately by VEP 116 (one process, no --fork, docker run timed
with /usr/bin/time -p as in the September series) and VEPYR (workers=1, timed with
time.perf_counter around vepyr.annotate as in measure_once.py), REPEATS times, so
both tools see the same machine state and the same core cache build.

Usage: interleaved_vep_vepyr.py OUT_DIR CHROM [REPEATS]
"""
import json, subprocess, sys, time
from pathlib import Path

D = Path("/Users/lukaszjezapkowicz/workspace/data_vepyr")
IMAGE = "ensemblorg/ensembl-vep:release_116.0"
INPUT = {"chr22": D / "runs/m5_completion_20261024/inputs/HG002_normalized_chr22.vcf.gz",
         "chr1": D / "runs/m5_completion_20261024/inputs/HG002_normalized_chr1.vcf.gz"}
PYN = D / ".venv-native/bin/python3"
CHILD = r'''
import json, sys, time, vepyr
t0 = time.perf_counter()
vepyr.annotate(vcf=sys.argv[1], cache_dir=sys.argv[2], everything=True, hgvs=True, reference_fasta=sys.argv[3],
               workers=1, output_vcf=sys.argv[4], show_progress=False)
print(json.dumps({"annotation_seconds": time.perf_counter() - t0}))
'''


def vep(chrom, out):
    inp = INPUT[chrom]
    t = out / "vep.time.txt"
    cmd = ["/usr/bin/time", "-p", "-o", str(t), "docker", "run", "--rm", "--user", "501:20", "--env", "HOME=/tmp",
           "-v", f"{D}/homo_sapiens_merged/116_GRCh38:/opt/vep/.vep/homo_sapiens_merged/116_GRCh38:ro",
           "-v", f"{D}/input:/ref:ro", "-v", f"{inp.parent}:/input:ro", "-v", f"{out}:/output", IMAGE, "vep",
           "--dir", "/opt/vep/.vep", "--cache", "--merged", "--offline", "--assembly", "GRCh38",
           "--input_file", f"/input/{inp.name}", "--output_file", "/output/vep.vcf", "--vcf", "--force_overwrite",
           "--no_stats", "--everything", "--hgvs", "--fasta", "/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa"]
    p = subprocess.run(cmd, capture_output=True, text=True)
    real = float(next(l.split()[1] for l in t.read_text().splitlines() if l.startswith("real")))
    (out / "vep.vcf").unlink(missing_ok=True)
    return {"tool": "vep", "returncode": p.returncode, "seconds": real}


def vepyr(chrom, out):
    p = subprocess.run([str(PYN), "-c", CHILD, str(INPUT[chrom]), str(D / "cache/116_GRCh38_merged"),
                        str(D / "input/Homo_sapiens.GRCh38.dna.primary_assembly.fa"), str(out / "vepyr.vcf")],
                       capture_output=True, text=True)
    (out / "vepyr.vcf").unlink(missing_ok=True)
    sec = json.loads(p.stdout.strip().splitlines()[-1])["annotation_seconds"] if p.returncode == 0 else None
    return {"tool": "vepyr", "returncode": p.returncode, "seconds": sec}


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
