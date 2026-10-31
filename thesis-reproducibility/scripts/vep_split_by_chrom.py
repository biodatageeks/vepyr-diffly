#!/usr/bin/env python3
"""VEP split by chromosome into independent processes (thesis item D10, 2026-10-31).

Production practice for VEP on a whole genome: one VEP process per chromosome, several at a time.
The HG002 chr1-22 input is cut into per-contig files with tabix -h (input lines unchanged), and
PARALLEL containers run at a time, largest chromosome first (record counts from bcftools index -s).
Every container runs the command of the interleaved series (one process, no --fork, merged cache,
--everything --hgvs --fasta) and reports its cgroup CPU time. Afterwards the output bodies are
hashed in input contig order; the result must equal the body of the single-process VEP run on the
whole input (5fee6224..., full-genome repeat 2).

Usage: vep_split_by_chrom.py OUT_DIR [PARALLEL, default 10]
"""
import concurrent.futures as cf
import hashlib, json, subprocess, sys, time
from pathlib import Path

D = Path("/Users/lukaszjezapkowicz/workspace/data_vepyr")
INPUT = D / "input/HG002_normalized.vcf.gz"
IMAGE = "ensemblorg/ensembl-vep:release_116.0"
EXPECTED = "5fee6224bf135d80075712141cc06194e0730fa253ff8b131a9d25d661647ee9"
WRAP = 'vep "$@"; rc=$?; echo CPUSTAT; cat /sys/fs/cgroup/cpu.stat; exit $rc'


def run_one(out: Path, contig: str, t0: float) -> dict:
    start = time.perf_counter() - t0
    cmd = ["docker", "run", "--rm", "--user", "501:20", "--env", "HOME=/tmp",
           "-v", f"{D}/homo_sapiens_merged/116_GRCh38:/opt/vep/.vep/homo_sapiens_merged/116_GRCh38:ro",
           "-v", f"{D}/input:/ref:ro", "-v", f"{out / 'shards'}:/input:ro", "-v", f"{out / 'vep'}:/output",
           IMAGE, "sh", "-c", WRAP, "vep",
           "--dir", "/opt/vep/.vep", "--cache", "--merged", "--offline", "--assembly", "GRCh38",
           "--input_file", f"/input/{contig}.vcf.gz", "--output_file", f"/output/{contig}.vcf", "--vcf",
           "--force_overwrite", "--no_stats", "--everything", "--hgvs", "--fasta", "/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa"]
    p = subprocess.run(cmd, capture_output=True, text=True)
    end = time.perf_counter() - t0
    stat = dict(l.split() for l in p.stdout.split("CPUSTAT\n", 1)[-1].splitlines() if len(l.split()) == 2)
    r = {"contig": contig, "returncode": p.returncode, "start_s": start, "end_s": end, "seconds": end - start,
         "cgroup_cpu_s": int(stat.get("usage_usec", 0)) / 1e6}
    if p.returncode:
        r["stderr_tail"] = p.stderr[-1500:]
    print(json.dumps(r), flush=True)
    return r


def main():
    out = Path(sys.argv[1]).resolve()
    parallel = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    (out / "shards").mkdir(parents=True, exist_ok=True)
    (out / "vep").mkdir(exist_ok=True)
    counts = [(l.split("\t")[0], int(l.split("\t")[2])) for l in
              subprocess.run(["bcftools", "index", "-s", str(INPUT)], capture_output=True, text=True, check=True).stdout.splitlines()]
    order = [c for c, _ in counts]
    for c in order:
        with open(out / "shards" / f"{c}.vcf.gz", "wb") as f:
            t = subprocess.Popen(["tabix", "-h", str(INPUT), c], stdout=subprocess.PIPE)
            subprocess.run(["bgzip", "-c"], stdin=t.stdout, stdout=f, check=True)
            t.wait()
    queue = [c for c, _ in sorted(counts, key=lambda x: -x[1])]
    t0 = time.perf_counter()
    started = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    with cf.ThreadPoolExecutor(parallel) as ex:
        runs = list(ex.map(lambda c: run_one(out, c, t0), queue))
    makespan = time.perf_counter() - t0
    h, n = hashlib.sha256(), 0
    for c in order:
        f = out / "vep" / f"{c}.vcf"
        if f.exists():
            with open(f, "rb") as fh:
                for line in fh:
                    if not line.startswith(b"#"):
                        h.update(line); n += 1
    body = h.hexdigest()
    res = {"parallel": parallel, "started_utc": started, "makespan_s": makespan, "records": n,
           "body_sha256": body, "equal_single_process_vep": body == EXPECTED,
           "sum_seconds": sum(r["seconds"] for r in runs), "sum_cgroup_cpu_s": sum(r["cgroup_cpu_s"] for r in runs),
           "failed": [r["contig"] for r in runs if r["returncode"]], "runs": runs, "record_counts": dict(counts)}
    (out / "result.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps({k: v for k, v in res.items() if k not in ("runs", "record_counts")}))
    for c in order:
        (out / "vep" / f"{c}.vcf").unlink(missing_ok=True)


if __name__ == "__main__":
    main()
