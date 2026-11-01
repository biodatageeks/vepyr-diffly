#!/usr/bin/env python3
"""Annotation split by contig as a Kubernetes Indexed Job (thesis items A1/A2, 2026-10-31).

The plan of distributed_annotate.py (LPT over contigs by record count) is run as one Indexed Job
with N pods. Pod i takes the contigs of node i of the plan:
  initContainer  copies its input shards, the cache fragments of its contigs (+ every entity's
                 chrom_manifest.json) and the FASTA with its index from the shared store to the
                 node's local disk (/work/<index>), and records bytes and seconds;
  main container one vepyr.annotate call per contig (bgzf output) on the local copies, then copies
                 the outputs back to the store (parts/) and records per-contig time, CPU, peak RSS,
                 records and body SHA-256.
The driver cuts the input with tabix -h, applies the Job, waits for it, merges the parts in input
order with bcftools concat --naive(-force) and checks every contig exactly once, the record count
and the body SHA-256 against the single-call output.

On the author's machine the cluster is k3s in the Colima VM (one node): the store is a macOS
directory reached over virtiofs (standing in for object storage). The local work directory is a
hostPath on the VM data disk (WORK_HOST), not an emptyDir: kubelet keeps emptyDir on the 19 GiB VM
system disk, too small for the cache of chr1-22 (~33 GB). The driver creates one directory per pod
before the Job and removes them afterwards.
All pods share one machine, so this tests orchestration, staging and correctness, not distribution.

Usage: k8s_annotate.py run INPUT_VCF_GZ N OUT_DIR [WORKERS, default 1]
"""
from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import distributed_annotate as da  # noqa: E402

IMAGE = "vepyr-vm:479d8f4"
FASTA = "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
WORK_HOST = "/var/lib/rancher/vepyr-work"
INIT = r'''
import json, os, shutil, sys, time
i = int(os.environ["JOB_COMPLETION_INDEX"]); p = json.load(open("/store/plan.json")); cs = p["nodes"][i]["contigs"]
W = f"/work/{i}"; t0 = time.time(); res = {}
def cp(src, dst):
    os.makedirs(os.path.dirname(dst), exist_ok=True); shutil.copyfile(src, dst); return os.path.getsize(dst)
for part, pairs in (
    ("input", [(f"/store/shards/{c}.vcf.gz{s}", f"{W}/in/{c}.vcf.gz{s}") for c in cs for s in ("", ".tbi")]),
    ("cache", [(f"/cache-src/{e}/{n}", f"{W}/cache/{e}/{n}") for e in sorted(os.listdir("/cache-src")) if os.path.isdir(f"/cache-src/{e}")
               for n in ["chrom_manifest.json"] + [f"{c}.parquet" for c in cs] if os.path.exists(f"/cache-src/{e}/{n}")]),
    ("fasta", [(f"/ref-src/{f}", f"{W}/ref/{f}") for f in ("%s", "%s.fai")])):
    t = time.time(); b = sum(cp(s, d) for s, d in pairs)
    res[part] = {"bytes": b, "seconds": time.time() - t, "files": len(pairs)}
os.makedirs(f"{W}/out", exist_ok=True)
json.dump({"index": i, "contigs": cs, "stage": res, "stage_s": time.time() - t0, "start_unix": t0, "end_unix": time.time()},
          open(f"/store/pods/{i}.stage.json", "w"), indent=1)
''' % (FASTA, FASTA)
MAIN = r'''
import gzip, hashlib, json, os, resource, shutil, time, vepyr
i = int(os.environ["JOB_COMPLETION_INDEX"]); p = json.load(open("/store/plan.json")); cs = p["nodes"][i]["contigs"]
w = int(os.environ["WORKERS"]); W = f"/work/{i}"; t00 = time.time(); rows = []
for c in cs:
    out = f"{W}/out/{c}.vcf.gz"; u0 = resource.getrusage(resource.RUSAGE_SELF); t0 = time.perf_counter()
    vepyr.annotate(vcf=f"{W}/in/{c}.vcf.gz", cache_dir=f"{W}/cache", everything=True, hgvs=True,
                   reference_fasta=f"{W}/ref/%s", workers=w, output_vcf=out, show_progress=False)
    dt = time.perf_counter() - t0; u1 = resource.getrusage(resource.RUSAGE_SELF)
    h = hashlib.sha256(); n = 0
    with gzip.open(out, "rb") as f:
        for line in f:
            if not line.startswith(b"#"):
                h.update(line); n += 1
    rows.append({"contig": c, "annotation_seconds": dt, "user_s": u1.ru_utime - u0.ru_utime, "sys_s": u1.ru_stime - u0.ru_stime,
                 "records": n, "body_sha256": h.hexdigest()})
t1 = time.time()
for c in cs:
    shutil.copyfile(f"{W}/out/{c}.vcf.gz", f"/store/parts/{c}.vcf.gz")
u = resource.getrusage(resource.RUSAGE_SELF)
json.dump({"index": i, "contigs": rows, "annotate_s": t1 - t00, "fetch_s": time.time() - t1, "max_rss_gib": u.ru_maxrss / 1024**2,
           "start_unix": t00, "end_unix": time.time(), "vepyr": vepyr.__version__}, open(f"/store/pods/{i}.json", "w"), indent=1)
''' % FASTA


def vm(cmd: str) -> None:
    da.sh(["colima", "ssh", "--", "sudo", "sh", "-c", cmd])


def job(name: str, n: int, store: Path, workers: int) -> dict:
    vol = lambda nm, path, ro: ({"name": nm, "hostPath": {"path": str(path), "type": "Directory"}}, {"name": nm, "mountPath": f"/{nm}", "readOnly": ro})
    vols = [vol("store", store, False), vol("cache-src", da.CACHE, True), vol("ref-src", da.FASTA.parent, True)]
    vols.append(vol("work", WORK_HOST, False))
    mounts = [m for _, m in vols]
    sec = {"runAsUser": 501, "runAsGroup": 20}
    return {"apiVersion": "batch/v1", "kind": "Job", "metadata": {"name": name},
            "spec": {"completionMode": "Indexed", "completions": n, "parallelism": n, "backoffLimit": 0,
                     "template": {"spec": {"restartPolicy": "Never", "securityContext": sec,
                                           "volumes": [v for v, _ in vols],
                                           "initContainers": [{"name": "stage", "image": IMAGE, "imagePullPolicy": "Never",
                                                               "command": ["python", "-c", INIT], "volumeMounts": mounts}],
                                           "containers": [{"name": "annotate", "image": IMAGE, "imagePullPolicy": "Never",
                                                           "command": ["python", "-c", MAIN], "env": [{"name": "WORKERS", "value": str(workers)}],
                                                           "volumeMounts": mounts}]}}}}


def run(inp: Path, n: int, out: Path, workers: int) -> None:
    out.mkdir(parents=True, exist_ok=True)
    nodes = [{"name": f"pod{i}", "host": None, "workdir": "", "python": [], "workers": workers, "speed": 1.0} for i in range(n)]
    (out / "nodes.json").write_text(json.dumps(nodes))
    da.plan(inp, out / "nodes.json", out)
    p = json.loads((out / "plan.json").read_text())
    for d in ("shards", "parts", "pods"):
        (out / d).mkdir(exist_ok=True)
    t0 = time.perf_counter()
    for c in p["contig_order"]:
        s = out / "shards" / f"{c}.vcf.gz"
        with open(s, "wb") as f:
            t = subprocess.Popen(["tabix", "-h", str(inp), c], stdout=subprocess.PIPE)
            subprocess.run(["bgzip", "-c"], stdin=t.stdout, stdout=f, check=True)
            t.wait()
        da.sh(["tabix", "-f", "-p", "vcf", str(s)])
    split_s = time.perf_counter() - t0
    name = f"vepyr-annotate-n{n}"
    subprocess.run(["kubectl", "delete", "job", name, "--ignore-not-found"], capture_output=True)
    (out / "job.json").write_text(json.dumps(job(name, n, out, workers), indent=1))
    vm(f"rm -rf {WORK_HOST} && mkdir -p " + " ".join(f"{WORK_HOST}/{i}" for i in range(n)) + f" && chown -R 501:20 {WORK_HOST}")
    t1 = time.perf_counter()
    da.sh(["kubectl", "apply", "-f", str(out / "job.json")])
    w = subprocess.run(["kubectl", "wait", f"job/{name}", "--for=condition=complete", "--timeout=3h"], capture_output=True, text=True)
    job_s = time.perf_counter() - t1
    pods = subprocess.run(["kubectl", "get", "pods", "-l", f"job-name={name}", "-o", "json"], capture_output=True, text=True).stdout
    (out / "pods_k8s.json").write_text(pods)
    if w.returncode:
        raise SystemExit(f"job did not complete: {w.stderr}")
    per = [json.loads((out / "pods" / f"{i}.json").read_text()) for i in range(n)]
    stage = [json.loads((out / "pods" / f"{i}.stage.json").read_text()) for i in range(n)]
    done = [c["contig"] for q in per for c in q["contigs"]]
    if sorted(done) != sorted(p["contig_order"]) or len(done) != len(set(done)):
        raise SystemExit("contigs not covered exactly once")
    t2 = time.perf_counter()
    parts = [str(out / "parts" / f"{c}.vcf.gz") for c in p["contig_order"]]
    heads = [da.sh(["bcftools", "view", "-h", "--no-version", x]).stdout.splitlines() for x in parts]
    differing = sorted({l for h in heads[1:] for l in set(h) ^ set(heads[0])})
    merged = out / "merged.vcf.gz"
    da.sh(["bcftools", "concat", "--naive-force" if differing else "--naive", "--no-version", "-o", str(merged), *parts])
    merge_s = time.perf_counter() - t2
    body, recs = da.body_sha256(merged)
    res = {"n": n, "workers": workers, "split_s": split_s, "job_s": job_s, "merge_s": merge_s, "makespan_s": time.perf_counter() - t0,
           "merged_records": recs, "records_match_parts": recs == sum(c["records"] for q in per for c in q["contigs"]),
           "merged_body_sha256": body, "pods": per, "stage": stage, "finished_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    (out / "run.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps({k: v for k, v in res.items() if k not in ("pods", "stage")}))
    for q, s in zip(per, stage):
        print(q["index"], [c["contig"] for c in q["contigs"]], "stage", round(s["stage_s"], 1),
              {k: (v["bytes"], round(v["seconds"], 1)) for k, v in s["stage"].items()}, "annotate", round(q["annotate_s"], 1),
              "fetch", round(q["fetch_s"], 1), "rss", round(q["max_rss_gib"], 2))
    subprocess.run(["kubectl", "delete", "job", name, "--wait"], capture_output=True)
    vm(f"rm -rf {WORK_HOST}")


if __name__ == "__main__":
    if len(sys.argv) < 5 or sys.argv[1] != "run":
        raise SystemExit(__doc__)
    run(Path(sys.argv[2]).resolve(), int(sys.argv[3]), Path(sys.argv[4]).resolve(), int(sys.argv[5]) if len(sys.argv) > 5 else 1)
