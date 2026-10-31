#!/usr/bin/env python3
"""Distributed annotation prototype: one input VCF split by contig across nodes (thesis items A1/A2, 2026-10-31).

The unit of work is one contig of the input. A node needs only the cache fragments of its own
contigs (<entity>/<contig>.parquet of every core entity), so the per-contig layout of the cache
decides how much each node reads and how much has to be copied to it.

  plan    per-contig record counts of the input (bcftools index -s), LPT assignment of contigs to
          nodes by record count divided by node speed, written to PLAN.json together with the
          input SHA-256.
  run     for every node in parallel: stage (per-contig input shards cut with tabix -h, i.e. the
          input lines unchanged; the node's cache fragments with every entity's chrom_manifest.json;
          the FASTA and its index), annotate (one process per node, one vepyr.annotate call per
          contig, bgzf output) and fetch the outputs. Per node: bytes and time of staging, per-contig
          annotation time, CPU, peak RSS, records, body SHA-256, fetch time. Then merge: the
          per-contig outputs concatenated in input order with bcftools concat --naive (bgzf blocks
          copied, header of the first output), every contig exactly once, record count equal to
          the sum over contigs, and the body SHA-256 of the result.
  single  the reference: one vepyr.annotate call over the whole input on one node with the full cache.

Nodes (JSON list): {"name": "m5", "host": null, "workdir": "/abs/dir", "python": ["/abs/python3"],
"workers": 4, "speed": 1.0}. host null = a local process (correctness and overhead only; several
local nodes share one machine and are not distribution in the thesis's sense); otherwise an ssh
destination, files copied with rsync. "python" is the command prefix that runs a script on the node
(it may be a docker run ... python3 with workdir mounted at the same path).

Usage: distributed_annotate.py plan INPUT_VCF_GZ NODES.json OUT_DIR
       distributed_annotate.py run OUT_DIR [CACHE_DIR FASTA]
       distributed_annotate.py single INPUT_VCF_GZ NODE_NAME NODES.json OUT_DIR [CACHE_DIR FASTA]
"""
from __future__ import annotations

import concurrent.futures as cf
import hashlib
import json
import os
import shlex
import subprocess
import sys
import time
from pathlib import Path

D = Path("/Users/lukaszjezapkowicz/workspace/data_vepyr")
CACHE = D / "cache/116_GRCh38_merged"
FASTA = D / "input/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
NODE = r'''
import gzip, hashlib, json, resource, sys, time, vepyr
a = json.loads(sys.argv[1])
rows = []
for c in a["contigs"]:
    inp, out = f"{a['workdir']}/in/{c}.vcf.gz", f"{a['workdir']}/out/{c}.vcf.gz"
    u0 = resource.getrusage(resource.RUSAGE_SELF); t0 = time.perf_counter()
    vepyr.annotate(vcf=inp, cache_dir=f"{a['workdir']}/cache", everything=True, hgvs=True,
                   reference_fasta=f"{a['workdir']}/ref/{a['fasta']}", workers=a["workers"],
                   output_vcf=out, show_progress=False)
    dt = time.perf_counter() - t0; u1 = resource.getrusage(resource.RUSAGE_SELF)
    h = hashlib.sha256(); n = 0
    with gzip.open(out, "rb") as f:
        for line in f:
            if not line.startswith(b"#"):
                h.update(line); n += 1
    rows.append({"contig": c, "annotation_seconds": dt, "user_s": u1.ru_utime - u0.ru_utime,
                 "sys_s": u1.ru_stime - u0.ru_stime, "records": n, "body_sha256": h.hexdigest()})
u = resource.getrusage(resource.RUSAGE_SELF)
rss = u.ru_maxrss / 1024**3 if sys.platform == "darwin" else u.ru_maxrss / 1024**2
print(json.dumps({"contigs": rows, "max_rss_gib": rss, "vepyr": vepyr.__version__}))
'''
SINGLE = r'''
import gzip, hashlib, json, resource, sys, time, vepyr
a = json.loads(sys.argv[1])
t0 = time.perf_counter()
vepyr.annotate(vcf=a["vcf"], cache_dir=a["cache"], everything=True, hgvs=True, reference_fasta=a["fasta"],
               workers=a["workers"], output_vcf=a["out"], show_progress=False)
dt = time.perf_counter() - t0
u = resource.getrusage(resource.RUSAGE_SELF)
h = hashlib.sha256(); n = 0
with gzip.open(a["out"], "rb") as f:
    for line in f:
        if not line.startswith(b"#"):
            h.update(line); n += 1
rss = u.ru_maxrss / 1024**3 if sys.platform == "darwin" else u.ru_maxrss / 1024**2
print(json.dumps({"annotation_seconds": dt, "user_s": u.ru_utime, "sys_s": u.ru_stime, "max_rss_gib": rss,
                  "records": n, "body_sha256": h.hexdigest(), "vepyr": vepyr.__version__}))
'''


def sh(cmd, **kw) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, check=True, capture_output=True, text=True, **kw)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def body_sha256(path: Path) -> tuple[str, int]:
    p = subprocess.Popen(["bgzip", "-dc", str(path)], stdout=subprocess.PIPE)
    h, n = hashlib.sha256(), 0
    for line in p.stdout:
        if not line.startswith(b"#"):
            h.update(line); n += 1
    if p.wait():
        raise SystemExit(f"bgzip -dc {path} failed")
    return h.hexdigest(), n


def plan(inp: Path, nodes_file: Path, out: Path) -> None:
    nodes = json.loads(nodes_file.read_text())
    counts = [(l.split("\t")[0], int(l.split("\t")[2])) for l in sh(["bcftools", "index", "-s", str(inp)]).stdout.splitlines()]
    load = {n["name"]: 0.0 for n in nodes}
    assign = {n["name"]: [] for n in nodes}
    speed = {n["name"]: n.get("speed", 1.0) for n in nodes}
    for contig, k in sorted(counts, key=lambda x: -x[1]):  # LPT: largest first, to the node finishing it earliest
        best = min(nodes, key=lambda n: (load[n["name"]] + k / speed[n["name"]], n["name"]))["name"]
        assign[best].append(contig); load[best] += k / speed[best]
    order = [c for c, _ in counts]  # bcftools index -s lists contigs in file order
    for n in nodes:
        n["contigs"] = sorted(assign[n["name"]], key=order.index)
    lb = max(max(k for _, k in counts) / max(speed.values()), sum(k for _, k in counts) / sum(speed.values()))
    out.mkdir(parents=True, exist_ok=True)
    p = {"input": str(inp), "input_sha256": sha256(inp), "contig_order": order, "records": dict(counts),
         "nodes": nodes, "predicted_load": load, "lower_bound": lb, "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    (out / "plan.json").write_text(json.dumps(p, indent=1) + "\n")
    print(json.dumps({n["name"]: (n["contigs"], round(load[n["name"]])) for n in nodes}), "LB", round(lb))


def on_node(node: dict, cmd: list[str], **kw) -> subprocess.CompletedProcess:
    if node.get("host"):
        cmd = ["ssh", node["host"], shlex.join(cmd)]
    return sh(cmd, **kw)


def stage(node: dict, shards: Path, cache: Path, fasta: Path) -> dict:
    """Input shards, cache fragments and FASTA of one node; returns bytes and seconds per part."""
    w, cs = Path(node["workdir"]), node["contigs"]
    files = {"input": [shards / f"{c}.vcf.gz{s}" for c in cs for s in ("", ".tbi")],
             "cache": [f for e in sorted(p for p in cache.iterdir() if p.is_dir())
                       for f in [e / "chrom_manifest.json", *(e / f"{c}.parquet" for c in cs)] if f.exists()],
             "fasta": [fasta, Path(str(fasta) + ".fai")]}
    dest = {"input": w / "in", "cache": w / "cache", "fasta": w / "ref"}
    res = {}
    on_node(node, ["mkdir", "-p", str(w / "in"), str(w / "cache"), str(w / "ref"), str(w / "out")])
    for part, fs in files.items():
        t0 = time.perf_counter()
        if node.get("host"):
            base = cache if part == "cache" else fs[0].parent
            lst = "\n".join(str(f.relative_to(base)) for f in fs) + "\n"
            sh(["rsync", "-a", "--files-from=-", str(base) + "/", f"{node['host']}:{dest[part]}/"], input=lst)
        else:
            for f in fs:
                t = dest[part] / (f.relative_to(cache) if part == "cache" else f.name)
                t.parent.mkdir(parents=True, exist_ok=True)
                if not t.exists():
                    os.link(f, t)  # same file system: hard link, nothing copied
        res[part] = {"bytes": sum(f.stat().st_size for f in fs), "seconds": time.perf_counter() - t0, "files": len(fs)}
    return res


def run_node(node: dict, shards: Path, cache: Path, fasta: Path, out: Path) -> dict:
    t0 = time.perf_counter()
    r = {"node": node["name"], "host": node.get("host"), "workers": node["workers"], "contigs": node["contigs"]}
    r["stage"] = stage(node, shards, cache, fasta)
    a = {"workdir": node["workdir"], "contigs": node["contigs"], "workers": node["workers"], "fasta": fasta.name}
    t1 = time.perf_counter()
    p = on_node(node, [*node["python"], "-c", NODE, json.dumps(a)])
    r["annotate_wall_s"] = time.perf_counter() - t1
    res = json.loads(p.stdout.strip().splitlines()[-1])
    r["per_contig"] = res.pop("contigs")
    r.update(res)
    t2 = time.perf_counter()
    got = out / "parts"
    got.mkdir(exist_ok=True)
    for c in node["contigs"]:
        src = Path(node["workdir"]) / "out" / f"{c}.vcf.gz"
        if node.get("host"):
            sh(["rsync", "-a", f"{node['host']}:{src}", str(got / src.name)])
        else:
            os.replace(src, got / src.name)
    r["fetch_s"] = time.perf_counter() - t2
    r["fetch_bytes"] = sum((got / f"{c}.vcf.gz").stat().st_size for c in node["contigs"])
    r["node_total_s"] = time.perf_counter() - t0
    return r


def run(out: Path, cache: Path, fasta: Path) -> None:
    p = json.loads((out / "plan.json").read_text())
    inp, shards = Path(p["input"]), out / "shards"
    if sha256(inp) != p["input_sha256"]:
        raise SystemExit("input changed since the plan")
    t0 = time.perf_counter()
    shards.mkdir(exist_ok=True)
    for c in p["contig_order"]:  # unchanged input lines: tabix -h, not bcftools view (which re-encodes records)
        s = shards / f"{c}.vcf.gz"
        with open(s, "wb") as f:
            t = subprocess.Popen(["tabix", "-h", str(inp), c], stdout=subprocess.PIPE)
            subprocess.run(["bgzip", "-c"], stdin=t.stdout, stdout=f, check=True)
            t.wait()
        sh(["tabix", "-f", "-p", "vcf", str(s)])
    split_s = time.perf_counter() - t0
    with cf.ThreadPoolExecutor(len(p["nodes"])) as ex:
        nodes = list(ex.map(lambda n: run_node(n, shards, cache, fasta, out), p["nodes"]))
    t1 = time.perf_counter()
    done = [c["contig"] for n in nodes for c in n["per_contig"]]
    if sorted(done) != sorted(p["contig_order"]) or len(done) != len(set(done)):
        raise SystemExit(f"contigs not covered exactly once: {sorted(done)}")
    parts = [str(out / "parts" / f"{c}.vcf.gz") for c in p["contig_order"]]
    merged = out / "merged.vcf.gz"
    heads = [sh(["bcftools", "view", "-h", "--no-version", p]).stdout.splitlines() for p in parts]
    differing = sorted({l for h in heads[1:] for l in set(h) ^ set(heads[0])})
    # --naive refuses parts whose headers differ; the body is what the check compares, so then force it and record why
    sh(["bcftools", "concat", "--naive-force" if differing else "--naive", "--no-version", "-o", str(merged), *parts])
    merge_s = time.perf_counter() - t1
    body, n = body_sha256(merged)
    per = {c["contig"]: c for nd in nodes for c in nd["per_contig"]}
    res = {"plan": p, "split_s": split_s, "merge_s": merge_s, "makespan_s": time.perf_counter() - t0,
           "nodes": nodes, "part_headers_identical": not differing, "part_header_differences": differing[:20], "merged_records": n, "records_match_parts": n == sum(c["records"] for c in per.values()),
           "merged_body_sha256": body, "finished_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    (out / "run.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps({k: v for k, v in res.items() if k not in ("plan", "nodes")}))
    for nd in nodes:
        print(nd["node"], {k: round(nd[k], 2) for k in ("annotate_wall_s", "fetch_s", "node_total_s", "max_rss_gib")},
              {k: (v["bytes"], round(v["seconds"], 2)) for k, v in nd["stage"].items()})


def single(inp: Path, name: str, nodes_file: Path, out: Path, cache: Path, fasta: Path) -> None:
    node = next(n for n in json.loads(nodes_file.read_text()) if n["name"] == name)
    out.mkdir(parents=True, exist_ok=True)
    a = {"vcf": str(inp), "cache": str(cache), "fasta": str(fasta), "workers": node["workers"],
         "out": str(out / "single.vcf.gz")}
    t0 = time.perf_counter()
    p = on_node(node, [*node["python"], "-c", SINGLE, json.dumps(a)])
    r = json.loads(p.stdout.strip().splitlines()[-1])
    r.update(node=name, workers=node["workers"], process_wall_s=time.perf_counter() - t0, input=str(inp),
             finished_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
    (out / "single.json").write_text(json.dumps(r, indent=1) + "\n")
    print(json.dumps(r))


def main() -> int:
    mode = sys.argv[1]
    if mode == "plan":
        plan(Path(sys.argv[2]).resolve(), Path(sys.argv[3]), Path(sys.argv[4]).resolve())
    elif mode == "run":
        a = sys.argv[2:]
        run(Path(a[0]).resolve(), Path(a[1]) if len(a) > 1 else CACHE, Path(a[2]) if len(a) > 2 else FASTA)
    elif mode == "single":
        a = sys.argv[2:]
        single(Path(a[0]).resolve(), a[1], Path(a[2]), Path(a[3]).resolve(),
               Path(a[4]) if len(a) > 4 else CACHE, Path(a[5]) if len(a) > 5 else FASTA)
    else:
        raise SystemExit(__doc__)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
