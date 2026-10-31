#!/usr/bin/env python3
"""Summary of the VEPYR-in-the-VM series (thesis item D7, vepyr_in_vm.py run).

Per chromosome and workers: mean, SD, min-max of the `docker run` wall time (measured as VEP's)
and of the annotation time inside the container, CPU (user + sys), busy cores, peak RSS and the
distinct output bodies. Then the ratio VEP / VEPYR with both tools in the same VM, both timed as
`docker run`: VEP from the interleaved series of the same day (results/interleaved-2026-10-31,
one process, no --fork) against VEPYR workers=1, with the 95% CI of ratio_ci.py (log scale,
delta method, t(0.975, 4) = 2.776), next to the native ratio of that series.

Usage: vm_analysis.py RUNS_JSONL INTERLEAVED_SUMMARY_JSON OUT_JSON
"""
import json, math, statistics as st, sys
from pathlib import Path

T = 2.776


def stats(xs):
    m = st.mean(xs); sd = st.stdev(xs) if len(xs) > 1 else 0.0
    return {"n": len(xs), "mean": m, "median": st.median(xs), "sd": sd, "cv": sd / m, "min": min(xs), "max": max(xs)}


def ratio(a, b):
    r = a["mean"] / b["mean"]
    se = math.sqrt(a["cv"] ** 2 / a["n"] + b["cv"] ** 2 / b["n"])
    return {"ratio": r, "ci95": [r * math.exp(-T * se), r * math.exp(T * se)]}


def main():
    rows = [json.loads(l) for l in Path(sys.argv[1]).read_text().splitlines()]
    inter = json.loads(Path(sys.argv[2]).read_text())
    if any(r["returncode"] for r in rows):
        raise SystemExit("failed run")
    res = {}
    for chrom in ("chr22", "chr1"):
        for w in sorted({r["workers"] for r in rows if r["chrom"] == chrom}):
            rs = [r for r in rows if r["chrom"] == chrom and r["workers"] == w]
            d = {"docker_wall": stats([r["docker_run_wall_s"] for r in rs]),
                 "annotation": stats([r["annotation_seconds"] for r in rs]),
                 "cpu": stats([r["user_s"] + r["sys_s"] for r in rs]),
                 "cores": stats([(r["user_s"] + r["sys_s"]) / r["annotation_seconds"] for r in rs]),
                 "max_rss_gib": stats([r["max_rss_gib"] for r in rs]),
                 "bodies": sorted({r["body_sha256"] for r in rs}), "records": sorted({r["records"] for r in rs})}
            if chrom in inter:
                d["equal_native_body"] = d["bodies"] == inter[chrom]["vepyr"]["bodies"]
            if w == 1 and chrom in inter:
                vep = inter[chrom]["vep"]["wall"]
                d["vep_over_vepyr_vm_docker_wall"] = ratio(vep, d["docker_wall"])
                d["vep_over_vepyr_vm_annotation"] = ratio(vep, d["annotation"])
                d["native_interleaved_ratio"] = inter[chrom]["ratio_wall"]["of_means"]
                d["vm_over_native_annotation"] = d["annotation"]["mean"] / inter[chrom]["vepyr"]["wall"]["mean"]
            res[f"{chrom}_w{w}"] = d
            s = d.get("vep_over_vepyr_vm_docker_wall")
            print(chrom, f"w{w}", "docker", round(d["docker_wall"]["mean"], 2), "annot", round(d["annotation"]["mean"], 2),
                  f"({d['annotation']['min']:.2f}-{d['annotation']['max']:.2f})", "cpu", round(d["cpu"]["mean"], 2),
                  "cores", round(d["cores"]["mean"], 2), "rss", round(d["max_rss_gib"]["mean"], 2),
                  "bodies", len(d["bodies"]), "eq-native", d.get("equal_native_body"),
                  "| S_vm", None if not s else (round(s["ratio"], 1), [round(x, 1) for x in s["ci95"]]),
                  "| vm/native", None if "vm_over_native_annotation" not in d else round(d["vm_over_native_annotation"], 2))
    Path(sys.argv[3]).write_text(json.dumps(res, indent=1) + "\n")


if __name__ == "__main__":
    main()
