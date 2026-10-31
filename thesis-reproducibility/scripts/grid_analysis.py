#!/usr/bin/env python3
"""Summary of the workers x CPUs grid (thesis item D13, vepyr_in_vm.py grid / grid-local).

chr1 baseline, workers p in {1, 2, 4, 8, 10} on c in {1, 2, 4, 10} VM CPUs (docker --cpuset-cpus).
Per cell: mean, SD, min-max of the annotation time, CPU time (user + sys), busy cores, peak RSS
and the distinct output bodies. Derived: speedup S(p, c) = T(1, 1) / T(p, c); single-CPU cost of
the partitioned path T(p, 1) / T(1, 1); CPU-time inflation CPU(p, c) / CPU(1, 1).

Usage: grid_analysis.py GRID_RUNS_JSONL OUT_JSON
"""
import json, statistics as st, sys
from pathlib import Path


def stats(xs):
    m = st.mean(xs); sd = st.stdev(xs) if len(xs) > 1 else 0.0
    return {"n": len(xs), "mean": m, "median": st.median(xs), "sd": sd, "cv": sd / m, "min": min(xs), "max": max(xs)}


def main():
    rows = [json.loads(l) for l in Path(sys.argv[1]).read_text().splitlines()]
    if any(r["returncode"] for r in rows):
        raise SystemExit("failed run")
    cells = {}
    for r in rows:
        cells.setdefault((r["workers"], r["cpus"]), []).append(r)
    res = {"cells": {}, "bodies": sorted({r["body_sha256"] for r in rows})}
    for (p, c), rs in sorted(cells.items()):
        res["cells"][f"p{p}_c{c}"] = {
            "workers": p, "cpus": c, "annotation": stats([r["annotation_seconds"] for r in rs]),
            "cpu": stats([r["user_s"] + r["sys_s"] for r in rs]),
            "cores": stats([(r["user_s"] + r["sys_s"]) / r["annotation_seconds"] for r in rs]),
            "max_rss_gib": stats([r["max_rss_gib"] for r in rs]), "bodies": sorted({r["body_sha256"] for r in rs})}
    base = res["cells"]["p1_c1"]
    for k, v in res["cells"].items():
        v["speedup_vs_p1_c1"] = base["annotation"]["mean"] / v["annotation"]["mean"]
        v["cpu_vs_p1_c1"] = v["cpu"]["mean"] / base["cpu"]["mean"]
    ps = sorted({p for p, _ in cells}); cs = sorted({c for _, c in cells})
    print("bodies:", len(res["bodies"]), res["bodies"][0][:12])
    for title, key in (("annotation time [s]", lambda v: v["annotation"]["mean"]), ("speedup vs p1 c1", lambda v: v["speedup_vs_p1_c1"]),
                       ("CPU time [s]", lambda v: v["cpu"]["mean"]), ("busy cores", lambda v: v["cores"]["mean"]),
                       ("CV of time [%]", lambda v: 100 * v["annotation"]["cv"])):
        print(f"\n{title}: rows workers, columns CPUs {cs}")
        for p in ps:
            print(f"  p={p:<3}", "  ".join(f"{key(res['cells'][f'p{p}_c{c}']):8.2f}" for c in cs))
    Path(sys.argv[2]).write_text(json.dumps(res, indent=1) + "\n")


if __name__ == "__main__":
    main()
