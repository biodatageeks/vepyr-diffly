#!/usr/bin/env python3
"""Summary of the interleaved VEP/VEPYR series (thesis items D5 and D6, interleaved_vep_vepyr.py).

Per chromosome and tool: mean, SD, min-max of the wall time and of the CPU time (VEP: cgroup
cpu.stat usage of the container; VEPYR: user + sys of the annotating process), cores = CPU/wall,
distinct output bodies. Ratios VEP/VEPYR: of the mean wall times with the 95% CI of ratio_ci.py
(log scale, delta method, t(0.975, 4) = 2.776, tools independent), the per-repeat ratios (both
tools ran back to back in every repeat), and the ratio of the mean CPU times.

Usage: interleaved_analysis.py RUN_DIR OUT_JSON   (RUN_DIR holds <chrom>/runs.jsonl)
"""
import json, math, statistics as st, sys
from pathlib import Path

T = 2.776


def stats(xs):
    m = st.mean(xs); sd = st.stdev(xs) if len(xs) > 1 else 0.0
    return {"n": len(xs), "mean": m, "sd": sd, "cv": sd / m, "min": min(xs), "max": max(xs)}


def main():
    run, out = Path(sys.argv[1]), Path(sys.argv[2])
    res = {}
    for f in sorted(run.glob("*/runs.jsonl")):
        rows = [json.loads(l) for l in f.read_text().splitlines()]
        if any(r["returncode"] for r in rows):
            raise SystemExit(f"failed run in {f}")
        chrom = f.parent.name
        cpu = {"vep": lambda r: r["cgroup_cpu_s"]["usage"], "vepyr": lambda r: r["user_s"] + r["sys_s"]}
        d = {}
        for tool in ("vep", "vepyr"):
            rs = sorted((r for r in rows if r["tool"] == tool), key=lambda r: r["repeat"])
            d[tool] = {"wall": stats([r["seconds"] for r in rs]), "cpu": stats([cpu[tool](r) for r in rs]),
                       "cores": stats([cpu[tool](r) / r["seconds"] for r in rs]),
                       "bodies": sorted({r["body_sha256"] for r in rs}), "repeats": [r["repeat"] for r in rs]}
        a, b = d["vep"]["wall"], d["vepyr"]["wall"]
        r = a["mean"] / b["mean"]
        se = math.sqrt(a["cv"] ** 2 / a["n"] + b["cv"] ** 2 / b["n"])
        by_rep = {}
        for x in rows:
            by_rep.setdefault(x["repeat"], {})[x["tool"]] = x
        paired = [v["vep"]["seconds"] / v["vepyr"]["seconds"] for _, v in sorted(by_rep.items()) if len(v) == 2]
        d["ratio_wall"] = {"of_means": r, "ci95": [r * math.exp(-T * se), r * math.exp(T * se)],
                           "per_repeat": paired, "per_repeat_min": min(paired), "per_repeat_max": max(paired)}
        d["ratio_cpu"] = d["vep"]["cpu"]["mean"] / d["vepyr"]["cpu"]["mean"]
        d["vep_body_equals_vepyr_body"] = d["vep"]["bodies"] == d["vepyr"]["bodies"]
        d["first_utc"], d["last_utc"] = rows[0]["utc"], rows[-1]["utc"]
        res[chrom] = d
        print(chrom, "VEP", round(a["mean"], 2), "VEPYR", round(b["mean"], 3), "S", round(r, 1),
              [round(x, 1) for x in d["ratio_wall"]["ci95"]], "paired", round(min(paired), 1), "-", round(max(paired), 1),
              "| CPU VEP", round(d["vep"]["cpu"]["mean"], 1), "VEPYR", round(d["vepyr"]["cpu"]["mean"], 2),
              "S_cpu", round(d["ratio_cpu"], 1), "| cores", round(d["vep"]["cores"]["mean"], 3), round(d["vepyr"]["cores"]["mean"], 2),
              "| bodies", len(d["vep"]["bodies"]), len(d["vepyr"]["bodies"]), "equal" if d["vep_body_equals_vepyr_body"] else "differ")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(res, indent=1) + "\n")


if __name__ == "__main__":
    main()
