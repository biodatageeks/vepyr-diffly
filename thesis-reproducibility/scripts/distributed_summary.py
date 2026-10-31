#!/usr/bin/env python3
"""Summary of distributed_annotate.py runs (thesis items A1/A2).

Per run directory (N<n>/run.json): makespan, split, staging, annotation and merge times, slowest
node, merged body SHA-256 against the single-node reference (single/single.json) and, if given,
against a known body SHA-256 (e.g. the VEP chr1-22 output). From the per-contig annotation times
of the N=1 run (all contigs on one node, one after another) it replays the LPT assignment for
N = 1..16 equal nodes: makespan of the annotation alone, its lower bound max(longest contig,
total / N), and the speedup over the sum of contig times.

Usage: distributed_summary.py RUN_ROOT OUT_JSON [EXPECTED_BODY_SHA256]
"""
import json, sys
from pathlib import Path


def lpt(times: dict[str, float], n: int) -> float:
    load = [0.0] * n
    for t in sorted(times.values(), reverse=True):
        load[load.index(min(load))] += t
    return max(load)


def main():
    root, out = Path(sys.argv[1]), Path(sys.argv[2])
    expect = sys.argv[3] if len(sys.argv) > 3 else None
    single = json.loads((root / "single/single.json").read_text())
    res = {"single": single, "single_equals_expected": None if expect is None else single["body_sha256"] == expect, "runs": {}}
    base = None
    for d in sorted(root.glob("N*/run.json"), key=lambda p: int(p.parent.name[1:])):
        r = json.loads(d.read_text())
        n = len(r["nodes"])
        per = {c["contig"]: c["annotation_seconds"] for nd in r["nodes"] for c in nd["per_contig"]}
        if n == 1:
            base = per
        res["runs"][n] = {
            "makespan_s": r["makespan_s"], "split_s": r["split_s"], "merge_s": r["merge_s"],
            "stage_s_max": max(sum(v["seconds"] for v in nd["stage"].values()) for nd in r["nodes"]),
            "annotate_wall_s_max": max(nd["annotate_wall_s"] for nd in r["nodes"]),
            "annotate_wall_s": {nd["node"]: nd["annotate_wall_s"] for nd in r["nodes"]},
            "max_rss_gib": max(nd["max_rss_gib"] for nd in r["nodes"]),
            "cache_bytes_per_node": {nd["node"]: nd["stage"]["cache"]["bytes"] for nd in r["nodes"]},
            "sum_contig_annotation_s": sum(per.values()),
            "merged_body_equals_single": r["merged_body_sha256"] == single["body_sha256"],
            "records_match_parts": r["records_match_parts"], "merged_records": r["merged_records"],
            "speedup_vs_single": single["annotation_seconds"] / r["makespan_s"],
        }
    if base:
        total, longest = sum(base.values()), max(base.values())
        res["lpt_from_n1"] = {n: {"makespan_s": lpt(base, n), "lower_bound_s": max(longest, total / n),
                                  "speedup": total / lpt(base, n)} for n in range(1, 17)}
        res["n1_contig_seconds"] = dict(sorted(base.items(), key=lambda x: -x[1]))
    out.write_text(json.dumps(res, indent=1) + "\n")
    print("single", round(single["annotation_seconds"], 1), "s", single["body_sha256"][:12], "expected:", res["single_equals_expected"])
    for n, r in res["runs"].items():
        print(f"N={n}", {k: (round(v, 2) if isinstance(v, float) else v) for k, v in r.items()
                         if k not in ("annotate_wall_s", "cache_bytes_per_node")})
    for n in (2, 4, 8, 16):
        if base:
            x = res["lpt_from_n1"][n]
            print(f"LPT N={n}: {x['makespan_s']:.1f} s (bound {x['lower_bound_s']:.1f}), speedup {x['speedup']:.2f}")


if __name__ == "__main__":
    main()
