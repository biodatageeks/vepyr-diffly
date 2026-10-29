#!/usr/bin/env python3
"""Scheduling-policy ablation of the cache-build queue by simulation (thesis item E3, 2026-10-30).

Replays the measured task times of each queue run under different dispatch orders with
greedy list scheduling on p identical processes (a task starts on the process that
becomes free first, in the order of the list). Task times are taken from the run as
measured, i.e. under that run's own contention; a different order could change them,
so the results show the effect of the order alone, not a measurement.

Policies:
  submitted     task-id order used by parallel_entity_queue.py (groups by creation, five
                entities per group); its simulated makespan is checked against the measured
                build time
  lpt_bytes     longest first by the queue's cost estimate (all_vars.gz bytes of the group);
                ties keep task-id order
  lpt_oracle    longest first by the measured task time (needs the times in advance)
  spt           shortest first by measured time (adversarial order)
  random        uniform random order, 1000 permutations: mean, 5th and 95th percentile
  static_rr     static round-robin split in task-id order (task i to process i mod p),
                no dynamic balancing; with p=10 and five entities per group this puts
                every variation task on two processes
  static_block  static split into p contiguous blocks of task ids (with p=10: two groups,
                i.e. two variation tasks, per process)
  static_lpt    static split by the queue's cost estimate: each task, in lpt_bytes order,
                to the process with the smallest estimated load so far
  static_oracle the same static LPT split, but by the measured task times

Lower bound LB = max(longest task, summed task time / p).

Usage: simulate_queue_policies.py OUT_JSON SUMMARY_JSON [SUMMARY_JSON ...]
  SUMMARY_JSON: queue_20_10_repeats/summary.json (key "repeats") or a single run summary
  (p<PP>_r<R>.summary.json of run_queue_scaling.py).
"""
from __future__ import annotations

import heapq
import json
import random
import statistics as st
import sys
from pathlib import Path

SEED = 20261030


def greedy(durations: list[float], order: list[int], p: int) -> float:
    free = [0.0] * p
    heapq.heapify(free)
    end = 0.0
    for i in order:
        t = heapq.heappop(free) + durations[i]
        end = max(end, t)
        heapq.heappush(free, t)
    return end


def static(durations: list[float], assign: list[int], p: int) -> float:
    load = [0.0] * p
    for i, proc in enumerate(assign):
        load[proc] += durations[i]
    return max(load)


def simulate(run: dict) -> dict:
    res = sorted(run["results"], key=lambda t: t["task"])
    p = run["processes"]
    d = [t["seconds"] for t in res]
    est = [t["bytes"] for t in res]
    n = len(d)
    ids = list(range(n))
    lb = max(max(d), sum(d) / p)
    out = {"processes": p, "tasks": n, "measured_s": run["seconds_wall"], "lower_bound_s": lb}
    pol = {
        "submitted": greedy(d, ids, p),
        "lpt_bytes": greedy(d, sorted(ids, key=lambda i: -est[i]), p),
        "lpt_oracle": greedy(d, sorted(ids, key=lambda i: -d[i]), p),
        "spt": greedy(d, sorted(ids, key=lambda i: d[i]), p),
        "static_rr": static(d, [i % p for i in ids], p),
        "static_block": static(d, [min(i * p // n, p - 1) for i in ids], p),
    }
    for name, cost in (("static_lpt", est), ("static_oracle", d)):
        load, assign = [0.0] * p, [0] * n
        for i in sorted(ids, key=lambda i: -cost[i]):
            k = min(range(p), key=lambda j: load[j])
            load[k] += cost[i]
            assign[i] = k
        pol[name] = static(d, assign, p)
    rng = random.Random(SEED)
    rnd = []
    for _ in range(1000):
        order = ids[:]
        rng.shuffle(order)
        rnd.append(greedy(d, order, p))
    rnd.sort()
    out["policies"] = {k: {"makespan_s": v, "over_lb": v / lb} for k, v in pol.items()}
    out["policies"]["random"] = {"makespan_s": st.mean(rnd), "over_lb": st.mean(rnd) / lb,
                                 "p5_s": rnd[49], "p95_s": rnd[949]}
    out["submitted_vs_measured"] = pol["submitted"] / run["seconds_wall"]
    return out


def main() -> int:
    out_path = Path(sys.argv[1])
    runs = []
    for src in sys.argv[2:]:
        d = json.loads(Path(src).read_text())
        for i, r in enumerate(d["repeats"] if "repeats" in d else [d], 1):
            item = simulate(r)
            item["source"] = f"{src}#{i}"
            runs.append(item)
    names = list(runs[0]["policies"])
    agg = {k: {"over_lb_min": min(r["policies"][k]["over_lb"] for r in runs),
               "over_lb_max": max(r["policies"][k]["over_lb"] for r in runs),
               "makespan_mean_s": st.mean(r["policies"][k]["makespan_s"] for r in runs)} for k in names}
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps({"runs": runs, "aggregate": agg}, indent=1) + "\n")
    print(f"{'policy':12} {'mean[s]':>9} {'T/LB min':>9} {'T/LB max':>9}")
    for k in names:
        a = agg[k]
        print(f"{k:12} {a['makespan_mean_s']:9.1f} {a['over_lb_min']:9.3f} {a['over_lb_max']:9.3f}")
    print("measured/LB:", [round(r["measured_s"] / r["lower_bound_s"], 3) for r in runs])
    print("simulated(submitted)/measured:", [round(r["submitted_vs_measured"], 4) for r in runs])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
