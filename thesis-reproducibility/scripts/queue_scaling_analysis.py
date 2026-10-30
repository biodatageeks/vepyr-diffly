#!/usr/bin/env python3
"""Scaling curve of the dynamic cache-build queue (thesis items E1/G1, 2026-10-29).

Reads the output of run_queue_scaling.py (vepyr branch thesis/cache-build-scaling):
runs.jsonl and p<PP>_r<R>.summary.json with per-task start/finish times. For every
process count p it reports the build time T_p (mean, median, SD n-1, min-max), the
speedup S_p = mean T_1 / mean T_p with the range mean T_1 / max T_p .. mean T_1 / min T_p,
the efficiency E_p = S_p / p, the summed task time, the longest task and the lower
bound LB_p = max(longest task, summed task time / p) with T_p / LB_p per run.
When LB_p is the summed task time W_p / p, E_p splits into W_1 / W_p (tasks running
longer under concurrency) times LB_p / T_p (idle time left by the schedule); both factors
are reported as work_ratio and lb_over_T (means over runs).
Only runs with exit code 0 and a summary are used. The power source before and after each
run is taken from the pmset snapshot; a run is marked on_battery if either says battery.
With --reruns DIR (run_queue_reruns.py output), the cells re-run on AC replace the
original ones, and the replaced runs are written to battery_vs_ac.json for comparison.
With --interference FILE (JSON list of {label, start_utc, end_utc}), every run whose
before..after interval overlaps a window is tagged with its label (e.g. a download that
ran during the series); counted per p in n_interference.

Usage: queue_scaling_analysis.py RUN_DIR OUT_DIR [FIG_PDF] [--reruns RERUN_DIR] [--interference FILE]
"""
from __future__ import annotations

import json
import re
import statistics as st
import sys
from datetime import datetime
from pathlib import Path


def load(run_dir: Path, windows: list[dict] = ()) -> list[dict]:
    runs = []
    for line in (run_dir / "runs.jsonl").read_text().splitlines():
        r = json.loads(line)
        summ = run_dir / f"p{r['processes']:02d}_r{r['repeat']}.summary.json"
        if r.get("exit_code") != 0 or not summ.exists():
            continue
        s = json.loads(summ.read_text())
        secs = [t["seconds"] for t in s["results"]]
        assert len(secs) == s["tasks"] == 100 and all(t["status"] == "done" for t in s["results"])
        p = r["processes"]
        lb = max(max(secs), sum(secs) / p)
        power = [("AC" if "'AC Power'" in r[k]["batt"] else "battery",
                  int(m.group(1)) if (m := re.search(r"(\d+)%", r[k]["batt"])) else None) for k in ("before", "after")]
        runs.append({"p": p, "repeat": r["repeat"], "T": s["seconds_wall"], "sum_task": sum(secs),
                     "power": power, "on_battery": any(src == "battery" for src, _ in power),
                     "rerun_on_ac": bool(r.get("rerun_on_ac")),
                     "start_utc": r["before"]["utc"], "end_utc": r["after"]["utc"],
                     "interference": [w["label"] for w in windows
                                      if datetime.fromisoformat(w["start_utc"]) < datetime.fromisoformat(r["after"]["utc"])
                                      and datetime.fromisoformat(r["before"]["utc"]) < datetime.fromisoformat(w["end_utc"])],
                     "max_task": max(secs), "LB": lb, "T_over_LB": s["seconds_wall"] / lb,
                     "files": r.get("files"), "therm_before": r["before"]["therm"].splitlines()[-1:],
                     "load_before": r["before"]["load"]})
    return runs


def summarize(runs: list[dict]) -> list[dict]:
    ps = sorted({r["p"] for r in runs})
    by = {p: [r for r in runs if r["p"] == p] for p in ps}
    t1 = st.mean(r["T"] for r in by[1]) if 1 in by else None
    out = []
    for p in ps:
        T = [r["T"] for r in by[p]]
        row = {"p": p, "n": len(T), "mean": st.mean(T), "median": st.median(T),
               "sd": st.stdev(T) if len(T) > 1 else None, "min": min(T), "max": max(T),
               "sum_task_mean": st.mean(r["sum_task"] for r in by[p]),
               "max_task_mean": st.mean(r["max_task"] for r in by[p]),
               "T_over_LB_min": min(r["T_over_LB"] for r in by[p]),
               "T_over_LB_max": max(r["T_over_LB"] for r in by[p]),
               "files": sorted({r["files"] for r in by[p]}),
               "n_on_battery": sum(r["on_battery"] for r in by[p]),
               "n_interference": sum(bool(r["interference"]) for r in by[p]),
               "lb_over_T": st.mean(r["LB"] / r["T"] for r in by[p])}
        if t1:
            row.update({"S": t1 / row["mean"], "S_lo": t1 / row["max"], "S_hi": t1 / row["min"],
                        "E": t1 / row["mean"] / p,
                        "work_ratio": st.mean(r["sum_task"] for r in by[1]) / row["sum_task_mean"]})
        out.append(row)
    return out


def figure(rows: list[dict], path: str) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = [r for r in rows if "S" in r]
    ps = [r["p"] for r in rows]
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.4, 3.0))
    a.plot([1, 10], [1, 10], color="#9a9a9a", lw=1, ls="--", label="skalowanie idealne")
    a.axvspan(4.5, 10.5, color="#eeeeee", zorder=0)
    a.text(7.5, 0.6, "$p>4$: więcej procesów\nniż rdzeni wydajnościowych", fontsize=6.5, color="#555555",
           ha="center", va="bottom")
    a.errorbar(ps, [r["S"] for r in rows],
               yerr=[[r["S"] - r["S_lo"] for r in rows], [r["S_hi"] - r["S"] for r in rows]],
               color="#3b6fb6", marker="o", ms=4, lw=1.5, capsize=3, label="dynamiczna kolejka")
    for r in rows:
        a.annotate(f"{r['S']:.2f}".replace(".", ","), (r["p"], r["S"]), textcoords="offset points",
                   xytext=(-6, 7), fontsize=7)
    a.set_xticks(ps); a.set_xlim(0.5, 10.5); a.set_ylim(0, 10.5)
    a.set_xlabel("liczba procesów $p$"); a.set_ylabel("przyspieszenie $S_p$")
    a.legend(fontsize=7, frameon=False, loc="upper left"); a.set_title("(a)", fontsize=9, loc="left")
    b.plot(ps, [r["sum_task_mean"] / 3600 for r in rows], color="#3b6fb6", marker="o", ms=4, lw=1.5)
    for r in rows:
        b.annotate(f"{r['sum_task_mean'] / 3600:.2f}".replace(".", ","), (r["p"], r["sum_task_mean"] / 3600),
                   textcoords="offset points", xytext=(-6, 7), fontsize=7)
    b.axvspan(4.5, 10.5, color="#eeeeee", zorder=0)
    b.set_xticks(ps); b.set_xlim(0.5, 10.5); b.set_ylim(0, 1.15 * max(r["sum_task_mean"] for r in rows) / 3600)
    b.set_xlabel("liczba procesów $p$"); b.set_ylabel("suma czasów 100 zadań [h]")
    b.set_title("(b)", fontsize=9, loc="left")
    for ax in (a, b):
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", color="#e3e3e3", lw=0.6)
    fig.tight_layout(); fig.savefig(path)


def main() -> int:
    args = sys.argv[1:]
    reruns, windows = None, []
    if "--reruns" in args:
        i = args.index("--reruns")
        reruns = Path(args[i + 1])
        del args[i:i + 2]
    if "--interference" in args:
        i = args.index("--interference")
        windows = json.loads(Path(args[i + 1]).read_text())
        del args[i:i + 2]
    run_dir, out = Path(args[0]), Path(args[1])
    out.mkdir(parents=True, exist_ok=True)
    runs = load(run_dir, windows)
    if reruns:
        new = {(r["p"], r["repeat"]): r for r in load(reruns, windows)}
        replaced = [r for r in runs if (r["p"], r["repeat"]) in new]
        runs = [r for r in runs if (r["p"], r["repeat"]) not in new] + list(new.values())
        runs.sort(key=lambda r: (r["p"], r["repeat"]))
        (out / "battery_vs_ac.json").write_text(json.dumps(
            [{"p": r["p"], "repeat": r["repeat"], "T_original": r["T"], "power_original": r["power"],
              "interference_original": r["interference"], "T_rerun_ac": new[(r["p"], r["repeat"])]["T"],
              "ratio": r["T"] / new[(r["p"], r["repeat"])]["T"]} for r in replaced], indent=1) + "\n")
    rows = summarize(runs)
    (out / "queue_scaling_runs.json").write_text(json.dumps(runs, indent=1) + "\n")
    (out / "queue_scaling.json").write_text(json.dumps(rows, indent=1) + "\n")
    cols = ["p", "n", "mean", "median", "sd", "min", "max", "S", "S_lo", "S_hi", "E",
            "sum_task_mean", "max_task_mean", "T_over_LB_min", "T_over_LB_max", "work_ratio", "lb_over_T", "n_on_battery", "n_interference"]
    lines = ["\t".join(cols)] + ["\t".join("" if r.get(c) is None else
                                           (f"{r[c]:.4f}" if isinstance(r.get(c), float) else str(r[c]))
                                           for c in cols) for r in rows]
    (out / "queue_scaling.tsv").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    if len(args) > 2:
        figure(rows, args[2])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
