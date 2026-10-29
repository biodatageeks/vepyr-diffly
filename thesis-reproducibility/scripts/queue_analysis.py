#!/usr/bin/env python3
"""Analysis of the dynamic cache-build queue runs (thesis, 2026-10-29).

From results/cache-build-2026-10/parallel_full_cache/queue_20_10_repeats/summary.json:
  * makespan vs the lower bound max(longest task, sum of task times / processes),
  * Pearson and Spearman correlation of the cost estimate (all_vars.gz bytes of the
    task's contig group) with the measured task time, per entity,
  * a Gantt chart of one repeat. Worker ids were not recorded, so tasks are placed
    on lanes greedily from their start/finish times (a lane is free once its last
    task finished); this reproduces the occupancy, not the actual process ids.
Usage (from thesis-reproducibility/): python3 scripts/queue_analysis.py [OUT_PDF]
Writes results/queue-analysis-2026-10-29/queue_analysis.json.
"""
import json, math, os, statistics as st, sys
from pathlib import Path
os.chdir(Path(__file__).resolve().parent.parent)
SRC = Path('results/cache-build-2026-10/parallel_full_cache/queue_20_10_repeats/summary.json')
OUT = Path('results/queue-analysis-2026-10-29'); OUT.mkdir(parents=True, exist_ok=True)
ENT = ['variation', 'transcript', 'exon', 'regulatory', 'motif']

def pearson(x, y):
    mx, my = st.mean(x), st.mean(y)
    return sum((a-mx)*(b-my) for a, b in zip(x, y)) / math.sqrt(sum((a-mx)**2 for a in x) * sum((b-my)**2 for b in y))

def ranks(v):
    order = sorted(range(len(v)), key=lambda i: v[i]); r = [0.0]*len(v); i = 0
    while i < len(v):
        j = i
        while j+1 < len(v) and v[order[j+1]] == v[order[i]]: j += 1
        for k in range(i, j+1): r[order[k]] = (i+j)/2
        i = j+1
    return r

d = json.loads(SRC.read_text())
rows = []
for rep in d['repeats']:
    res, p = rep['results'], rep['processes']
    tsum = sum(t['seconds'] for t in res); tmax = max(t['seconds'] for t in res)
    longest = max(res, key=lambda t: t['seconds'])
    corr = {}
    for e in ENT:
        ts = [t for t in res if t['entity'] == e]
        x = [t['bytes'] for t in ts]; y = [t['seconds'] for t in ts]
        corr[e] = {'pearson': pearson(x, y), 'spearman': pearson(ranks(x), ranks(y)), 'sum_s': sum(y)}
    rows.append({'repeat': rep['repeat'], 'makespan_s': rep['seconds_wall'], 'sum_task_s': tsum,
                 'sum_over_p_s': tsum/p, 'max_task_s': tmax, 'max_task': [longest['entity'], longest['group']],
                 'lower_bound_s': max(tmax, tsum/p), 'makespan_over_lb': rep['seconds_wall']/max(tmax, tsum/p),
                 'variation_share': corr['variation']['sum_s']/tsum, 'corr': corr})
json.dump(rows, open(OUT/'queue_analysis.json', 'w'), indent=1)
for r in rows:
    print(r['repeat'], round(r['makespan_s'], 1), round(r['lower_bound_s'], 1), round(r['makespan_over_lb'], 3),
          round(r['variation_share'], 3), {e: round(c['pearson'], 2) for e, c in r['corr'].items()})

if len(sys.argv) > 1:
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    rep = sorted(d['repeats'], key=lambda r: r['seconds_wall'])[len(d['repeats'])//2]  # median repeat
    res = sorted(rep['results'], key=lambda t: t['started'])
    t0 = min(t['started'] for t in res)
    lanes = []  # last finish per lane
    place = []
    for t in res:
        free = [i for i, f in enumerate(lanes) if f <= t['started'] + 1e-6]
        i = free[0] if free else len(lanes)
        if i == len(lanes): lanes.append(0)
        lanes[i] = t['finished']; place.append((i, t))
    colors = {'variation': '#3b6fb6', 'transcript': '#d08a2e', 'exon': '#5a9e5a', 'regulatory': '#a05aa0', 'motif': '#888888'}
    fig, ax = plt.subplots(figsize=(7.2, 3.4))
    for i, t in place:
        ax.barh(i, t['seconds'], left=t['started']-t0, color=colors[t['entity']], edgecolor='white', linewidth=0.4, height=0.8)
    lb = max(max(t['seconds'] for t in res), sum(t['seconds'] for t in res)/rep['processes'])
    ax.axvline(lb, color='black', ls='--', lw=0.8)
    ax.text(lb, len(lanes)-0.3, ' dolna granica', fontsize=7, va='bottom')
    ax.set_xlabel('czas od uruchomienia kolejki [s]'); ax.set_ylabel('proces (odtworzony tor)')
    ax.set_yticks(range(len(lanes))); ax.set_yticklabels([str(i+1) for i in range(len(lanes))], fontsize=7)
    ax.invert_yaxis(); ax.set_xlim(0, rep['seconds_wall']*1.02)
    ax.legend(handles=[Patch(color=c, label=e) for e, c in colors.items()], fontsize=7, ncol=5, loc='upper center', bbox_to_anchor=(0.5, -0.22), frameon=False)
    fig.tight_layout(); fig.savefig(sys.argv[1]); print('figure', sys.argv[1], 'repeat', rep['repeat'], 'lanes', len(lanes))
