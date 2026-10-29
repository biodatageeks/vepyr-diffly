#!/usr/bin/env python3
"""Figure for chapter 5 (thesis, 2026-10-29): (a) annotation speedup S_p vs number
of processes/workers with the ideal line, whiskers = T_1 mean / T_p max .. T_p min;
(b) time vs peak RSS for chr1-22 (VEP single process, VEPYR workers 1-10).
Reads results/dispersion-2026-10-29/dispersion.json (scripts/dispersion_stats.py) and
results/completion-2026-10/m5_completion_20261024/S5/result.json.
Usage (from thesis-reproducibility/): python3 scripts/fig_annotation_scaling.py OUT_PDF"""
import json, os, sys
from pathlib import Path
os.chdir(Path(__file__).resolve().parent.parent)
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
D = {(o['group'], o['label']): o for o in json.load(open('results/dispersion-2026-10-29/dispersion.json'))}
S5 = json.load(open('results/completion-2026-10/m5_completion_20261024/S5/result.json'))
fig, (a, b) = plt.subplots(1, 2, figsize=(7.4, 3.2), gridspec_kw={'width_ratios': [1.25, 1]})
series = [
    ('VEPYR chr1', 'vepyr_scaling', 'chr1 w{}', [1, 2, 4, 8, 10], lambda p: p, '#3b6fb6', 'o', '-'),
    ('VEPYR chr22', 'vepyr_scaling', 'chr22 w{}', [1, 2, 4, 8, 10], lambda p: p, '#3b6fb6', 's', '--'),
    ('VEPYR chr1–22', 'full', 'VEPYR w{}', [1, 2, 4, 8, 10], lambda p: p, '#1f3f73', 'D', '-'),
    ('VEPYR chr1, wszystkie wtyczki', 'allplugins', 'chr1 w{}', [1, 2, 4, 8, 10], lambda p: p, '#5a9e5a', '^', '-'),
    ('VEP chr1 (--fork)', 'vep_fork', 'chr1 fork {}', ['none', '1', '2', '4', '8', '10'], lambda k: 1 if k == 'none' else int(k) + 1, '#d08a2e', 'o', '-'),
    ('VEP chr22 (--fork)', 'vep_fork', 'chr22 fork {}', ['none', '1', '2', '4', '8', '10'], lambda k: 1 if k == 'none' else int(k) + 1, '#d08a2e', 's', '--'),
]
for name, grp, lab, keys, xf, col, mk, ls in series:
    t1 = D[(grp, lab.format(keys[0]))]['mean']
    xs, ys, lo, hi = [], [], [], []
    for k in keys:
        o = D[(grp, lab.format(k))]
        xs.append(xf(k)); s = t1 / o['mean']; ys.append(s)
        lo.append(s - t1 / o['max']); hi.append(t1 / o['min'] - s)
    a.errorbar(xs, ys, yerr=[lo, hi], color=col, marker=mk, ls=ls, ms=3.5, lw=1, capsize=2, label=name)
a.plot([1, 11], [1, 11], color='black', lw=0.6, ls=':', label='skalowanie idealne')
a.axvspan(4.5, 11.5, color='0.92', zorder=0)
a.text(4.7, 10.2, 'powyżej 4 rdzeni\nwydajnościowych', fontsize=6.5, color='0.35', va='top')
a.set_xlabel('liczba procesów VEP / wykonawców VEPYR'); a.set_ylabel('przyspieszenie $S_p$')
a.set_xlim(0.5, 11.5); a.set_ylim(0, 11); a.set_xticks([1, 2, 3, 4, 5, 8, 9, 10, 11])
a.legend(fontsize=6.2, loc='upper left', frameon=False); a.set_title('(a)', fontsize=9, loc='left')
vep_t = D[('full', 'VEP none')]['mean']
b.scatter([vep_t], [2.73], color='#d08a2e', marker='o', zorder=3)
b.annotate('VEP, 1 proces', (vep_t, 2.73), textcoords='offset points', xytext=(-60, 8), fontsize=7)
ws = [1, 2, 4, 8, 10]
tv = [D[('full', f'VEPYR w{w}')]['mean'] for w in ws]; rv = [S5[str(w)]['mean_rss_gib'] for w in ws]
b.plot(tv, rv, color='#1f3f73', marker='D', ms=3.5, lw=1)
for w, t, r in zip(ws, tv, rv):
    b.annotate(f'w={w}', (t, r), textcoords='offset points', xytext=(4, -3), fontsize=6.5)
b.set_xscale('log'); b.set_xlabel('czas adnotacji chr1–22 [s] (skala log.)'); b.set_ylabel('szczytowe RSS [GiB]')
b.set_ylim(0, 10); b.set_title('(b)', fontsize=9, loc='left')
fig.tight_layout(); fig.savefig(sys.argv[1]); print('ok', sys.argv[1])
