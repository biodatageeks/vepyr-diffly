#!/usr/bin/env python3
"""95% confidence intervals of the main time ratios (appendix A, tab:ratio-ci):
ratio of means on the log scale, delta method, t(0.975, 4) = 2.776, the two tools
treated as independent samples (n = 5 each). Input: results/dispersion-2026-10-29/
dispersion.json (scripts/dispersion_stats.py). Output: results/dispersion-2026-10-29/ratio_ci.json."""
import json, math, os
from pathlib import Path
os.chdir(Path(__file__).resolve().parent.parent)
D = {(o['group'], o['label']): o for o in json.load(open('results/dispersion-2026-10-29/dispersion.json'))}
T = 2.776
def ci(a, b):
    A, B = D[a], D[b]; r = A['mean'] / B['mean']
    se = math.sqrt((A['sd'] / A['mean']) ** 2 / A['n'] + (B['sd'] / B['mean']) ** 2 / B['n'])
    return [r, r * math.exp(-T * se), r * math.exp(T * se)]
out = {}
for c in ('chr1', 'chr22'):
    for k in ('baseline', 'clinvar', 'alphamissense', 'dbnsfp', 'spliceai', 'cadd'):
        vep = ('vep_fork', c + ' fork none') if k == 'baseline' else ('vep_plugins', f'{c} {k}')
        out[f'{c} {k} 1proc'] = ci(vep, ('vepyr_plugins', f'{c} {k}'))
    out[f'{c} best'] = ci(('vep_fork', c + ' fork 10'), ('vepyr_scaling', c + ' w10'))
out['full w10'] = ci(('full', 'VEP none'), ('full', 'VEPYR w10'))
out['full w1'] = ci(('full', 'VEP none'), ('full', 'VEPYR w1'))
out['queue vs native'] = ci(('cache', 'native seq'), ('cache', 'queue 10'))
json.dump(out, open('results/dispersion-2026-10-29/ratio_ci.json', 'w'), indent=1)
for k, v in out.items(): print(k, [round(x, 2) for x in v])
