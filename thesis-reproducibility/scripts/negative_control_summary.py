#!/usr/bin/env python3
"""Summarise negative_control.py output into one table (thesis, 2026-10-29).
Usage: negative_control_summary.py OUT_DIR  -> OUT_DIR/summary.json and summary.tsv
For each class: what was injected and what each comparator reported. E2E counters are
those of compare_vcfs (VEP side = original); vepyr-diffly counters are left_only /
right_only / unequal rows of the variant and consequence tiers (left = original)."""
import json, sys
from pathlib import Path
out = Path(sys.argv[1])
rows = []
for d in sorted(p for p in out.iterdir() if (p / 'result.json').exists()):
    r = json.loads((d / 'result.json').read_text())
    e = r['e2e']
    strict = {k: v for k, v in e['field_mismatch_counts'].items() if v}
    fmt = {k: v for k, v in e['field_format_mismatch_counts'].items() if v}
    order = {k: v for k, v in e['field_order_mismatch_counts'].items() if v}
    def tiers(x):
        if 'tiers' not in x:
            return {'error': x}
        return {t: {k: x['tiers'][t][k] for k in ('left_only_rows', 'right_only_rows', 'unequal_rows')} for t in ('variant', 'consequence')}
    rows.append({
        'class': r['class'], 'expected': r['expected'],
        'e2e': {'variants_only_in_vep': e['variants_only_in_vep'], 'variants_only_in_vepyr': e['variants_only_in_vepyr'],
                'csq_count_mismatch': e['csq_entry_count_mismatch'], 'entries_only_in_vep': e['csq_entries_only_in_vep'],
                'entries_only_in_vepyr': e['csq_entries_only_in_vepyr'], 'csq_order_mismatch': e['csq_order_mismatch'],
                'field_mismatch': strict, 'format_only': fmt, 'order_only': order,
                'ledger_rows': e.get('mismatch_ledger', {}).get('rows') if isinstance(e.get('mismatch_ledger'), dict) else None},
        'diffly_full': tiers(r['diffly_full']), 'diffly_plugins': tiers(r['diffly_plugins']), 'seconds': r['seconds']})
(out / 'summary.json').write_text(json.dumps(rows, indent=1) + '\n')
with (out / 'summary.tsv').open('w') as f:
    f.write('class\texpected\te2e\tdiffly_full\tdiffly_plugins\n')
    for x in rows:
        f.write('\t'.join([x['class'], json.dumps(x['expected']), json.dumps(x['e2e']), json.dumps(x['diffly_full']), json.dumps(x['diffly_plugins'])]) + '\n')
for x in rows:
    e = x['e2e']
    e_short = {k: v for k, v in e.items() if v not in (0, {}, None)}
    print(x['class'], '| exp', x['expected'], '| E2E', e_short, '| DF', x['diffly_full'], '| DP', x['diffly_plugins'])
