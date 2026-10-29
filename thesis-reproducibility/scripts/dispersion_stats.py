#!/usr/bin/env python3
"""Dispersion of every timing configuration reported in the thesis (appendix A,
tab:repro-dispersion): n, mean, median, sample standard deviation (n-1), min, max and
coefficient of variation, computed from the same raw files as the means.

Usage (from thesis-reproducibility/):  python3 scripts/dispersion_stats.py
Writes results/dispersion-2026-10-29/dispersion.json and dispersion.tsv.

Repeat 1 of the chr1/chr22 VEPYR worker-scaling series is taken unrounded from
chapter5-m5-.../benchmarks/vepyr_chr*/workers*.metrics.json (see README)."""
import os
from pathlib import Path
os.chdir(Path(__file__).resolve().parent.parent)
import json, statistics as st
R='results/'
A9=R+'annotation-2026-09/'
CH5=A9+'chapter5-m5-20260924T061203Z/benchmarks/'
rows=[]  # (group, label, values)
def add(g,l,v): rows.append((g,l,[float(x) for x in v]))

s=json.load(open(A9+'5repeats_m5/summary_5repeats.json'))
for chrom in ('chr1','chr22'):
    for p in ('1','2','4','8','10'):
        vals=[r['annotation_seconds'] for r in s[f'{chrom}_scaling'][p]]
        m=json.load(open(CH5+f'vepyr_{chrom}/workers{p}.metrics.json'))
        vals[0]=m['annotation_seconds']
        add('vepyr_scaling',f'{chrom} w{p}',vals)
v=json.load(open(A9+'vep_5repeats_m5/summary_vep_5repeats.json'))
for chrom in ('chr1','chr22'):
    for f in ('none','1','2','4','8','10'):
        add('vep_fork',f'{chrom} fork {f}',v[f'{chrom}_forks'][f])
    for pl in ('clinvar','alphamissense','dbnsfp','spliceai','cadd'):
        add('vep_plugins',f'{chrom} {pl}',v[f'{chrom}_plugins'][pl])
full=json.load(open(R+'vep-full-genome-2026-10/vep_full_repeats_m5_20261023/summary.json'))
add('full','VEP none',[r['wall_s'] for r in full['runs'] if r['repeat']>=2])
s5=[json.loads(l) for l in open(R+'completion-2026-10/m5_completion_20261024/S5/raw_runs.jsonl')]
for p in (1,2,4,8,10):
    add('full',f'VEPYR w{p}',[r['annotation_seconds'] for r in s5 if r['workers']==p])
ps=[json.loads(l) for l in open(R+'plugin-series-2026-10-25/raw_runs.jsonl')]
for chrom in ('chr1','chr22'):
    for c in ('baseline','clinvar','alphamissense','dbnsfp','spliceai','cadd'):
        add('vepyr_plugins',f'{chrom} {c}',[r['annotation_seconds'] for r in ps if r['chrom']==chrom and r['config']==c])
sc=[json.loads(l) for l in open(R+'annotation-2026-10/plugin_scaling_m5_5repeats_20261013/runs.jsonl')]
for chrom in ('1','22'):
    for p in (1,2,4,8,10):
        vals=[r['annotation_s'] for r in sc if r['chrom']==chrom and r['workers']==p and r['returncode']==0]
        if vals: add('allplugins',f'chr{chrom} w{p}',vals)
s6=[json.loads(l) for l in open(R+'completion-2026-10/m5_completion_20261024/S6/raw_runs.jsonl')]
for sc_ in ('baseline','all_plugins'):
    for n in (1,2,4,8):
        add('multiproc',f'{sc_} n{n}',[r['group_wall_s'] for r in s6 if r['scenario']==sc_ and r['n']==n])
C=R+'cache-build-2026-10/'
add('cache','baseline 9cpu24gb',[x['seconds_wall'] for x in json.load(open(C+'full_cache_20/baseline_9cpu24gb_5_summary.json'))])
add('cache','baseline 10cpu',[x['seconds_wall'] for x in json.load(open(C+'full_cache_20/baseline_10cpu_5_summary.json'))])
add('cache','native seq',[x['wall_seconds'] for x in json.load(open(C+'full_core_5x/summary.json'))])
add('cache','10 containers',[x['wall_seconds'] for x in json.load(open(C+'full_cache_20/distributed_5_summary.json'))])
add('cache','queue 10',[x['seconds_wall'] for x in json.load(open(C+'parallel_full_cache/queue_20_10_repeats/summary.json'))['repeats']])
add('cache','full build',[json.load(open(R+f'completion-2026-10/m5_completion_20261024/S11/repeat_{i}/summary.json'))['wall_seconds'] for i in range(1,6)])
out=[]
for g,l,v in rows:
    m=st.mean(v); sd=st.stdev(v)
    out.append(dict(group=g,label=l,n=len(v),mean=m,median=st.median(v),sd=sd,min=min(v),max=max(v),cv=100*sd/m,cv_pop=100*st.pstdev(v)/m))
for o in out:
    print(f"{o['group']:14s} {o['label']:22s} n={o['n']} mean={o['mean']:10.3f} med={o['median']:10.3f} sd={o['sd']:8.3f} min={o['min']:10.3f} max={o['max']:10.3f} cv={o['cv']:5.2f} cvp={o['cv_pop']:5.2f}")

out_dir=Path('results/dispersion-2026-10-29'); out_dir.mkdir(parents=True, exist_ok=True)
json.dump(out,open(out_dir/'dispersion.json','w'),indent=1)
with open(out_dir/'dispersion.tsv','w') as t:
    t.write('group\tconfig\tn\tmean_s\tmedian_s\tsd_s\tmin_s\tmax_s\tcv_pct\n')
    for o in out:
        t.write(f"{o['group']}\t{o['label']}\t{o['n']}\t{o['mean']:.4f}\t{o['median']:.4f}\t{o['sd']:.4f}\t{o['min']:.4f}\t{o['max']:.4f}\t{o['cv']:.2f}\n")
