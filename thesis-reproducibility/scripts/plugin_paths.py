#!/usr/bin/env python3
"""Which plugin code paths and semantic limits the validation data exercised
(thesis, C4/C6). Over an annotated VCF with all five plugins counts: variants and
CSQ entries with each plugin's values; dbNSFP values holding '&' lists (several
transcripts); variants whose SpliceAI values name more than one gene; variants with
AlphaMissense values on more than one transcript and on non-MANE transcripts; SNVs
and indels without CADD; indels with SpliceAI values (the SNV-only source gives none).
Usage: plugin_paths.py ANNOTATED_VCF"""
import re, sys, json
from collections import Counter
path = sys.argv[1]
c = Counter()
for line in open(path):
    if line.startswith('##INFO=<ID=CSQ'):
        F = line.split('Format: ')[1].rstrip('">\n').split('|'); ix = {k: i for i, k in enumerate(F)}
        dbn = [k for k in ("SIFT4G_score,SIFT4G_pred,Polyphen2_HDIV_score,Polyphen2_HVAR_score,MutationTaster_score,"
                           "MutationTaster_pred,PROVEAN_score,PROVEAN_pred,VEST4_score,MetaSVM_score,MetaSVM_pred,"
                           "MetaLR_score,MetaLR_pred,REVEL_score,GERP++_RS,phyloP100way_vertebrate,"
                           "phastCons100way_vertebrate,CADD_raw,CADD_phred").split(",") if k in ix]
        continue
    if line[0] == '#': continue
    f = line.split('\t', 8); ref, alt = f[3], f[4]
    snv = len(ref) == 1 and len(alt) == 1
    c['variants'] += 1; c['snv' if snv else 'indel'] += 1
    m = re.search(r'CSQ=([^;\t]*)', f[7])
    if not m: continue
    ents = [e.split('|') for e in m.group(1).split(',')]
    g = lambda e, k: e[ix[k]] if k in ix and ix[k] < len(e) else ''
    cadd = any(g(e, 'CADD_PHRED') for e in ents)
    c['snv_with_cadd' if snv and cadd else 'snv_without_cadd' if snv else 'indel_with_cadd' if cadd else 'indel_without_cadd'] += 1
    sa_syms = {g(e, 'SpliceAI_pred_SYMBOL') for e in ents if g(e, 'SpliceAI_pred_SYMBOL')}
    if sa_syms: c['variants_with_spliceai'] += 1; c['indel_with_spliceai'] += (not snv)
    if len(sa_syms) > 1: c['variants_spliceai_multi_gene'] += 1
    am = [e for e in ents if g(e, 'am_pathogenicity')]
    if am:
        c['variants_with_am'] += 1
        if len({g(e, 'Feature') for e in am}) > 1: c['variants_am_multi_transcript'] += 1
        if any(not g(e, 'MANE_SELECT') and not g(e, 'Feature').startswith('NM_') for e in am): pass
        c['am_entries'] += len(am); c['am_entries_non_mane'] += sum(1 for e in am if not g(e, 'MANE_SELECT'))
    dbv = [e for e in ents if any(g(e, k) for k in dbn)]
    if dbv:
        c['variants_with_dbnsfp'] += 1
        if any('&' in g(e, k) for e in dbv for k in dbn): c['variants_dbnsfp_list_values'] += 1
    if any(g(e, 'ClinVar_CLNSIG') for e in ents): c['variants_with_clinvar'] += 1
print(json.dumps(dict(c), indent=1))
