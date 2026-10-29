#!/usr/bin/env python3
"""What the validation data actually exercised (thesis, B6): counts over an annotated
VCF of variant classes, indel lengths, SO consequence terms, IMPACT, 3'-shifted HGVS
(HGVS_OFFSET), RefSeq alignment edits (BAM_EDIT, REFSEQ_MATCH, GIVEN_REF/USED_REF),
incomplete CDS flags, empty Existing_variation and the number of non-empty values per
CSQ field. The VEPYR chr1-22 baseline output used here is byte-identical (body) to the
VEP output, so the counts describe both.

Usage: csq_coverage.py ANNOTATED_VCF OUT_JSON
"""
import json, sys
from collections import Counter


def main():
    path, out = sys.argv[1], sys.argv[2]
    fields = None
    vclass, indel_len, cons, impact, flags = Counter(), Counter(), Counter(), Counter(), Counter()
    nonempty = Counter(); n_rec = n_ent = 0
    special = Counter()
    with open(path) as f:
        for line in f:
            if line.startswith("##INFO=<ID=CSQ"):
                fields = line.split("Format: ")[1].rstrip('">\n').split("|")
                ix = {k: i for i, k in enumerate(fields)}
                continue
            if line[0] == "#":
                continue
            c = line.split("\t", 8)
            ref, alt = c[3], c[4]
            n_rec += 1
            if len(ref) != len(alt):
                indel_len[min(abs(len(alt) - len(ref)), 51)] += 1
            info = c[7]
            p = info.find("CSQ=")
            if p < 0:
                special["records_without_csq"] += 1
                continue
            csq = info[p + 4:].split(";", 1)[0]
            first = True
            for ent in csq.split(","):
                v = ent.split("|")
                n_ent += 1
                if first and "VARIANT_CLASS" in ix:
                    vclass[v[ix["VARIANT_CLASS"]]] += 1
                first = False
                for t in v[ix["Consequence"]].split("&"):
                    cons[t] += 1
                impact[v[ix["IMPACT"]]] += 1
                for k, i in ix.items():
                    if i < len(v) and v[i] not in ("", "."):
                        nonempty[k] += 1
                for fl in (v[ix["FLAGS"]].split("&") if "FLAGS" in ix and v[ix["FLAGS"]] else []):
                    flags[fl] += 1
                for k in ("HGVS_OFFSET", "BAM_EDIT", "GIVEN_REF", "USED_REF", "REFSEQ_MATCH"):
                    if k in ix and v[ix[k]] not in ("", "."):
                        special[f"entries_with_{k}"] += 1
                if "GIVEN_REF" in ix and "USED_REF" in ix and v[ix["GIVEN_REF"]] != v[ix["USED_REF"]]:
                    special["entries_given_ref_ne_used_ref"] += 1
                if "Existing_variation" in ix and v[ix["Existing_variation"]] in ("", "."):
                    special["entries_without_existing_variation"] += 1
    res = {"records": n_rec, "csq_entries": n_ent, "variant_class": dict(vclass.most_common()),
           "indel_length": {str(k): v for k, v in sorted(indel_len.items())}, "consequence_terms": dict(cons.most_common()),
           "impact": dict(impact.most_common()), "flags": dict(flags.most_common()), "special": dict(special),
           "nonempty_per_field": {k: nonempty.get(k, 0) for k in fields}}
    json.dump(res, open(out, "w"), indent=1)
    print(json.dumps({k: res[k] for k in ("records", "csq_entries", "variant_class", "impact", "special")}))
    print(len(cons), "SO terms;", sorted(cons.items(), key=lambda kv: kv[1])[:12])


if __name__ == "__main__":
    main()
