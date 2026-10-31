#!/usr/bin/env python3
"""Record-by-record difference of two annotated VCFs with the same CSQ format (thesis item B17, 2026-10-31).

Used for VEP 116 on chr22 without --fork and with --fork 10: counts differing records and CSQ
entries, the CSQ fields (by name) in which they differ, and the (old value, new value, SYMBOL,
Feature) combinations.

Usage: csq_record_diff.py A.vcf B.vcf OUT_JSON
"""
import json, re, sys
from collections import Counter


def body(path):
    fmt = None
    for line in open(path):
        if line.startswith("##INFO=<ID=CSQ"):
            fmt = re.search(r"Format: ([^\"]+)", line).group(1).split("|")
        elif not line.startswith("#"):
            yield fmt, line.rstrip("\n").split("\t")


def main():
    a, b, out = sys.argv[1:4]
    n = 0; recs = []; fields = Counter(); combos = Counter(); entries = 0; other = 0
    for (fmt, x), (_, y) in zip(body(a), body(b)):
        n += 1
        if x == y:
            continue
        recs.append(f"{x[0]}:{x[1]}")
        if x[:7] != y[:7] or x[8:] != y[8:]:
            other += 1; continue
        ca = next(t for t in x[7].split(";") if t.startswith("CSQ="))[4:].split(",")
        cb = next(t for t in y[7].split(";") if t.startswith("CSQ="))[4:].split(",")
        if len(ca) != len(cb):
            other += 1; continue
        for e1, e2 in zip(ca, cb):
            if e1 != e2:
                entries += 1
                s1, s2 = e1.split("|"), e2.split("|")
                for j, (v1, v2) in enumerate(zip(s1, s2)):
                    if v1 != v2:
                        fields[fmt[j]] += 1
                        combos[f"{v1} -> {v2} | {s1[fmt.index('SYMBOL')]} | {s1[fmt.index('Feature')]}"] += 1
    res = {"a": a, "b": b, "records": n, "differing_records": len(recs), "differing_csq_entries": entries,
           "records_differing_outside_csq_values": other, "fields": dict(fields), "combinations": dict(combos.most_common()),
           "first_last_differing": recs[:1] + recs[-1:]}
    open(out, "w").write(json.dumps(res, indent=1) + "\n")
    print(json.dumps({k: v for k, v in res.items() if k != "combinations"}), "\n", json.dumps(res["combinations"]))


if __name__ == "__main__":
    main()
