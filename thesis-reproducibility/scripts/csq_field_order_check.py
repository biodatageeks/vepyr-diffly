#!/usr/bin/env python3
"""Do two annotated VCFs differ only in the order of the CSQ fields? (thesis item C1, 2026-10-31)

With plugins, VEPYR writes the plugin fields of CSQ in a different order than VEP, so the
output bodies differ byte by byte. This reads each file's CSQ Format from its header, maps every
CSQ entry to {field name: value} and compares record by record: all non-INFO columns, the other
INFO keys in order, and the CSQ entries in order as name -> value maps.

Usage: csq_field_order_check.py A.vcf[.gz] B.vcf[.gz] OUT_JSON
"""
import gzip, itertools, json, re, sys


def records(path):
    fmt = None
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as f:
        for line in f:
            if line.startswith("##INFO=<ID=CSQ"):
                fmt = re.search(r"Format: ([^\"]+)", line).group(1).split("|")
            if line.startswith("#"):
                continue
            col = line.rstrip("\n").split("\t")
            info = col[7].split(";")
            csq = next(x[4:] for x in info if x.startswith("CSQ="))
            yield fmt, col[:7] + col[8:], [x for x in info if not x.startswith("CSQ=")], \
                [dict(zip(fmt, e.split("|"))) for e in csq.split(",")]


def main():
    a, b, out = sys.argv[1:4]
    n = differ = 0
    first = None
    fa = fb = None
    for ra, rb in itertools.zip_longest(records(a), records(b)):
        n += 1
        if ra is None or rb is None:  # different record counts
            differ += 1
            first = first or {"record": n, "pos": None}
            continue
        (fa, xa, ia, ca), (fb, xb, ib, cb) = ra, rb
        if (xa, ia, ca) != (xb, ib, cb):
            differ += 1
            first = first or {"record": n, "pos": xa[:2]}
    res = {"a": a, "b": b, "records_compared": n, "records_differing_after_name_mapping": differ,
           "first_difference": first, "fields_a": len(fa), "fields_b": len(fb),
           "same_field_set": set(fa) == set(fb), "same_field_order": fa == fb,
           "only_a": sorted(set(fa) - set(fb)), "only_b": sorted(set(fb) - set(fa))}
    with open(out, "w") as f:
        f.write(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res))


if __name__ == "__main__":
    main()
