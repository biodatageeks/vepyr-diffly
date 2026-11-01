#!/usr/bin/env python3
"""VEP and VEPYR on the edge-case records of edge_variants.py, entry by entry (thesis item B8, 2026-11-01).

For every case of cases.json: whether each tool wrote the record, the number of CSQ entries and of
fields per entry (against the header's Format), and, for every (Feature, Allele) entry, Consequence and
HGVSc of both tools. A case is "equal" when both tools wrote the record and every entry has the same
(Allele, Feature, Consequence, HGVSc).

Usage: b8_vep_vs_vepyr.py CASES_JSON VEP_VCF VEPYR_VCF OUT_JSON
"""
import json
import sys
from pathlib import Path


def read(path: Path):
    fields, out = None, {}
    for line in path.read_text().splitlines():
        if line.startswith("##INFO=<ID=CSQ"):
            fields = line.split("Format: ")[1].rstrip('">').split("|")
        if line.startswith("#"):
            continue
        c = line.split("\t")
        csq = next((x[4:] for x in c[7].split(";") if x.startswith("CSQ=")), "")
        ents = [e.split("|") for e in csq.split(",")] if csq else []
        rows = [dict(zip(fields, e)) for e in ents]
        out.setdefault(f"{c[0]}:{c[1]}:{c[3]}:{c[4]}", []).append({
            "entries": len(ents), "fields_per_entry": sorted({len(e) for e in ents}),
            "by_feature": {f"{r.get('Feature', '')}|{r.get('Allele', '')}": [r.get("Consequence", ""), r.get("HGVSc", "")]
                           for r in rows}})
    return fields, out


def main() -> int:
    cases = json.loads(Path(sys.argv[1]).read_text())
    fv, vep = read(Path(sys.argv[2]))
    fy, vepyr = read(Path(sys.argv[3]))
    rows, seen = [], {}
    for c in cases:
        k = f"chr22:{c['pos']}:{c['ref']}:{c['alt']}"
        i = seen[k] = seen.get(k, -1) + 1
        v = vep.get(k, [None] * (i + 1))[i] if i < len(vep.get(k, [])) else None
        y = vepyr.get(k, [None] * (i + 1))[i] if i < len(vepyr.get(k, [])) else None
        row = {"case": c["case"], "key": k, "vep_written": v is not None, "vepyr_written": y is not None}
        if v and y:
            feats = sorted(set(v["by_feature"]) | set(y["by_feature"]))
            diff = {f: {"vep": v["by_feature"].get(f), "vepyr": y["by_feature"].get(f)} for f in feats
                    if v["by_feature"].get(f) != y["by_feature"].get(f)}
            row.update(vep_entries=v["entries"], vepyr_entries=y["entries"], vep_fields=v["fields_per_entry"],
                       vepyr_fields=y["fields_per_entry"], equal=not diff, differing_entries=len(diff),
                       differences=dict(list(diff.items())[:4]))
        else:
            row["equal"] = v is None and y is None
        rows.append(row)
    res = {"header_fields_vep": len(fv), "header_fields_vepyr": len(fy), "cases": len(rows),
           "equal": sum(r["equal"] for r in rows), "rows": rows}
    Path(sys.argv[4]).write_text(json.dumps(res, indent=1) + "\n")
    for r in rows:
        print(f"{r['case']:24} equal={r['equal']!s:5} vep={r['vep_written']!s:5} vepyr={r['vepyr_written']!s:5} "
              f"ent={r.get('vep_entries')}/{r.get('vepyr_entries')} fields={r.get('vep_fields')}/{r.get('vepyr_fields')} "
              f"diff={r.get('differing_entries')}")
    print(f"header fields VEP/VEPYR {len(fv)}/{len(fy)}; equal {res['equal']}/{len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
