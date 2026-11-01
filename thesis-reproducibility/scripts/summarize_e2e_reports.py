#!/usr/bin/env python3
"""Per-chromosome summary of run_comparison.py reports (thesis item B12, 2026-11-01).

Reads fast_chr<N>_<PROFILE>_116_report.json (and the matching mismatches.jsonl) from REPORT_DIR
for chr1-22 and prints one JSON document: per chromosome the compared units, one-sided units,
CSQ entry count mismatches, one-sided CSQ entries, CSQ order mismatches, number of compared fields,
fields at 100%, summed field value / order / format mismatches, mismatch ledger lines and the
annotation time, plus totals.

Usage: summarize_e2e_reports.py REPORT_DIR PROFILE > summary.json
"""
import json
import sys
from pathlib import Path


def main() -> None:
    d, profile = Path(sys.argv[1]), sys.argv[2]
    rows = []
    for n in range(1, 23):
        r = json.loads((d / f"fast_chr{n}_{profile}_116_report.json").read_text())
        c = r["comparison"]
        rates = c["field_match_rates"]
        led = d / f"fast_chr{n}_{profile}_116_mismatches.jsonl"
        rows.append({
            "chrom": r["chrom"],
            "units": c["variants_compared"],
            "one_sided": c["variants_only_in_vepyr"] + c["variants_only_in_vep"],
            "csq_count_mismatch": c["csq_entry_count_mismatch"],
            "csq_entries_one_sided": c["csq_entries_only_in_vepyr"] + c["csq_entries_only_in_vep"],
            "csq_order_mismatch": c["csq_order_mismatch"],
            "fields": len(rates),
            "fields_100": sum(1 for v in rates.values() if v == 100.0),
            "field_value_mismatches": sum(c["field_mismatch_counts"].values()),
            "field_order_mismatches": sum(c["field_order_mismatch_counts"].values()),
            "field_format_mismatches": sum(c["field_format_mismatch_counts"].values()),
            "fields_only_one_side": len(c["fields_only_in_vepyr"]) + len(c["fields_only_in_vep"]),
            "ledger_lines": sum(1 for _ in led.open()) if led.exists() else None,
            "annotation_s": r["annotation"]["time_s"],
            "vepyr_rev": r["build"]["vepyr_rev"],
            "vepyr_dirty": r["build"]["vepyr_dirty"],
            "plugin_cache_path": r.get("plugin_cache_path"),
            "reference_header_time": r["reference_identity"]["header"].split("time=\"")[1].split("\"")[0],
        })
    keys = ["units", "one_sided", "csq_count_mismatch", "csq_entries_one_sided", "csq_order_mismatch",
            "field_value_mismatches", "field_order_mismatches", "field_format_mismatches", "fields_only_one_side",
            "ledger_lines", "annotation_s"]
    total = {k: sum(x[k] or 0 for x in rows) for k in keys}
    print(json.dumps({"profile": profile, "report_dir": str(d), "chromosomes": rows, "total": total}, indent=1))


if __name__ == "__main__":
    main()
