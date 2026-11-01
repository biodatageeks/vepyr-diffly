#!/usr/bin/env python3
"""Exon files of the queue cache that differ from the reference on common columns (thesis B15, 2026-11-01).

For every exon file marked not identical by `compare_cache_content.py --common-columns`, recomputes the
order-independent row digests without the source_file column (path of the VEP raw-cache file a record
was read from) and counts the positions at which the sorted source_file columns differ (a rough measure, not a row count).
Run with the .venv-native interpreter (polars).

Usage: b15_exon_source_file.py CONTENT_COMMON_JSON QUEUE_RUN_DIR REFERENCE_CACHE OUT_JSON
"""
import json, sys
from pathlib import Path
import polars as pl


def digest(path: Path, cols):
    h = pl.struct([pl.col(c) for c in cols]).hash(seed=0)
    return pl.scan_parquet(path).select(pl.len().alias("n"), (h % (1 << 32)).cast(pl.Int64).sum().alias("s32"),
                                        h.bitwise_xor().alias("x64")).collect(engine="streaming").row(0)


def main():
    res, run, ref, out = json.loads(Path(sys.argv[1]).read_text()), Path(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4])
    rows = []
    for f in [x["file"] for x in res if x["entity"] == "exon" and x.get("identical") is False]:
        qp, rp = next(run.glob(f"task_*/cache/exon/{f}")), ref / "exon" / f
        cols = sorted((set(pl.scan_parquet(qp).collect_schema().names()) & set(pl.scan_parquet(rp).collect_schema().names()))
                      - {"source_file"})
        a, b = digest(qp, cols), digest(rp, cols)
        sq = pl.read_parquet(qp, columns=["source_file"])["source_file"].sort()
        sr = pl.read_parquet(rp, columns=["source_file"])["source_file"].sort()
        rows.append({"file": f, "rows": a[0], "identical_without_source_file": a == b,
                     "source_file_sorted_positions_differing": int((sq != sr).sum())})
    summary = {"files": len(rows), "identical_without_source_file": sum(r["identical_without_source_file"] for r in rows),
               "rows": sum(r["rows"] for r in rows), "source_file_sorted_positions_differing": sum(r["source_file_sorted_positions_differing"] for r in rows)}
    out.write_text(json.dumps({"summary": summary, "files": rows}, indent=1) + "\n")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
