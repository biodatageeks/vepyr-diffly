#!/usr/bin/env python3
"""Full key-uniqueness check of the CADD and SpliceAI plugin caches (thesis item B18).

The thesis checked the uniqueness of the lookup keys on a sample. Here every row of every
chromosome file of the plugin cache used in the plugin runs (data_vepyr/plugin_cache,
built 2026-08-27; or another root given as argument) is checked. The key is the
manifest's key_columns plus match_columns (CADD: start, end, allele_string; SpliceAI: the
same plus symbol); chrom is constant within a file. The files are read in batches of
key columns; rows whose start equals the last start of a batch are carried over to the
next one, so every group of equal start is checked whole while memory stays bounded.
Also checked: start is non-decreasing within each file (the order the lookup relies on).

Usage: plugin_key_uniqueness.py OUT_JSON [PLUGIN_CACHE_ROOT] [PLUGIN ...]
Run with the .venv-native interpreter (polars, pyarrow).
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import polars as pl
import pyarrow.parquet as pq

ROOT = Path("/Users/lukaszjezapkowicz/workspace/data_vepyr/plugin_cache")
BATCH = 4_000_000


def check_file(path: Path, keys: list[str]) -> dict:
    f = pq.ParquetFile(path)
    rows = dups = unsorted = 0
    carry, last = None, None
    for rb in f.iter_batches(batch_size=BATCH, columns=keys):
        df = pl.from_arrow(rb).with_columns(pl.col("start").cast(pl.Int64))  # diff on unsigned would wrap
        rows += df.height
        if last is not None and df["start"][0] < last:
            unsorted += 1
        unsorted += int((df["start"].diff().drop_nulls() < 0).sum())
        if carry is not None:
            df = pl.concat([carry, df])
        last = df["start"][-1]
        carry, df = df.filter(pl.col("start") == last), df.filter(pl.col("start") != last)
        dups += df.height - df.unique(subset=keys).height
    if carry is not None:
        dups += carry.height - carry.unique(subset=keys).height
    return {"rows": rows, "duplicate_rows": dups, "order_violations": unsorted}


def main() -> int:
    out = Path(sys.argv[1])
    root = Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT
    plugins = sys.argv[3:] or ["cadd", "spliceai"]
    res = {"root": str(root)}
    for plugin in plugins:
        d = root / "plugin" / plugin
        m = json.loads((d / "manifest.json").read_text())
        keys = [k for k in m["key_columns"] if k != "chrom"] + [c["column"] for c in m["match_columns"]]
        per, t0 = {}, time.time()
        for c in m["chroms"]:
            per[c["chrom"]] = check_file(d / c["file"], keys)
            per[c["chrom"]]["manifest_rows"] = c["rows"]
            print(plugin, c["chrom"], json.dumps(per[c["chrom"]]), round(time.time() - t0), "s", flush=True)
        res[plugin] = {"keys": keys, "files": per, "seconds": time.time() - t0,
                       "rows": sum(v["rows"] for v in per.values()),
                       "duplicate_rows": sum(v["duplicate_rows"] for v in per.values()),
                       "order_violations": sum(v["order_violations"] for v in per.values()),
                       "rows_match_manifest": all(v["rows"] == v["manifest_rows"] for v in per.values())}
        out.write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps({p: {k: v for k, v in res[p].items() if k != "files"} for p in plugins}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
