#!/usr/bin/env python3
"""Full key-uniqueness check of the CADD and SpliceAI plugin caches (thesis item B18).

The thesis checked the uniqueness of the lookup keys on a sample. Here every row of every
chromosome file of the plugin cache used in the plugin runs (data_vepyr/plugin_cache,
built 2026-08-27; or another root given as argument) is checked. The key is the
manifest's key_columns plus match_columns (CADD: start, end, allele_string; SpliceAI: the
same plus symbol); chrom is constant within a file. The files are read in batches of
key columns; rows whose start equals the last start of a batch are carried over to the
next one, so every group of equal start is checked whole while memory stays bounded.
The rows of a file are stored by tier (warm = 0, then cold = 1; ch. 4), each tier sorted by
start, so the order is checked within a tier and every key of a later tier is also looked up
among the keys of the earlier tiers (the warm tier is small enough to keep in memory).
The first run (2026-11-02, results/plugin-key-uniqueness-2026-11-02/first_run.json) checked
the order over the whole file and found exactly one decrease per file: the tier boundary.

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
KEEP_MAX = 20_000_000  # rows of a tier kept for the cross-tier check


def check_file(path: Path, keys: list[str]) -> dict:
    """Rows are stored by tier (warm = 0 first, then cold = 1), each tier sorted by start."""
    f = pq.ParquetFile(path)
    rows = dups = unsorted = cross = 0
    tiers: list[int] = []        # tier of each run of equal tier, in file order
    carry = last = None          # equal-start group carried to the next batch; last start of the run
    earlier = None               # unique keys of the finished runs (warm is small), for the cross-tier check
    current: list | None = []    # key frames of the current run; None once it exceeds KEEP_MAX rows
    current_rows, complete = 0, True
    for rb in f.iter_batches(batch_size=BATCH, columns=keys + ["tier"]):
        df = pl.from_arrow(rb).with_columns(pl.col("start").cast(pl.Int64))  # diff on unsigned would wrap
        rows += df.height
        df = df.with_columns(run=(pl.col("tier") != pl.col("tier").shift()).fill_null(False).cum_sum())
        for _, part in df.group_by("run", maintain_order=True):
            t = int(part["tier"][0])
            part = part.drop("tier", "run")
            if not tiers or tiers[-1] != t:  # a new run starts: close the previous one
                if carry is not None:
                    dups += carry.height - carry.unique(subset=keys).height
                if tiers:
                    if current is None:
                        complete = False
                    else:
                        done = pl.concat(current).unique(subset=keys)
                        earlier = done if earlier is None else pl.concat([earlier, done]).unique(subset=keys)
                tiers.append(t)
                carry = last = None
                current, current_rows = [], 0
            if last is not None and part["start"][0] < last:
                unsorted += 1
            unsorted += int((part["start"].diff().drop_nulls() < 0).sum())
            if earlier is not None:
                cross += part.join(earlier, on=keys, how="semi").height
            if current is not None:
                current.append(part.select(keys))
                current_rows += part.height
                if current_rows > KEEP_MAX:
                    current = None
            if carry is not None:
                part = pl.concat([carry, part])
            last = part["start"][-1]
            carry, part = part.filter(pl.col("start") == last), part.filter(pl.col("start") != last)
            dups += part.height - part.unique(subset=keys).height
    if carry is not None:
        dups += carry.height - carry.unique(subset=keys).height
    return {"rows": rows, "duplicate_rows_within_tier": dups, "order_violations_within_tier": unsorted,
            "tiers_in_order": tiers, "keys_also_in_earlier_tier": cross, "cross_tier_check_complete": complete}


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
                       "duplicate_rows_within_tier": sum(v["duplicate_rows_within_tier"] for v in per.values()),
                       "order_violations_within_tier": sum(v["order_violations_within_tier"] for v in per.values()),
                       "keys_also_in_earlier_tier": sum(v["keys_also_in_earlier_tier"] for v in per.values()),
                       "tiers_in_order": sorted({tuple(v["tiers_in_order"]) for v in per.values()}),
                       "cross_tier_check_complete": all(v["cross_tier_check_complete"] for v in per.values()),
                       "rows_match_manifest": all(v["rows"] == v["manifest_rows"] for v in per.values())}
        out.write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps({p: {k: v for k, v in res[p].items() if k != "files"} for p in plugins}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
