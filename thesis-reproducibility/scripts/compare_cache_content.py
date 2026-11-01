#!/usr/bin/env python3
"""Content check of a queue-built cache against the reference build (thesis, B15).

For every (entity, contig) Parquet file of the queue entities (five in the scaling runs;
also translation_core and translation_sift when the run built them, as the E5 queue does), compares the row
count and two order-independent digests of the rows (sum of 32-bit row hashes and XOR
of 64-bit row hashes, polars hash over all columns) between the per-task outputs of
one dynamic-queue run and the core cache built sequentially by VEPYR 0.7.0 (S2).
Row order inside a file is not compared. Column sets must be equal.
With --common-columns the digests cover only the columns present on both sides (in sorted order);
the columns found on one side only are recorded, and one-sided files get their row count. Added
2026-11-01 (B15): the queue's vepyr build writes object_hash and raw_object_json (exon also _rn),
which VEPYR 0.7.0 does not.

Usage: compare_cache_content.py QUEUE_RUN_DIR REFERENCE_CACHE OUT_JSON [--common-columns]
"""
import json, sys, time
from pathlib import Path
import polars as pl

ENTITIES = ["variation", "transcript", "exon", "regulatory", "motif", "translation_core", "translation_sift"]


def digest(path: Path, only=None):
    lf = pl.scan_parquet(path)
    if only is not None:
        lf = lf.select(only)
    cols = lf.collect_schema().names()
    h = pl.struct(pl.all()).hash(seed=0)
    r = lf.select(pl.len().alias("n"), (h % (1 << 32)).cast(pl.Int64).sum().alias("s32"),
                  h.bitwise_xor().alias("x64")).collect(engine="streaming").row(0)
    return {"rows": int(r[0]), "s32": int(r[1] or 0), "x64": int(r[2] or 0), "columns": sorted(cols)}


def main():
    qdir, ref, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    common = "--common-columns" in sys.argv[4:]
    names = lambda p: pl.scan_parquet(p).collect_schema().names()
    nrows = lambda p: int(pl.scan_parquet(p).select(pl.len()).collect().item())
    res, t0 = [], time.time()
    built = {p.name for p in qdir.glob("task_*/cache/*") if p.is_dir()}
    for e in [x for x in ENTITIES if x in built]:
        qfiles = {p.name: p for p in qdir.glob(f"task_*/cache/{e}/*.parquet")}
        rfiles = {p.name: p for p in (ref / e).glob("*.parquet")}
        for name in sorted(set(qfiles) | set(rfiles)):
            row = {"entity": e, "file": name, "in_queue": name in qfiles, "in_reference": name in rfiles}
            if name in qfiles and name in rfiles:
                only = None
                if common:
                    qc, rc = names(qfiles[name]), names(rfiles[name])
                    only = sorted(set(qc) & set(rc))
                    row.update(queue_only_columns=sorted(set(qc) - set(rc)), reference_only_columns=sorted(set(rc) - set(qc)))
                a, b = digest(qfiles[name], only), digest(rfiles[name], only)
                row.update(queue=a, reference=b, identical=(a == b))
            elif common:
                row.update(rows=nrows(qfiles.get(name) or rfiles[name]))
            res.append(row)
        done = [r for r in res if r["entity"] == e]
        print(e, len(done), "files;", sum(r.get("identical", False) for r in done), "identical;",
              sum(not r["in_queue"] for r in done), "only in reference;", sum(not r["in_reference"] for r in done),
              "only in queue", round(time.time() - t0), "s", flush=True)
    out.write_text(json.dumps(res, indent=1) + "\n")


if __name__ == "__main__":
    main()
