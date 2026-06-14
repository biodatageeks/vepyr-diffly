#!/usr/bin/env python3
"""Verify two VEPYR built caches are logically identical.

Correctness gate for the cache-parallelization plan: a parallel build MUST
produce the same logical cache content as the serial baseline (only then is VEP
parity preserved). Parquet bytes legitimately differ (row-group boundaries,
compression), so we compare LOGICAL content, not bytes:

  1. identical set of output files (relative paths under the version dir)
  2. identical row count per file
  3. identical order-independent multiset content hash per file
     (polars hash_rows -> sort hashes -> sha256), robust to row ordering

Exit code 0 == identical, 1 == differences found.

Usage:
  compare_cache.py --baseline <cacheA> --candidate <cacheB>
where each path is the cache_dir passed to build_cache (we descend into
<cache_dir>/parquet/<version>/ automatically, or accept the version dir directly).
"""

from __future__ import annotations

import argparse
import hashlib
import io
import sys
from pathlib import Path

import polars as pl


def find_version_dir(cache_dir: Path) -> Path:
    """Resolve to the directory holding the entity subdirs."""
    parquet_root = cache_dir / "parquet"
    search_roots = [parquet_root] if parquet_root.is_dir() else [cache_dir]
    for root in search_roots:
        # version dir like 115_GRCh38_ensembl, or cache_dir itself may be it
        candidates = [p for p in root.iterdir() if p.is_dir()] if root.is_dir() else []
        for cand in candidates:
            if (cand / "variation").is_dir() or (cand / "transcript").is_dir():
                return cand
        if (root / "variation").is_dir() or (root / "transcript").is_dir():
            return root
    raise FileNotFoundError(f"could not locate entity dirs under {cache_dir}")


def list_parquet(version_dir: Path) -> dict[str, Path]:
    """Map relative path -> absolute path for every .parquet under version_dir."""
    return {
        str(p.relative_to(version_dir)): p
        for p in sorted(version_dir.rglob("*.parquet"))
    }


def content_fingerprint(path: Path) -> tuple[int, str]:
    """Return (row_count, order-independent sha256 of sorted row hashes + schema)."""
    df = pl.read_parquet(path)
    schema_repr = ";".join(f"{n}:{t}" for n, t in zip(df.columns, df.dtypes))
    # u64 hash per row (content-deterministic), sorted -> order-independent.
    # Serialize via Arrow IPC (no numpy dependency) and digest the bytes.
    row_hashes = df.hash_rows().sort()
    buf = io.BytesIO()
    row_hashes.to_frame("h").write_ipc(buf, compression="uncompressed")
    digest = hashlib.sha256()
    digest.update(schema_repr.encode())
    digest.update(b"\x00")
    digest.update(buf.getvalue())
    return df.height, digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    args = parser.parse_args()

    base_dir = find_version_dir(args.baseline)
    cand_dir = find_version_dir(args.candidate)
    base_files = list_parquet(base_dir)
    cand_files = list_parquet(cand_dir)

    print(f"baseline:  {base_dir}  ({len(base_files)} parquet)")
    print(f"candidate: {cand_dir}  ({len(cand_files)} parquet)")

    problems: list[str] = []

    only_base = sorted(set(base_files) - set(cand_files))
    only_cand = sorted(set(cand_files) - set(base_files))
    for rel in only_base:
        problems.append(f"missing in candidate: {rel}")
    for rel in only_cand:
        problems.append(f"extra in candidate:   {rel}")

    shared = sorted(set(base_files) & set(cand_files))
    for rel in shared:
        b_rows, b_hash = content_fingerprint(base_files[rel])
        c_rows, c_hash = content_fingerprint(cand_files[rel])
        if b_rows != c_rows:
            problems.append(f"row count differs {rel}: base={b_rows} cand={c_rows}")
        elif b_hash != c_hash:
            problems.append(f"content differs {rel}: rows={b_rows} (hash mismatch)")
        else:
            print(f"  OK  {rel}  rows={b_rows}")

    if problems:
        print("\nDIFFERENCES FOUND:", file=sys.stderr)
        for p in problems:
            print(f"  - {p}", file=sys.stderr)
        return 1
    print(f"\nIDENTICAL: {len(shared)} files match logically.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
