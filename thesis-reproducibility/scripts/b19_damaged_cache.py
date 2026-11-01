#!/usr/bin/env python3
"""VEPYR on a damaged core cache (thesis item B19, second part; 2026-11-01).

A chr22-only copy of the core cache (hard links to cache/116_GRCh38_merged, chr22.parquet of every
entity, with the entity's chrom_manifest.json copied unchanged) is annotated with the HG002 chr22 slice (baseline: everything, hgvs, FASTA, workers=1) as it is
and after one change at a time, made on a private copy of the file:
  intact               no change (control: the body must equal the chr22 baseline body 19efe3c7...)
  transcript_truncated transcript/chr22.parquet cut to half its size
  exon_missing         exon/chr22.parquet removed
  variation_bytes      variation/chr22.parquet with 4 KiB of zero bytes written in the middle
  transcript_bytes     the same for transcript/chr22.parquet
  translation_missing  translation_core/chr22.parquet removed
Reported per case: return code, last lines of stderr, number of output records and whether the body
equals the baseline body. The question is whether a damaged cache stops the run with an error or
changes the output silently.

Usage: b19_damaged_cache.py OUT_DIR   (.venv-native interpreter)
"""
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

D = Path("/Users/lukaszjezapkowicz/workspace/data_vepyr")
SRC = D / "cache/116_GRCh38_merged"
INPUT = D / "runs/m5_completion_20261024/inputs/HG002_normalized_chr22.vcf.gz"
FASTA = D / "input/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
BASELINE = "19efe3c74723b858"
CHILD = ("import sys, vepyr; vepyr.annotate(vcf=sys.argv[1], cache_dir=sys.argv[2], everything=True, hgvs=True, "
         "reference_fasta=sys.argv[3], workers=1, output_vcf=sys.argv[4], show_progress=False)")


def make_cache(root: Path) -> None:
    for ent in sorted(p.name for p in SRC.iterdir() if p.is_dir()):
        (root / ent).mkdir(parents=True)
        os.link(SRC / ent / "chr22.parquet", root / ent / "chr22.parquet")
        shutil.copyfile(SRC / ent / "chrom_manifest.json", root / ent / "chrom_manifest.json")


def private(path: Path) -> None:  # replace the hard link by a copy before changing it
    tmp = path.with_suffix(".tmp")
    shutil.copyfile(path, tmp)
    path.unlink()
    tmp.rename(path)


def damage(root: Path, case: str) -> str:
    if case == "transcript_truncated":
        p = root / "transcript/chr22.parquet"; private(p)
        n = p.stat().st_size
        with open(p, "r+b") as f:
            f.truncate(n // 2)
        return f"{n} -> {n // 2} bytes"
    if case in ("exon_missing", "translation_missing"):
        p = root / ("exon" if case == "exon_missing" else "translation_core") / "chr22.parquet"
        p.unlink()
        return f"removed {p.relative_to(root)}"
    if case in ("variation_bytes", "transcript_bytes"):
        p = root / ("variation" if case == "variation_bytes" else "transcript") / "chr22.parquet"; private(p)
        n = p.stat().st_size
        with open(p, "r+b") as f:
            f.seek(n // 2); f.write(b"\0" * 4096)
        return f"4096 zero bytes at offset {n // 2} of {n}"
    return "none"


def body(path: Path) -> tuple[str, int]:
    h, n = hashlib.sha256(), 0
    with open(path, "rb") as f:
        for line in f:
            if not line.startswith(b"#"):
                h.update(line); n += 1
    return h.hexdigest(), n


def main() -> int:
    out = Path(sys.argv[1]).resolve()
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for case in ("intact", "transcript_truncated", "exon_missing", "variation_bytes", "transcript_bytes",
                 "translation_missing"):
        root = out / f"cache_{case}"
        shutil.rmtree(root, ignore_errors=True)
        make_cache(root)
        change = damage(root, case)
        vcf = out / f"{case}.vcf"
        vcf.unlink(missing_ok=True)
        p = subprocess.run([sys.executable, "-c", CHILD, str(INPUT), str(root), str(FASTA), str(vcf)],
                           capture_output=True, text=True)
        row = {"case": case, "change": change, "returncode": p.returncode,
               "stderr_tail": p.stderr.strip().splitlines()[-6:]}
        if vcf.exists():
            sha, n = body(vcf)
            row.update(records=n, body_sha256=sha, body_equal_baseline=sha.startswith(BASELINE))
            vcf.unlink()
        rows.append(row)
        shutil.rmtree(root)
        print(json.dumps(row)[:600], flush=True)
    (out / "summary.json").write_text(json.dumps(rows, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
