#!/usr/bin/env python3
"""A second, harder validation set (thesis item B7, 2026-10-30).

The E2E validation in the thesis uses GIAB HG002, mostly common SNVs. Two more inputs:
  clinvar  ClinVar GRCh38 20260706 (plugin_input/clinvar/clinvar.vcf.gz, the plugin source),
           chromosomes 1-22: every record with one ALT of A/C/G/T bases and a REF of A/C/G/T
           bases (clinical variants, many indels and splice-region variants);
  panel    synthetic exon-boundary panel on chr22: for every exon of every protein_coding
           transcript in Ensembl 116 GFF3 (raw/d20/Homo_sapiens.GRCh38.116.gff3.gz), positions
           -3..+3 around both exon ends; at each position an SNV, a 1-bp and a 3-bp deletion
           and a 1-bp insertion, REF from the primary-assembly FASTA (left-aligned VCF form with
           the anchor base); duplicates removed.
Both are annotated by VEP 116 (merged, --everything --hgvs --fasta, --fork 10) and VEPYR
0.7.0 (everything, hgvs, FASTA, workers=10), then compared with the E2E comparator of the
vepyr repository (comparison.compare.compare_vcfs) and the body SHA-256. VEP --fork 10 is
used for time; B17/D9 check that it gives the same body as a single process.

Usage: second_dataset.py OUT_DIR [clinvar|panel ...]   (Colima running)
"""
from __future__ import annotations

import contextlib
import gzip
import hashlib
import io
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

D = Path("/Users/lukaszjezapkowicz/workspace/data_vepyr")
PYN = os.environ.get("PYN", str(D / ".venv-native/bin/python3"))
E2E_DIR = D / "tools/vepyr/e2e-testing/scripts"
VEP_CACHE = D / "homo_sapiens_merged/116_GRCh38"
VEPYR_CACHE = D / "cache/116_GRCh38_merged"
FASTA = D / "input/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
CLINVAR = D / "plugin_input/clinvar/clinvar.vcf.gz"
GFF3 = D / "raw/d20/Homo_sapiens.GRCh38.116.gff3.gz"
IMAGE = "ensemblorg/ensembl-vep:release_116.0"
BASES = re.compile(r"[ACGT]+")


class Fasta:
    def __init__(self, path: Path):
        self.f = open(path, "rb")
        self.idx = {}
        for line in open(str(path) + ".fai"):
            name, length, offset, lb, lw = line.split("\t")[:5]
            self.idx[name] = (int(length), int(offset), int(lb), int(lw))

    def seq(self, chrom: str, start: int, end: int) -> str:
        """1-based inclusive."""
        length, offset, lb, lw = self.idx[chrom]
        out = []
        for p in range(start - 1, min(end, length)):
            self.f.seek(offset + (p // lb) * lw + p % lb)
            out.append(self.f.read(1))
        return b"".join(out).decode().upper()


def write_vcf(path: Path, contigs: list[str], rows) -> int:
    n = 0
    with path.open("w") as f:
        f.write("##fileformat=VCFv4.2\n" + "".join(f"##contig=<ID={c}>\n" for c in contigs))
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, vid, ref, alt in rows:
            f.write(f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t.\tPASS\t.\n")
            n += 1
    gz = Path(str(path) + ".gz")
    subprocess.run(f"bgzip -f {path} && tabix -f -p vcf {gz}", shell=True, check=True)
    return n


def clinvar_input(out: Path) -> tuple[Path, dict]:
    chroms = [str(i) for i in range(1, 23)]
    rows, skipped = {c: [] for c in chroms}, 0
    with gzip.open(CLINVAR, "rt") as f:
        for line in f:
            if line[0] == "#":
                continue
            c = line.split("\t", 6)
            if c[0] not in rows:
                continue
            if BASES.fullmatch(c[3]) and BASES.fullmatch(c[4]):
                rows[c[0]].append((f"chr{c[0]}", int(c[1]), c[2], c[3], c[4]))
            else:
                skipped += 1
    allrows = [r for ch in chroms for r in sorted(set(rows[ch]), key=lambda r: (r[1], r[3], r[4]))]
    n = write_vcf(out / "clinvar_chr1_22.vcf", [f"chr{c}" for c in chroms], allrows)
    snv = sum(len(r[3]) == 1 == len(r[4]) for r in allrows)
    return out / "clinvar_chr1_22.vcf.gz", {"records": n, "snv": snv, "indel_or_mnv": n - snv,
                                            "skipped_non_acgt_or_multi_alt": skipped}


def panel_input(out: Path) -> tuple[Path, dict]:
    pc, exons = set(), []
    with gzip.open(GFF3, "rt") as f:
        for line in f:
            if not line.startswith("22\t"):
                continue
            c = line.rstrip("\n").split("\t")
            attrs = dict(kv.split("=", 1) for kv in c[8].split(";") if "=" in kv)
            if attrs.get("biotype") == "protein_coding" and attrs.get("ID", "").startswith("transcript:"):
                pc.add(attrs["ID"])
            elif c[2] == "exon":
                exons.append((attrs.get("Parent"), int(c[3]), int(c[4])))
    fa = Fasta(FASTA)
    ends = sorted({p for parent, s, e in exons if parent in pc for p in (s, e)})
    vars_ = set()
    for b in ends:
        for p in range(b - 3, b + 4):
            ctx = fa.seq("22", p - 1, p + 3)  # p-1 .. p+3
            if len(ctx) < 5 or "N" in ctx:
                continue
            anchor, base = ctx[0], ctx[1]
            vars_.add((p, base, "ACGT"[("ACGT".index(base) + 1) % 4]))   # SNV
            vars_.add((p - 1, anchor + base, anchor))                    # 1-bp deletion of p
            vars_.add((p - 1, anchor + ctx[1:4], anchor))                # 3-bp deletion p..p+2
            vars_.add((p, base, base + "T"))                             # 1-bp insertion after p
    rows = [("chr22", p, ".", r, a) for p, r, a in sorted(vars_)]
    n = write_vcf(out / "panel_chr22.vcf", ["chr22"], rows)
    return out / "panel_chr22.vcf.gz", {"records": n, "protein_coding_transcripts": len(pc),
                                        "exon_boundaries": len(ends)}


def body(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for line in f:
            if not line.startswith(b"#"):
                h.update(line)
    return h.hexdigest()


def annotate_both(inp: Path, out: Path, tag: str) -> dict:
    vep_vcf, vepyr_vcf = out / f"{tag}_vep.vcf", out / f"{tag}_vepyr.vcf"
    cmd = ["docker", "run", "--rm", "--user", "501:20", "--env", "HOME=/tmp",
           "-v", f"{VEP_CACHE}:/opt/vep/.vep/homo_sapiens_merged/116_GRCh38:ro",
           "-v", f"{inp.parent}:/input:ro", "-v", f"{FASTA.parent}:/refdir:ro", "-v", f"{out}:/output",
           IMAGE, "vep", "--dir", "/opt/vep/.vep", "--cache", "--merged", "--offline", "--assembly", "GRCh38",
           "--input_file", f"/input/{inp.name}", "--output_file", f"/output/{vep_vcf.name}", "--vcf",
           "--force_overwrite", "--no_stats", "--everything", "--hgvs", "--fasta", f"/refdir/{FASTA.name}",
           "--fork", "10"]
    t0 = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True)
    (out / f"{tag}_vep.log").write_text(p.stdout + p.stderr)
    t1 = time.time()
    child = ("import sys, vepyr; vepyr.annotate(vcf=sys.argv[1], cache_dir=sys.argv[3], everything=True, hgvs=True, "
             "reference_fasta=sys.argv[4], workers=10, output_vcf=sys.argv[2], show_progress=False)")
    q = subprocess.run([PYN, "-c", child, str(inp), str(vepyr_vcf), str(VEPYR_CACHE), str(FASTA)],
                       capture_output=True, text=True)
    (out / f"{tag}_vepyr.log").write_text(q.stdout + q.stderr)
    res = {"vep_returncode": p.returncode, "vepyr_returncode": q.returncode, "vep_seconds": t1 - t0,
           "vepyr_seconds": time.time() - t1}
    if p.returncode or q.returncode:
        return res
    res["vep_body_sha256"], res["vepyr_body_sha256"] = body(vep_vcf), body(vepyr_vcf)
    res["identical"] = res["vep_body_sha256"] == res["vepyr_body_sha256"]
    sys.path.insert(0, str(E2E_DIR))
    from comparison import compare  # noqa: E402
    with contextlib.redirect_stdout(io.StringIO()):
        res["e2e"] = compare.compare_vcfs(str(vepyr_vcf), str(vep_vcf), tag,
                                          mismatch_ledger_path=str(out / f"{tag}_e2e_ledger.jsonl"))
    for f in (vep_vcf, vepyr_vcf):
        subprocess.run(["bgzip", "-f", str(f)], check=True)
    return res


def main() -> int:
    out = Path(sys.argv[1]).resolve()
    out.mkdir(parents=True, exist_ok=True)
    which = sys.argv[2:] or ["panel", "clinvar"]
    res = {}
    for tag in which:
        inp, counts = {"clinvar": clinvar_input, "panel": panel_input}[tag](out)
        res[tag] = {"input": counts, **annotate_both(inp, out, tag)}
        print(tag, json.dumps({k: v for k, v in res[tag].items() if k != "e2e"}, default=str), flush=True)
        (out / "summary.json").write_text(json.dumps(res, indent=1, default=str) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
