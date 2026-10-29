#!/usr/bin/env python3
"""Behaviour of VEPYR 0.7.0 on unusual contigs (thesis item B19, 2026-10-29).

Annotates two tiny VCFs with the configuration of the timing series (everything, hgvs,
reference FASTA, workers=1):
  * edge.vcf: one SNV each on chr22, chrX, chrM, chrUn_KI270742v1 and the made-up
    contig chrFAKE;
  * mt.vcf: the chrM SNV (m.3243A>G) again, with the contig spelled MT.
For every input record the summary lists whether it reached the output, the number of
CSQ entries and the set of consequence terms, so silently dropped records and
contig-name dependent annotations are visible. No VEP run is involved.

Usage: edge_contigs.py OUT_DIR
Environment: PYN (.venv-native python with vepyr 0.7.0), CACHE, FASTA.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

D = Path("/Users/lukaszjezapkowicz/workspace/data_vepyr")
PYN = os.environ.get("PYN", str(D / ".venv-native/bin/python3"))
CACHE = os.environ.get("CACHE", str(D / "cache/116_GRCh38_merged"))
FASTA = os.environ.get("FASTA", str(D / "input/Homo_sapiens.GRCh38.dna.primary_assembly.fa"))

HEAD = "##fileformat=VCFv4.2\n"
COLS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
INPUTS = {
    "edge": [("chr22", 17951272, "C", "G"), ("chrX", 153866000, "A", "G"), ("chrM", 3243, "A", "G"),
             ("chrUn_KI270742v1", 10000, "A", "G"), ("chrFAKE", 100, "A", "G")],
    "mt": [("MT", 3243, "A", "G")],
}

CHILD = r'''
import sys, vepyr
vepyr.annotate(vcf=sys.argv[1], cache_dir=sys.argv[3], everything=True, hgvs=True, reference_fasta=sys.argv[4],
               workers=1, output_vcf=sys.argv[2], show_progress=False)
'''


def write_input(path: Path, records) -> None:
    contigs = "".join(f"##contig=<ID={c}>\n" for c in dict.fromkeys(r[0] for r in records))
    body = "".join(f"{c}\t{p}\t.\t{r}\t{a}\t.\tPASS\t.\n" for c, p, r, a in records)
    path.write_text(HEAD + contigs + COLS + body)


def summarize(path: Path) -> dict:
    out = {}
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            continue
        c = line.split("\t")
        csq = next((x[4:] for x in c[7].split(";") if x.startswith("CSQ=")), "")
        entries = [e.split("|") for e in csq.split(",")] if csq else []
        out[f"{c[0]}:{c[1]}:{c[3]}:{c[4]}"] = {
            "csq_entries": len(entries),
            "consequences": sorted({t for e in entries for t in e[1].split("&")}),
            "existing_variation": sorted({x for e in entries for x in e[17].split("&") if x}) if entries else [],
        }
    return out


def main() -> int:
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    summary = {}
    for name, records in INPUTS.items():
        vin, vout = out / f"{name}.vcf", out / f"{name}.annotated.vcf"
        write_input(vin, records)
        p = subprocess.run([PYN, "-c", CHILD, str(vin), str(vout), CACHE, FASTA], capture_output=True, text=True)
        (out / f"{name}.stderr.txt").write_text(p.stderr)
        got = summarize(vout) if p.returncode == 0 else {}
        summary[name] = {
            "returncode": p.returncode,
            "stderr_lines": len(p.stderr.splitlines()),
            "records_in": len(records),
            "records_out": len(got),
            "records": {f"{c}:{pos}:{r}:{a}": got.get(f"{c}:{pos}:{r}:{a}", "missing from output")
                        for c, pos, r, a in records},
        }
    (out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
