#!/usr/bin/env python3
"""Edge-case variant representations (thesis item B8, 2026-10-29): VEPYR side.

Writes a small chr22 VCF with representation edge cases in a coding exon (REF taken
from the reference FASTA unless a case is deliberately wrong), annotates it with VEPYR
0.7.0 (timing-series configuration, workers=1 on the plain VCF and workers=2 on the
bgzip+tabix copy) and lists per input record: present in the output, CSQ entries,
consequence terms and HGVSc/HGVSp values. The VEP side and the comparison are run by
post_e1_vep_checks.py b8 on the same input.

Usage: edge_variants.py OUT_DIR
Environment: PYN, CACHE, FASTA (as in edge_contigs.py).
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


def ref(pos: int, n: int = 1) -> str:
    out = subprocess.run(["samtools", "faidx", FASTA, f"22:{pos}-{pos + n - 1}"], capture_output=True,
                         text=True, check=True).stdout
    return "".join(out.splitlines()[1:]).upper()


def other(base: str, k: int = 0) -> str:
    return [b for b in "ACGT" if b != base][k]


def cases() -> list[tuple[str, int, str, str, str, str]]:
    """(case, pos, REF, ALT, ID, INFO) on chr22, sorted by position."""
    c = []
    p = 16591593                                     # missense SNV site in HG002
    r = ref(p)
    c.append(("multiallelic_snv", p, r, f"{other(r, 0)},{other(r, 1)}", ".", "."))
    p = 16591600
    r = ref(p)
    c.append(("spanning_deletion_star", p, r, f"{other(r)},*", ".", "."))
    p = 16591605
    r = ref(p)
    c.append(("gvcf_non_ref", p, r, f"{other(r)},<NON_REF>", ".", "."))
    p = 16591610
    r = ref(p, 2)
    c.append(("mnv_2bp", p, r, other(r[0]) + other(r[1]), ".", "."))
    p = 16591620
    r = ref(p, 3)
    c.append(("complex_delins", p, r, other(r[0], 1), ".", "."))
    p = 16591630
    r = ref(p)
    c.append(("ref_mismatch", p, other(r), r, ".", "."))
    p = 16591640
    r = ref(p)
    c.append(("no_alt_dot", p, r, ".", ".", "."))
    p = 16591650
    r = ref(p)
    c.append(("duplicate_record", p, r, other(r), ".", "."))
    c.append(("duplicate_record", p, r, other(r), ".", "."))
    # Homopolymer GGGGG at 19029751-19029755 (C at 19029750).
    s = ref(19029750, 6)
    assert s == "CGGGGG", s
    c.append(("del_left_aligned", 19029750, "CG", "C", ".", "."))
    c.append(("ins_left_aligned", 19029750, "C", "CG", ".", "."))
    c.append(("del_right_shifted", 19029754, "GG", "G", ".", "."))
    c.append(("ins_right_shifted", 19029755, "G", "GG", ".", "."))
    p = 19029760
    r = ref(p, 60)
    c.append(("long_deletion_60bp", p, r, r[0], ".", "."))
    p = 19029830
    r = ref(p)
    c.append(("symbolic_del", p, r, "<DEL>", ".", f"SVTYPE=DEL;END={p + 500};SVLEN=-500"))
    p = 19029840
    r = ref(p)
    c.append(("symbolic_dup", p, r, "<DUP>", ".", f"SVTYPE=DUP;END={p + 300};SVLEN=300"))
    p = 19029850
    r = ref(p)
    c.append(("symbolic_ins", p, r, "<INS>", ".", "SVTYPE=INS;SVLEN=120"))
    return c


HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr22,length=50818468>
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=NON_REF,Description="Any other allele">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="SV length">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""

CHILD = r'''
import sys, vepyr
vepyr.annotate(vcf=sys.argv[1], cache_dir=sys.argv[3], everything=True, hgvs=True, reference_fasta=sys.argv[4],
               workers=int(sys.argv[5]), output_vcf=sys.argv[2], show_progress=False)
'''


def summarize(path: Path) -> dict:
    fields, out = None, {}
    for line in path.read_text().splitlines():
        if line.startswith("##INFO=<ID=CSQ"):
            fields = line.split("Format: ")[1].rstrip('">').split("|")
        if line.startswith("#"):
            continue
        c = line.split("\t")
        csq = next((x[4:] for x in c[7].split(";") if x.startswith("CSQ=")), "")
        ents = [dict(zip(fields, e.split("|"))) for e in csq.split(",")] if csq else []
        key = f"{c[0]}:{c[1]}:{c[3]}:{c[4]}"
        out.setdefault(key, []).append({
            "csq_entries": len(ents),
            "alleles": sorted({e["Allele"] for e in ents}),
            "consequences": sorted({t for e in ents for t in e["Consequence"].split("&")}),
            "hgvsc": sorted({e["HGVSc"] for e in ents if e.get("HGVSc")})[:6],
            "hgvsp": sorted({e["HGVSp"] for e in ents if e.get("HGVSp")})[:6],
        })
    return out


def main() -> int:
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    cs = cases()
    vcf = out / "edge_variants.vcf"
    vcf.write_text(HEADER + "".join(f"chr22\t{p}\t{i}\t{r}\t{a}\t.\tPASS\t{info}\n" for _, p, r, a, i, info in cs))
    (out / "cases.json").write_text(json.dumps([{"case": n, "pos": p, "ref": r, "alt": a, "info": info}
                                                for n, p, r, a, _, info in cs], indent=1) + "\n")
    gz = out / "edge_variants.vcf.gz"
    subprocess.run(["bgzip", "-f", "-k", str(vcf)], check=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", str(gz)], check=True)
    summary = {}
    for workers, inp in ((1, vcf), (2, gz)):
        o = out / f"vepyr_w{workers}.vcf"
        p = subprocess.run([PYN, "-c", CHILD, str(inp), str(o), CACHE, FASTA, str(workers)],
                           capture_output=True, text=True)
        (out / f"vepyr_w{workers}.stderr.txt").write_text(p.stderr)
        got = summarize(o) if p.returncode == 0 else {}
        summary[f"workers{workers}"] = {
            "returncode": p.returncode, "stderr_tail": p.stderr[-600:],
            "records_in": len(cs), "records_out": sum(len(v) for v in got.values()),
            "cases": [{"case": n, "key": f"chr22:{pos}:{r}:{a}",
                       "out": got.get(f"chr22:{pos}:{r}:{a}", "missing from output")} for n, pos, r, a, _, _ in cs],
        }
    (out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    for w, s in summary.items():
        print(w, "rc", s["returncode"], "in", s["records_in"], "out", s["records_out"])
        for x in s["cases"]:
            o = x["out"]
            print(f"  {x['case']:24} {x['key']:38}",
                  o if isinstance(o, str) else [(e["csq_entries"], e["consequences"][:3], e["hgvsc"][:1]) for e in o])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
