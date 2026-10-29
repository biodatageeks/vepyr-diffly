#!/usr/bin/env python3
"""Official vs corrected VEP CADD plugin (thesis item B5).

Runs VEP 116 (merged cache, --everything --hgvs --fasta, as in the plugin series) on
HG002 chr22 or chr1 with only the CADD plugin, using the official CADD.pm of
Ensembl/VEP_plugins release/116, and compares CADD_RAW/CADD_PHRED per CSQ entry
(variant, Allele, Feature) with the five-plugin reference built with the corrected
plugin (mwiewior/VEP_plugins 7a1f6450fd12). Both runs read the same CADD v1.7 files.

Reported per chromosome: CSQ entries and variants with a CADD value only in the
corrected reference (missed by the official plugin), only in the official run, equal,
or different; split by SNV/indel; first examples.

Usage: cadd_official_vs_fixed.py OUT_DIR CHROM   (CHROM = chr22 or chr1)
"""
from __future__ import annotations

import gzip
import hashlib
import json
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

D = Path("/Users/lukaszjezapkowicz/workspace/data_vepyr")
OFFICIAL_URL = "https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/116/CADD.pm"
INPUT = D / "runs/m5_completion_20261024/inputs/HG002_normalized_{c}.vcf.gz"
REFERENCE = D / "output/116/plugins/HG002_{c}_5plugins_vep116_caddfix.vcf.gz"
IMAGE = "ensemblorg/ensembl-vep:release_116.0"


def entries(line: str, fields: list[str]):
    c = line.rstrip("\n").split("\t")
    csq = next((x[4:] for x in c[7].split(";") if x.startswith("CSQ=")), "")
    key = f"{c[0]}:{c[1]}:{c[3]}:{c[4]}"
    snv = len(c[3]) == 1 and len(c[4]) == 1
    out = {}
    for e in csq.split(",") if csq else []:
        d = dict(zip(fields, e.split("|")))
        out[(d["Allele"], d["Feature"])] = (d.get("CADD_RAW", ""), d.get("CADD_PHRED", ""))
    return key, snv, out


def fields_of(header_line: str) -> list[str]:
    return header_line.split("Format: ")[1].rstrip('">\n').split("|")


def read(path: Path):
    op = gzip.open if path.suffix == ".gz" else open
    fields = None
    with op(path, "rt") as f:
        for line in f:
            if line.startswith("##INFO=<ID=CSQ"):
                fields = fields_of(line)
            if line.startswith("#"):
                continue
            yield entries(line, fields)


def main() -> int:
    out, chrom = Path(sys.argv[1]).resolve(), sys.argv[2]
    out.mkdir(parents=True, exist_ok=True)
    plugins = out / "plugins"
    plugins.mkdir(exist_ok=True)
    pm = plugins / "CADD.pm"
    if not pm.exists():
        pm.write_bytes(urllib.request.urlopen(OFFICIAL_URL).read())
    sha = hashlib.sha256(pm.read_bytes()).hexdigest()
    inp = Path(str(INPUT).format(c=chrom))
    vcf = out / f"{chrom}_cadd_official.vcf"
    cmd = ["docker", "run", "--rm", "--user", "501:20", "--env", "HOME=/tmp",
           "-v", f"{D}/homo_sapiens_merged/116_GRCh38:/opt/vep/.vep/homo_sapiens_merged/116_GRCh38:ro",
           "-v", f"{inp}:/input.vcf.gz:ro", "-v", f"{D}/input:/refdir:ro", "-v", f"{out}:/output",
           "-v", f"{D}/plugin_input:/pluginwork:ro", "-v", f"{plugins}:/plugins:ro",
           IMAGE, "vep", "--dir", "/opt/vep/.vep", "--cache", "--merged", "--offline", "--assembly", "GRCh38",
           "--input_file", "/input.vcf.gz", "--output_file", f"/output/{vcf.name}", "--vcf", "--force_overwrite",
           "--no_stats", "--everything", "--hgvs", "--fasta", "/refdir/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
           "--dir_plugins", "/plugins",
           "--plugin", "CADD,snv=/pluginwork/cadd/whole_genome_SNVs.tsv.gz,indels=/pluginwork/cadd/gnomad.genomes.r4.0.indel.tsv.gz"]
    (out / f"{chrom}_command.txt").write_text(" ".join(cmd) + "\n")
    if not vcf.exists():
        t0 = time.time()
        p = subprocess.run(cmd, capture_output=True, text=True)
        (out / f"{chrom}_vep.log").write_text(f"seconds={time.time() - t0:.1f}\n" + p.stdout + p.stderr)
        if p.returncode != 0:
            print(p.stderr[-2000:])
            return 1
    cnt = {k: {"snv": 0, "indel": 0} for k in ("only_fixed", "only_official", "equal", "different")}
    var = {k: {"snv": 0, "indel": 0} for k in ("variants_any_missed", "variants_with_cadd_fixed")}
    examples, n = [], 0
    for (k1, snv, a), (k2, _, b) in zip(read(REFERENCE.parent / REFERENCE.name.format(c=chrom)), read(vcf)):
        assert k1 == k2, (k1, k2)
        n += 1
        t = "snv" if snv else "indel"
        missed = False
        has = False
        for ek in set(a) | set(b):
            fa, fb = a.get(ek, ("", "")), b.get(ek, ("", ""))
            if fa[1] or fa[0]:
                has = True
            if (fa[0] or fa[1]) and not (fb[0] or fb[1]):
                cnt["only_fixed"][t] += 1
                missed = True
                if len(examples) < 15:
                    examples.append({"variant": k1, "allele_feature": list(ek), "fixed": fa, "official": fb})
            elif (fb[0] or fb[1]) and not (fa[0] or fa[1]):
                cnt["only_official"][t] += 1
            elif fa == fb:
                cnt["equal"][t] += 1
            else:
                cnt["different"][t] += 1
        var["variants_any_missed"][t] += missed
        var["variants_with_cadd_fixed"][t] += has
    res = {"chrom": chrom, "records": n, "official_cadd_pm_sha256": sha, "entries": cnt, "variants": var,
           "examples": examples}
    (out / f"{chrom}_summary.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps({k: v for k, v in res.items() if k != "examples"}, indent=1))
    vcf.unlink()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
