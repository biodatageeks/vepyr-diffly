#!/usr/bin/env python3
"""All five plugins on chr1-22 (thesis item B12, 2026-10-30).

The plugin validation in the thesis covers chr1 and chr22. Here every chromosome 1-22 of
HG002 is annotated with all five plugins by VEP 116 and by VEPYR 0.7.0 and the outputs are
compared per chromosome:
  slices   cut from the local plugin sources exactly as build_vep_plugin_reference.sh cuts
           them over HTTP (`tabix -h SRC REGION | bgzip`, tabix -p vcf / -s1 -b2 -e2), same
           names; input = header + chromosome of input/HG002_normalized.vcf.gz;
  VEP      the reference command of the builder (pinned plugin code incl. corrected CADD.pm)
           plus --fork 10 (B17/D9 check that --fork leaves the body unchanged);
  VEPYR    everything, hgvs, FASTA, workers=10, plugin cache of the plugin series
           (cache/plugin_cache_116), all five plugins in VEP's CSQ order.
Per chromosome: records, body SHA-256 of both, identical or not; if not, the E2E comparator
(comparison.compare.compare_vcfs) with its mismatch ledger and both outputs kept (bgzip).
Slices of a finished chromosome are deleted (CADD alone is ~83 GB in total).

Usage: plugins_genome.py OUT_DIR [CHROM ...]   (OUT_DIR under data_vepyr; Colima running)
"""
from __future__ import annotations

import contextlib
import hashlib
import io
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

D = Path("/Users/lukaszjezapkowicz/workspace/data_vepyr")
PI = D / "plugin_input"
PYN = D / ".venv-native/bin/python3"
E2E_DIR = D / "tools/vepyr/e2e-testing/scripts"
BUILDER = E2E_DIR / "build_vep_plugin_reference.sh"
PLUGINS = D / "output/116/plugins/plugin_code"
NORM = D / "input/HG002_normalized.vcf.gz"
FASTA = D / "input/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
CACHE = D / "cache/116_GRCh38_merged"
PLUGIN_ROOT = D / "cache/plugin_cache_116"
IMAGE = "ensemblorg/ensembl-vep:release_116.0"
SLICES = [  # source, region template, slice name template, kind (as in the builder)
    ("clinvar/clinvar.vcf.gz", "{n}", "clinvar_chr{n}.vcf.gz", "vcf"),
    ("spliceai/spliceai_scores.masked.snv.ensembl_mane.grch38.110.vcf.gz", "{n}", "spliceai_chr{n}.vcf.gz", "vcf"),
    ("alphamissense/AlphaMissense_hg38.bgz.tsv.gz", "chr{n}", "alphamissense_chr{n}.tsv.gz", "tsv"),
    ("dbnsfp/dbNSFP5.3.1a_grch38.gz", "{n}", "dbNSFP5.3.1a_grch38_chr{n}.gz", "tsv"),
    ("cadd/whole_genome_SNVs.tsv.gz", "{n}", "cadd_snv_chr{n}.tsv.gz", "tsv"),
    ("cadd/gnomad.genomes.r4.0.indel.tsv.gz", "{n}", "cadd_indel_chr{n}.tsv.gz", "tsv"),
]
VEPYR = r'''
import json, sys, time, vepyr
t0 = time.perf_counter()
vepyr.annotate(vcf=sys.argv[1], cache_dir=sys.argv[2], everything=True, hgvs=True, reference_fasta=sys.argv[3],
               plugin_cache_root=sys.argv[4], workers=10, output_vcf=sys.argv[5], show_progress=False)
print(json.dumps({"seconds": time.perf_counter() - t0}))
'''


def sh(cmd: str) -> None:
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")


def body(path: Path) -> tuple[str, int]:
    h, n = hashlib.sha256(), 0
    with open(path, "rb") as f:
        for line in f:
            if not line.startswith(b"#"):
                h.update(line); n += 1
    return h.hexdigest(), n


def dbnsfp_cols() -> str:
    txt = subprocess.run(["sed", "-n", '/^DBNSFP_COLS="/,/"$/p', str(BUILDER)], capture_output=True, text=True).stdout
    return txt.replace("\\\n", "").strip().removeprefix('DBNSFP_COLS="').rstrip('"')


def one(out: Path, n: str) -> dict:
    work = out / f"work_chr{n}"
    sl = work / "slices"
    sl.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    for src, region, name, kind in SLICES:
        dst = sl / name.format(n=n)
        if not (dst.exists() and Path(f"{dst}.tbi").exists()):
            sh(f"tabix -h {PI / src} {region.format(n=n)} | bgzip -c > {dst}.partial.gz && mv {dst}.partial.gz {dst}")
            sh(f"tabix -f -p vcf {dst}" if kind == "vcf" else f"tabix -f -s 1 -b 2 -e 2 {dst}")
    inp = work / f"HG002_norm_chr{n}.vcf.gz"
    sh(f"{{ tabix -H {NORM}; tabix {NORM} chr{n}; }} | bgzip -c > {inp} && tabix -f -p vcf {inp}")
    r = {"chrom": f"chr{n}", "slicing_s": time.time() - t0}
    vep_vcf, vepyr_vcf = work / f"chr{n}_vep_fork10.vcf", work / f"chr{n}_vepyr_w10.vcf"
    s = f"/data/{sl.relative_to(D)}"
    cmd = ["docker", "run", "--rm", "--user", "501:20", "-v", f"{D}:/data", "-v", f"{PLUGINS}:/plugins:ro", IMAGE,
           "vep", "--cache", "--cache_version", "116", "--dir_cache", "/data", "--offline", "--merged",
           "--everything", "--no_stats", "--force_overwrite", "--vcf", "--fasta", f"/data/{FASTA.relative_to(D)}",
           "--input_file", f"/data/{inp.relative_to(D)}", "--output_file", f"/data/{vep_vcf.relative_to(D)}",
           "--dir_plugins", "/plugins", "--fork", "10",
           "--custom", f"{s}/clinvar_chr{n}.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN,CLNVC,CLNVI",
           "--plugin", f"SpliceAI,snv={s}/spliceai_chr{n}.vcf.gz,indel={s}/spliceai_chr{n}.vcf.gz",
           "--plugin", f"CADD,snv={s}/cadd_snv_chr{n}.tsv.gz,indels={s}/cadd_indel_chr{n}.tsv.gz",
           "--plugin", f"AlphaMissense,file={s}/alphamissense_chr{n}.tsv.gz",
           "--plugin", f"dbNSFP,{s}/dbNSFP5.3.1a_grch38_chr{n}.gz,{dbnsfp_cols()}"]
    t1 = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True)
    r["vep_seconds"], r["vep_returncode"] = time.time() - t1, p.returncode
    (out / f"chr{n}_vep.log").write_text(" ".join(cmd) + "\n" + p.stdout + p.stderr)
    q = subprocess.run([str(PYN), "-c", VEPYR, str(inp), str(CACHE), str(FASTA), str(PLUGIN_ROOT), str(vepyr_vcf)],
                       capture_output=True, text=True)
    r["vepyr_returncode"] = q.returncode
    (out / f"chr{n}_vepyr.log").write_text(q.stdout + q.stderr)
    if p.returncode or q.returncode:
        return r
    r["vepyr_seconds"] = json.loads(q.stdout.strip().splitlines()[-1])["seconds"]
    (r["vep_body_sha256"], r["vep_records"]), (r["vepyr_body_sha256"], r["vepyr_records"]) = body(vep_vcf), body(vepyr_vcf)
    r["identical"] = r["vep_body_sha256"] == r["vepyr_body_sha256"]
    if not r["identical"]:
        sys.path.insert(0, str(E2E_DIR))
        from comparison import compare  # noqa: E402
        with contextlib.redirect_stdout(io.StringIO()):
            r["e2e"] = compare.compare_vcfs(str(vepyr_vcf), str(vep_vcf), f"chr{n}_plugins",
                                            mismatch_ledger_path=str(out / f"chr{n}_e2e_ledger.jsonl"))
        for f in (vep_vcf, vepyr_vcf):
            sh(f"bgzip -f {f} && mv {f}.gz {out}/")
    shutil.rmtree(work)
    return r


def main() -> int:
    out = Path(sys.argv[1]).resolve()
    if not str(out).startswith(f"{D}/"):
        raise SystemExit(f"OUT_DIR must be under {D} (the only directory mounted into the VEP container)")
    out.mkdir(parents=True, exist_ok=True)
    chroms = [c.removeprefix("chr") for c in sys.argv[2:]] or [str(i) for i in range(1, 23)]
    res_path = out / "summary.json"
    res = json.loads(res_path.read_text()) if res_path.exists() else {}
    for n in chroms:
        if res.get(f"chr{n}", {}).get("vepyr_body_sha256"):
            continue
        res[f"chr{n}"] = one(out, n)
        res_path.write_text(json.dumps(res, indent=1, default=str) + "\n")
        print(json.dumps({k: v for k, v in res[f"chr{n}"].items() if k != "e2e"}, default=str), flush=True)
    done = [v for v in res.values() if "identical" in v]
    print(json.dumps({"chromosomes": len(done), "identical": sum(v["identical"] for v in done),
                      "records": sum(v["vep_records"] for v in done)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
