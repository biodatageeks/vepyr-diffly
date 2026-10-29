#!/usr/bin/env python3
"""Plugin cache of known provenance and its build cost (thesis items C1 and C2, 2026-10-30).

All plugin runs in the thesis used one plugin cache built on 2026-08-27, before the
manifests of vepyr-plugins v0.1.1/v0.2.0, with no code version or source digests recorded.

  c1  Rebuild the five plugin caches for chr1 and chr22 with the studied VEPYR (0.7.0 @
      479d8f4, vepyr.build_plugin_cache) from the manifests of vepyr-plugins v0.2.0 (local
      clone), MD5 of every source checked against the manifest (strict; "warn" for the
      re-compressed AlphaMissense file, as its manifest documents). The cache manifests then
      record version@commit and the verified digests. Checks: row digests of every
      plugin-chromosome file against the 2026-08-27 copy (compare_cache_content.digest), and
      chr22/chr1 annotation with all five plugins on the new cache, whose body SHA-256 must
      equal the VEP 116 reference (results/references-2026-10-25).
  c2  Build cost of the plugin-chromosome units (5 plugins x chr1, chr22): every unit in its
      own process and its own cache root (no shared manifest), sources not re-hashed
      (verify_source=False), first P=1 (sequential), then P in PARALLEL concurrent units,
      longest first by the P=1 times. Per unit: wall time, peak RSS, CPU; per pass: makespan.

Usage: plugin_cache_rebuild.py c1 OUT_DIR
       plugin_cache_rebuild.py c2 OUT_DIR [PARALLEL, default 5]
Run with the .venv-native interpreter (vepyr 0.7.0 @ 479d8f4, polars).
"""
from __future__ import annotations

import concurrent.futures as cf
import json
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_cache_content import digest  # noqa: E402

D = Path("/Users/lukaszjezapkowicz/workspace/data_vepyr")
PKG = Path(__file__).resolve().parents[1]
PI = D / "plugin_input"
REPO = D / "tools/vepyr-plugins"
VERSION = "v0.2.0"
CORE = D / "cache/116_GRCh38_merged"
OLD = D / "plugin_cache"
FASTA = D / "input/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
INPUT = {"chr22": D / "runs/m5_completion_20261024/inputs/HG002_normalized_chr22.vcf.gz",
         "chr1": D / "runs/m5_completion_20261024/inputs/HG002_normalized_chr1.vcf.gz"}
CHROMS = ["chr1", "chr22"]
SOURCES = {
    "clinvar": str(PI / "clinvar/clinvar.vcf.gz"),
    "alphamissense": str(PI / "alphamissense/AlphaMissense_hg38.bgz.tsv.gz"),
    "spliceai": str(PI / "spliceai/spliceai_scores.masked.snv.ensembl_mane.grch38.110.vcf.gz"),
    "dbnsfp": str(PI / "dbnsfp/dbNSFP5.3.1a_grch38.gz"),
    "cadd": {"snv": str(PI / "cadd/whole_genome_SNVs.tsv.gz"), "indel": str(PI / "cadd/gnomad.genomes.r4.0.indel.tsv.gz")},
}
BUILD = r'''
import json, resource, sys, time, vepyr
a = json.loads(sys.argv[1])
t0 = time.perf_counter()
res = vepyr.build_plugin_cache(a["plugin"], a["version"], source_path=a["source"], cache_dir=a["cache_dir"],
                               plugin_cache_root=a["root"], chroms=a["chroms"], plugins_repo=a["repo"],
                               overwrite=True, verify_source=a["verify"])
dt = time.perf_counter() - t0
u = resource.getrusage(resource.RUSAGE_SELF)
print(json.dumps({"seconds": dt, "max_rss_gib": u.ru_maxrss / 1024**3, "user_s": u.ru_utime, "sys_s": u.ru_stime,
                  "chroms": [list(x) for x in res]}))
'''
ANNOTATE = r'''
import hashlib, json, sys, time, vepyr
vcf, out, cache, fasta, plugins = sys.argv[1:6]
t0 = time.perf_counter()
vepyr.annotate(vcf=vcf, cache_dir=cache, everything=True, hgvs=True, reference_fasta=fasta,
               plugin_cache_root=plugins, workers=1, output_vcf=out, show_progress=False)
dt = time.perf_counter() - t0
h = hashlib.sha256(); n = 0
with open(out, "rb") as f:
    for line in f:
        if not line.startswith(b"#"):
            h.update(line); n += 1
print(json.dumps({"annotation_seconds": dt, "body_sha256": h.hexdigest(), "records": n}))
'''


def child(code: str, *args: str, log: Path | None = None) -> dict:
    p = subprocess.run([sys.executable, "-c", code, *args], capture_output=True, text=True)
    if log:
        log.write_text(p.stdout + p.stderr)
    if p.returncode:
        raise SystemExit(f"child failed:\n{p.stderr[-2000:]}")
    return json.loads(p.stdout.strip().splitlines()[-1])


def build(plugin: str, root: Path, chroms: list[str], verify, log: Path) -> dict:
    a = {"plugin": plugin, "version": VERSION, "source": SOURCES[plugin], "cache_dir": str(CORE), "root": str(root),
         "chroms": chroms, "repo": str(REPO), "verify": verify}
    t0 = time.perf_counter()
    r = child(BUILD, json.dumps(a), log=log)
    r["process_wall_s"] = time.perf_counter() - t0
    return r


def c1(out: Path) -> dict:
    root = out / "plugin_cache"
    res = {"version": VERSION, "repo_commit": subprocess.run(["git", "-C", str(REPO), "rev-parse", VERSION + "^{commit}"],
                                                            capture_output=True, text=True).stdout.strip(),
           "builds": {}, "digests": {}, "annotation": {}}
    for plugin in SOURCES:
        res["builds"][plugin] = build(plugin, root, CHROMS, "warn" if plugin == "alphamissense" else True,
                                      out / f"c1_build_{plugin}.log")
        print(plugin, json.dumps(res["builds"][plugin]), flush=True)
    for m in sorted(root.rglob("manifest.json")):
        (out / f"c1_manifest_{m.parent.name}.json").write_text(m.read_text())
    for plugin in SOURCES:
        for chrom in CHROMS:
            new = next(root.rglob(f"{plugin}/{chrom}.parquet"))
            old = OLD / "plugin" / plugin / f"{chrom}.parquet"
            a, b = digest(old), digest(new)
            res["digests"][f"{plugin}:{chrom}"] = {"old": a, "new": b, "identical_rows": a == b}
            print(plugin, chrom, "identical" if a == b else "DIFFERENT", a["rows"], b["rows"], flush=True)
    m = next(root.rglob("manifest.json"))  # <root>/plugin/<name>/manifest.json, as in data_vepyr/plugin_cache
    plugin_root = m.parents[2] if m.parents[1].name == "plugin" else m.parents[1]
    for chrom in ("chr22", "chr1"):
        ref = (PKG / f"results/references-2026-10-25/{chrom}_reference_body.sha256").read_text().split()[0]
        vcf = out / f"{chrom}_all_plugins_new_cache.vcf"
        r = child(ANNOTATE, str(INPUT[chrom]), str(vcf), str(CORE), str(FASTA), str(plugin_root))
        vcf.unlink(missing_ok=True)
        r.update(vep_reference_body_sha256=ref, equal_vep_reference=r["body_sha256"] == ref)
        res["annotation"][chrom] = r
        print(chrom, json.dumps(r), flush=True)
    return res


def c2(out: Path, parallel: int) -> dict:
    units = [(p, c) for p in SOURCES for c in CHROMS]
    res = {"units": [f"{p}:{c}" for p, c in units], "passes": {}}

    def one(tag: str, unit: tuple[str, str]) -> dict:
        p, c = unit
        r = build(p, out / tag / f"{p}_{c}", [c], False, out / tag / f"{p}_{c}.log")
        r.update(unit=f"{p}:{c}", finished=time.perf_counter() - t0)
        return r

    for tag, workers in (("seq", 1), (f"par{parallel}", parallel)):
        (out / tag).mkdir(parents=True, exist_ok=True)
        if tag != "seq":  # longest first, by the sequential pass
            dur = {r["unit"]: r["seconds"] for r in res["passes"]["seq"]["units"]}
            units.sort(key=lambda u: -dur[f"{u[0]}:{u[1]}"])
        t0 = time.perf_counter()
        with cf.ThreadPoolExecutor(workers) as ex:
            done = list(ex.map(lambda u: one(tag, u), units))
        res["passes"][tag] = {"workers": workers, "makespan_s": time.perf_counter() - t0,
                              "sum_unit_s": sum(r["seconds"] for r in done),
                              "max_unit_rss_gib": max(r["max_rss_gib"] for r in done), "units": done}
        print(tag, json.dumps({k: v for k, v in res["passes"][tag].items() if k != "units"}), flush=True)
    return res


def main() -> int:
    mode, out = sys.argv[1], Path(sys.argv[2]).resolve()
    out.mkdir(parents=True, exist_ok=True)
    import vepyr
    if mode == "c1":
        res = c1(out)
    elif mode == "c2":
        res = c2(out, int(sys.argv[3]) if len(sys.argv) > 3 else 5)
    else:
        raise SystemExit(__doc__)
    res["vepyr"] = vepyr.__version__
    (out / f"{mode}_summary.json").write_text(json.dumps(res, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
