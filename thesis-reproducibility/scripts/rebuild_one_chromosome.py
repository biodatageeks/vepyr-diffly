#!/usr/bin/env python3
"""Rebuild of one chromosome in a finished core cache (thesis item E9, 2026-10-30).

Checks the claim of chapter 6 that chromosome manifests allow rebuilding a chosen part of
the cache without removing the entries of the other chromosomes. On an APFS clone of the
core cache built by VEPYR 0.7.0 (S2, 2026-10-29; `cp -c`, the original is not modified),
chr22 is rebuilt with the public API, vepyr.build_cache(..., chroms=["chr22"],
overwrite=True). Then:
  * files: SHA-256 of every file of the clone against the original; files added or removed;
  * chrom_manifest.json of every entity: entries of the other contigs unchanged, chr22 entry;
  * chr22 shards: row count and order-independent row digests (as in compare_cache_content.py);
  * annotation of the HG002 chr22 slice with the rebuilt clone: with all five plugins the body
    SHA-256 must equal the VEP 116 reference (4d8ff253..., results/references-2026-10-25),
    without plugins it must equal the same annotation with the original cache.
Run with the .venv-native interpreter (vepyr 0.7.0 @ 479d8f4, polars).

Usage: rebuild_one_chromosome.py OUT_DIR
"""
import hashlib, json, os, shutil, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_cache_content import digest  # noqa: E402

D = Path("/Users/lukaszjezapkowicz/workspace/data_vepyr")
PKG = Path(__file__).resolve().parents[1]
ORIGINAL = D / "cache/116_GRCh38_merged"
VEP_CACHE = D / "homo_sapiens_merged/116_GRCh38"
PLUGIN_ROOT_ALL = D / "cache/plugin_cache_116"
FASTA = D / "input/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
INPUT = D / "runs/m5_completion_20261024/inputs/HG002_normalized_chr22.vcf.gz"
REFERENCE_BODY = (PKG / "results/references-2026-10-25/chr22_reference_body.sha256").read_text().split()[0]

ANNOTATE = r'''
import hashlib, json, sys, time, vepyr
vcf, out, cache, fasta, plugins = sys.argv[1:6]
t0 = time.perf_counter()
vepyr.annotate(vcf=vcf, cache_dir=cache, everything=True, hgvs=True, reference_fasta=fasta,
               plugin_cache_root=plugins or None, workers=1, output_vcf=out, show_progress=False)
dt = time.perf_counter() - t0
h = hashlib.sha256(); n = 0
with open(out, "rb") as f:
    for line in f:
        if not line.startswith(b"#"):
            h.update(line); n += 1
print(json.dumps({"annotation_seconds": dt, "body_sha256": h.hexdigest(), "records": n}))
'''


def sha(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def tree(root: Path) -> dict[str, str]:
    files = sorted(p for p in root.rglob("*") if p.is_file())
    with ThreadPoolExecutor(8) as ex:
        return dict(zip((str(p.relative_to(root)) for p in files), ex.map(sha, files)))


def annotate(out: Path, cache: Path, plugins: str) -> dict:
    p = subprocess.run([sys.executable, "-c", ANNOTATE, str(INPUT), str(out), str(cache), str(FASTA), plugins],
                       capture_output=True, text=True)
    if p.returncode:
        raise SystemExit(p.stderr[-2000:])
    out.unlink(missing_ok=True)
    return json.loads(p.stdout.strip().splitlines()[-1])


def main() -> int:
    import vepyr
    out = Path(sys.argv[1]).resolve()
    clone_root = out / "cache"
    clone = clone_root / ORIGINAL.name
    if clone.exists():
        raise SystemExit(f"{clone} exists; remove it first")
    clone_root.mkdir(parents=True)
    res = {"vepyr": vepyr.__version__, "original": str(ORIGINAL)}

    t0 = time.time()
    subprocess.run(["cp", "-c", "-R", str(ORIGINAL), str(clone)], check=True)
    res["clone_seconds"] = time.time() - t0
    t0 = time.time()
    before = tree(ORIGINAL)
    res["hash_original_seconds"] = time.time() - t0

    t0 = time.perf_counter()
    built = vepyr.build_cache(116, str(clone_root), cache_type="merged", chroms=["chr22"],
                              local_cache=str(VEP_CACHE), overwrite=True, show_progress=False)
    res["rebuild_seconds"] = time.perf_counter() - t0
    res["rebuilt_files"] = [[str(Path(p).relative_to(clone)), n] for p, n in built]

    after = tree(clone)
    changed = sorted(k for k in before.keys() & after.keys() if before[k] != after[k])
    res["files"] = {"original": len(before), "after": len(after), "added": sorted(after.keys() - before.keys()),
                    "removed": sorted(before.keys() - after.keys()), "changed": changed,
                    "changed_outside_chr22_and_manifests": [k for k in changed if not (
                        k.endswith("/chr22.parquet") or k.endswith("/chrom_manifest.json"))]}

    manifests = {}
    for m in sorted(ORIGINAL.glob("*/chrom_manifest.json")):
        a = json.loads(m.read_text())
        b = json.loads((clone / m.relative_to(ORIGINAL)).read_text())
        other = lambda xs: sorted((x for x in xs if x["chrom"] != "chr22"), key=lambda x: x["chrom"])  # noqa: E731
        c22 = lambda xs: [x for x in xs if x["chrom"] == "chr22"]  # noqa: E731
        manifests[m.parent.name] = {"entries": [len(a), len(b)], "others_unchanged": other(a) == other(b),
                                    "chr22": [c22(a), c22(b)]}
    res["manifests"] = manifests

    shards = {}
    for p in sorted(ORIGINAL.glob("*/chr22.parquet")):
        x, y = digest(p), digest(clone / p.relative_to(ORIGINAL))
        shards[p.parent.name] = {"original": x, "rebuilt": y, "identical_rows": x == y}
    res["chr22_shards"] = shards

    res["annotation"] = {
        "plugins_rebuilt": annotate(out / "plugins_rebuilt.vcf", clone, str(PLUGIN_ROOT_ALL)),
        "baseline_rebuilt": annotate(out / "baseline_rebuilt.vcf", clone, ""),
        "baseline_original": annotate(out / "baseline_original.vcf", ORIGINAL, ""),
        "vep_reference_body_sha256": REFERENCE_BODY,
    }
    a = res["annotation"]
    a["plugins_equal_vep_reference"] = a["plugins_rebuilt"]["body_sha256"] == REFERENCE_BODY
    a["baseline_equal_original_cache"] = a["baseline_rebuilt"]["body_sha256"] == a["baseline_original"]["body_sha256"]

    (out / "summary.json").write_text(json.dumps(res, indent=1) + "\n")
    shutil.rmtree(clone_root)
    print(json.dumps({k: v for k, v in res.items() if k not in ("rebuilt_files", "chr22_shards")}, indent=1)[:4000])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
