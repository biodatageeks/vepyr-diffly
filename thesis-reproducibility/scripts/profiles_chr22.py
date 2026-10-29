#!/usr/bin/env python3
"""Ensembl and RefSeq cache profiles on chr22 (thesis item B13, 2026-10-30).

The validation in the thesis uses the merged profile only. For the ensembl and refseq
profiles of release 116 this script, on the HG002 chr22 slice:
  1. builds the VEPYR 0.7.0 cache for chr22 only (vepyr.build_cache(..., chroms=["chr22"]))
     from the raw Ensembl caches: homo_sapiens/116_GRCh38 (present) and
     homo_sapiens_refseq/116_GRCh38 (info/synonym files and chromosome 22 extracted from
     homo_sapiens_refseq_vep_116_GRCh38.tar.gz, BSD `sum` checked against Ensembl CHECKSUMS);
  2. annotates with VEP 116 (image and options of the timing series: --cache --offline
     --everything --hgvs --fasta --vcf --no_stats; --refseq for refseq) and with VEPYR
     (workers=1, everything, hgvs, reference FASTA);
  3. compares the SHA-256 of the record bodies and, if they differ, runs vepyr-diffly
     compare-existing (preset ensembl_everything; DIFFLY_SRC selects the checkout) to classify
     the differences.
One run of each tool per profile (times are recorded, not a timing series).
Run with the .venv-native interpreter (vepyr 0.7.0 @ 479d8f4); Colima running.

Usage: profiles_chr22.py OUT_DIR [ensembl|refseq ...]
"""
import hashlib, json, os, subprocess, sys, time
from pathlib import Path

D = Path("/Users/lukaszjezapkowicz/workspace/data_vepyr")
FASTA = D / "input/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
CHR22 = D / "runs/m5_completion_20261024/inputs/HG002_normalized_chr22.vcf.gz"
IMAGE = "ensemblorg/ensembl-vep:release_116.0"
TARBALL = D / "raw/vep_cache/homo_sapiens_refseq_vep_116_GRCh38.tar.gz"
CHECKSUMS = D / "raw/vep_cache/CHECKSUMS_release116"
# vepyr-diffly with the allele fix (branch fix/csq-allele-common-prefix, 1d657ce) via DIFFLY_SRC
DIFFLY_SRC = Path(os.environ.get("DIFFLY_SRC", str(D / "tools/vepyr-diffly/src")))
RAW = {"ensembl": D / "homo_sapiens/116_GRCh38", "refseq": D / "homo_sapiens_refseq/116_GRCh38"}
SPECIES_DIR = {"ensembl": "homo_sapiens", "refseq": "homo_sapiens_refseq"}

ANNOTATE = r'''
import json, sys, time, vepyr
vcf, out, cache, fasta = sys.argv[1:5]
t0 = time.perf_counter()
vepyr.annotate(vcf=vcf, cache_dir=cache, everything=True, hgvs=True, reference_fasta=fasta,
               workers=1, output_vcf=out, show_progress=False)
print(json.dumps({"annotation_seconds": time.perf_counter() - t0}))
'''


def body(path: Path) -> tuple[str, int]:
    h, n = hashlib.sha256(), 0
    with open(path, "rb") as f:
        for line in f:
            if not line.startswith(b"#"):
                h.update(line); n += 1
    return h.hexdigest(), n


def ensure_refseq_raw(res: dict) -> None:
    raw = RAW["refseq"]
    if (raw / "info.txt").exists() and (raw / "22").is_dir():
        return
    want = next(l.split()[:2] for l in CHECKSUMS.read_text().splitlines() if l.endswith(TARBALL.name))
    got = subprocess.run(["sum", str(TARBALL)], capture_output=True, text=True, check=True).stdout.split()[:2]
    res["refseq_tarball_sum"] = {"expected": want, "got": got, "bytes": TARBALL.stat().st_size}
    if got != want:
        raise SystemExit(f"checksum mismatch for {TARBALL}: {got} != {want}")
    t0 = time.time()
    prefix = "homo_sapiens_refseq/116_GRCh38"
    subprocess.run(["tar", "-xzf", str(TARBALL), "-C", str(D), "--include", f"{prefix}/*.txt",
                    "--include", f"{prefix}/22"], check=True)
    res["refseq_extract_seconds"] = time.time() - t0
    if not ((raw / "info.txt").exists() and (raw / "22").is_dir()):
        raise SystemExit(f"extraction did not produce {raw}/info.txt and {raw}/22")


def vep(profile: str, out_dir: Path, name: str) -> float:
    cmd = ["docker", "run", "--rm", "--user", "501:20", "--env", "HOME=/tmp",
           "-v", f"{RAW[profile]}:/opt/vep/.vep/{SPECIES_DIR[profile]}/116_GRCh38:ro",
           "-v", f"{CHR22}:/input/{CHR22.name}:ro", "-v", f"{FASTA.parent}:/refdir:ro", "-v", f"{out_dir}:/output",
           IMAGE, "vep", "--dir", "/opt/vep/.vep", "--cache", "--offline", "--assembly", "GRCh38",
           "--input_file", f"/input/{CHR22.name}", "--output_file", f"/output/{name}",
           "--vcf", "--force_overwrite", "--no_stats", "--everything", "--hgvs", "--fasta", f"/refdir/{FASTA.name}"]
    if profile == "refseq":
        cmd.append("--refseq")
    (out_dir / f"{name}.command.txt").write_text(" ".join(cmd) + "\n")
    t0 = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True)
    (out_dir / f"{name}.log").write_text(p.stdout + p.stderr)
    if p.returncode:
        raise SystemExit(f"VEP {profile} failed:\n{p.stderr[-2000:]}")
    return time.time() - t0


def main() -> int:
    import vepyr
    out = Path(sys.argv[1]).resolve()
    out.mkdir(parents=True, exist_ok=True)
    profiles = sys.argv[2:] or ["ensembl", "refseq"]
    res = {"vepyr": vepyr.__version__, "input": str(CHR22),
           "image_id": subprocess.run(["docker", "image", "inspect", "--format", "{{.Id}} {{.Architecture}}", IMAGE],
                                      capture_output=True, text=True).stdout.strip()}
    if "refseq" in profiles:
        ensure_refseq_raw(res)
    for prof in profiles:
        r = res[prof] = {}
        t0 = time.perf_counter()
        built = vepyr.build_cache(116, str(out / "cache"), cache_type=prof, chroms=["chr22"],
                                  local_cache=str(RAW[prof]), overwrite=True, show_progress=False)
        r["vepyr_build_seconds"] = time.perf_counter() - t0
        r["vepyr_build_files"] = len(built)
        cache = out / "cache" / f"116_GRCh38_{prof}"

        vep_vcf, vepyr_vcf = out / f"chr22_{prof}_vep.vcf", out / f"chr22_{prof}_vepyr.vcf"
        r["vep_seconds"] = vep(prof, out, vep_vcf.name)
        p = subprocess.run([sys.executable, "-c", ANNOTATE, str(CHR22), str(vepyr_vcf), str(cache), str(FASTA)],
                           capture_output=True, text=True)
        if p.returncode:
            raise SystemExit(f"VEPYR {prof} failed:\n{p.stderr[-2000:]}")
        r["vepyr_seconds"] = json.loads(p.stdout.strip().splitlines()[-1])["annotation_seconds"]
        (r["vep_body_sha256"], r["vep_records"]), (r["vepyr_body_sha256"], r["vepyr_records"]) = body(vep_vcf), body(vepyr_vcf)
        r["identical"] = r["vep_body_sha256"] == r["vepyr_body_sha256"]
        if not r["identical"]:
            d = out / f"diffly_{prof}"
            cmd = [sys.executable, "-m", "vepyr_diffly.cli", "compare-existing", "--preset", "ensembl_everything",
                   "--left-vcf", str(vep_vcf), "--right-vcf", str(vepyr_vcf), "--output-dir", str(d)]
            p = subprocess.run(cmd, capture_output=True, text=True, env={**os.environ, "PYTHONPATH": str(DIFFLY_SRC)})
            (out / f"diffly_{prof}.log").write_text(p.stdout + p.stderr)
            r["diffly_returncode"] = p.returncode
            if (d / "summary.json").exists():
                r["diffly_summary"] = json.loads((d / "summary.json").read_text())
        for f in (vep_vcf, vepyr_vcf):
            subprocess.run(["bgzip", "-f", str(f)], check=True)
        (out / "summary.json").write_text(json.dumps(res, indent=1) + "\n")
        print(prof, json.dumps({k: v for k, v in r.items() if k != "diffly_summary"}), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
