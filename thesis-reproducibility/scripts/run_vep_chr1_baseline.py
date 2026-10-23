#!/usr/bin/env python3
"""Run Ensembl VEP 116 baseline on chr1 (HG002, 323,430 variants) and measure warm execution time."""

import os
import shutil
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path("/Users/lukaszjezapkowicz/workspace/data_vepyr")
OUT_DIR = ROOT / "runs/vep_chr1_warmup"
OUT_DIR.mkdir(parents=True, exist_ok=True)
(OUT_DIR / ".metadata_never_index").touch()

OUT_VCF = OUT_DIR / "HG002_vep_chr1_warm.vcf"
TIME_LOG = OUT_DIR / "chr1_none.time.txt"
RUN_LOG = OUT_DIR / "vep_chr1.log"

CMD = [
    "/usr/bin/caffeinate", "-d", "-i", "-s", "-m",
    "/usr/bin/time", "-l", "-p",
    "docker", "run", "--rm",
    "--user", "501:20",
    "--env", "HOME=/tmp",
    "-v", f"{ROOT}/homo_sapiens_merged/116_GRCh38:/opt/vep/.vep/homo_sapiens_merged/116_GRCh38:ro",
    "-v", f"{ROOT}/runs/chapter5-m5-20260924T061203Z/correctness/input/HG002_normalized_chr1.vcf.gz:/input.vcf.gz:ro",
    "-v", f"{ROOT}/input:/refdir:ro",
    "-v", f"{OUT_DIR}:/output",
    "ensemblorg/ensembl-vep:release_116.0",
    "vep",
    "--dir", "/opt/vep/.vep",
    "--cache", "--merged", "--offline",
    "--assembly", "GRCh38",
    "--input_file", "/input.vcf.gz",
    "--output_file", f"/output/{OUT_VCF.name}",
    "--vcf", "--force_overwrite", "--no_stats",
    "--everything", "--hgvs",
    "--fasta", "/refdir/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
]

def main():
    print(f"[{datetime.now(timezone.utc).isoformat()}] Starting Ensembl VEP 116 baseline on chr1...", flush=True)
    if OUT_VCF.exists():
        OUT_VCF.unlink()
    
    t0 = time.perf_counter()
    with open(RUN_LOG, "w") as log_f, open(TIME_LOG, "w") as time_f:
        proc = subprocess.Popen(
            CMD,
            stdout=log_f,
            stderr=subprocess.PIPE,
            text=True,
        )
        _, stderr_text = proc.communicate()
        time_f.write(stderr_text)
    
    duration = time.perf_counter() - t0
    print(f"[{datetime.now(timezone.utc).isoformat()}] Finished VEP in {duration:.2f} s ({duration/60:.2f} min) with return code {proc.returncode}", flush=True)
    
    if OUT_VCF.exists():
        vcf_size_mb = OUT_VCF.stat().st_size / (1024 * 1024)
        print(f"Generated VCF size: {vcf_size_mb:.1f} MB", flush=True)
        # Unlink to preserve disk space
        OUT_VCF.unlink()
        print("Cleaned up temporary VCF.", flush=True)

if __name__ == "__main__":
    main()
