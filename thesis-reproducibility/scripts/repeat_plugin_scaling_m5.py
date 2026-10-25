#!/usr/bin/env python3
import json, os, shutil, subprocess, time
from pathlib import Path

ROOT = Path(__file__).resolve().parent / "plugin_scaling_m5_5repeats_20261013"
DATA = Path("/Users/lukaszjezapkowicz/workspace/data_vepyr")
RUNNER = Path(__file__).resolve().parent / "full_correctness_m5_20261012/scripts/run_comparison.py"
PYTHON = DATA / ".venv-native/bin/python"
VCF_ROOT = DATA / "runs/chapter5-m5-20260924T061203Z/correctness/input"
OUT = ROOT / "runs.jsonl"
ROOT.mkdir(parents=True, exist_ok=True)

def main():
    with OUT.open("a") as log:
        for chrom in os.environ.get("PLUGIN_REPEAT_CHROMS", "1,22").split(","):
            for workers in (1, 2, 4, 8, 10):
                for repeat in range(1, 6):
                    work = ROOT / f"chr{chrom}/workers{workers}/repeat{repeat}"
                    work.mkdir(parents=True, exist_ok=True)
                    cmd = [str(PYTHON), str(RUNNER), "--release", "116", "--profile", "merged_plugins",
                           "--chroms", chrom, "--workers", str(workers), "--skip-compare",
                           "--vcf", str(VCF_ROOT / f"HG002_normalized_chr{chrom}.vcf.gz"),
                           "--cache-dir", str(DATA / "cache/116_GRCh38_merged"),
                           "--plugin-cache", str(DATA / "plugin_cache"),
                           "--fasta", str(DATA / "input/Homo_sapiens.GRCh38.dna.primary_assembly.fa")]
                    env = os.environ.copy(); env["DATA_VEPYR_DIR"] = str(DATA); env["VEPYR_E2E_ROOT"] = str(work)
                    (work / "command.txt").write_text(" ".join(cmd) + "\n")
                    started = time.perf_counter()
                    p = subprocess.run(cmd, env=env, stdout=(work/"stdout.log").open("w"), stderr=(work/"stderr.log").open("w"))
                    elapsed = time.perf_counter() - started
                    reports = list((work / "reports").glob("*.json"))
                    report = json.loads(reports[0].read_text()) if reports else {}
                    row = {"chrom": chrom, "workers": workers, "repeat": repeat,
                           "returncode": p.returncode, "elapsed_wall_s": round(elapsed, 3),
                           "annotation_s": report.get("annotation", {}).get("time_s"),
                           "input_variants": report.get("input_variants"),
                           "output_variants": report.get("annotation", {}).get("output_variants")}
                    log.write(json.dumps(row) + "\n"); log.flush()
                    print(json.dumps(row), flush=True)
                    # Keep logs/reports, remove generated bulky result payloads.
                    for child in (work / "results", work / "output"):
                        if child.exists(): shutil.rmtree(child)

if __name__ == "__main__": main()
