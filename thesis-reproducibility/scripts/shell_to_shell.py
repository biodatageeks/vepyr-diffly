#!/usr/bin/env python3
"""Shell-to-shell timing of VEPYR (thesis item D18, 2026-10-30).

The VEPYR series time vepyr.annotate inside a running interpreter (perf_counter), while VEP
is timed as the whole `docker run` command. Here VEPYR is timed the same way as VEP, as
one command from the shell: the `vepyr annotate` CLI of the native venv (interpreter
start, import, annotation, output), with hyperfine (--warmup 1, RUNS runs), on the chr22
and chr1 baseline with workers 1 and 10. The input preparation that VEPYR needs for
workers > 1 and VEP does not (bgzip + tabix of a plain VCF) is timed separately.
One output of each setting is checked against the body SHA-256 of the native series
(chr1: e9703e5e...; chr22: the workers=1 output of this script is the reference for
workers=10).

Usage: shell_to_shell.py OUT_DIR [RUNS]   (needs hyperfine; Colima can stay stopped)
"""
from __future__ import annotations

import hashlib
import json
import shlex
import subprocess
import sys
from pathlib import Path

D = Path("/Users/lukaszjezapkowicz/workspace/data_vepyr")
VEPYR = D / ".venv-native/bin/vepyr"
CACHE = D / "cache/116_GRCh38_merged"
FASTA = D / "input/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
INPUT = {"chr22": D / "runs/m5_completion_20261024/inputs/HG002_normalized_chr22.vcf.gz",
         "chr1": D / "runs/m5_completion_20261024/inputs/HG002_normalized_chr1.vcf.gz"}
NATIVE_BODY = {"chr1": "e9703e5e71f9e85438575d0382a6ea797b668b09db799a9e7825e602c1f81ffa"}


def body(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for line in f:
            if not line.startswith(b"#"):
                h.update(line)
    return h.hexdigest()


def cli(inp: Path, out: Path, workers: int) -> str:
    return shlex.join([str(VEPYR), "annotate", "-i", str(inp), "-o", str(out), "--dir_cache", str(CACHE),
                       "--fasta", str(FASTA), "--everything", "--hgvsc", "--fork", str(workers), "--no_progress"])


def main() -> int:
    out = Path(sys.argv[1]).resolve()
    runs = sys.argv[2] if len(sys.argv) > 2 else "5"
    out.mkdir(parents=True, exist_ok=True)
    res = {"vepyr_version": subprocess.run([str(VEPYR), "--version"], capture_output=True, text=True).stdout.strip(),
           "hyperfine": subprocess.run(["hyperfine", "--version"], capture_output=True, text=True).stdout.strip()}
    for chrom, inp in INPUT.items():
        plain = out / f"{chrom}.vcf"
        subprocess.run(f"bgzip -dc {shlex.quote(str(inp))} > {shlex.quote(str(plain))}", shell=True, check=True)
        o = out / f"{chrom}_out.vcf"
        cmds = {f"{chrom}_w1": cli(inp, o, 1), f"{chrom}_w10": cli(inp, o, 10),
                f"{chrom}_bgzip_tabix": f"bgzip -f -k {shlex.quote(str(plain))} && tabix -f -p vcf {shlex.quote(str(plain))}.gz"}
        for name, c in cmds.items():
            js = out / f"{name}.json"
            subprocess.run(["hyperfine", "--warmup", "1", "--runs", runs, "--export-json", str(js), "-n", name, c],
                           check=True, stdout=subprocess.DEVNULL)
            r = json.loads(js.read_text())["results"][0]
            res[name] = {k: r[k] for k in ("mean", "stddev", "median", "min", "max", "user", "system")} | {"command": c}
            if name.endswith(("_w1", "_w10")):
                res[name]["body_sha256"] = body(o)
            print(name, json.dumps({k: res[name][k] for k in ("mean", "stddev", "min", "max")}), flush=True)
        ref = NATIVE_BODY.get(chrom, res[f"{chrom}_w1"]["body_sha256"])
        res[f"{chrom}_identical"] = {w: res[f"{chrom}_{w}"]["body_sha256"] == ref for w in ("w1", "w10")}
        for f in (o, plain, Path(f"{plain}.gz"), Path(f"{plain}.gz.tbi")):
            f.unlink(missing_ok=True)
    (out / "summary.json").write_text(json.dumps(res, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
