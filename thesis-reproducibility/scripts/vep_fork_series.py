#!/usr/bin/env python3
"""VEP 116 at its fastest measured setting, --fork 10, on the full input and with all five
plugins (thesis items D9 and D16, 2026-10-30).

Modes
  full     chr1-22 HG002 (input/HG002_normalized.vcf.gz) with the command of the full-input
           series (run_vep_full_repeats.sh) plus --fork 10; the body SHA-256 must equal the
           single-process run (5fee6224...).
  plugins  chr1 or chr22 with all five plugins: slices, input and pinned plugin code
           (corrected CADD.pm) of build_vep_plugin_reference.sh, the reference command plus
           --fork 10; the body SHA-256 must equal
           results/references-2026-10-25/chrN_reference_body.sha256. The slices are made once
           per chromosome by running the builder (rclone serves plugin_input), which also
           writes one more single-process reference (its time includes slicing). A work
           directory with slices can be reused through PLUGIN_WORK.
Per run: wall time of `docker run` (/usr/bin/time -p, as in the other VEP series); CPU time
of all processes in the container (cgroup cpu.stat after vep exits; the command runs under
`sh -c` for that); memory: cgroup memory.peak (includes page cache charged to the
container) and, sampled every INTERVAL seconds with `docker exec`, the sums of Rss and Pss
(/proc/PID/smaps_rollup) over the VEP parent and its fork children. Pss divides pages shared
between the children, so its sum does not count them several times.

Usage: vep_fork_series.py full OUT_DIR [REPEATS]
       vep_fork_series.py plugins OUT_DIR CHROM [REPEATS]
Env: INTERVAL (seconds between memory samples, default 5), PLUGIN_WORK, KEEP_FIRST (path: the
repeat-1 output is bgzipped there instead of deleted, e.g. as the chr1-22 VEP reference for B2).
"""
from __future__ import annotations

import hashlib
import json
import os
import shlex
import subprocess
import sys
import time
from pathlib import Path

D = Path("/Users/lukaszjezapkowicz/workspace/data_vepyr")
PKG = Path(__file__).resolve().parents[1]
IMAGE = "ensemblorg/ensembl-vep:release_116.0"
FORK = 10
INTERVAL = float(os.environ.get("INTERVAL", "5"))
BUILDER = D / "tools/vepyr/e2e-testing/scripts/build_vep_plugin_reference.sh"
PLUGINS = D / "output/116/plugins/plugin_code"
FULL_BODY = "5fee6224bf135d80075712141cc06194e0730fa253ff8b131a9d25d661647ee9"
FASTA = "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
WRAP = ('vep "$@"; rc=$?; echo CGROUP; cat /sys/fs/cgroup/cpu.stat; '
        'echo "memory_peak $(cat /sys/fs/cgroup/memory.peak 2>/dev/null)"; exit $rc')
SAMPLER = ('for d in /proc/[0-9]*; do c=$(cat $d/comm 2>/dev/null); case "$c" in perl|vep) '
           'echo "pid ${d#/proc/}"; grep -E "^(Rss|Pss):" $d/smaps_rollup 2>/dev/null;; esac; done')


def body(path: Path) -> tuple[str, int]:
    h, n = hashlib.sha256(), 0
    with open(path, "rb") as f:
        for line in f:
            if not line.startswith(b"#"):
                h.update(line); n += 1
    return h.hexdigest(), n


def sample(name: str) -> dict | None:
    p = subprocess.run(["docker", "exec", name, "sh", "-c", SAMPLER], capture_output=True, text=True, timeout=60)
    if p.returncode or "pid" not in p.stdout:
        return None
    s = {"t": time.time(), "procs": 0, "rss_kb": 0, "pss_kb": 0}
    for line in p.stdout.splitlines():
        k, v = line.split()[:2]
        if k == "pid":
            s["procs"] += 1
        elif k in ("Rss:", "Pss:"):
            s[k[:3].lower() + "_kb"] += int(v)
    return s


def run_vep(run: Path, name: str, mounts: list[str], args: list[str], output: Path, expected_body: str) -> dict:
    run.mkdir(parents=True, exist_ok=True)
    t = run / "time.txt"
    cmd = ["/usr/bin/time", "-p", "-o", str(t), "docker", "run", "--rm", "--name", name, "--user", "501:20",
           "--env", "HOME=/tmp", *mounts, IMAGE, "sh", "-c", WRAP, "vep", *args, "--fork", str(FORK)]
    (run / "command.txt").write_text(shlex.join(cmd) + "\n")
    subprocess.run(["docker", "rm", "-f", name], capture_output=True)
    samples = []
    with open(run / "stdout.txt", "w") as so, open(run / "stderr.txt", "w") as se:
        t0 = time.time()
        p = subprocess.Popen(cmd, stdout=so, stderr=se)
        time.sleep(min(INTERVAL, 3))
        while p.poll() is None:
            s = sample(name)
            if s:
                s["t"] -= t0
                samples.append(s)
            time.sleep(INTERVAL)
        rc = p.wait()
    (run / "memory_samples.jsonl").write_text("".join(json.dumps(s) + "\n" for s in samples))
    out = (run / "stdout.txt").read_text()
    stat = dict(l.split()[:2] for l in out.split("CGROUP\n", 1)[-1].splitlines() if len(l.split()) == 2)
    r = {"repeat_dir": run.name, "exit_code": rc,
         "wall_s": float(next(l.split()[1] for l in t.read_text().splitlines() if l.startswith("real"))),
         "cpu_s": {k: int(stat[f"{k}_usec"]) / 1e6 for k in ("usage", "user", "system") if f"{k}_usec" in stat},
         "cgroup_memory_peak_gib": int(stat["memory_peak"]) / 1024**3 if stat.get("memory_peak", "").isdigit() else None,
         "samples": len(samples),
         "peak_sum_rss_gib": max((s["rss_kb"] for s in samples), default=0) / 1024**2,
         "peak_sum_pss_gib": max((s["pss_kb"] for s in samples), default=0) / 1024**2,
         "max_procs": max((s["procs"] for s in samples), default=0)}
    if r["cpu_s"].get("usage"):
        r["avg_busy_cores"] = r["cpu_s"]["usage"] / r["wall_s"]
    if output.exists():
        r["body_sha256"], r["records"] = body(output)
        r["body_identical"] = r["body_sha256"] == expected_body
        keep = os.environ.get("KEEP_FIRST")
        if keep and run.name.endswith("_r1") and not Path(keep).exists():
            subprocess.run(f"bgzip -c {shlex.quote(str(output))} > {shlex.quote(keep)} && tabix -p vcf {shlex.quote(keep)}",
                           shell=True, check=True)
            r["kept_output"] = keep
        output.unlink()
    (run / "result.json").write_text(json.dumps(r, indent=1) + "\n")
    return r


def plugin_work(out: Path, chrom: str) -> tuple[Path, dict]:
    work = Path(os.environ.get("PLUGIN_WORK", out / f"work_chr{chrom}")).resolve()
    if not str(work).startswith(f"{D}/"):
        raise SystemExit(f"plugin work dir must be under {D}")
    slices = work / "slices"
    need = [slices / f"{p}_chr{chrom}.{e}" for p, e in (("clinvar", "vcf.gz"), ("spliceai", "vcf.gz"),
            ("cadd_snv", "tsv.gz"), ("cadd_indel", "tsv.gz"), ("alphamissense", "tsv.gz"))]
    need.append(slices / f"dbNSFP5.3.1a_grch38_chr{chrom}.gz")
    info = {}
    if not all(p.exists() for p in need):
        work.mkdir(parents=True, exist_ok=True)
        rcl = subprocess.Popen(["rclone", "serve", "http", str(D / "plugin_input"), "--addr", "127.0.0.1:18080",
                                "--read-only"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        try:
            for _ in range(30):
                if subprocess.run(["curl", "-sf", "--head", "http://127.0.0.1:18080/clinvar/clinvar.vcf.gz"],
                                  capture_output=True).returncode == 0:
                    break
                time.sleep(1)
            ref = work / "output" / f"HG002_chr{chrom}_5plugins_builder.vcf"
            env = {**os.environ, "DATA_VEPYR_DIR": str(D), "VEP_PLUGIN_DIR": str(PLUGINS),
                   "VEP_PLUGIN_SOURCE_URL": "http://127.0.0.1:18080",
                   "VEP_NORMALIZED_VCF": str(D / "input/HG002_normalized.vcf.gz"),
                   "VEP_OUTPUT_VCF": str(ref), "VEP_KEEP_PLAIN": "0"}
            t0 = time.time()
            with open(out / f"builder_chr{chrom}.log", "w") as log:
                rc = subprocess.run([str(BUILDER), chrom, str(work)], env=env, stdout=log, stderr=log).returncode
            info["builder_seconds_incl_slicing"] = time.time() - t0
        finally:
            rcl.terminate()
        if rc:
            raise SystemExit(f"builder failed for chr{chrom}, see builder_chr{chrom}.log")
        gz = Path(str(ref) + ".gz")
        info["builder_body_sha256"] = subprocess.run(f"bgzip -dc {shlex.quote(str(gz))} | grep -v '^#' | shasum -a 256",
                                                     shell=True, capture_output=True, text=True).stdout.split()[0]
    return work, info


def main() -> int:
    mode, out = sys.argv[1], Path(sys.argv[2]).resolve()
    out.mkdir(parents=True, exist_ok=True)
    env = {"utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "fork": FORK, "interval_s": INTERVAL,
           "image": subprocess.run(["docker", "image", "inspect", "--format", "{{.Id}} {{.Architecture}}", IMAGE],
                                   capture_output=True, text=True).stdout.strip(),
           "colima": [l for l in (Path.home() / ".colima/default/colima.yaml").read_text().splitlines()
                      if l.split(":")[0] in ("cpu", "memory", "vmType", "mountType")]}
    if mode == "full":
        reps = int(sys.argv[3]) if len(sys.argv) > 3 else 5
        tag, expected = "full", FULL_BODY
    elif mode == "plugins":
        chrom = sys.argv[3].removeprefix("chr")
        reps = int(sys.argv[4]) if len(sys.argv) > 4 else 5
        tag = f"plugins_chr{chrom}"
        expected = (PKG / f"results/references-2026-10-25/chr{chrom}_reference_body.sha256").read_text().split()[0]
        work, info = plugin_work(out, chrom)
        env.update(info, plugin_work=str(work), expected_body_sha256=expected,
                   builder_body_identical=info.get("builder_body_sha256", expected) == expected)
    else:
        raise SystemExit(__doc__)
    (out / f"{tag}_environment.json").write_text(json.dumps(env, indent=1) + "\n")

    for rep in range(1, reps + 1):
        run = out / f"{tag}_r{rep}"
        if (run / "result.json").exists():
            continue
        if mode == "full":
            name = f"/output/HG002_vep_fork{FORK}.vcf"
            mounts = ["-v", f"{D}/homo_sapiens_merged/116_GRCh38:/opt/vep/.vep/homo_sapiens_merged/116_GRCh38:ro",
                      "-v", f"{D}/input:/input:ro", "-v", f"{run}:/output"]
            args = ["--dir", "/opt/vep/.vep", "--cache", "--merged", "--offline", "--assembly", "GRCh38",
                    "--input_file", "/input/HG002_normalized.vcf.gz", "--output_file", name,
                    "--vcf", "--force_overwrite", "--no_stats", "--everything", "--hgvs", "--fasta", f"/input/{FASTA}"]
            output = run / Path(name).name
        else:
            if not str(run).startswith(f"{D}/"):
                raise SystemExit(f"OUT_DIR must be under {D} in plugins mode (only {D} is mounted)")
            s = f"/data/{work.relative_to(D)}/slices"
            output = run / f"HG002_chr{chrom}_5plugins_fork{FORK}.vcf"
            cols = subprocess.run(["sed", "-n", '/^DBNSFP_COLS="/,/"$/p', str(BUILDER)], capture_output=True,
                                  text=True).stdout.replace("\\\n", "").strip().removeprefix('DBNSFP_COLS="').rstrip('"')
            mounts = ["-v", f"{D}:/data", "-v", f"{PLUGINS}:/plugins:ro"]
            args = ["--cache", "--cache_version", "116", "--dir_cache", "/data", "--offline", "--merged",
                    "--everything", "--no_stats", "--force_overwrite", "--vcf", "--fasta", f"/data/input/{FASTA}",
                    "--input_file", f"/data/{work.relative_to(D)}/input/HG002_norm_chr{chrom}.vcf.gz",
                    "--output_file", f"/data/{output.relative_to(D)}", "--dir_plugins", "/plugins",
                    "--custom", f"{s}/clinvar_chr{chrom}.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN,CLNVC,CLNVI",
                    "--plugin", f"SpliceAI,snv={s}/spliceai_chr{chrom}.vcf.gz,indel={s}/spliceai_chr{chrom}.vcf.gz",
                    "--plugin", f"CADD,snv={s}/cadd_snv_chr{chrom}.tsv.gz,indels={s}/cadd_indel_chr{chrom}.tsv.gz",
                    "--plugin", f"AlphaMissense,file={s}/alphamissense_chr{chrom}.tsv.gz",
                    "--plugin", f"dbNSFP,{s}/dbNSFP5.3.1a_grch38_chr{chrom}.gz,{cols}"]
        r = run_vep(run, f"vep_{tag}_r{rep}", mounts, args, output, expected)
        r.update(mode=mode, repeat=rep)
        with (out / f"{tag}_runs.jsonl").open("a") as f:
            f.write(json.dumps(r) + "\n")
        print(json.dumps(r), flush=True)
        if r["exit_code"] != 0:
            return 1
        time.sleep(10)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
