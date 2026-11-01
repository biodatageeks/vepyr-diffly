#!/usr/bin/env python3
"""One memory measurement method for VEP and VEPYR (thesis item D16, 2026-11-01).

In the timing series the peak RSS of VEP and of VEPYR was taken in different ways (VEP: docker stats of
the container; VEPYR: ru_maxrss of the process on macOS), and for VEP --fork the RSS of the forked
children counts the pages they share with the parent several times. Here both tools run in containers
in the same Colima VM (vz, virtiofs mounts as in the timing series and in D7: image vepyr-vm:479d8f4 for
VEPYR) and one sampler inside the VM reads, every 0.1 s, for all processes of the container's cgroup:
  sum of Rss and sum of Pss (/proc/PID/smaps_rollup; Pss divides each shared page among the processes
  that map it, so the forked children's shared pages are counted once), and
  memory.peak of the container cgroup (/sys/fs/cgroup/docker/ID, cgroupfs driver, cgroup v2: anonymous and page-cache pages charged to it).
Reported per configuration: max over samples of sum Rss and of sum Pss, the last memory.peak read before
the container exited, the number of samples, the docker run wall time (one run: not a timing series)
and the body SHA-256 of the output.

Configurations: chr22 baseline (everything, hgvs, FASTA): VEP without --fork, --fork 4, --fork 10;
VEPYR workers 1 and 10. chr22 with the five plugins: VEP --fork 10 (command and plugin slices of
vep_plugin_stability.sh, B17) and VEPYR workers 10 (plugin cache of the plugin series). chr1 baseline:
VEP --fork 10, VEPYR workers 1 and 10. REPEATS rounds in shuffled order.

Usage: memory_uniform.py OUT_DIR B17_OUT_DIR [REPEATS, default 2]   (Colima running)
"""
from __future__ import annotations

import hashlib
import json
import random
import subprocess
import sys
import time
from pathlib import Path

D = Path("/Users/lukaszjezapkowicz/workspace/data_vepyr")
VEP_IMAGE = "ensemblorg/ensembl-vep:release_116.0"
VEPYR_IMAGE = "vepyr-vm:479d8f4"
FASTA = "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
INPUT = {"chr22": D / "runs/m5_completion_20261024/inputs/HG002_normalized_chr22.vcf.gz",
         "chr1": D / "runs/m5_completion_20261024/inputs/HG002_normalized_chr1.vcf.gz"}
PLUGIN_CACHE = D / "plugin_cache"  # cache/plugin_cache_116 is a symlink to it (plugin cache of the plugin series)
VEPYR_CHILD = r'''
import sys, vepyr
kw = dict(vcf=sys.argv[1], cache_dir="/cache", everything=True, hgvs=True, reference_fasta="/ref/%s",
          workers=int(sys.argv[3]), output_vcf=sys.argv[2], show_progress=False)
if len(sys.argv) > 4:
    kw["plugin_cache_root"] = sys.argv[4]
vepyr.annotate(**kw)
''' % FASTA
SAMPLER = r'''
cg=/sys/fs/cgroup/docker/$1
i=0; while [ ! -d "$cg" ] && [ $i -lt 100 ]; do sleep 0.05; i=$((i+1)); done
maxr=0; maxp=0; peak=0; n=0
while [ -d "$cg" ]; do
  set -- $(for p in $(cat "$cg/cgroup.procs" 2>/dev/null); do cat /proc/$p/smaps_rollup 2>/dev/null; done |
           awk '/^Rss:/{r+=$2} /^Pss:/{p+=$2} END{print r+0, p+0}') "$1"
  [ "$1" -gt "$maxr" ] && maxr=$1
  [ "$2" -gt "$maxp" ] && maxp=$2
  v=$(cat "$cg/memory.peak" 2>/dev/null) && peak=$v
  n=$((n+1)); sleep 0.1
done
echo "$maxr $maxp $peak $n"
'''


def vep_cmd(chrom: str, out: Path, fork: int | None, b17: Path | None) -> list[str]:
    base = ["--rm=false", "--user", "501:20", "--env", "HOME=/tmp"]
    if b17 is None:
        inp = INPUT[chrom]
        cmd = ["docker", "run", "-d", *base,
               "-v", f"{D}/homo_sapiens_merged/116_GRCh38:/opt/vep/.vep/homo_sapiens_merged/116_GRCh38:ro",
               "-v", f"{inp}:/input/{inp.name}:ro", "-v", f"{D}/input:/refdir:ro", "-v", f"{out}:/output",
               VEP_IMAGE, "vep", "--dir", "/opt/vep/.vep", "--cache", "--merged", "--offline", "--assembly", "GRCh38",
               "--input_file", f"/input/{inp.name}", "--output_file", "/output/out.vcf", "--vcf", "--force_overwrite",
               "--no_stats", "--everything", "--hgvs", "--fasta", f"/refdir/{FASTA}"]
    else:  # the --fork 10 command of vep_plugin_stability.sh (B17)
        builder = D / "tools/vepyr/e2e-testing/scripts/build_vep_plugin_reference.sh"
        text = builder.read_text()
        cols = text.split('DBNSFP_COLS="', 1)[1].split('"', 1)[0].replace("\\\n", "").replace("\n", "").replace(" ", "")
        s = "/data/" + str(b17 / "work/slices").removeprefix(str(D) + "/")
        cmd = ["docker", "run", "-d", *base, "-v", f"{D}:/data", "-v", f"{D}/output/116/plugins/plugin_code:/plugins:ro",
               "-v", f"{out}:/output", VEP_IMAGE,
               "vep", "--cache", "--cache_version", "116", "--dir_cache", "/data", "--offline", "--merged",
               "--everything", "--no_stats", "--force_overwrite", "--vcf", "--fasta", f"/data/input/{FASTA}",
               "--input_file", "/data/" + str(b17 / "work/input/HG002_norm_chr22.vcf.gz").removeprefix(str(D) + "/"),
               "--output_file", "/output/out.vcf", "--dir_plugins", "/plugins",
               "--custom", f"{s}/clinvar_chr22.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN,CLNVC,CLNVI",
               "--plugin", f"SpliceAI,snv={s}/spliceai_chr22.vcf.gz,indel={s}/spliceai_chr22.vcf.gz",
               "--plugin", f"CADD,snv={s}/cadd_snv_chr22.tsv.gz,indels={s}/cadd_indel_chr22.tsv.gz",
               "--plugin", f"AlphaMissense,file={s}/alphamissense_chr22.tsv.gz",
               "--plugin", f"dbNSFP,{s}/dbNSFP5.3.1a_grch38_chr22.gz,{cols}"]
    if fork:
        cmd += ["--fork", str(fork)]
    return cmd


def vepyr_cmd(chrom: str, out: Path, workers: int, plugins: bool) -> list[str]:
    inp = INPUT[chrom]
    cmd = ["docker", "run", "-d", "--rm=false", "--user", "501:20", "--env", "HOME=/tmp",
           "-v", f"{D}/cache/116_GRCh38_merged:/cache:ro", "-v", f"{D}/input:/ref:ro",
           "-v", f"{inp.parent}:/input:ro", "-v", f"{out}:/output"]
    if plugins:
        cmd += ["-v", f"{PLUGIN_CACHE}:/plugins:ro"]
    cmd += [VEPYR_IMAGE, "python", "-c", VEPYR_CHILD, f"/input/{inp.name}", "/output/out.vcf", str(workers)]
    return cmd + (["/plugins"] if plugins else [])


def body_sha(path: Path) -> tuple[str, int]:
    h, n = hashlib.sha256(), 0
    with open(path, "rb") as f:
        for line in f:
            if not line.startswith(b"#"):
                h.update(line); n += 1
    return h.hexdigest(), n


def run(cfg: dict, out: Path, b17: Path) -> dict:
    work = out / "work"
    work.mkdir(exist_ok=True)
    (work / "out.vcf").unlink(missing_ok=True)
    if cfg["tool"] == "VEP":
        cmd = vep_cmd(cfg["chrom"], work, cfg["fork"], b17 if cfg["plugins"] else None)
    else:
        cmd = vepyr_cmd(cfg["chrom"], work, cfg["workers"], cfg["plugins"])
    t0 = time.perf_counter()
    cid = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout.strip()
    script = out / "sampler.sh"
    script.write_text(SAMPLER)
    sampler = subprocess.Popen(["colima", "ssh", "--", "sudo", "sh", str(script), cid],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    rc = int(subprocess.run(["docker", "wait", cid], capture_output=True, text=True).stdout.strip())
    wall = time.perf_counter() - t0
    s_out, s_err = sampler.communicate(timeout=120)
    logs = subprocess.run(["docker", "logs", cid], capture_output=True, text=True)
    subprocess.run(["docker", "rm", cid], capture_output=True)
    maxr, maxp, peak, n = (int(x) for x in s_out.split()) if s_out.strip() else (0, 0, 0, 0)
    row = {**cfg, "returncode": rc, "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "wall_s": wall,
           "max_sum_rss_gib": maxr / 2**20, "max_sum_pss_gib": maxp / 2**20, "cgroup_memory_peak_gib": peak / 2**30,
           "samples": n}
    if rc == 0 and (work / "out.vcf").exists():
        row["body_sha256"], row["records"] = body_sha(work / "out.vcf")
    else:
        row["log_tail"] = (logs.stdout + logs.stderr)[-1500:] + s_err[-500:]
    (work / "out.vcf").unlink(missing_ok=True)
    return row


def main() -> int:
    out, b17 = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()
    reps = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    out.mkdir(parents=True, exist_ok=True)
    cfgs = ([{"tool": "VEP", "chrom": "chr22", "plugins": False, "fork": f, "workers": None} for f in (None, 4, 10)]
            + [{"tool": "VEPYR", "chrom": "chr22", "plugins": False, "fork": None, "workers": w} for w in (1, 10)]
            + [{"tool": "VEP", "chrom": "chr22", "plugins": True, "fork": 10, "workers": None},
               {"tool": "VEPYR", "chrom": "chr22", "plugins": True, "fork": None, "workers": 10},
               {"tool": "VEP", "chrom": "chr1", "plugins": False, "fork": 10, "workers": None}]
            + [{"tool": "VEPYR", "chrom": "chr1", "plugins": False, "fork": None, "workers": w} for w in (1, 10)])
    rng = random.Random(20261101)
    rows = []
    with open(out / "runs.jsonl", "a") as log:
        for r in range(1, reps + 1):
            block = list(cfgs)
            rng.shuffle(block)
            for cfg in block:
                row = run({**cfg, "repeat": r}, out, b17)
                rows.append(row)
                log.write(json.dumps(row) + "\n"); log.flush()
                print(json.dumps({k: row.get(k) for k in ("tool", "chrom", "plugins", "fork", "workers", "returncode",
                                                           "wall_s", "max_sum_rss_gib", "max_sum_pss_gib",
                                                           "cgroup_memory_peak_gib", "samples")}), flush=True)
    summary = {}
    for row in rows:
        k = f"{row['tool']} {row['chrom']}{' plugins' if row['plugins'] else ''} " + (
            f"fork={row['fork']}" if row["tool"] == "VEP" else f"workers={row['workers']}")
        s = summary.setdefault(k, {"runs": 0, "max_sum_rss_gib": [], "max_sum_pss_gib": [], "cgroup_memory_peak_gib": [],
                                   "wall_s": [], "bodies": set()})
        s["runs"] += 1
        for m in ("max_sum_rss_gib", "max_sum_pss_gib", "cgroup_memory_peak_gib", "wall_s"):
            s[m].append(round(row[m], 3))
        s["bodies"].add(row.get("body_sha256", "failed"))
    for s in summary.values():
        s["bodies"] = sorted(s["bodies"])
    (out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
