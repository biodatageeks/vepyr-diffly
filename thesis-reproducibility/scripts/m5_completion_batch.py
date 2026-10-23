#!/usr/bin/env python3
"""M5 completion batch for the thesis (decided 2026-10-23).

Waits until the full-genome VEP repeats (runs/vep_full_repeats_m5_20261023) have
finished, then runs, strictly sequentially, the measurements that are missing or
lack raw artifacts. Every step writes its own directory under OUT and the shared
status.json read by tools/m5_progress.py. A failed step is recorded and the batch
continues with steps that do not depend on it.

Usage:
  nohup caffeinate -is python3 tools/m5_completion_batch.py > runs/m5_completion_20261024/batch.out 2>&1 &
  python3 tools/m5_completion_batch.py --dry-run     # path/environment checks only
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import statistics
import subprocess
import sys
import threading
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

D = Path("/Users/lukaszjezapkowicz/workspace/data_vepyr")
THESIS = Path("/Users/lukaszjezapkowicz/Desktop/magisterka/praca/magisterska_jezapkowicz")
OUT = D / "runs/m5_completion_20261024"
VEP_RERUN = D / "runs/vep_full_repeats_m5_20261023"
KEPT_VCF = VEP_RERUN / "kept_reference/HG002_vep_forknone_repeat2.vcf"
REF_WGS = D / "output/116/HG002_annotated_wgs_everything_hgvs_merged.vcf.gz"
PYN = D / ".venv-native/bin/python3"      # vepyr 0.7.0 @ 479d8f4 (annotation series)
PYL = D / ".venv-lpt/bin/python3"         # vepyr-lpt (cache-build-opt-b-native) used on 2026-10-19
CACHE_DIR = D / "cache"
CORE_CACHE = CACHE_DIR / "116_GRCh38_merged"
PLUGIN_ROOT_ALL = D / "cache/plugin_cache_116"
VEP_CACHE = D / "homo_sapiens_merged/116_GRCh38"
FASTA = D / "input/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
IN_WGS = D / "input/HG002_normalized.vcf.gz"
# Recreated by step S0 with the same command as on 2026-09-24 (correctness/commands.log):
# bcftools view -r chrN -Oz input/HG002_normalized.vcf.gz (the original slices were deleted).
IN_CHR = {
    "chr1": OUT / "inputs/HG002_normalized_chr1.vcf.gz",
    "chr22": OUT / "inputs/HG002_normalized_chr22.vcf.gz",
}
MEASURE = D / "tools/measure_once.py"
IMAGE = "ensemblorg/ensembl-vep:release_116.0"
RUN_COMPARISON = D / "tools/vepyr/e2e-testing/scripts/run_comparison.py"
GATE = D / "tools/vepyr/e2e-testing/scripts/verify_parity_gate.py"
E2E_REPORTS = D / "tools/vepyr/e2e-testing/reports"
FULL_BUILD_SCRIPT = THESIS / "runs/full_e2e_all.py"
FULL_BUILD_RUN1 = THESIS / "runs/full_e2e_all"
STATUS = OUT / "status.json"
LOG = OUT / "batch.log"

DBNSFP_COLS = ("SIFT4G_score,SIFT4G_pred,Polyphen2_HDIV_score,Polyphen2_HVAR_score,MutationTaster_score,"
               "MutationTaster_pred,PROVEAN_score,PROVEAN_pred,VEST4_score,MetaSVM_score,MetaSVM_pred,"
               "MetaLR_score,MetaLR_pred,REVEL_score,GERP++_RS,phyloP100way_vertebrate,"
               "phastCons100way_vertebrate,CADD_raw,CADD_phred")
PLUGIN_ARGS = {
    "clinvar": ["--custom", "/pluginwork/clinvar/clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN,CLNVC,CLNVI"],
    "spliceai": ["--plugin", "SpliceAI,snv=/pluginwork/spliceai/spliceai_scores.masked.snv.ensembl_mane.grch38.110.vcf.gz,"
                 "indel=/pluginwork/spliceai/spliceai_scores.masked.snv.ensembl_mane.grch38.110.vcf.gz"],
    "cadd": ["--plugin", "CADD,snv=/pluginwork/cadd/whole_genome_SNVs.tsv.gz,indels=/pluginwork/cadd/gnomad.genomes.r4.0.indel.tsv.gz"],
    "alphamissense": ["--plugin", "AlphaMissense,file=/pluginwork/alphamissense/AlphaMissense_hg38.bgz.tsv.gz"],
    "dbnsfp": ["--plugin", f"dbNSFP,/pluginwork/dbnsfp/dbNSFP5.3.1a_grch38.gz,{DBNSFP_COLS}"],
}

# step id, name, estimated minutes, dependencies
STEPS = [
    ("S0", "Odtworzenie wejść chr1 i chr22 (bcftools view -r, jak 24.09)", 3, []),
    ("S1", "VEP pełny genom -> referencja E2E (bgzip+tabix, kontrola sumy, usunięcie kopii)", 25, []),
    ("S2", "Odbudowa rdzeniowej pamięci podręcznej VEPYR 0.7.0 (cache/116_GRCh38_merged)", 90, []),
    ("S3", "Narzut uruchomienia kontenera Docker, 5x", 2, []),
    ("S4", "Deterministyczność VEPYR workers=1 vs 8 (chr22, chr1)", 6, ["S0", "S2"]),
    ("S5", "VEPYR pełne wejście, workers 10/8/4/2/1, 5x", 65, ["S2"]),
    ("S6", "Test niezależnych procesów na M5 (chr22, n=1,2,4,8; bazowy i komplet wtyczek), 5x", 20, ["S0", "S2"]),
    ("S7", "Nagłówki CSQ VEP dla każdej wtyczki osobno (mały plik wejściowy)", 20, ["S0"]),
    ("S8", "Rozmiary pamięci podręcznej chr1+chr22 (VEP, Parquet, źródła i Parquet wtyczek)", 35, ["S2"]),
    ("S10", "Walidacja E2E pełnego genomu (chr1-22), VEPYR 0.7.0 vs referencja z S1 + bramka", 120, ["S1", "S2"]),
    ("S9", "VEP RSS: scenariusz bazowy chr22 5x i chr1 5x (próbkowanie VmHWM)", 245, ["S0"]),
    ("S11", "Pełna budowa z encjami translacji (środowisko z 19.10), powtórzenia 2-5", 290, []),
]


def now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def log(msg: str) -> None:
    line = f"{now()} {msg}"
    print(line, flush=True)
    with LOG.open("a") as f:
        f.write(line + "\n")


def load_status() -> dict:
    if STATUS.exists():
        return json.loads(STATUS.read_text())
    return {
        "created": now(),
        "phase": "init",
        "steps": [{"id": s, "name": n, "eta_min": e, "deps": d, "state": "pending"} for s, n, e, d in STEPS],
    }


def save_status(st: dict) -> None:
    st["updated"] = now()
    tmp = STATUS.with_suffix(".tmp")
    tmp.write_text(json.dumps(st, indent=2, ensure_ascii=False))
    tmp.replace(STATUS)


def set_step(st: dict, sid: str, **kw) -> None:
    for s in st["steps"]:
        if s["id"] == sid:
            s.update(kw)
    save_status(st)


def run(cmd, log_path: Path | None = None, env=None, timeout=None, check=True) -> subprocess.CompletedProcess:
    cmd = [str(c) for c in cmd]
    if log_path:
        with log_path.open("a") as f:
            f.write(f"$ {' '.join(cmd)}\n")
    p = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=timeout)
    if log_path:
        with log_path.open("a") as f:
            f.write(p.stdout[-20000:] + p.stderr[-20000:] + f"[exit {p.returncode}]\n")
    if check and p.returncode != 0:
        raise RuntimeError(f"command failed ({p.returncode}): {' '.join(cmd)[:300]}\n{p.stderr[-2000:]}")
    return p


def real_time(path: Path) -> float | None:
    m = re.search(r"^real\s+([0-9.]+)", path.read_text(), re.M) if path.exists() else None
    return float(m.group(1)) if m else None


def body_sha256(path: Path) -> tuple[str, int]:
    h = hashlib.sha256()
    n = 0
    opener = subprocess.Popen(["gzip", "-cd", str(path)], stdout=subprocess.PIPE) if str(path).endswith(".gz") else None
    stream = opener.stdout if opener else path.open("rb")
    for line in stream:
        if not line.startswith(b"#"):
            h.update(line)
            n += 1
    if opener:
        opener.wait()
    else:
        stream.close()
    return h.hexdigest(), n


def disk_free_gb() -> float:
    return shutil.disk_usage(D).free / 1e9


# ---------------------------------------------------------------- steps

def s0(d: Path) -> dict:
    res = {}
    for chrom, path in IN_CHR.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        run(["bcftools", "view", "-r", chrom, "-Oz", "-o", path, IN_WGS], d / "cmd.log")
        run(["bcftools", "index", "-f", "-t", path], d / "cmd.log")
        n = int(subprocess.run(f"bcftools view -H {path} | wc -l", shell=True, capture_output=True, text=True).stdout)
        res[chrom] = {"path": str(path), "records": n}
    expected = {"chr1": 323430, "chr22": 50861}
    if any(res[c]["records"] != expected[c] for c in expected):
        raise RuntimeError(f"unexpected record counts {res}")
    return res


def s1(d: Path) -> dict:
    rec_expected = (VEP_RERUN / "repeat_2/output_records.txt").read_text().strip()
    sha_expected = (VEP_RERUN / "repeat_2/output_body.sha256").read_text().strip()
    if REF_WGS.exists() and not KEPT_VCF.exists():
        sha, n = body_sha256(REF_WGS)
        return {"skipped_compress": True, "body_sha256": sha, "records": n, "matches_repeat2": sha == sha_expected}
    REF_WGS.parent.mkdir(parents=True, exist_ok=True)
    tmp = REF_WGS.with_suffix(".tmp.gz")
    with tmp.open("wb") as f:
        p = subprocess.run(["bgzip", "-@", "8", "-c", str(KEPT_VCF)], stdout=f)
    if p.returncode != 0:
        raise RuntimeError("bgzip failed")
    tmp.replace(REF_WGS)
    run(["tabix", "-f", "-p", "vcf", REF_WGS], d / "cmd.log")
    sha, n = body_sha256(REF_WGS)
    ok = sha == sha_expected and str(n) == rec_expected
    res = {"reference": str(REF_WGS), "bytes": REF_WGS.stat().st_size, "body_sha256": sha, "records": n,
           "expected_sha256": sha_expected, "expected_records": rec_expected, "matches_repeat2": ok}
    if not ok:
        raise RuntimeError(f"compressed reference does not match repeat 2: {res}")
    KEPT_VCF.unlink()
    res["kept_vcf_removed"] = True
    return res


def s2(d: Path) -> dict:
    code = (
        "import json,time,vepyr,importlib.metadata as m\n"
        "t0=time.perf_counter()\n"
        f"r=vepyr.build_cache(116, {str(CACHE_DIR)!r}, cache_type='merged', local_cache={str(VEP_CACHE)!r}, overwrite=True, show_progress=False)\n"
        "print(json.dumps({'wall_seconds': time.perf_counter()-t0, 'files': len(r), 'vepyr': m.version('vepyr')}))\n"
    )
    p = run([PYN, "-c", code], d / "cmd.log", timeout=6 * 3600)
    res = json.loads(p.stdout.strip().splitlines()[-1])
    res["bytes"] = sum(f.stat().st_size for f in CORE_CACHE.rglob("*") if f.is_file())
    res["entities"] = sorted(x.name for x in CORE_CACHE.iterdir() if x.is_dir())
    return res


def s3(d: Path) -> dict:
    times = []
    for i in range(1, 6):
        t = d / f"time_{i}.txt"
        run(["/usr/bin/time", "-p", "-o", t, "docker", "run", "--rm", IMAGE, "true"], d / "cmd.log")
        times.append(real_time(t))
        time.sleep(2)
    return {"runs_s": times, "mean_s": statistics.mean(times), "median_s": statistics.median(times),
            "min_s": min(times), "max_s": max(times)}


ANNOTATE_KEEP = r'''
import hashlib, json, sys, time, vepyr
vcf, out, workers = sys.argv[1], sys.argv[2], int(sys.argv[3])
t0 = time.perf_counter()
vepyr.annotate(vcf=vcf, cache_dir=sys.argv[4], everything=True, hgvs=True, reference_fasta=sys.argv[5],
               workers=workers, output_vcf=out, show_progress=False)
dt = time.perf_counter() - t0
h = hashlib.sha256(); n = 0
with open(out, "rb") as f:
    for line in f:
        if not line.startswith(b"#"):
            h.update(line); n += 1
print(json.dumps({"annotation_seconds": dt, "body_sha256": h.hexdigest(), "records": n}))
'''


def s4(d: Path) -> dict:
    res = {}
    for chrom in ("chr22", "chr1"):
        for w in (1, 8):
            out = d / f"{chrom}_w{w}.vcf"
            p = run([PYN, "-c", ANNOTATE_KEEP, IN_CHR[chrom], out, w, CORE_CACHE, FASTA], d / "cmd.log", timeout=3600)
            res[f"{chrom}_w{w}"] = json.loads(p.stdout.strip().splitlines()[-1])
            out.unlink(missing_ok=True)
        res[f"{chrom}_identical"] = res[f"{chrom}_w1"]["body_sha256"] == res[f"{chrom}_w8"]["body_sha256"]
    return res


def measure(input_vcf: Path, workers: int, out_vcf: Path, plugin_root: Path | None = None) -> dict:
    cmd = [PYN, MEASURE, "--input-vcf", input_vcf, "--cache-dir", CORE_CACHE, "--reference-fasta", FASTA,
           "--workers", workers, "--output-vcf", out_vcf]
    if plugin_root:
        cmd += ["--plugin-cache-root", plugin_root]
    t0 = time.perf_counter()
    p = run(cmd, timeout=4 * 3600)
    r = json.loads(p.stdout.strip().splitlines()[-1])
    r["wall_seconds"] = time.perf_counter() - t0
    return r


def s5(d: Path) -> dict:
    raw = d / "raw_runs.jsonl"
    done = {(r["repeat"], r["workers"]) for r in map(json.loads, raw.read_text().splitlines())} if raw.exists() else set()
    for rep in range(1, 6):
        for w in (10, 8, 4, 2, 1):
            if (rep, w) in done:
                continue
            r = measure(IN_WGS, w, D / "output/5repeats_tmp/wgs_completion.vcf")
            r.update(repeat=rep, workers=w, timestamp=now())
            with raw.open("a") as f:
                f.write(json.dumps(r) + "\n")
            log(f"S5 repeat {rep} workers {w}: {r['annotation_seconds']:.2f}s RSS {r['max_rss_gb']} GiB")
            time.sleep(3)
    runs = list(map(json.loads, raw.read_text().splitlines()))
    summ = {}
    for w in (1, 2, 4, 8, 10):
        xs = [r for r in runs if r["workers"] == w]
        summ[str(w)] = {"n": len(xs), "mean_s": statistics.mean(r["annotation_seconds"] for r in xs),
                        "mean_rss_gib": statistics.mean(r["max_rss_gb"] for r in xs),
                        "rsd_pct": 100 * statistics.stdev(r["annotation_seconds"] for r in xs) / statistics.mean(r["annotation_seconds"] for r in xs)}
    return summ


def s6(d: Path) -> dict:
    raw = d / "raw_runs.jsonl"
    done = {(r["scenario"], r["n"], r["repeat"]) for r in map(json.loads, raw.read_text().splitlines())} if raw.exists() else set()
    for scenario, root in (("baseline", None), ("all_plugins", PLUGIN_ROOT_ALL)):
        for n in (1, 2, 4, 8):
            for rep in range(1, 6):
                if (scenario, n, rep) in done:
                    continue
                procs = []
                t0 = time.perf_counter()
                for i in range(n):
                    cmd = [str(c) for c in [PYN, MEASURE, "--input-vcf", IN_CHR["chr22"], "--cache-dir", CORE_CACHE,
                                            "--reference-fasta", FASTA, "--workers", 1,
                                            "--output-vcf", d / f"mp_{scenario}_{n}_{i}.vcf"]]
                    if root:
                        cmd += ["--plugin-cache-root", str(root)]
                    procs.append(subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True))
                outs = [p.communicate() for p in procs]
                group = time.perf_counter() - t0
                if any(p.returncode for p in procs):
                    raise RuntimeError(f"S6 process failed: {[o[1][-500:] for o in outs]}")
                per = [json.loads(o[0].strip().splitlines()[-1]) for o in outs]
                r = {"scenario": scenario, "n": n, "repeat": rep, "group_wall_s": group,
                     "annotation_s": [x["annotation_seconds"] for x in per], "timestamp": now()}
                with raw.open("a") as f:
                    f.write(json.dumps(r) + "\n")
                time.sleep(3)
    runs = list(map(json.loads, raw.read_text().splitlines()))
    summ = {}
    for scenario in ("baseline", "all_plugins"):
        t1 = statistics.mean(r["group_wall_s"] for r in runs if r["scenario"] == scenario and r["n"] == 1)
        for n in (1, 2, 4, 8):
            xs = [r["group_wall_s"] for r in runs if r["scenario"] == scenario and r["n"] == n]
            summ[f"{scenario}_n{n}"] = {"mean_group_s": statistics.mean(xs), "median_group_s": statistics.median(xs),
                                         "ratio_mean": statistics.mean(xs) / t1, "n_runs": len(xs)}
    return summ


def vep_docker(name: str, mounts: list[str], vep_args: list[str], d: Path, tag: str, sample_rss: bool):
    t = d / f"{tag}.time.txt"
    cmd = ["/usr/bin/time", "-p", "-o", str(t), "docker", "run", "--rm", "--name", name, "--user", "501:20",
           "--env", "HOME=/tmp", *mounts, IMAGE, "vep", *vep_args]
    (d / f"{tag}.command.txt").write_text(" ".join(cmd) + "\n")
    subprocess.run(["docker", "rm", "-f", name], capture_output=True)
    proc = subprocess.Popen(cmd, stdout=(d / f"{tag}.stdout.txt").open("w"), stderr=(d / f"{tag}.stderr.txt").open("w"))
    peak = [0]
    sampler = ('for p in /proc/[0-9]*; do c=$(cat $p/comm 2>/dev/null); case "$c" in perl|vep) '
               'awk "/^VmHWM/{print \\$2}" $p/status;; esac; done')

    def sample():
        while proc.poll() is None:
            s = subprocess.run(["docker", "exec", name, "sh", "-c", sampler], capture_output=True, text=True)
            vals = [int(x) for x in s.stdout.split() if x.isdigit()]
            if vals:
                peak[0] = max(peak[0], max(vals))
            time.sleep(5)
    th = threading.Thread(target=sample, daemon=True)
    if sample_rss:
        th.start()
    rc = proc.wait()
    if sample_rss:
        th.join(timeout=30)
    return rc, real_time(t), peak[0]


def s7(d: Path) -> dict:
    tiny = d / "tiny_chr22.vcf.gz"
    with tiny.open("wb") as f:
        hdr = subprocess.run(["bcftools", "view", "-h", str(IN_CHR["chr22"])], capture_output=True).stdout
        body = subprocess.run(f"bcftools view -H {IN_CHR['chr22']} | head -300", shell=True, capture_output=True).stdout
        f.write(subprocess.run(["bgzip", "-c"], input=hdr + body, capture_output=True).stdout)
    plugdir = d / "plugins"
    plugdir.mkdir(exist_ok=True)
    for name, src in (("AlphaMissense.pm", "alphamissense/AlphaMissense.pm"), ("CADD.pm", "cadd/CADD.pm"),
                      ("dbNSFP.pm", "dbnsfp/dbNSFP.pm"), ("SpliceAI.pm", "spliceai/SpliceAI.pm")):
        shutil.copy2(D / "plugin_input" / src, plugdir / name)
    fixed_cadd = D / "tmp/CADD.pm.7a1f6450fd12"
    if fixed_cadd.exists():
        shutil.copy2(fixed_cadd, plugdir / "CADD.pm")
    res = {}
    for scen in ["baseline", *PLUGIN_ARGS]:
        out = d / f"{scen}.vcf"
        mounts = ["-v", f"{VEP_CACHE}:/opt/vep/.vep/homo_sapiens_merged/116_GRCh38:ro", "-v", f"{tiny}:/input.vcf.gz:ro",
                  "-v", f"{D}/input:/refdir:ro", "-v", f"{d}:/output"]
        args = ["--dir", "/opt/vep/.vep", "--cache", "--merged", "--offline", "--assembly", "GRCh38",
                "--input_file", "/input.vcf.gz", "--output_file", f"/output/{out.name}", "--vcf", "--force_overwrite",
                "--no_stats", "--everything", "--hgvs", "--fasta", "/refdir/Homo_sapiens.GRCh38.dna.primary_assembly.fa"]
        if scen != "baseline":
            mounts += ["-v", f"{D}/plugin_input:/pluginwork:ro", "-v", f"{plugdir}:/plugins:ro"]
            args += ["--dir_plugins", "/plugins", *PLUGIN_ARGS[scen]]
        rc, wall, _ = vep_docker(f"vep_hdr_{scen}", mounts, args, d, scen, False)
        fields = []
        if out.exists():
            for line in out.open():
                if line.startswith("##INFO=<ID=CSQ"):
                    fields = line.split("Format: ")[1].rstrip('">\n').split("|")
                    break
        res[scen] = {"exit": rc, "wall_s": wall, "n_fields": len(fields), "fields": fields}
        out.unlink(missing_ok=True)
    base = set(res["baseline"]["fields"])
    for scen in PLUGIN_ARGS:
        res[scen]["added_fields"] = [f for f in res[scen]["fields"] if f not in base]
    return res


def tabix_names(path: Path) -> list[str]:
    return subprocess.run(["tabix", "-l", str(path)], capture_output=True, text=True).stdout.split()


def s8(d: Path) -> dict:
    def du(paths):
        return sum(f.stat().st_size for p in paths for f in ([p] if p.is_file() else p.rglob("*")) if f.is_file())
    res = {"vep_native_chr1_chr22_bytes": du([VEP_CACHE / "1", VEP_CACHE / "22"])}
    core = [f for f in CORE_CACHE.rglob("*") if f.is_file() and re.fullmatch(r"chr(1|22)\.parquet", f.name)]
    res["vepyr_core_files"] = sorted(str(f.relative_to(CORE_CACHE)) for f in core)
    res["vepyr_core_chr1_chr22_bytes"] = du(core)
    plug_root = D / "plugin_cache/plugin"
    plug = [f for f in plug_root.rglob("*") if f.is_file() and re.search(r"(^|[^0-9])(chr)?(1|22)\.parquet$", f.name)]
    res["plugin_parquet_files"] = sorted(str(f.relative_to(plug_root)) for f in plug)
    res["plugin_parquet_chr1_chr22_bytes"] = du(plug)
    sources = {
        "clinvar": D / "plugin_input/clinvar/clinvar.vcf.gz",
        "spliceai": D / "plugin_input/spliceai/spliceai_scores.masked.snv.ensembl_mane.grch38.110.vcf.gz",
        "cadd_snv": D / "plugin_input/cadd/whole_genome_SNVs.tsv.gz",
        "cadd_indel": D / "plugin_input/cadd/gnomad.genomes.r4.0.indel.tsv.gz",
        "alphamissense": D / "plugin_input/alphamissense/AlphaMissense_hg38.bgz.tsv.gz",
        "dbnsfp": D / "plugin_input/dbnsfp/dbNSFP5.3.1a_grch38.gz",
    }
    slices = {}
    for key, src in sources.items():
        names = tabix_names(src)
        for want in ("1", "22"):
            seq = want if want in names else (f"chr{want}" if f"chr{want}" in names else None)
            if seq is None:
                slices[f"{key}_{want}"] = None
                continue
            tmp = d / f"slice_{key}_{want}.gz"
            with tmp.open("wb") as f:
                p1 = subprocess.Popen(["tabix", "-h", str(src), seq], stdout=subprocess.PIPE)
                subprocess.run(["bgzip", "-@", "4", "-c"], stdin=p1.stdout, stdout=f)
                p1.wait()
            slices[f"{key}_{want}"] = tmp.stat().st_size
            tmp.unlink()
    res["plugin_source_slices_bytes"] = slices
    res["plugin_source_chr1_chr22_bytes"] = sum(v for v in slices.values() if v)
    return res


def s10(d: Path) -> dict:
    backup = d / "reports_backup_before"
    if not backup.exists():
        shutil.copytree(E2E_REPORTS, backup)
    env = dict(os.environ, DATA_VEPYR_DIR=str(D))
    p = run([PYN, RUN_COMPARISON, "--release", "116", "--profile", "merged", "--chroms", "1-22", "--workers", "4",
             "--bgzf", "--force", "--vep", REF_WGS], d / "run_comparison.log", env=env, timeout=10 * 3600, check=False)
    newdir = d / "reports"
    newdir.mkdir(exist_ok=True)
    for f in E2E_REPORTS.glob("fast_chr*_merged_116_*"):
        shutil.copy2(f, newdir / f.name)
    g = run([PYN, GATE, "--release", "116", "--profile", "merged", "--chroms", "1-22", "--report-dir", newdir],
            d / "gate.log", env=env, check=False)
    # restore the 2026-09-24 evidence in the shared reports directory
    for f in backup.iterdir():
        shutil.copy2(f, E2E_REPORTS / f.name)
    summ = {"run_comparison_exit": p.returncode, "gate_exit": g.returncode, "gate_tail": g.stdout[-3000:]}
    tot = {"variants": 0, "only_vepyr": 0, "only_vep": 0, "csq_count_mismatch": 0, "csq_order_mismatch": 0}
    for f in sorted(newdir.glob("fast_chr*_merged_116_report.json")):
        c = json.loads(f.read_text())["comparison"]
        tot["variants"] += c["variants_compared"]
        tot["only_vepyr"] += c["variants_only_in_vepyr"]
        tot["only_vep"] += c["variants_only_in_vep"]
        tot["csq_count_mismatch"] += c["csq_entry_count_mismatch"]
        tot["csq_order_mismatch"] += c["csq_order_mismatch"]
    summ["totals"] = tot
    if p.returncode != 0:
        raise RuntimeError(f"run_comparison exit {p.returncode}; see run_comparison.log; {summ}")
    return summ


def s9(d: Path) -> dict:
    raw = d / "raw_runs.jsonl"
    done = {(r["chrom"], r["repeat"]) for r in map(json.loads, raw.read_text().splitlines())} if raw.exists() else set()
    for chrom in ("chr22", "chr1"):
        for rep in range(1, 6):
            if (chrom, rep) in done:
                continue
            out = d / f"{chrom}_r{rep}.vcf"
            mounts = ["-v", f"{VEP_CACHE}:/opt/vep/.vep/homo_sapiens_merged/116_GRCh38:ro",
                      "-v", f"{IN_CHR[chrom]}:/input.vcf.gz:ro", "-v", f"{D}/input:/refdir:ro", "-v", f"{d}:/output"]
            args = ["--dir", "/opt/vep/.vep", "--cache", "--merged", "--offline", "--assembly", "GRCh38",
                    "--input_file", "/input.vcf.gz", "--output_file", f"/output/{out.name}", "--vcf",
                    "--force_overwrite", "--no_stats", "--everything", "--hgvs",
                    "--fasta", "/refdir/Homo_sapiens.GRCh38.dna.primary_assembly.fa"]
            rc, wall, hwm_kb = vep_docker(f"vep_rss_{chrom}_{rep}", mounts, args, d, f"{chrom}_r{rep}", True)
            n = sum(1 for line in out.open() if not line.startswith("#")) if out.exists() else None
            out.unlink(missing_ok=True)
            r = {"chrom": chrom, "repeat": rep, "exit": rc, "wall_s": wall, "peak_rss_kb": hwm_kb,
                 "peak_rss_gb": round(hwm_kb * 1024 / 1e9, 3), "peak_rss_gib": round(hwm_kb / 1024**2, 3),
                 "records": n, "timestamp": now()}
            with raw.open("a") as f:
                f.write(json.dumps(r) + "\n")
            log(f"S9 {chrom} repeat {rep}: exit {rc}, {wall}s, VmHWM {r['peak_rss_gib']} GiB")
            if rc != 0:
                raise RuntimeError(f"VEP failed for {chrom} repeat {rep}")
            time.sleep(5)
    runs = list(map(json.loads, raw.read_text().splitlines()))
    return {c: {"mean_wall_s": statistics.mean(r["wall_s"] for r in runs if r["chrom"] == c),
                "mean_peak_rss_gib": statistics.mean(r["peak_rss_gib"] for r in runs if r["chrom"] == c),
                "mean_peak_rss_gb": statistics.mean(r["peak_rss_gb"] for r in runs if r["chrom"] == c)}
            for c in ("chr22", "chr1")}


def s11(d: Path) -> dict:
    r1 = d / "repeat_1"
    if not r1.exists():
        r1.mkdir(parents=True)
        for name in ("summary.json", "progress.json", "command.json"):
            shutil.copy2(FULL_BUILD_RUN1 / name, r1 / name)
    old_cache = FULL_BUILD_RUN1 / "cache"
    if old_cache.exists():
        log(f"S11 removing output of the 2026-10-19 run ({old_cache}); its JSON artifacts are kept")
        shutil.rmtree(old_cache)
    src = FULL_BUILD_SCRIPT.read_text()
    for rep in range(2, 6):
        rd = d / f"repeat_{rep}"
        if (rd / "summary.json").exists():
            continue
        rd.mkdir(exist_ok=True)
        script = rd / "full_build.py"
        script.write_text(src.replace('ROOT = Path(__file__).resolve().parent / "full_e2e_all"',
                                      f'ROOT = Path({str(rd)!r})'))
        if f"ROOT = Path({str(rd)!r})" not in script.read_text():
            raise RuntimeError("could not retarget ROOT in full_e2e_all.py copy")
        run([PYL, script], rd / "run.log", timeout=4 * 3600)
        if (rd / "cache").exists():
            shutil.rmtree(rd / "cache")
        log(f"S11 repeat {rep}: {json.loads((rd / 'summary.json').read_text())['wall_seconds']:.1f}s")
    walls = [json.loads((d / f"repeat_{i}/summary.json").read_text())["wall_seconds"] for i in range(1, 6)]
    return {"wall_seconds": walls, "mean_s": statistics.mean(walls), "rsd_pct": 100 * statistics.stdev(walls) / statistics.mean(walls)}


FUNCS = {"S0": s0, "S1": s1, "S2": s2, "S3": s3, "S4": s4, "S5": s5, "S6": s6, "S7": s7, "S8": s8, "S9": s9, "S10": s10, "S11": s11}


def vep_reruns_finished() -> bool:
    if not (VEP_RERUN / "summary.json").exists():
        return False
    p = subprocess.run(["docker", "ps", "-q", "--filter", "name=vep_full_repeat_"], capture_output=True, text=True)
    return p.stdout.strip() == ""


def dry_run() -> int:
    ok = True
    for p in [D, VEP_RERUN, KEPT_VCF, PYN, PYL, VEP_CACHE, FASTA, IN_WGS, MEASURE, RUN_COMPARISON,
              GATE, E2E_REPORTS, FULL_BUILD_SCRIPT, FULL_BUILD_RUN1 / "summary.json", PLUGIN_ROOT_ALL]:
        e = Path(p).exists()
        ok &= e
        print(f"{'OK ' if e else 'MISSING'} {p}")
    for py, mod in ((PYN, "import vepyr,importlib.metadata as m; print('vepyr', m.version('vepyr'))"),
                    (PYL, "import vepyr; print('vepyr-lpt build_cache_entity', hasattr(vepyr,'build_cache_entity'))")):
        p = subprocess.run([str(py), "-c", mod], capture_output=True, text=True)
        print(p.stdout.strip() or p.stderr.strip()[-300:])
    src = FULL_BUILD_SCRIPT.read_text()
    print("full_e2e_all ROOT line present:", 'ROOT = Path(__file__).resolve().parent / "full_e2e_all"' in src)
    print(f"disk free: {disk_free_gb():.1f} GB")
    return 0 if ok else 1


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--only", nargs="*", help="run only these step ids")
    a = ap.parse_args()
    if a.dry_run:
        return dry_run()
    OUT.mkdir(parents=True, exist_ok=True)
    st = load_status()
    st["pid"] = os.getpid()
    st["phase"] = "waiting_for_vep_reruns"
    save_status(st)
    log("batch started; waiting for VEP full-genome repeats to finish")
    while not vep_reruns_finished():
        time.sleep(60)
    st["phase"] = "running"
    save_status(st)
    log(f"VEP repeats finished; disk free {disk_free_gb():.1f} GB")
    failed = set()
    for sid, name, _, deps in STEPS:
        if a.only and sid not in a.only:
            continue
        cur = next(s for s in st["steps"] if s["id"] == sid)
        if cur["state"] == "done":
            continue
        if any(x in failed for x in deps):
            set_step(st, sid, state="skipped", detail=f"dependency failed: {deps}")
            failed.add(sid)
            continue
        d = OUT / sid
        d.mkdir(exist_ok=True)
        set_step(st, sid, state="running", started=now(), disk_free_gb=round(disk_free_gb(), 1))
        log(f"{sid} start: {name}")
        t0 = time.time()
        try:
            res = FUNCS[sid](d)
            (d / "result.json").write_text(json.dumps(res, indent=2, ensure_ascii=False, default=str))
            set_step(st, sid, state="done", finished=now(), minutes=round((time.time() - t0) / 60, 1))
            log(f"{sid} done in {(time.time() - t0) / 60:.1f} min")
        except Exception as e:  # noqa: BLE001 - record and continue
            (d / "error.txt").write_text(traceback.format_exc())
            set_step(st, sid, state="failed", finished=now(), minutes=round((time.time() - t0) / 60, 1), detail=str(e)[:500])
            log(f"{sid} FAILED: {e}")
            failed.add(sid)
    st["phase"] = "finished"
    save_status(st)
    log("batch finished")
    return 0


if __name__ == "__main__":
    sys.exit(main())
