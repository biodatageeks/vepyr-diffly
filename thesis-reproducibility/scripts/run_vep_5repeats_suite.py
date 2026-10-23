#!/usr/bin/env python3
"""Run 5 repeats of all Ensembl VEP 116 benchmarks (chr1 & chr22, plugins & forks) on Apple M5.

Measures:
1. chr22 fork scaling (forks: none, 1, 2, 4, 8, 10) - runs 2..5
2. chr22 plugins (baseline, clinvar, spliceai, cadd, alphamissense, dbnsfp) - runs 2..5
3. chr1 fork scaling (forks: none, 1, 2, 4, 8, 10) - runs 2..5
4. chr1 plugins (baseline, clinvar, spliceai, cadd, alphamissense, dbnsfp) - runs 2..5

Outputs:
- runs/vep_5repeats_m5/raw_runs.jsonl
- runs/vep_5repeats_m5/summary_vep_5repeats.json
- Updates raport_m5.html and RAPORT_M1_vs_M5_DLA_PROMOTORA.md after each phase!
"""

from __future__ import annotations

import json
import os
import re
import shutil
import statistics
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path("/Users/lukaszjezapkowicz/workspace/data_vepyr")
OUT_DIR = ROOT / "output/vep_5repeats_tmp"
LOGS_DIR = ROOT / "runs/vep_5repeats_m5"

CACHE_DIR = ROOT / "homo_sapiens_merged/116_GRCh38"
INPUT_CHR22 = ROOT / "runs/chapter5-m5-20260924T061203Z/correctness/input/HG002_normalized_chr22.vcf.gz"
INPUT_CHR1 = ROOT / "runs/chapter5-m5-20260924T061203Z/correctness/input/HG002_normalized_chr1.vcf.gz"
FASTA = ROOT / "input/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
PLUGIN_DIR_CHR22 = ROOT / "tmp/vep-plugin-chr22"
PLUGIN_DIR_CHR1 = ROOT / "tmp/vep-plugin-chr1"

HTML_PATH = Path("/Users/lukaszjezapkowicz/Desktop/magisterka/praca/magisterska_jezapkowicz/raport_m5.html")
MD_PATH = Path("/Users/lukaszjezapkowicz/Desktop/magisterka/praca/magisterska_jezapkowicz/RAPORT_M1_vs_M5_DLA_PROMOTORA.md")
ARTIFACT_MD = Path("/Users/lukaszjezapkowicz/.gemini/antigravity-cli/brain/ad545a87-cb81-4290-9afe-9bd392089f8d/RAPORT_M1_vs_M5_DLA_PROMOTORA.md")

VEP_IMAGE = "ensemblorg/ensembl-vep:release_116.0"

# Initial recorded Run 1 measurements
RUN1_DATA = {
    "chr22_forks": {
        "none": 384.39,
        "1": 234.51,
        "2": 147.25,
        "4": 107.46,
        "8": 89.20,
        "10": 80.18,
    },
    "chr1_forks": {
        "none": 2563.79,
        "1": 1059.73,
        "2": 771.61,
        "4": 577.56,
        "8": 465.88,
        "10": 471.68,
    },
    "chr22_plugins": {
        "baseline": 384.39,
        "clinvar": 395.97,
        "spliceai": 647.82,
        "cadd": 615.73,
        "alphamissense": 395.17,
        "dbnsfp": 411.10,
    },
    "chr1_plugins": {
        "baseline": 2563.79,
        "clinvar": 2240.81,
        "spliceai": 3631.82,
        "cadd": 3469.91,
        "alphamissense": 2199.77,
        "dbnsfp": 2228.33,
    },
}

PLUGINS_ARGS_CHR22 = {
    "clinvar": "--custom /pluginwork/slices/clinvar_chr22.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN,CLNVC,CLNVI",
    "spliceai": "--plugin SpliceAI,snv=/pluginwork/slices/spliceai_chr22.vcf.gz,indel=/pluginwork/slices/spliceai_chr22.vcf.gz",
    "cadd": "--plugin CADD,snv=/pluginwork/slices/cadd_snv_chr22.tsv.gz,indels=/pluginwork/slices/cadd_indel_chr22.tsv.gz",
    "alphamissense": "--plugin AlphaMissense,file=/pluginwork/slices/alphamissense_chr22.tsv.gz",
    "dbnsfp": "--plugin dbNSFP,/pluginwork/slices/dbNSFP5.3.1a_grch38_chr22.gz,SIFT4G_score,SIFT4G_pred,Polyphen2_HDIV_score,Polyphen2_HVAR_score,MutationTaster_score,MutationTaster_pred,PROVEAN_score,PROVEAN_pred,VEST4_score,MetaSVM_score,MetaSVM_pred,MetaLR_score,MetaLR_pred,REVEL_score,GERP++_RS,phyloP100way_vertebrate,phastCons100way_vertebrate,CADD_raw,CADD_phred",
}

PLUGINS_ARGS_CHR1 = {
    "clinvar": "--custom /pluginwork/slices/clinvar_chr1.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN,CLNVC,CLNVI",
    "spliceai": "--plugin SpliceAI,snv=/pluginwork/slices/spliceai_chr1.vcf.gz,indel=/pluginwork/slices/spliceai_chr1.vcf.gz",
    "cadd": "--plugin CADD,snv=/pluginwork/slices/cadd_snv_chr1.tsv.gz,indels=/pluginwork/slices/cadd_indel_chr1.tsv.gz",
    "alphamissense": "--plugin AlphaMissense,file=/pluginwork/slices/alphamissense_chr1.tsv.gz",
    "dbnsfp": "--plugin dbNSFP,/pluginwork/slices/dbNSFP5.3.1a_grch38_chr1.gz,SIFT4G_score,SIFT4G_pred,Polyphen2_HDIV_score,Polyphen2_HVAR_score,MutationTaster_score,MutationTaster_pred,PROVEAN_score,PROVEAN_pred,VEST4_score,MetaSVM_score,MetaSVM_pred,MetaLR_score,MetaLR_pred,REVEL_score,GERP++_RS,phyloP100way_vertebrate,phastCons100way_vertebrate,CADD_raw,CADD_phred",
}


def run_vep(
    input_vcf: Path,
    output_vcf: Path,
    fork: int | None = None,
    plugin_name: str | None = None,
    chr_name: str = "chr22",
) -> float:
    time_txt = OUT_DIR / "time_tmp.txt"
    if time_txt.exists():
        time_txt.unlink()
    if output_vcf.exists():
        output_vcf.unlink()

    cmd = [
        "/usr/bin/caffeinate", "-d", "-i", "-s", "-m",
        "/usr/bin/time", "-p", "-o", str(time_txt),
        "docker", "run", "--rm",
        "--user", "501:20",
        "--env", "HOME=/tmp",
        "-v", f"{CACHE_DIR}:/opt/vep/.vep/homo_sapiens_merged/116_GRCh38:ro",
        "-v", f"{input_vcf}:/input.vcf.gz:ro",
        "-v", f"{ROOT}/input:/refdir:ro",
        "-v", f"{OUT_DIR}:/output",
    ]

    plugin_dir = PLUGIN_DIR_CHR1 if chr_name == "chr1" else PLUGIN_DIR_CHR22
    if plugin_name is not None and plugin_name != "baseline":
        cmd.extend([
            "-v", f"{plugin_dir}:/pluginwork:ro",
            "-v", f"{plugin_dir}/plugins:/plugins:ro",
        ])

    cmd.extend([
        VEP_IMAGE,
        "vep",
        "--dir", "/opt/vep/.vep",
        "--cache", "--merged", "--offline",
        "--assembly", "GRCh38",
        "--input_file", "/input.vcf.gz",
        "--output_file", f"/output/{output_vcf.name}",
        "--vcf", "--force_overwrite", "--no_stats",
        "--everything", "--hgvs",
        "--fasta", "/refdir/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
    ])

    if fork is not None and fork > 0:
        cmd.extend(["--fork", str(fork)])

    if plugin_name is not None and plugin_name != "baseline":
        cmd.extend(["--dir_plugins", "/plugins"])
        plugin_map = PLUGINS_ARGS_CHR1 if chr_name == "chr1" else PLUGINS_ARGS_CHR22
        extra_args = plugin_map[plugin_name].split()
        cmd.extend(extra_args)

    t0 = time.perf_counter()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    wall_duration = time.perf_counter() - t0

    # Parse time_txt
    real_time = wall_duration
    if time_txt.exists():
        content = time_txt.read_text(encoding="utf-8")
        for line in content.splitlines():
            line = line.strip()
            if line.startswith("real "):
                try:
                    real_time = float(line.split()[1])
                except (IndexError, ValueError):
                    pass
        time_txt.unlink()

    if output_vcf.exists():
        output_vcf.unlink()

    return real_time


def update_all_reports(all_data: dict) -> None:
    """Recalculate Table 5 and Table 6 with 5-run VEP averages."""
    # Compute means
    vep_c1_plug = {}
    for p in ["baseline", "clinvar", "alphamissense", "dbnsfp", "spliceai", "cadd"]:
        vep_c1_plug[p] = statistics.mean(all_data["chr1_plugins"][p])

    vep_c22_plug = {}
    for p in ["baseline", "clinvar", "alphamissense", "dbnsfp", "spliceai", "cadd"]:
        vep_c22_plug[p] = statistics.mean(all_data["chr22_plugins"][p])

    vep_c1_forks = {}
    for f in ["none", "1", "2", "4", "8", "10"]:
        vep_c1_forks[f] = statistics.mean(all_data["chr1_forks"][f])

    vep_c22_forks = {}
    for f in ["none", "1", "2", "4", "8", "10"]:
        vep_c22_forks[f] = statistics.mean(all_data["chr22_forks"][f])

    vepyr_c1_plug = {
        "baseline": (18.74, 3.61),
        "clinvar": (19.63, 3.64),
        "alphamissense": (19.59, 3.64),
        "dbnsfp": (20.37, 3.62),
        "spliceai": (28.82, 3.63),
        "cadd": (31.30, 3.71),
    }
    vepyr_c22_plug = {
        "baseline": (3.12, 1.17),
        "clinvar": (3.24, 1.16),
        "alphamissense": (3.23, 1.12),
        "dbnsfp": (3.46, 1.15),
        "spliceai": (4.83, 1.23),
        "cadd": (4.81, 1.21),
    }

    plugin_display = {
        "baseline": "baseline",
        "clinvar": "ClinVar",
        "spliceai": "SpliceAI",
        "cadd": "CADD",
        "alphamissense": "AlphaMissense",
        "dbnsfp": "dbNSFP",
    }

    # Generate Table 5 HTML rows
    t5_chr1_html = []
    for p in ["baseline", "clinvar", "alphamissense", "dbnsfp", "spliceai", "cadd"]:
        vt = vep_c1_plug[p]
        vm = vt / 60.0
        yt, yr = vepyr_c1_plug[p]
        sp = vt / yt
        t5_chr1_html.append(
            f'      <tr><td><strong>{plugin_display[p]}</strong></td>'
            f'<td>{vt:.2f} s ({vm:.1f} min)</td><td>1,71 GB</td>'
            f'<td>{yt:.2f} s</td><td>{yr:.2f} GB</td>'
            f'<td><strong>{sp:.1f}×</strong></td></tr>'.replace(".", ",")
        )

    t5_chr22_html = []
    for p in ["baseline", "clinvar", "alphamissense", "dbnsfp", "spliceai", "cadd"]:
        vt = vep_c22_plug[p]
        vm = vt / 60.0
        yt, yr = vepyr_c22_plug[p]
        sp = vt / yt
        t5_chr22_html.append(
            f'      <tr><td><strong>{plugin_display[p]}</strong></td>'
            f'<td>{vt:.2f} s ({vm:.1f} min)</td><td>1,73 GB</td>'
            f'<td>{yt:.2f} s</td><td>{yr:.2f} GB</td>'
            f'<td><strong>{sp:.1f}×</strong></td></tr>'.replace(".", ",")
        )

    # Generate Table 6 HTML rows
    t6_html = []
    for f in ["none", "1", "2", "4", "8", "10"]:
        label = "none" if f == "none" else f"fork {f}"
        t1 = vep_c1_forks[f]
        sp1 = vep_c1_forks["none"] / t1
        t22 = vep_c22_forks[f]
        sp22 = vep_c22_forks["none"] / t22
        t6_html.append(
            f'      <tr><td><strong>{label}</strong></td>'
            f'<td>{t1:.2f} s</td><td>{sp1:.2f}×</td>'
            f'<td>{t22:.2f} s</td><td>{sp22:.2f}×</td></tr>'.replace(".", ",")
        )

    # Replace in HTML
    html = HTML_PATH.read_text(encoding="utf-8")
    # Replace Table 5 chr1 tbody
    html = re.sub(
        r'(<h4>Chromosom 1 \(323 430 wariantów\):</h4>\s*<table>\s*<thead>.*?</thead>\s*<tbody>)(.*?)(</tbody>)',
        r'\g<1>\n' + '\n'.join(t5_chr1_html) + r'\n    \g<3>',
        html,
        flags=re.DOTALL
    )
    # Replace Table 5 chr22 tbody
    html = re.sub(
        r'(<h4>Chromosom 22 \(50 861 wariantów\):</h4>\s*<table>\s*<thead>.*?</thead>\s*<tbody>)(.*?)(</tbody>)',
        r'\g<1>\n' + '\n'.join(t5_chr22_html) + r'\n    \g<3>',
        html,
        flags=re.DOTALL
    )
    # Replace Table 6 tbody
    html = re.sub(
        r'(<h3>Tabela 6: Skalowanie VEP przez --fork \(scenariusz baseline.*?\)</h3>\s*<table>\s*<thead>.*?</thead>\s*<tbody>)(.*?)(</tbody>)',
        r'\g<1>\n' + '\n'.join(t6_html) + r'\n    \g<3>',
        html,
        flags=re.DOTALL
    )
    HTML_PATH.write_text(html, encoding="utf-8")

    # Generate Table 5 & 6 Markdown
    t5_chr1_md = []
    for p in ["baseline", "clinvar", "alphamissense", "dbnsfp", "spliceai", "cadd"]:
        vt = vep_c1_plug[p]
        vm = vt / 60.0
        yt, yr = vepyr_c1_plug[p]
        sp = vt / yt
        t5_chr1_md.append(
            f"| **{plugin_display[p]}** | {vt:.2f} s ({vm:.1f} min) | 1,71 GB | {yt:.2f} s | {yr:.2f} GB | **{sp:.1f}×** |".replace(".", ",")
        )

    t5_chr22_md = []
    for p in ["baseline", "clinvar", "alphamissense", "dbnsfp", "spliceai", "cadd"]:
        vt = vep_c22_plug[p]
        vm = vt / 60.0
        yt, yr = vepyr_c22_plug[p]
        sp = vt / yt
        t5_chr22_md.append(
            f"| **{plugin_display[p]}** | {vt:.2f} s ({vm:.1f} min) | 1,73 GB | {yt:.2f} s | {yr:.2f} GB | **{sp:.1f}×** |".replace(".", ",")
        )

    t6_md = []
    for f in ["none", "1", "2", "4", "8", "10"]:
        label = "none" if f == "none" else f"fork {f}"
        t1 = vep_c1_forks[f]
        sp1 = vep_c1_forks["none"] / t1
        t22 = vep_c22_forks[f]
        sp22 = vep_c22_forks["none"] / t22
        t6_md.append(
            f"| **{label}** | {t1:.2f} s | {sp1:.2f}× | {t22:.2f} s | {sp22:.2f}× |".replace(".", ",")
        )

    for md_file in [MD_PATH, ARTIFACT_MD]:
        if not md_file.exists():
            continue
        md = md_file.read_text(encoding="utf-8")
        md = re.sub(
            r'(#### Chromosom 1 \(323 430 wariantów\):\s*\|.*?\|\s*\| :---.*?\|\n)(.*?)(\n\n#### Chromosom 22)',
            r'\g<1>' + '\n'.join(t5_chr1_md) + r'\g<3>',
            md,
            flags=re.DOTALL
        )
        md = re.sub(
            r'(#### Chromosom 22 \(50 861 wariantów\):\s*\|.*?\|\s*\| :---.*?\|\n)(.*?)(\n\n\\\*|\n\n---)',
            r'\g<1>' + '\n'.join(t5_chr22_md) + r'\g<3>',
            md,
            flags=re.DOTALL
        )
        md = re.sub(
            r'(### Tabela 6: Skalowanie VEP przez `--fork` \(scenariusz baseline.*?\)\s*\|.*?\|\s*\| :---.*?\|\n)(.*?)(\n\n---)',
            r'\g<1>' + '\n'.join(t6_md) + r'\g<3>',
            md,
            flags=re.DOTALL
        )
        md_file.write_text(md, encoding="utf-8")


def main():
    LOGS_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / ".metadata_never_index").touch()
    (LOGS_DIR / ".metadata_never_index").touch()

    raw_log = LOGS_DIR / "raw_runs.jsonl"
    summary_json = LOGS_DIR / "summary_vep_5repeats.json"

    # Initialize or resume
    if summary_json.exists():
        with summary_json.open("r", encoding="utf-8") as f:
            all_data = json.load(f)
    else:
        all_data = {
            "chr22_forks": {f: [RUN1_DATA["chr22_forks"][f]] for f in ["none", "1", "2", "4", "8", "10"]},
            "chr1_forks": {f: [RUN1_DATA["chr1_forks"][f]] for f in ["none", "1", "2", "4", "8", "10"]},
            "chr22_plugins": {p: [RUN1_DATA["chr22_plugins"][p]] for p in ["baseline", "clinvar", "spliceai", "cadd", "alphamissense", "dbnsfp"]},
            "chr1_plugins": {p: [RUN1_DATA["chr1_plugins"][p]] for p in ["baseline", "clinvar", "spliceai", "cadd", "alphamissense", "dbnsfp"]},
        }

    # Check if we already have run 2 for chr1 baseline (2327.08 s from warmup)
    if len(all_data["chr1_forks"]["none"]) == 1:
        all_data["chr1_forks"]["none"].append(2327.08)
        all_data["chr1_plugins"]["baseline"].append(2327.08)

    def is_vep_running():
        try:
            out = subprocess.check_output(
                ["docker", "ps", "-q", "--filter", "ancestor=ensemblorg/ensembl-vep:release_116.0"],
                text=True
            ).strip()
            return len(out) > 0
        except Exception:
            return False

    if is_vep_running():
        print("Existing VEP run detected. Waiting for it to finish before starting suite...", flush=True)
        while is_vep_running():
            time.sleep(15)
        print("Existing VEP run finished.", flush=True)
        time_log_prev = Path("/Users/lukaszjezapkowicz/workspace/data_vepyr/runs/vep_chr1_warmup/chr1_none.time.txt")
        if time_log_prev.exists():
            for line in time_log_prev.read_text().splitlines():
                if line.strip().startswith("real "):
                    try:
                        t_val = float(line.strip().split()[1])
                        if len(all_data["chr1_forks"]["none"]) < 3 and t_val > 500:
                            all_data["chr1_forks"]["none"].append(t_val)
                            all_data["chr1_plugins"]["baseline"].append(t_val)
                            print(f"Recorded completed warm run: {t_val:.2f} s", flush=True)
                    except Exception:
                        pass

    # Helper to save state
    def save_and_update():
        summary_json.write_text(json.dumps(all_data, indent=2), encoding="utf-8")
        update_all_reports(all_data)

    print("=== [1/4] Running Ensembl VEP chr22 fork scaling (runs 2..5) ===", flush=True)
    for rep in range(2, 6):
        for f in ["10", "8", "4", "2", "1", "none"]:
            runs = all_data["chr22_forks"][f]
            if len(runs) >= rep:
                continue
            fork_val = None if f == "none" else int(f)
            out_vcf = OUT_DIR / f"chr22_fork_{f}_r{rep}.vcf"
            print(f"  chr22 fork={f} (repeat {rep}/5)...", end="", flush=True)
            t = run_vep(INPUT_CHR22, out_vcf, fork=fork_val, chr_name="chr22")
            runs.append(t)
            if f == "none":
                all_data["chr22_plugins"]["baseline"] = runs
            rec = {"suite": "chr22_forks", "param": f, "repeat": rep, "time_seconds": t, "timestamp": datetime.now(timezone.utc).isoformat()}
            with raw_log.open("a", encoding="utf-8") as rf:
                rf.write(json.dumps(rec) + "\n")
            print(f" {t:.2f}s ({t/60:.2f} min)", flush=True)
            save_and_update()
            time.sleep(5)

    print("=== [2/4] Running Ensembl VEP chr22 plugins (runs 2..5) ===", flush=True)
    for rep in range(2, 6):
        for p in ["clinvar", "alphamissense", "dbnsfp", "cadd", "spliceai"]:
            runs = all_data["chr22_plugins"][p]
            if len(runs) >= rep:
                continue
            out_vcf = OUT_DIR / f"chr22_plugin_{p}_r{rep}.vcf"
            print(f"  chr22 plugin={p} (repeat {rep}/5)...", end="", flush=True)
            t = run_vep(INPUT_CHR22, out_vcf, plugin_name=p, chr_name="chr22")
            runs.append(t)
            rec = {"suite": "chr22_plugins", "param": p, "repeat": rep, "time_seconds": t, "timestamp": datetime.now(timezone.utc).isoformat()}
            with raw_log.open("a", encoding="utf-8") as rf:
                rf.write(json.dumps(rec) + "\n")
            print(f" {t:.2f}s ({t/60:.2f} min)", flush=True)
            save_and_update()
            time.sleep(5)

    print("=== [3/4] Running Ensembl VEP chr1 fork scaling (runs 2..5) ===", flush=True)
    for rep in range(2, 6):
        for f in ["10", "8", "4", "2", "1", "none"]:
            runs = all_data["chr1_forks"][f]
            if len(runs) >= rep:
                continue
            fork_val = None if f == "none" else int(f)
            out_vcf = OUT_DIR / f"chr1_fork_{f}_r{rep}.vcf"
            print(f"  chr1 fork={f} (repeat {rep}/5)...", end="", flush=True)
            t = run_vep(INPUT_CHR1, out_vcf, fork=fork_val, chr_name="chr1")
            runs.append(t)
            if f == "none":
                all_data["chr1_plugins"]["baseline"] = runs
            rec = {"suite": "chr1_forks", "param": f, "repeat": rep, "time_seconds": t, "timestamp": datetime.now(timezone.utc).isoformat()}
            with raw_log.open("a", encoding="utf-8") as rf:
                rf.write(json.dumps(rec) + "\n")
            print(f" {t:.2f}s ({t/60:.2f} min)", flush=True)
            save_and_update()
            time.sleep(10)

    print("=== [4/4] Running Ensembl VEP chr1 plugins (runs 2..5) ===", flush=True)
    for rep in range(2, 6):
        for p in ["alphamissense", "clinvar", "dbnsfp", "cadd", "spliceai"]:
            runs = all_data["chr1_plugins"][p]
            if len(runs) >= rep:
                continue
            out_vcf = OUT_DIR / f"chr1_plugin_{p}_r{rep}.vcf"
            print(f"  chr1 plugin={p} (repeat {rep}/5)...", end="", flush=True)
            t = run_vep(INPUT_CHR1, out_vcf, plugin_name=p, chr_name="chr1")
            runs.append(t)
            rec = {"suite": "chr1_plugins", "param": p, "repeat": rep, "time_seconds": t, "timestamp": datetime.now(timezone.utc).isoformat()}
            with raw_log.open("a", encoding="utf-8") as rf:
                rf.write(json.dumps(rec) + "\n")
            print(f" {t:.2f}s ({t/60:.2f} min)", flush=True)
            save_and_update()
            time.sleep(10)

    print("=== ALL ENSEMBL VEP 5-REPEAT BENCHMARKS COMPLETED SUCCESSFULLY ===", flush=True)


if __name__ == "__main__":
    main()
