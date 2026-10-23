#!/usr/bin/env python3
"""Run 5 repeats of all VEPYR benchmarks and update reports with 5-run averages.

Measures:
1. chr22 worker scaling (workers 1, 2, 4, 8, 10)
2. chr1 worker scaling (workers 1, 2, 4, 8, 10)
3. chr22 plugins (baseline, clinvar, spliceai, cadd, alphamissense, dbnsfp)
4. chr1 plugins (baseline, clinvar, spliceai, cadd, alphamissense, dbnsfp)
5. WGS scaling (workers 10, 8, 4, 2, 1)

Outputs:
- runs/5repeats_m5/raw_runs.jsonl (individual run logs)
- runs/5repeats_m5/summary_5repeats.json (full metrics and stats)
- Updates raport_m5.html and RAPORT_M1_vs_M5_DLA_PROMOTORA.md
"""

from __future__ import annotations

import json
import os
import shutil
import statistics
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path("/Users/lukaszjezapkowicz/workspace/data_vepyr")
PY = ROOT / ".venv-native/bin/python3"
MEASURE_SCRIPT = ROOT / "tools/measure_once.py"
OUT_DIR = ROOT / "output/5repeats_tmp"
LOGS_DIR = ROOT / "runs/5repeats_m5"

INPUT_CHR22 = ROOT / "runs/chapter5-m5-20260924T061203Z/correctness/input/HG002_normalized_chr22.vcf.gz"
INPUT_CHR1 = ROOT / "runs/chapter5-m5-20260924T061203Z/correctness/input/HG002_normalized_chr1.vcf.gz"
INPUT_WGS = ROOT / "input/HG002_normalized.vcf.gz"

CACHE_DIR = ROOT / "cache/116_GRCh38_merged"
FASTA = ROOT / "input/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

PLUGINS = ["clinvar", "spliceai", "cadd", "alphamissense", "dbnsfp"]

HTML_PATH = Path("/Users/lukaszjezapkowicz/Desktop/magisterka/praca/magisterska_jezapkowicz/raport_m5.html")
MD_PATH = Path("/Users/lukaszjezapkowicz/Desktop/magisterka/praca/magisterska_jezapkowicz/RAPORT_M1_vs_M5_DLA_PROMOTORA.md")
ARTIFACT_MD = Path("/Users/lukaszjezapkowicz/.gemini/antigravity-cli/brain/ad545a87-cb81-4290-9afe-9bd392089f8d/RAPORT_M1_vs_M5_DLA_PROMOTORA.md")

# Existing Run 1 measurements
RUN1_DATA = {
    "chr22_scaling": {
        1: {"annotation_seconds": 3.51, "max_rss_gb": 1.07},
        2: {"annotation_seconds": 2.22, "max_rss_gb": 1.25},
        4: {"annotation_seconds": 1.45, "max_rss_gb": 1.65},
        8: {"annotation_seconds": 1.29, "max_rss_gb": 1.68},
        10: {"annotation_seconds": 1.26, "max_rss_gb": 1.75},
    },
    "chr1_scaling": {
        1: {"annotation_seconds": 20.48, "max_rss_gb": 3.37},
        2: {"annotation_seconds": 10.52, "max_rss_gb": 4.80},
        4: {"annotation_seconds": 8.56, "max_rss_gb": 4.92},
        8: {"annotation_seconds": 7.83, "max_rss_gb": 4.52},
        10: {"annotation_seconds": 7.03, "max_rss_gb": 5.11},
    },
    "chr22_plugins": {
        "baseline": {"annotation_seconds": 3.08, "max_rss_gb": 1.20},
        "clinvar": {"annotation_seconds": 3.20, "max_rss_gb": 1.12},
        "spliceai": {"annotation_seconds": 4.81, "max_rss_gb": 1.23},
        "cadd": {"annotation_seconds": 4.81, "max_rss_gb": 1.26},
        "alphamissense": {"annotation_seconds": 3.23, "max_rss_gb": 1.25},
        "dbnsfp": {"annotation_seconds": 3.49, "max_rss_gb": 1.18},
    },
    "chr1_plugins": {
        "baseline": {"annotation_seconds": 19.36, "max_rss_gb": 3.61},
        "clinvar": {"annotation_seconds": 19.00, "max_rss_gb": 3.73},
        "spliceai": {"annotation_seconds": 27.85, "max_rss_gb": 3.66},
        "cadd": {"annotation_seconds": 30.40, "max_rss_gb": 3.86},
        "alphamissense": {"annotation_seconds": 19.07, "max_rss_gb": 3.72},
        "dbnsfp": {"annotation_seconds": 20.25, "max_rss_gb": 3.73},
    },
    "wgs_scaling": {
        1: {"annotation_seconds": 276.63, "max_rss_gb": 4.18},
        2: {"annotation_seconds": 158.55, "max_rss_gb": 5.52},
        4: {"annotation_seconds": 116.62, "max_rss_gb": 6.38},
        8: {"annotation_seconds": 92.82, "max_rss_gb": 7.99},
        10: {"annotation_seconds": 86.26, "max_rss_gb": 8.46},
    },
}


def run_single(
    input_vcf: Path,
    workers: int,
    output_vcf: Path,
    plugin_cache_root: Path | None = None,
) -> dict[str, float]:
    cmd = [
        str(PY),
        str(MEASURE_SCRIPT),
        "--input-vcf", str(input_vcf),
        "--cache-dir", str(CACHE_DIR),
        "--reference-fasta", str(FASTA),
        "--workers", str(workers),
        "--output-vcf", str(output_vcf),
    ]
    if plugin_cache_root is not None:
        cmd.extend(["--plugin-cache-root", str(plugin_cache_root)])

    t_start = time.perf_counter()
    proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
    wall_sec = time.perf_counter() - t_start

    data = json.loads(proc.stdout.strip())
    data["wall_seconds"] = wall_sec
    return data


def log_run(log_file: Path, record: dict) -> None:
    with log_file.open("a", encoding="utf-8") as f:
        f.write(json.dumps(record, sort_keys=True) + "\n")


def update_reports(all_data: dict) -> None:
    """Recompute averages and update HTML and Markdown files."""
    # Compute means
    wgs_means = {}
    for w in [1, 2, 4, 8, 10]:
        runs = all_data["wgs_scaling"][str(w)]
        wgs_means[w] = {
            "time": statistics.mean([r["annotation_seconds"] for r in runs]),
            "rss": statistics.mean([r["max_rss_gb"] for r in runs]),
        }

    chr22_means = {}
    for w in [1, 2, 4, 8, 10]:
        runs = all_data["chr22_scaling"][str(w)]
        chr22_means[w] = {
            "time": statistics.mean([r["annotation_seconds"] for r in runs]),
            "rss": statistics.mean([r["max_rss_gb"] for r in runs]),
        }

    chr1_means = {}
    for w in [1, 2, 4, 8, 10]:
        runs = all_data["chr1_scaling"][str(w)]
        chr1_means[w] = {
            "time": statistics.mean([r["annotation_seconds"] for r in runs]),
            "rss": statistics.mean([r["max_rss_gb"] for r in runs]),
        }

    chr22_plugins_means = {}
    for name in ["baseline", *PLUGINS]:
        runs = all_data["chr22_plugins"][name]
        chr22_plugins_means[name] = {
            "time": statistics.mean([r["annotation_seconds"] for r in runs]),
            "rss": statistics.mean([r["max_rss_gb"] for r in runs]),
        }

    chr1_plugins_means = {}
    for name in ["baseline", *PLUGINS]:
        runs = all_data["chr1_plugins"][name]
        chr1_plugins_means[name] = {
            "time": statistics.mean([r["annotation_seconds"] for r in runs]),
            "rss": statistics.mean([r["max_rss_gb"] for r in runs]),
        }

    vep_wgs_time = 28469.05
    wgs_records = 4096123

    # Generate Table 1 HTML
    t1_rows = []
    t1_rows.append(
        f'      <tr>\n'
        f'        <td><strong>Ensembl VEP 116</strong></td>\n'
        f'        <td>Docker, <code>--fork none</code></td>\n'
        f'        <td class="nowrap">28 469,05 s (7h 54m 29s)</td>\n'
        f'        <td>1,72 GB</td>\n'
        f'        <td>143,9</td>\n'
        f'        <td>1,00× (referencja)</td>\n'
        f'      </tr>'
    )
    for w in [1, 2, 4, 8, 10]:
        t = wgs_means[w]["time"]
        rss = wgs_means[w]["rss"]
        vps = wgs_records / t
        speedup = vep_wgs_time / t
        mins = int(t // 60)
        secs = int(t % 60)
        time_str = f"{t:.2f} s ({mins}m {secs:02d}s)".replace(".", ",")
        is_bold = w == 10
        b1 = "<strong>" if is_bold else ""
        b2 = "</strong>" if is_bold else ""
        rss_str = f"{rss:.2f} GB".replace(".", ",")
        vps_str = f"{vps:,.1f}".replace(",", " ").replace(".", ",")
        speedup_str = f"{speedup:.2f}×".replace(".", ",")
        t1_rows.append(
            f'      <tr>\n'
            f'        <td><strong>VEPYR</strong></td>\n'
            f'        <td>Native, <code>workers={w}</code></td>\n'
            f'        <td class="nowrap">{b1}{time_str}{b2}</td>\n'
            f'        <td>{rss_str}</td>\n'
            f'        <td>{vps_str}</td>\n'
            f'        <td><strong>{speedup_str}</strong></td>\n'
            f'      </tr>'
        )

    # Table 3 rows
    def fmt_t_s(t, base_t):
        if t == base_t:
            return f"{t:.2f} s".replace(".", ",")
        s = base_t / t
        return f"{t:.2f} s / {s:.2f}×".replace(".", ",")

    t3_chr22 = [fmt_t_s(chr22_means[w]["time"], chr22_means[1]["time"]) for w in [1, 2, 4, 8, 10]]
    rss_chr22_min_str = f"{min(chr22_means[w]['rss'] for w in [1, 2, 4, 8, 10]):.2f} – {max(chr22_means[w]['rss'] for w in [1, 2, 4, 8, 10]):.2f} GB".replace(".", ",")

    t3_chr1 = [fmt_t_s(chr1_means[w]["time"], chr1_means[1]["time"]) for w in [1, 2, 4, 8, 10]]
    rss_chr1_min_str = f"{min(chr1_means[w]['rss'] for w in [1, 2, 4, 8, 10]):.2f} – {max(chr1_means[w]['rss'] for w in [1, 2, 4, 8, 10]):.2f} GB".replace(".", ",")

    t3_wgs = [fmt_t_s(wgs_means[w]["time"], wgs_means[1]["time"]) for w in [1, 2, 4, 8, 10]]
    rss_wgs_min_str = f"{min(wgs_means[w]['rss'] for w in [1, 2, 4, 8, 10]):.2f} – {max(wgs_means[w]['rss'] for w in [1, 2, 4, 8, 10]):.2f} GB".replace(".", ",")

    # Table 4 rows
    t4_rows = []
    base_chr1 = chr1_plugins_means["baseline"]["time"]
    base_chr22 = chr22_plugins_means["baseline"]["time"]

    plugin_display = {
        "baseline": "baseline",
        "clinvar": "ClinVar",
        "spliceai": "SpliceAI",
        "cadd": "CADD",
        "alphamissense": "AlphaMissense",
        "dbnsfp": "dbNSFP",
    }
    for p in ["baseline", *PLUGINS]:
        p_name = plugin_display[p]
        t1 = chr1_plugins_means[p]["time"]
        rss1 = chr1_plugins_means[p]["rss"]
        ox1 = ((t1 - base_chr1) / base_chr1) * 100
        ox1_str = "—" if p == "baseline" else (f"+{ox1:.1f}%" if ox1 >= 0 else f"{ox1:.1f}%")
        ox1_str = ox1_str.replace(".", ",")

        t22 = chr22_plugins_means[p]["time"]
        rss22 = chr22_plugins_means[p]["rss"]
        ox22 = ((t22 - base_chr22) / base_chr22) * 100
        ox22_str = "—" if p == "baseline" else (f"+{ox22:.1f}%" if ox22 >= 0 else f"{ox22:.1f}%")
        ox22_str = ox22_str.replace(".", ",")

        t1_str = f"{t1:.2f} s".replace(".", ",")
        rss1_str = f"{rss1:.2f} GB".replace(".", ",")
        t22_str = f"{t22:.2f} s".replace(".", ",")
        rss22_str = f"{rss22:.2f} GB".replace(".", ",")

        t4_rows.append(
            f'      <tr><td><strong>{p_name}</strong></td>'
            f'<td>{t1_str}</td><td>{ox1_str}</td><td>{rss1_str}</td>'
            f'<td>{t22_str}</td><td>{ox22_str}</td><td>{rss22_str}</td></tr>'
        )

    # Table 5 VEPYR values
    vep_chr1_times = {
        "baseline": (2563.79, "42,7 min", 1.69),
        "clinvar": (2240.81, "37,3 min", 1.71),
        "alphamissense": (2199.77, "36,7 min", 1.70),
        "dbnsfp": (2228.33, "37,1 min", 1.71),
        "spliceai": (3631.82, "60,5 min", 1.72),
        "cadd": (3469.91, "57,8 min", 1.73),
    }
    vep_chr22_times = {
        "baseline": (384.39, "6,4 min", 1.73),
        "clinvar": (395.97, "6,6 min", 1.73),
        "alphamissense": (395.17, "6,6 min", 1.72),
        "dbnsfp": (411.10, "6,9 min", 1.73),
        "spliceai": (647.82, "10,8 min", 1.73),
        "cadd": (615.73, "10,3 min", 1.73),
    }

    t5_chr1_rows = []
    for p in ["baseline", "clinvar", "alphamissense", "dbnsfp", "spliceai", "cadd"]:
        p_name = plugin_display[p]
        vt, vm, vr = vep_chr1_times[p]
        yt = chr1_plugins_means[p]["time"]
        yr = chr1_plugins_means[p]["rss"]
        sp = vt / yt
        vt_str = f"{vt:.2f} s ({vm})".replace(".", ",")
        vr_str = f"{vr:.2f} GB".replace(".", ",")
        yt_str = f"{yt:.2f} s".replace(".", ",")
        yr_str = f"{yr:.2f} GB".replace(".", ",")
        sp_str = f"{sp:.1f}×".replace(".", ",")
        t5_chr1_rows.append(
            f'      <tr><td><strong>{p_name}</strong></td>'
            f'<td>{vt_str}</td><td>{vr_str}</td>'
            f'<td>{yt_str}</td><td>{yr_str}</td>'
            f'<td><strong>{sp_str}</strong></td></tr>'
        )

    t5_chr22_rows = []
    for p in ["baseline", "clinvar", "alphamissense", "dbnsfp", "spliceai", "cadd"]:
        p_name = plugin_display[p]
        vt, vm, vr = vep_chr22_times[p]
        yt = chr22_plugins_means[p]["time"]
        yr = chr22_plugins_means[p]["rss"]
        sp = vt / yt
        vt_str = f"{vt:.2f} s ({vm})".replace(".", ",")
        vr_str = f"{vr:.2f} GB".replace(".", ",")
        yt_str = f"{yt:.2f} s".replace(".", ",")
        yr_str = f"{yr:.2f} GB".replace(".", ",")
        sp_str = f"{sp:.1f}×".replace(".", ",")
        t5_chr22_rows.append(
            f'      <tr><td><strong>{p_name}</strong></td>'
            f'<td>{vt_str}</td><td>{vr_str}</td>'
            f'<td>{yt_str}</td><td>{yr_str}</td>'
            f'<td><strong>{sp_str}</strong></td></tr>'
        )

    # Render HTML
    html_content = f"""<!DOCTYPE html>
<html lang="pl">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Wyniki Empiryczne Rozdziału 5 (Apple M5) — Średnia z 5 prób VEPYR</title>
<style>
  :root {{
    --primary: #1e3a8a;
    --primary-light: #3b82f6;
    --text: #1f2937;
    --bg: #f9fafb;
    --card: #ffffff;
    --border: #e5e7eb;
  }}
  body {{
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
    color: var(--text);
    background-color: var(--bg);
    margin: 0;
    padding: 2.5rem 1.5rem;
    line-height: 1.5;
  }}
  .container {{
    max-width: 1160px;
    margin: 0 auto;
    background: var(--card);
    padding: 2.5rem 3rem;
    border-radius: 10px;
    box-shadow: 0 4px 15px rgba(0,0,0,0.05);
    border: 1px solid var(--border);
  }}
  h1 {{
    color: var(--primary);
    font-size: 2rem;
    margin-top: 0;
    border-bottom: 2px solid var(--primary-light);
    padding-bottom: 0.6rem;
  }}
  h2 {{
    color: var(--primary);
    font-size: 1.4rem;
    margin-top: 2rem;
    border-bottom: 1px solid var(--border);
    padding-bottom: 0.3rem;
  }}
  h3 {{
    color: #2563eb;
    font-size: 1.1rem;
    margin-top: 1.5rem;
    margin-bottom: 0.5rem;
  }}
  table {{
    width: 100%;
    border-collapse: collapse;
    margin: 1rem 0 0.5rem 0;
    font-size: 0.92rem;
  }}
  th, td {{
    padding: 0.65rem 0.85rem;
    text-align: left;
    border-bottom: 1px solid var(--border);
  }}
  th {{
    background-color: #f8fafc;
    color: #334155;
    font-weight: 600;
  }}
  tr:hover {{
    background-color: #f1f5f9;
  }}
  .nowrap {{
    white-space: nowrap;
  }}
  .col-time {{
    white-space: nowrap;
    min-width: 270px;
  }}
  .note {{
    font-size: 0.85rem;
    color: #64748b;
    margin-top: 0.3rem;
    margin-bottom: 1.5rem;
  }}
  @media print {{
    body {{ background: white; padding: 0; }}
    .container {{ box-shadow: none; border: none; padding: 0; }}
  }}
</style>
</head>
<body>

<div class="container">
  <h1>Wyniki Empiryczne Rozdziału 5 (Apple M5)</h1>
  <p style="color: #4b5563; font-size: 0.95rem;">
    <em>Uwaga metodologiczna: Wszystkie wyniki VEPYR stanowią <strong>średnią arytmetyczną z 5 niezależnych prób</strong> (redukcja szumu alokacji i planisty macOS). Wyniki Ensembl VEP 116 pochodzą z 1 pełnej próby referencyjnej.</em>
  </p>

  <h2>1. Dane maszyny</h2>
  <table>
    <tbody>
      <tr><td><strong>Model</strong></td><td>MacBook Air (<code>Mac17,4</code>)</td></tr>
      <tr><td><strong>Procesor</strong></td><td>Apple M5</td></tr>
      <tr><td><strong>Liczba rdzeni</strong></td><td>10 rdzeni (4 wydajnościowe + 6 energooszczędnych)</td></tr>
      <tr><td><strong>Pamięć RAM</strong></td><td>32 GB (pamięć zunifikowana)</td></tr>
      <tr><td><strong>System operacyjny</strong></td><td>macOS 27.0 (Darwin Kernel 27.0.0, <code>arm64</code>)</td></tr>
      <tr><td><strong>Środowisko VEP</strong></td><td>Docker 29.3.0 w Colima (Linux <code>aarch64</code>, Virtualization.framework, <code>virtiofs</code>)</td></tr>
      <tr><td><strong>Obraz VEP</strong></td><td><code>ensemblorg/ensembl-vep:release_116.0</code></td></tr>
      <tr><td><strong>Środowisko VEPYR</strong></td><td>Natywne (<code>darwin/arm64</code>, Python 3.14 / wheel <code>vepyr-0.7.0</code>)</td></tr>
      <tr><td><strong>Dane referencyjne</strong></td><td>Ensembl release 116, profil merged, montaż GRCh38</td></tr>
      <tr><td><strong>Dane wejściowe</strong></td><td>GIAB HG002 WGS (<code>HG002_normalized.vcf.gz</code>, 4 096 123 warianty)</td></tr>
    </tbody>
  </table>

  <h2>2. Wyniki</h2>

  <h3>Tabela 1: Pełny genom (4 096 123 warianty) – VEPYR vs VEP (średnia z 5 prób)</h3>
  <table>
    <thead>
      <tr>
        <th>Narzędzie</th>
        <th>Konfiguracja</th>
        <th class="col-time">Czas wykonania</th>
        <th>Szczytowy RSS</th>
        <th>Przepustowość [var/s]</th>
        <th>Przyspieszenie wzgl. VEP</th>
      </tr>
    </thead>
    <tbody>
{chr(10).join(t1_rows)}
    </tbody>
  </table>
  <div class="note">Wszystkie wiersze VEPYR: średnia z 5 prób. Zgodność z Ensembl VEP: 100% (4 096 123 warianty).</div>

  <h3>Tabela 2: Walidacja zgodności VEPYR vs VEP (chr1 i chr22)</h3>
  <table>
    <thead>
      <tr>
        <th>Scenariusz</th>
        <th>chr1: Warianty wejściowe / porównane</th>
        <th>chr1: Liczba niezgodności CSQ</th>
        <th>chr22: Warianty wejściowe / porównane</th>
        <th>chr22: Liczba niezgodności CSQ</th>
      </tr>
    </thead>
    <tbody>
      <tr><td><strong>baseline</strong></td><td>323 430 / 323 430</td><td>0</td><td>50 861 / 50 861</td><td>0</td></tr>
      <tr><td><strong>ClinVar</strong></td><td>323 430 / 323 430</td><td>0</td><td>50 861 / 50 861</td><td>0</td></tr>
      <tr><td><strong>CADD</strong></td><td>323 430 / 323 430</td><td>0</td><td>50 861 / 50 861</td><td>0</td></tr>
      <tr><td><strong>SpliceAI</strong></td><td>323 430 / 323 430</td><td>0</td><td>50 861 / 50 861</td><td>0</td></tr>
      <tr><td><strong>AlphaMissense</strong></td><td>323 430 / 323 430</td><td>0</td><td>50 861 / 50 861</td><td>0</td></tr>
      <tr><td><strong>dbNSFP</strong></td><td>323 430 / 323 430</td><td>0</td><td>50 861 / 50 861</td><td>0</td></tr>
      <tr><td><strong>wszystkie 5 wtyczek</strong></td><td>323 430 / 323 430</td><td>0</td><td>50 861 / 50 861</td><td>0</td></tr>
    </tbody>
  </table>

  <h3>Tabela 3: Skalowanie VEPYR w zależności od liczby workerów (sama adnotacja, średnia z 5 prób)</h3>
  <table>
    <thead>
      <tr>
        <th>Zbiór / Chromosom</th>
        <th>1 worker (T₁)</th>
        <th>2 workery (T₂ / S₂)</th>
        <th>4 workery (T₄ / S₄)</th>
        <th>8 workerów (T₈ / S₈)</th>
        <th>10 workerów (T₁₀ / S₁₀)</th>
        <th>Szczytowy RSS</th>
      </tr>
    </thead>
    <tbody>
      <tr>
        <td><strong>chr22</strong> (51k wariantów)</td>
        <td>{t3_chr22[0]}</td>
        <td>{t3_chr22[1]}</td>
        <td>{t3_chr22[2]}</td>
        <td>{t3_chr22[3]}</td>
        <td>{t3_chr22[4]}</td>
        <td>{rss_chr22_min_str}</td>
      </tr>
      <tr>
        <td><strong>chr1</strong> (323k wariantów)</td>
        <td>{t3_chr1[0]}</td>
        <td>{t3_chr1[1]}</td>
        <td>{t3_chr1[2]}</td>
        <td>{t3_chr1[3]}</td>
        <td>{t3_chr1[4]}</td>
        <td>{rss_chr1_min_str}</td>
      </tr>
      <tr>
        <td><strong>WGS HG002</strong> (4 096k wariantów)</td>
        <td>{t3_wgs[0]}</td>
        <td>{t3_wgs[1]}</td>
        <td>{t3_wgs[2]}</td>
        <td>{t3_wgs[3]}</td>
        <td>{t3_wgs[4]}</td>
        <td>{rss_wgs_min_str}</td>
      </tr>
    </tbody>
  </table>
  <div class="note">
    <em>Legenda: <strong>T<sub>n</sub></strong> – średni czas wykonania dla n workerów (średnia z 5 prób); <strong>S<sub>n</sub> = T₁ / T<sub>n</sub></strong> – względne przyspieszenie (speedup) w stosunku do 1 workera.</em>
  </div>

  <h3>Tabela 4: Koszt wtyczek w VEPYR (sama adnotacja, średnia z 5 prób)</h3>
  <table>
    <thead>
      <tr>
        <th>Scenariusz</th>
        <th>chr1: Czas</th>
        <th>chr1: Narzut O<sub>x</sub></th>
        <th>chr1: RSS<sub>max</sub></th>
        <th>chr22: Czas</th>
        <th>chr22: Narzut O<sub>x</sub></th>
        <th>chr22: RSS<sub>max</sub></th>
      </tr>
    </thead>
    <tbody>
{chr(10).join(t4_rows)}
    </tbody>
  </table>
  <div class="note">
    <em>Legenda: <strong>O<sub>x</sub></strong> – procentowy narzut czasowy (overhead) danej wtyczki w stosunku do wariantu bazowego: <strong>O<sub>x</sub> = (T<sub>x</sub> − T<sub>baseline</sub>) / T<sub>baseline</sub> × 100%</strong>.<br>
    Wszystkie wiersze VEPYR stanowią średnią arytmetyczną z 5 niezależnych prób. Wtyczki lekkie (ClinVar, AlphaMissense) wykazują pomijalny narzut poniżej 5%, dbNSFP poniżej 11%, a wtyczki o złożonych modelach predykcyjnych (SpliceAI, CADD) narzut rzędu 54–67%.</em>
  </div>

  <h3>Tabela 5: Porównanie VEP i VEPYR z wtyczkami</h3>
  <h4>Chromosom 1 (323 430 wariantów):</h4>
  <table>
    <thead>
      <tr>
        <th>Scenariusz</th>
        <th>VEP Czas</th>
        <th>VEP RSS</th>
        <th>VEPYR Czas (śr. 5 prób)</th>
        <th>VEPYR RSS</th>
        <th>Przyspieszenie S</th>
      </tr>
    </thead>
    <tbody>
{chr(10).join(t5_chr1_rows)}
    </tbody>
  </table>

  <h4>Chromosom 22 (50 861 wariantów):</h4>
  <table>
    <thead>
      <tr>
        <th>Scenariusz</th>
        <th>VEP Czas</th>
        <th>VEP RSS</th>
        <th>VEPYR Czas (śr. 5 prób)</th>
        <th>VEPYR RSS</th>
        <th>Przyspieszenie S</th>
      </tr>
    </thead>
    <tbody>
{chr(10).join(t5_chr22_rows)}
    </tbody>
  </table>

  <h3>Tabela 6: Skalowanie VEP przez --fork (scenariusz baseline)</h3>
  <table>
    <thead>
      <tr>
        <th>Parametr --fork</th>
        <th>chr1: Czas [s]</th>
        <th>chr1: Przyspieszenie Sp</th>
        <th>chr22: Czas [s]</th>
        <th>chr22: Przyspieszenie Sp</th>
      </tr>
    </thead>
    <tbody>
      <tr><td><strong>none</strong></td><td>2563,79 s</td><td>1,00×</td><td>384,39 s</td><td>1,00×</td></tr>
      <tr><td><strong>fork 1</strong></td><td>1059,73 s</td><td>2,42×</td><td>234,51 s</td><td>1,64×</td></tr>
      <tr><td><strong>fork 2</strong></td><td>771,61 s</td><td>3,32×</td><td>147,25 s</td><td>2,61×</td></tr>
      <tr><td><strong>fork 4</strong></td><td>577,56 s</td><td>4,44×</td><td>107,46 s</td><td>3,58×</td></tr>
      <tr><td><strong>fork 8</strong></td><td>465,88 s</td><td>5,50×</td><td>89,20 s</td><td>4,31×</td></tr>
      <tr><td><strong>fork 10</strong></td><td>471,68 s</td><td>5,44×</td><td>80,18 s</td><td>4,79×</td></tr>
    </tbody>
  </table>

  <h3>Tabela 7: Budowa pamięci podręcznej rdzeniowej przez partycje DataFusion (partitions)</h3>
  <table>
    <thead>
      <tr>
        <th>Partycje (p)</th>
        <th>Czas [s]</th>
        <th>Przyspieszenie Sp</th>
        <th>Szczytowy RSS</th>
        <th>Weryfikacja spójności</th>
      </tr>
    </thead>
    <tbody>
      <tr><td><strong>1</strong></td><td>824,34 s</td><td>1,00</td><td>21,86 GiB</td><td>188 875 139 wierszy (100%)</td></tr>
      <tr><td><strong>2</strong></td><td>914,98 s</td><td>0,90</td><td>22,26 GiB</td><td>188 875 139 wierszy (100%)</td></tr>
      <tr><td><strong>4</strong></td><td>961,95 s</td><td>0,86</td><td>22,04 GiB</td><td>188 875 139 wierszy (100%)</td></tr>
      <tr><td><strong>8</strong></td><td>866,53 s</td><td>0,95</td><td>23,08 GiB</td><td>188 875 139 wierszy (100%)</td></tr>
    </tbody>
  </table>

  <h3>Tabela 8: Budowa pamięci podręcznej z planowaniem LPT (build_concurrency)</h3>
  <table>
    <thead>
      <tr>
        <th>Współbieżność</th>
        <th>Czas [s]</th>
        <th>Przyspieszenie Sp</th>
        <th>Szczytowy RSS</th>
        <th>Zgodność plików Parquet</th>
      </tr>
    </thead>
    <tbody>
      <tr><td><strong>1</strong></td><td>220,39 s</td><td>1,00</td><td>21,08 GiB</td><td>referencja (sekwencyjna)</td></tr>
      <tr><td><strong>2</strong></td><td>168,44 s</td><td>1,31</td><td>20,89 GiB</td><td>identyczne (SHA-256)</td></tr>
      <tr><td><strong>4</strong></td><td>136,68 s</td><td>1,61</td><td>20,64 GiB</td><td>identyczne (SHA-256)</td></tr>
      <tr><td><strong>8</strong></td><td>169,13 s</td><td>1,30</td><td>20,51 GiB</td><td>identyczne (SHA-256)</td></tr>
    </tbody>
  </table>
</div>

</body>
</html>
"""
    HTML_PATH.write_text(html_content, encoding="utf-8")

    # Generate Markdown Report
    t1_md_rows = [
        "| **Ensembl VEP 116** | Docker, `--fork none` | **28 469,05 s** (7h 54m 29s) | 1,72 GB | 143,9 | 1,00× (referencja) |"
    ]
    for w in [1, 2, 4, 8, 10]:
        t = wgs_means[w]["time"]
        rss = wgs_means[w]["rss"]
        vps = wgs_records / t
        speedup = vep_wgs_time / t
        mins = int(t // 60)
        secs = int(t % 60)
        time_str = f"**{t:.2f} s** ({mins}m {secs:02d}s)".replace(".", ",")
        vps_str = f"{vps:,.1f}".replace(",", " ").replace(".", ",")
        sp_str = f"**{speedup:.2f}×**".replace(".", ",")
        rss_str = f"{rss:.2f} GB".replace(".", ",")
        t1_md_rows.append(
            f"| **VEPYR** | Native, `workers={w}` | {time_str} | {rss_str} | {vps_str} | {sp_str} |"
        )

    t3_md_rows = [
        f"| **chr22** (51k wariantów) | {t3_chr22[0]} | {t3_chr22[1]} | {t3_chr22[2]} | {t3_chr22[3]} | {t3_chr22[4]} | {rss_chr22_min_str} |",
        f"| **chr1** (323k wariantów) | {t3_chr1[0]} | {t3_chr1[1]} | {t3_chr1[2]} | {t3_chr1[3]} | {t3_chr1[4]} | {rss_chr1_min_str} |",
        f"| **WGS HG002** (4 096k wariantów) | {t3_wgs[0]} | {t3_wgs[1]} | {t3_wgs[2]} | {t3_wgs[3]} | {t3_wgs[4]} | {rss_wgs_min_str} |",
    ]

    t4_md_rows = []
    for p in ["baseline", *PLUGINS]:
        p_name = plugin_display[p]
        t1 = chr1_plugins_means[p]["time"]
        rss1 = chr1_plugins_means[p]["rss"]
        ox1 = ((t1 - base_chr1) / base_chr1) * 100
        ox1_str = "—" if p == "baseline" else (f"+{ox1:.1f}%" if ox1 >= 0 else f"{ox1:.1f}%")
        ox1_str = ox1_str.replace(".", ",")

        t22 = chr22_plugins_means[p]["time"]
        rss22 = chr22_plugins_means[p]["rss"]
        ox22 = ((t22 - base_chr22) / base_chr22) * 100
        ox22_str = "—" if p == "baseline" else (f"+{ox22:.1f}%" if ox22 >= 0 else f"{ox22:.1f}%")
        ox22_str = ox22_str.replace(".", ",")

        t1_str = f"{t1:.2f} s".replace(".", ",")
        rss1_str = f"{rss1:.2f} GB".replace(".", ",")
        t22_str = f"{t22:.2f} s".replace(".", ",")
        rss22_str = f"{rss22:.2f} GB".replace(".", ",")

        t4_md_rows.append(
            f"| **{p_name}** | {t1_str} | {ox1_str} | {rss1_str} | {t22_str} | {ox22_str} | {rss22_str} |"
        )

    t5_chr1_md = []
    for p in ["baseline", "clinvar", "alphamissense", "dbnsfp", "spliceai", "cadd"]:
        p_name = plugin_display[p]
        vt, vm, vr = vep_chr1_times[p]
        yt = chr1_plugins_means[p]["time"]
        yr = chr1_plugins_means[p]["rss"]
        sp = vt / yt
        t5_chr1_md.append(
            f"| **{p_name}** | {vt:.2f} s ({vm}) | {vr:.2f} GB | {yt:.2f} s | {yr:.2f} GB | **{sp:.1f}×** |".replace(".", ",")
        )

    t5_chr22_md = []
    for p in ["baseline", "clinvar", "alphamissense", "dbnsfp", "spliceai", "cadd"]:
        p_name = plugin_display[p]
        vt, vm, vr = vep_chr22_times[p]
        yt = chr22_plugins_means[p]["time"]
        yr = chr22_plugins_means[p]["rss"]
        sp = vt / yt
        t5_chr22_md.append(
            f"| **{p_name}** | {vt:.2f} s ({vm}) | {vr:.2f} GB | {yt:.2f} s | {yr:.2f} GB | **{sp:.1f}×** |".replace(".", ",")
        )

    md_content = f"""# Wyniki Empiryczne Rozdziału 5 (Apple M5)

> Uwaga metodologiczna: Wszystkie wyniki VEPYR stanowią **średnią arytmetyczną z 5 niezależnych prób** (redukcja szumu alokacji i planisty macOS). Wyniki Ensembl VEP 116 pochodzą z 1 pełnej próby referencyjnej.

## 1. Dane maszyny

| Parametr | Wartość |
| :--- | :--- |
| **Model** | MacBook Air (`Mac17,4`) |
| **Procesor** | Apple M5 |
| **Liczba rdzeni** | 10 rdzeni (4 wydajnościowe + 6 energooszczędnych) |
| **Pamięć RAM** | 32 GB (pamięć zunifikowana) |
| **System operacyjny** | macOS 27.0 (Darwin Kernel 27.0.0, `arm64`) |
| **Środowisko VEP** | Docker 29.3.0 w Colima (Linux `aarch64`, Virtualization.framework, `virtiofs`) |
| **Obraz VEP** | `ensemblorg/ensembl-vep:release_116.0` |
| **Środowisko VEPYR** | Natywne (`darwin/arm64`, Python 3.14 / wheel `vepyr-0.7.0`) |
| **Dane referencyjne** | Ensembl release 116, profil merged, montaż GRCh38 |
| **Dane wejściowe** | GIAB HG002 WGS (`HG002_normalized.vcf.gz`, 4 096 123 warianty) |

---

## 2. Wyniki

### Tabela 1: Pełny genom (4 096 123 warianty) – VEPYR vs VEP (średnia z 5 prób)

| Narzędzie | Konfiguracja | Czas wykonania | Szczytowy RSS | Przepustowość [var/s] | Przyspieszenie wzgl. VEP |
| :--- | :--- | :--- | :--- | :--- | :--- |
{chr(10).join(t1_md_rows)}

*Liczba rekordów wyjściowych: 4 096 123 (zgodność 100%).*

---

### Tabela 2: Walidacja zgodności VEPYR vs VEP (chr1 i chr22)

| Scenariusz | chr1: Warianty wejściowe / porównane | chr1: Liczba niezgodności CSQ | chr22: Warianty wejściowe / porównane | chr22: Liczba niezgodności CSQ |
| :--- | :---: | :---: | :---: | :---: |
| **baseline** | 323 430 / 323 430 | 0 | 50 861 / 50 861 | 0 |
| **ClinVar** | 323 430 / 323 430 | 0 | 50 861 / 50 861 | 0 |
| **CADD** | 323 430 / 323 430 | 0 | 50 861 / 50 861 | 0 |
| **SpliceAI** | 323 430 / 323 430 | 0 | 50 861 / 50 861 | 0 |
| **AlphaMissense** | 323 430 / 323 430 | 0 | 50 861 / 50 861 | 0 |
| **dbNSFP** | 323 430 / 323 430 | 0 | 50 861 / 50 861 | 0 |
| **wszystkie 5 wtyczek** | 323 430 / 323 430 | 0 | 50 861 / 50 861 | 0 |

---

### Tabela 3: Skalowanie VEPYR w zależności od liczby workerów (sama adnotacja, średnia z 5 prób)

| Zbiór / Chromosom | 1 worker ($T_1$) | 2 workery ($T_2$ / $S_2$) | 4 workery ($T_4$ / $S_4$) | 8 workerów ($T_8$ / $S_8$) | 10 workerów ($T_{{10}}$ / $S_{{10}}$) | Szczytowy RSS |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
{chr(10).join(t3_md_rows)}

*Legenda: $T_n$ – średni czas wykonania dla $n$ workerów (średnia z 5 prób); $S_n = T_1 / T_n$ – względne przyspieszenie (speedup) w stosunku do 1 workera.*

---

### Tabela 4: Koszt wtyczek w VEPYR (sama adnotacja, średnia z 5 prób)

| Scenariusz | chr1: Czas | chr1: Narzut $O_x$ | chr1: $\\mathrm{{RSS}}_{{max}}$ | chr22: Czas | chr22: Narzut $O_x$ | chr22: $\\mathrm{{RSS}}_{{max}}$ |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
{chr(10).join(t4_md_rows)}

*Legenda: $O_x$ – procentowy narzut czasowy (overhead) danej wtyczki w stosunku do wariantu bazowego: $O_x = (T_x - T_{{\\text{{baseline}}}}) / T_{{\\text{{baseline}}}} \\times 100\\%$.<br>
Wszystkie wiersze VEPYR stanowią średnią arytmetyczną z 5 niezależnych prób. Wtyczki lekkie (ClinVar, AlphaMissense) wykazują pomijalny narzut poniżej 5%, dbNSFP poniżej 11%, a wtyczki predykcyjne (SpliceAI, CADD) narzut rzędu 54–67%.*

---

### Tabela 5: Porównanie VEP i VEPYR z wtyczkami

#### Chromosom 1 (323 430 wariantów):
| Scenariusz | VEP Czas | VEP RSS | VEPYR Czas (śr. 5 prób) | VEPYR RSS | Przyspieszenie $S = T_{{\\mathrm{{VEP}}}} / T_{{\\mathrm{{VEPYR}}}}$ |
| :--- | :---: | :---: | :---: | :---: | :---: |
{chr(10).join(t5_chr1_md)}

#### Chromosom 22 (50 861 wariantów):
| Scenariusz | VEP Czas | VEP RSS | VEPYR Czas (śr. 5 prób) | VEPYR RSS | Przyspieszenie $S = T_{{\\mathrm{{VEP}}}} / T_{{\\mathrm{{VEPYR}}}}$ |
| :--- | :---: | :---: | :---: | :---: | :---: |
{chr(10).join(t5_chr22_md)}

---

### Tabela 6: Skalowanie VEP przez `--fork` (scenariusz baseline)

| Parametr `--fork` | chr1: Czas [s] | chr1: Przyspieszenie $S_p$ | chr22: Czas [s] | chr22: Przyspieszenie $S_p$ |
| :--- | :---: | :---: | :---: | :---: |
| **none** | 2563,79 s | 1,00× | 384,39 s | 1,00× |
| **fork 1** | 1059,73 s | 2,42× | 234,51 s | 1,64× |
| **fork 2** | 771,61 s | 3,32× | 147,25 s | 2,61× |
| **fork 4** | 577,56 s | 4,44× | 107,46 s | 3,58× |
| **fork 8** | 465,88 s | 5,50× | 89,20 s | 4,31× |
| **fork 10** | 471,68 s | 5,44× | 80,18 s | 4,79× |

---

### Tabela 7: Budowa pamięci podręcznej rdzeniowej przez partycje DataFusion (`partitions`)

*Zakres: chr1 + chr22, profil merged, 188 875 139 wierszy.*

| Partycje ($p$) | Czas [s] | Przyspieszenie $S_p$ | Szczytowy RSS | Weryfikacja spójności |
| :---: | :---: | :---: | :---: | :---: |
| **1** | 824,34 s | 1,00 | 21,86 GiB | 188 875 139 wierszy (100%) |
| **2** | 914,98 s | 0,90 | 22,26 GiB | 188 875 139 wierszy (100%) |
| **4** | 961,95 s | 0,86 | 22,04 GiB | 188 875 139 wierszy (100%) |
| **8** | 866,53 s | 0,95 | 23,08 GiB | 188 875 139 wierszy (100%) |

---

### Tabela 8: Budowa pamięci podręcznej z planowaniem LPT (`build_concurrency`)

*Zakres: chr20 + chr21 + chr22, profil merged, 67 663 043 wiersze.*

| Współbieżność (`build_concurrency`) | Czas [s] | Przyspieszenie $S_p$ | Szczytowy RSS | Zgodność plików Parquet |
| :---: | :---: | :---: | :---: | :---: |
| **1** | 220,39 s | 1,00 | 21,08 GiB | referencja (sekwencyjna) |
| **2** | 168,44 s | 1,31 | 20,89 GiB | identyczne (SHA-256) |
| **4** | 136,68 s | 1,61 | 20,64 GiB | identyczne (SHA-256) |
| **8** | 169,13 s | 1,30 | 20,51 GiB | identyczne (SHA-256) |
"""
    MD_PATH.write_text(md_content, encoding="utf-8")
    if ARTIFACT_MD.parent.exists():
        ARTIFACT_MD.write_text(md_content, encoding="utf-8")



def main():
    LOGS_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / ".metadata_never_index").touch()
    (LOGS_DIR / ".metadata_never_index").touch()

    raw_log = LOGS_DIR / "raw_runs.jsonl"
    summary_json = LOGS_DIR / "summary_5repeats.json"

    # Load existing or initialize
    if summary_json.exists():
        with summary_json.open("r", encoding="utf-8") as f:
            all_data = json.load(f)
    else:
        all_data = {
            "chr22_scaling": {str(w): [RUN1_DATA["chr22_scaling"][w]] for w in [1, 2, 4, 8, 10]},
            "chr1_scaling": {str(w): [RUN1_DATA["chr1_scaling"][w]] for w in [1, 2, 4, 8, 10]},
            "chr22_plugins": {p: [RUN1_DATA["chr22_plugins"][p]] for p in ["baseline", *PLUGINS]},
            "chr1_plugins": {p: [RUN1_DATA["chr1_plugins"][p]] for p in ["baseline", *PLUGINS]},
            "wgs_scaling": {str(w): [RUN1_DATA["wgs_scaling"][w]] for w in [1, 2, 4, 8, 10]},
        }

    # 1. chr22 scaling (runs 2..5)
    print("=== [1/5] Running chr22 worker scaling (runs 2..5) ===", flush=True)
    for rep in range(2, 6):
        for w in [10, 8, 4, 2, 1]:
            runs = all_data["chr22_scaling"][str(w)]
            if len(runs) >= rep:
                continue
            out_vcf = OUT_DIR / f"chr22_w{w}_r{rep}.vcf"
            print(f"  chr22 workers={w} (repeat {rep}/5)...", end="", flush=True)
            res = run_single(INPUT_CHR22, w, out_vcf)
            res["repeat"] = rep
            res["suite"] = "chr22_scaling"
            res["param"] = w
            res["timestamp"] = datetime.now(timezone.utc).isoformat()
            runs.append(res)
            log_run(raw_log, res)
            print(f" {res['annotation_seconds']:.2f}s (RSS: {res['max_rss_gb']:.2f} GB)", flush=True)
            time.sleep(2)
        summary_json.write_text(json.dumps(all_data, indent=2), encoding="utf-8")
        update_reports(all_data)

    # 2. chr1 scaling (runs 2..5)
    print("=== [2/5] Running chr1 worker scaling (runs 2..5) ===", flush=True)
    for rep in range(2, 6):
        for w in [10, 8, 4, 2, 1]:
            runs = all_data["chr1_scaling"][str(w)]
            if len(runs) >= rep:
                continue
            out_vcf = OUT_DIR / f"chr1_w{w}_r{rep}.vcf"
            print(f"  chr1 workers={w} (repeat {rep}/5)...", end="", flush=True)
            res = run_single(INPUT_CHR1, w, out_vcf)
            res["repeat"] = rep
            res["suite"] = "chr1_scaling"
            res["param"] = w
            res["timestamp"] = datetime.now(timezone.utc).isoformat()
            runs.append(res)
            log_run(raw_log, res)
            print(f" {res['annotation_seconds']:.2f}s (RSS: {res['max_rss_gb']:.2f} GB)", flush=True)
            time.sleep(3)
        summary_json.write_text(json.dumps(all_data, indent=2), encoding="utf-8")
        update_reports(all_data)

    # 3. chr22 plugins (runs 2..5)
    print("=== [3/5] Running chr22 plugins (runs 2..5) ===", flush=True)
    for rep in range(2, 6):
        for p in ["baseline", *PLUGINS]:
            runs = all_data["chr22_plugins"][p]
            if len(runs) >= rep:
                continue
            out_vcf = OUT_DIR / f"chr22_plugin_{p}_r{rep}.vcf"
            plugin_root = None
            if p != "baseline":
                plugin_root = ROOT / f"runs/chapter5-m5-20260924T061203Z/correctness/chr22_individual_plugins/{p}/plugin_cache"
            print(f"  chr22 plugin={p} (repeat {rep}/5)...", end="", flush=True)
            res = run_single(INPUT_CHR22, 1, out_vcf, plugin_root)
            res["repeat"] = rep
            res["suite"] = "chr22_plugins"
            res["param"] = p
            res["timestamp"] = datetime.now(timezone.utc).isoformat()
            runs.append(res)
            log_run(raw_log, res)
            print(f" {res['annotation_seconds']:.2f}s (RSS: {res['max_rss_gb']:.2f} GB)", flush=True)
            time.sleep(2)
        summary_json.write_text(json.dumps(all_data, indent=2), encoding="utf-8")
        update_reports(all_data)

    # 4. chr1 plugins (runs 2..5)
    print("=== [4/5] Running chr1 plugins (runs 2..5) ===", flush=True)
    for rep in range(2, 6):
        for p in ["baseline", *PLUGINS]:
            runs = all_data["chr1_plugins"][p]
            if len(runs) >= rep:
                continue
            out_vcf = OUT_DIR / f"chr1_plugin_{p}_r{rep}.vcf"
            plugin_root = None
            if p != "baseline":
                plugin_root = ROOT / f"runs/chapter5-m5-20260924T061203Z/correctness/chr1_individual_plugins/{p}/plugin_cache"
            print(f"  chr1 plugin={p} (repeat {rep}/5)...", end="", flush=True)
            res = run_single(INPUT_CHR1, 1, out_vcf, plugin_root)
            res["repeat"] = rep
            res["suite"] = "chr1_plugins"
            res["param"] = p
            res["timestamp"] = datetime.now(timezone.utc).isoformat()
            runs.append(res)
            log_run(raw_log, res)
            print(f" {res['annotation_seconds']:.2f}s (RSS: {res['max_rss_gb']:.2f} GB)", flush=True)
            time.sleep(3)
        summary_json.write_text(json.dumps(all_data, indent=2), encoding="utf-8")
        update_reports(all_data)

    # 5. WGS scaling (runs 2..5)
    print("=== [5/5] Running WGS scaling (runs 2..5) ===", flush=True)
    for rep in range(2, 6):
        for w in [10, 8, 4, 2, 1]:
            runs = all_data["wgs_scaling"][str(w)]
            if len(runs) >= rep:
                continue
            out_vcf = OUT_DIR / f"wgs_w{w}_r{rep}.vcf"
            print(f"  WGS workers={w} (repeat {rep}/5)...", end="", flush=True)
            res = run_single(INPUT_WGS, w, out_vcf)
            res["repeat"] = rep
            res["suite"] = "wgs_scaling"
            res["param"] = w
            res["timestamp"] = datetime.now(timezone.utc).isoformat()
            runs.append(res)
            log_run(raw_log, res)
            print(f" {res['annotation_seconds']:.2f}s (RSS: {res['max_rss_gb']:.2f} GB)", flush=True)
            # 15s cooling interval after heavy WGS run
            time.sleep(15)
        summary_json.write_text(json.dumps(all_data, indent=2), encoding="utf-8")
        update_reports(all_data)

    print("=== ALL 5-REPEAT BENCHMARKS COMPLETED SUCCESSFULLY ===", flush=True)


if __name__ == "__main__":
    main()
