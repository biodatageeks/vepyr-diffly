#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

from dotenv import load_dotenv


REPO = Path(__file__).resolve().parents[1]
SRC = REPO / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import build_chr_cache as build_script


PLUGINS = ["clinvar", "spliceai", "cadd", "alphamissense", "dbnsfp"]
BASE_COLUMNS = ["chrom", "start", "ref", "alt", "gene_symbol", "consequence_terms"]
PLUGIN_COLUMNS = BASE_COLUMNS + [
    "clnrevstat",
    "ds_dg",
    "raw_score",
    "am_class",
    "polyphen2_hdiv_score",
]

BENCHMARK_MODES = [
    ("vep", [], False),
    ("vep_with_plugins", PLUGINS, False),
    ("vep_fjall", [], True),
    ("vep_with_plugins_fjall", PLUGINS, True),
]

SMOKE_VARIANTS = [
    ("chr1", 66926, "AG", "A"),
    ("chr1", 65420, "C", "A"),
    ("chr1", 10001, "T", "A"),
    ("chr1", 69094, "G", "A"),
    ("chr10", 47057, "C", "A"),
]


def parse_args() -> argparse.Namespace:
    load_dotenv(REPO / ".env")
    parser = argparse.ArgumentParser()
    parser.add_argument("--sizes", default="100,1000,10000,100000")
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=REPO / ".cache/vepyr_cache/115_GRCh38_vep",
    )
    parser.add_argument(
        "--reference-fasta",
        type=Path,
        default=Path(os.environ["VEPYR_DIFFLY_REFERENCE_FASTA"]),
    )
    parser.add_argument(
        "--vepyr-path",
        type=Path,
        default=Path(os.environ.get("VEPYR_DIFFLY_VEPYR_PATH", REPO / "vepyr")),
    )
    parser.add_argument("--skip-install", action="store_true")
    parser.add_argument(
        "--input-mode",
        choices=("synthetic",),
        default="synthetic",
    )
    return parser.parse_args()


def write_synthetic_vcf(path: Path, rows: int) -> None:
    lines = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=chr1>",
        "##contig=<ID=chr10>",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    for idx in range(rows):
        chrom, pos, ref, alt = SMOKE_VARIANTS[idx % len(SMOKE_VARIANTS)]
        lines.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_case(
    input_vcf: Path,
    cache_dir: Path,
    reference_fasta: Path,
    plugins: list[str],
    *,
    use_fjall: bool,
) -> dict[str, object]:
    import vepyr

    columns = PLUGIN_COLUMNS if plugins else BASE_COLUMNS
    started = time.perf_counter()
    lazy = vepyr.annotate(
        str(input_vcf),
        str(cache_dir),
        reference_fasta=str(reference_fasta),
        use_fjall=use_fjall,
        plugins=plugins,
        everything=True,
    )
    df = lazy.select(columns).collect()
    elapsed = time.perf_counter() - started
    return {
        "elapsed_s": round(elapsed, 3),
        "rows": df.height,
        "cols": df.width,
    }


def main() -> int:
    args = parse_args()
    sizes = [int(part) for part in args.sizes.split(",") if part.strip()]

    if not args.skip_install:
        build_script.ensure_vepyr_installed(args.vepyr_path)
    build_script.ensure_vepyr_import_path(args.vepyr_path)

    bench_dir = Path("/tmp/vepyr-annotation-bench")
    bench_dir.mkdir(parents=True, exist_ok=True)

    print(
        json.dumps(
            {
                "cache_dir": str(args.cache_dir),
                "reference_fasta": str(args.reference_fasta),
                "sizes": sizes,
                "input_mode": args.input_mode,
                "plugins": PLUGINS,
                "modes": [mode for mode, _, _ in BENCHMARK_MODES],
                "note": "Benchmark materializes selected columns only.",
            },
            indent=2,
        ),
        flush=True,
    )

    results: list[dict[str, object]] = []
    for size in sizes:
        input_vcf = bench_dir / f"synthetic_{size}.vcf"
        write_synthetic_vcf(input_vcf, size)
        for mode, plugins, use_fjall in BENCHMARK_MODES:
            payload = run_case(
                input_vcf,
                args.cache_dir,
                args.reference_fasta,
                plugins,
                use_fjall=use_fjall,
            )
            result = {
                "variants": size,
                "mode": mode,
                "plugins_enabled": bool(plugins),
                "use_fjall": use_fjall,
                **payload,
            }
            results.append(result)
            print(json.dumps(result), flush=True)

    print("FINAL_RESULTS_JSON=" + json.dumps(results), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
