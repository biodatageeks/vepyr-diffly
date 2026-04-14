#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import pyarrow.parquet as pq

import build_chr_cache as build_script


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


@dataclass(frozen=True)
class PluginExpectation:
    key_columns: tuple[str, ...]
    value_columns: tuple[str, ...]


PLUGIN_EXPECTATIONS: dict[str, PluginExpectation] = {
    "clinvar": PluginExpectation(
        key_columns=("chrom", "pos", "ref", "alt"),
        value_columns=("clnsig", "clnrevstat", "clndn", "clnvc", "clnvi"),
    ),
    "spliceai": PluginExpectation(
        key_columns=("chrom", "pos", "ref", "alt"),
        value_columns=("symbol", "ds_ag", "ds_al", "ds_dg", "ds_dl"),
    ),
    "cadd": PluginExpectation(
        key_columns=("chrom", "pos", "ref", "alt"),
        value_columns=("raw_score", "phred_score"),
    ),
    "alphamissense": PluginExpectation(
        key_columns=("chrom", "pos", "ref", "alt"),
        value_columns=("genome", "uniprot_id", "transcript_id", "protein_variant", "am_pathogenicity", "am_class"),
    ),
    "dbnsfp": PluginExpectation(
        key_columns=("chrom", "pos", "ref", "alt"),
        value_columns=("sift4g_score", "polyphen2_hdiv_score", "revel_score", "cadd_phred"),
    ),
}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    from vepyr_diffly.settings import env_path, load_repo_env

    load_repo_env()
    parser = argparse.ArgumentParser(
        description=(
            "Create a mini VCF from existing plugin cache rows, run vepyr.annotate(), "
            "and verify that plugin-specific annotation fields are populated."
        )
    )
    parser.add_argument(
        "--plugins",
        default="clinvar,spliceai,cadd,alphamissense,dbnsfp",
        help="Comma-separated plugin list.",
    )
    parser.add_argument("--release", type=int, default=115)
    parser.add_argument("--assembly", default="GRCh38")
    parser.add_argument(
        "--cache-flavor",
        choices=("vep", "merged", "refseq"),
        default="vep",
    )
    parser.add_argument(
        "--cache-root",
        type=Path,
        default=Path(".cache/vepyr_cache"),
        help="Root cache directory containing <release>_<assembly>_<flavor>/.",
    )
    parser.add_argument(
        "--reference-fasta",
        type=Path,
        default=env_path("VEPYR_DIFFLY_REFERENCE_FASTA"),
        help="Reference FASTA path used by vepyr.annotate().",
    )
    parser.add_argument(
        "--vepyr-path",
        type=Path,
        default=env_path("VEPYR_DIFFLY_VEPYR_PATH"),
        help="Local vepyr checkout to install/import.",
    )
    parser.add_argument(
        "--skip-install",
        action="store_true",
        help="Skip pip install -e for the local vepyr checkout.",
    )
    parser.add_argument(
        "--keep-vcf",
        action="store_true",
        help="Keep the generated mini input VCF in /tmp for inspection.",
    )
    return parser.parse_args(argv)


def _plugin_parquet_files(partitioned_cache_dir: Path, plugin_name: str) -> list[Path]:
    plugin_dir = partitioned_cache_dir / plugin_name
    return sorted(plugin_dir.glob("*.parquet"))


def _read_first_row(path: Path) -> dict[str, object]:
    table = pq.read_table(path)
    if table.num_rows == 0:
        raise SystemExit(f"empty parquet file: {path}")
    row = table.slice(0, 1).to_pylist()[0]
    return row


def _normalize_chrom(value: object) -> str:
    chrom = str(value)
    return chrom if chrom.startswith("chr") else f"chr{chrom}"


def _build_variant_row(plugin_name: str, row: dict[str, object]) -> dict[str, object]:
    return {
        "plugin": plugin_name,
        "chrom": _normalize_chrom(row["chrom"]),
        "pos": int(row["pos"]),
        "ref": str(row["ref"]),
        "alt": str(row["alt"]),
    }


def _write_vcf(path: Path, variants: list[dict[str, object]]) -> None:
    contigs = sorted({variant["chrom"] for variant in variants})
    lines = ["##fileformat=VCFv4.2"]
    lines.extend(f"##contig=<ID={chrom}>" for chrom in contigs)
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    for variant in variants:
        lines.append(
            f"{variant['chrom']}\t{variant['pos']}\t.\t{variant['ref']}\t{variant['alt']}\t.\t.\t."
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _find_matching_row(df, variant: dict[str, object]):
    rows = df.filter(
        (df["chrom"] == variant["chrom"])
        & (df["start"] == variant["pos"])
        & (df["ref"] == variant["ref"])
        & (df["alt"] == variant["alt"])
    )
    if rows.height == 0:
        return None
    return rows.to_dicts()[0]


def _has_non_null_value(row: dict[str, object], fields: tuple[str, ...]) -> bool:
    for field in fields:
        value = row.get(field)
        if value is not None:
            return True
    return False


def _log(message: str) -> None:
    print(message, flush=True)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    plugins = build_script.parse_plugins(args.plugins)
    if args.reference_fasta is None:
        raise SystemExit("missing --reference-fasta or VEPYR_DIFFLY_REFERENCE_FASTA")

    partitioned_cache_dir = build_script.resolve_partitioned_cache_dir(
        cache_dir=args.cache_root.resolve(),
        release=args.release,
        assembly=args.assembly,
        cache_flavor=args.cache_flavor,
    )
    if not partitioned_cache_dir.exists():
        raise SystemExit(f"partitioned cache directory not found: {partitioned_cache_dir}")

    if not args.skip_install:
        if args.vepyr_path is None:
            raise SystemExit("missing --vepyr-path or VEPYR_DIFFLY_VEPYR_PATH")
        build_script.ensure_vepyr_installed(args.vepyr_path)
    if args.vepyr_path is not None:
        build_script.ensure_vepyr_import_path(args.vepyr_path)

    import vepyr

    variants: list[dict[str, object]] = []
    for plugin_name in plugins:
        parquet_files = _plugin_parquet_files(partitioned_cache_dir, plugin_name)
        if not parquet_files:
            raise SystemExit(
                f"missing plugin cache for {plugin_name}: no parquet files in {partitioned_cache_dir / plugin_name}"
            )
        row = _read_first_row(parquet_files[0])
        variant = _build_variant_row(plugin_name, row)
        variants.append(variant)
        _log(
            f"[{plugin_name}] cache variant: {variant['chrom']}:{variant['pos']} {variant['ref']}>{variant['alt']}"
        )

    temp_vcf = Path(tempfile.mkstemp(prefix="plugin-annot-smoke-", suffix=".vcf")[1])
    try:
        _write_vcf(temp_vcf, variants)
        _log(f"[input] VCF: {temp_vcf}")
        _log("[annotate] running vepyr.annotate()")

        lazy = vepyr.annotate(
            str(temp_vcf),
            str(partitioned_cache_dir),
            reference_fasta=str(args.reference_fasta),
            use_fjall=False,
            plugins=plugins,
            everything=True,
        )

        selected_columns = [
            "chrom",
            "start",
            "ref",
            "alt",
        ]
        for plugin_name in plugins:
            selected_columns.extend(PLUGIN_EXPECTATIONS[plugin_name].value_columns)
        selected_columns = list(dict.fromkeys(selected_columns))
        df = lazy.select(selected_columns).collect()

        failures: list[str] = []
        report: list[dict[str, object]] = []
        for variant in variants:
            plugin_name = str(variant["plugin"])
            row = _find_matching_row(df, variant)
            if row is None:
                failures.append(
                    f"{plugin_name}: missing annotated row for {variant['chrom']}:{variant['pos']} {variant['ref']}>{variant['alt']}"
                )
                continue
            fields = PLUGIN_EXPECTATIONS[plugin_name].value_columns
            populated = _has_non_null_value(row, fields)
            report.append(
                {
                    "plugin": plugin_name,
                    "variant": f"{variant['chrom']}:{variant['pos']} {variant['ref']}>{variant['alt']}",
                    "populated_fields": [field for field in fields if row.get(field) is not None],
                }
            )
            _log(
                f"[{plugin_name}] annotate: variant={variant['chrom']}:{variant['pos']} {variant['ref']}>{variant['alt']} "
                f"populated={','.join(field for field in fields if row.get(field) is not None) or 'NONE'}"
            )
            if not populated:
                failures.append(
                    f"{plugin_name}: annotated row present but all plugin fields are NULL for "
                    f"{variant['chrom']}:{variant['pos']} {variant['ref']}>{variant['alt']}"
                )

        _log("[summary]")
        print(json.dumps(report, indent=2), flush=True)

        if failures:
            _log("annotation smoke test FAILED:")
            for failure in failures:
                _log(f"  - {failure}")
            return 1

        _log("annotation smoke test PASSED")
        return 0
    finally:
        if not args.keep_vcf:
            try:
                temp_vcf.unlink()
            except FileNotFoundError:
                pass


if __name__ == "__main__":
    raise SystemExit(main())
