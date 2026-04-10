#!/usr/bin/env python3
from __future__ import annotations

import argparse
import shutil
import sys
import tempfile
from collections.abc import Iterable
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

import build_chr_cache as build_script


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


EXPECTED_SCHEMAS: dict[str, list[tuple[str, pa.DataType]]] = {
    "clinvar": [
        ("chrom", pa.string()),
        ("pos", pa.uint32()),
        ("ref", pa.string()),
        ("alt", pa.string()),
        ("clnsig", pa.string()),
        ("clnrevstat", pa.string()),
        ("clndn", pa.string()),
        ("clnvc", pa.string()),
        ("clnvi", pa.string()),
        ("af_esp", pa.float32()),
        ("af_exac", pa.float32()),
        ("af_tgp", pa.float32()),
    ],
    "spliceai": [
        ("chrom", pa.string()),
        ("pos", pa.uint32()),
        ("ref", pa.string()),
        ("alt", pa.string()),
        ("symbol", pa.string()),
        ("ds_ag", pa.float32()),
        ("ds_al", pa.float32()),
        ("ds_dg", pa.float32()),
        ("ds_dl", pa.float32()),
        ("dp_ag", pa.int32()),
        ("dp_al", pa.int32()),
        ("dp_dg", pa.int32()),
        ("dp_dl", pa.int32()),
    ],
    "cadd": [
        ("chrom", pa.string()),
        ("pos", pa.uint32()),
        ("ref", pa.string()),
        ("alt", pa.string()),
        ("raw_score", pa.float32()),
        ("phred_score", pa.float32()),
    ],
    "alphamissense": [
        ("chrom", pa.string()),
        ("pos", pa.uint32()),
        ("ref", pa.string()),
        ("alt", pa.string()),
        ("genome", pa.string()),
        ("uniprot_id", pa.string()),
        ("transcript_id", pa.string()),
        ("protein_variant", pa.string()),
        ("am_pathogenicity", pa.float32()),
        ("am_class", pa.string()),
    ],
    "dbnsfp": [
        ("chrom", pa.string()),
        ("pos", pa.uint32()),
        ("ref", pa.string()),
        ("alt", pa.string()),
        ("sift4g_score", pa.string()),
        ("sift4g_pred", pa.string()),
        ("polyphen2_hdiv_score", pa.string()),
        ("polyphen2_hvar_score", pa.string()),
        ("lrt_score", pa.string()),
        ("lrt_pred", pa.string()),
        ("mutationtaster_score", pa.string()),
        ("mutationtaster_pred", pa.string()),
        ("fathmm_score", pa.string()),
        ("fathmm_pred", pa.string()),
        ("provean_score", pa.string()),
        ("provean_pred", pa.string()),
        ("vest4_score", pa.string()),
        ("metasvm_score", pa.string()),
        ("metasvm_pred", pa.string()),
        ("metalr_score", pa.string()),
        ("metalr_pred", pa.string()),
        ("revel_score", pa.float32()),
        ("gerp_rs", pa.float32()),
        ("phylop100way", pa.float32()),
        ("phylop30way", pa.float32()),
        ("phastcons100way", pa.float32()),
        ("phastcons30way", pa.float32()),
        ("siphy_29way", pa.float32()),
        ("cadd_raw", pa.float32()),
        ("cadd_phred", pa.float32()),
    ],
}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    from vepyr_diffly.settings import env_int, env_path, load_repo_env

    load_repo_env()
    parser = argparse.ArgumentParser(
        description="Run a plugin round-trip check: convert -> read back -> verify row counts and schema."
    )
    parser.add_argument(
        "--plugins",
        default="clinvar,spliceai,cadd,alphamissense",
        help="Comma-separated plugin list.",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        help="Optional output cache directory. Defaults to a temporary directory.",
    )
    parser.add_argument(
        "--keep-cache",
        action="store_true",
        help="Keep the output cache directory after the run.",
    )
    parser.add_argument(
        "--preview-rows",
        type=int,
        default=env_int("VEPYR_DIFFLY_PLUGIN_PREVIEW_ROWS") or 1000,
        help="Limit conversion to the first N data rows from each source.",
    )
    parser.add_argument("--partitions", type=int, default=8)
    parser.add_argument(
        "--vepyr-path",
        type=Path,
        default=env_path("VEPYR_DIFFLY_VEPYR_PATH"),
        help="Local vepyr checkout to install editable before building.",
    )
    parser.add_argument(
        "--skip-install",
        action="store_true",
        help="Skip `pip install -e <vepyr_path> --no-build-isolation` before import.",
    )
    parser.add_argument(
        "--clinvar-source",
        type=Path,
        default=env_path(build_script.PLUGIN_SOURCE_ENV_VARS["clinvar"]),
    )
    parser.add_argument(
        "--spliceai-source",
        type=Path,
        default=env_path(build_script.PLUGIN_SOURCE_ENV_VARS["spliceai"]),
    )
    parser.add_argument(
        "--cadd-snv-source",
        type=Path,
        default=env_path(build_script.PLUGIN_SOURCE_ENV_VARS["cadd_snv"]),
    )
    parser.add_argument(
        "--cadd-indel-source",
        type=Path,
        default=env_path(build_script.PLUGIN_SOURCE_ENV_VARS["cadd_indel"]),
    )
    parser.add_argument(
        "--alphamissense-source",
        type=Path,
        default=env_path(build_script.PLUGIN_SOURCE_ENV_VARS["alphamissense"]),
    )
    parser.add_argument(
        "--dbnsfp-source",
        type=Path,
        default=env_path(build_script.PLUGIN_SOURCE_ENV_VARS["dbnsfp"]),
    )
    return parser.parse_args(argv)


def _count_table_rows(path: Path) -> int:
    parquet = pq.ParquetFile(path)
    return parquet.metadata.num_rows


def _count_rows_in_parquet_dir(plugin_dir: Path) -> tuple[int, list[Path]]:
    files = sorted(plugin_dir.glob("*.parquet"))
    total_rows = sum(_count_table_rows(path) for path in files)
    return total_rows, files


def _validate_schema(plugin_name: str, parquet_files: Iterable[Path]) -> list[str]:
    expected = EXPECTED_SCHEMAS[plugin_name]
    expected_names = [name for name, _ in expected]
    errors: list[str] = []
    for path in parquet_files:
        schema = pq.ParquetFile(path).schema_arrow
        actual_names = schema.names
        if actual_names != expected_names:
            errors.append(
                f"{path.name}: expected columns {expected_names}, got {actual_names}"
            )
            continue
        for name, expected_type in expected:
            actual_type = schema.field(name).type
            if actual_type != expected_type:
                errors.append(
                    f"{path.name}:{name} expected {expected_type} got {actual_type}"
                )
    return errors


def _count_data_rows(path: Path) -> int:
    rows = 0
    with build_script._open_text_reader(path) as reader:
        for line in reader:
            if line.startswith("#") or line.isspace():
                continue
            rows += 1
    return rows


def _count_vcf_allele_rows(path: Path) -> int:
    rows = 0
    with build_script._open_text_reader(path) as reader:
        for line in reader:
            if line.startswith("#") or line.isspace():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            alts = [alt for alt in parts[4].split(",") if alt]
            rows += max(1, len(alts))
    return rows


def _expected_source_rows(plugin_name: str, source: object) -> int:
    if plugin_name == "cadd":
        assert isinstance(source, dict)
        return _count_data_rows(source["snv"]) + _count_data_rows(source["indel"])
    if plugin_name in {"clinvar", "spliceai"}:
        assert isinstance(source, Path)
        return _count_vcf_allele_rows(source)
    assert isinstance(source, Path)
    return _count_data_rows(source)


def _print_plugin_header(plugin_name: str, preview_rows: int) -> None:
    print(f"[{plugin_name}] round-trip start (preview_rows={preview_rows})")


def _print_plugin_step(plugin_name: str, step: str, detail: str) -> None:
    print(f"[{plugin_name}] {step}: {detail}")


def _materialize_plugin(
    *,
    plugin_name: str,
    source: object,
    preview_rows: int,
    cache_dir: Path,
    partitions: int,
    vepyr_module,
    temp_paths: list[Path],
) -> tuple[int, Path]:
    if plugin_name == "cadd":
        assert isinstance(source, dict)
        snv_preview, _ = build_script._prepare_preview_source(
            "cadd_snv",
            source["snv"],
            preview_rows,
            temp_paths,
        )
        indel_preview, _ = build_script._prepare_preview_source(
            "cadd_indel",
            source["indel"],
            preview_rows,
            temp_paths,
        )
        expected_rows = _expected_source_rows(
            plugin_name,
            {"snv": snv_preview, "indel": indel_preview},
        )
        _print_plugin_step(
            "cadd",
            "convert",
            f"source_snv={snv_preview} source_indel={indel_preview}",
        )
        vepyr_module.build_plugin(
            "cadd",
            {"snv": str(snv_preview), "indel": str(indel_preview)},
            str(cache_dir),
            partitions=partitions,
        )
        return expected_rows, cache_dir / "cadd"

    assert isinstance(source, Path)
    preview_source, _ = build_script._prepare_preview_source(
        plugin_name,
        source,
        preview_rows,
        temp_paths,
    )
    expected_rows = _expected_source_rows(plugin_name, preview_source)
    _print_plugin_step(plugin_name, "convert", f"source={preview_source}")
    vepyr_module.build_plugin(
        plugin_name,
        str(preview_source),
        str(cache_dir),
        partitions=partitions,
    )
    return expected_rows, cache_dir / plugin_name


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    plugins = build_script.parse_plugins(args.plugins)
    plugin_sources = build_script.resolve_plugin_sources(args, plugins)

    temp_cache_dir: Path | None = None
    if args.cache_dir is None:
        temp_cache_dir = Path(tempfile.mkdtemp(prefix="plugin-round-trip-"))
        cache_dir = temp_cache_dir
    else:
        cache_dir = args.cache_dir.resolve()
        cache_dir.mkdir(parents=True, exist_ok=True)

    if not args.skip_install:
        if args.vepyr_path is None:
            raise SystemExit("missing --vepyr-path or VEPYR_DIFFLY_VEPYR_PATH")
        build_script.ensure_vepyr_installed(args.vepyr_path)
    if args.vepyr_path is not None:
        build_script.ensure_vepyr_import_path(args.vepyr_path)

    import vepyr

    temp_paths: list[Path] = []
    failures: list[str] = []
    try:
        for plugin_name in plugins:
            _print_plugin_header(plugin_name, args.preview_rows)
            expected_rows, parquet_dir = _materialize_plugin(
                plugin_name=plugin_name,
                source=plugin_sources[plugin_name],
                preview_rows=args.preview_rows,
                cache_dir=cache_dir,
                partitions=args.partitions,
                vepyr_module=vepyr,
                temp_paths=temp_paths,
            )
            actual_rows, parquet_files = _count_rows_in_parquet_dir(parquet_dir)
            if not parquet_files:
                failures.append(f"{plugin_name}: no parquet files written to {parquet_dir}")
                _print_plugin_step(plugin_name, "read-back", f"parquet_dir={parquet_dir} files=0")
                _print_plugin_step(plugin_name, "round-trip", "FAILED")
                continue
            _print_plugin_step(
                plugin_name,
                "read-back",
                f"parquet_dir={parquet_dir} files={len(parquet_files)}",
            )
            row_ok = actual_rows == expected_rows
            _print_plugin_step(
                plugin_name,
                "verify rows",
                f"expected={expected_rows} actual={actual_rows} result={'OK' if row_ok else 'FAILED'}",
            )
            if not row_ok:
                failures.append(
                    f"{plugin_name}: row count mismatch, expected {expected_rows}, got {actual_rows}"
                )

            schema_errors = _validate_schema(plugin_name, parquet_files)
            if schema_errors:
                _print_plugin_step(
                    plugin_name,
                    "verify schema",
                    f"FAILED ({len(schema_errors)} issue(s))",
                )
                failures.extend(f"{plugin_name}: {error}" for error in schema_errors)
            else:
                _print_plugin_step(plugin_name, "verify schema", "OK")

            _print_plugin_step(
                plugin_name,
                "round-trip",
                "PASSED" if row_ok and not schema_errors else "FAILED",
            )
    finally:
        for temp_path in temp_paths:
            try:
                temp_path.unlink()
            except FileNotFoundError:
                pass
        if temp_cache_dir is not None and not args.keep_cache:
            shutil.rmtree(temp_cache_dir, ignore_errors=True)

    if failures:
        print("round-trip check FAILED:")
        for failure in failures:
            print(f"  - {failure}")
        return 1

    print("round-trip check PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
