#!/usr/bin/env python3
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path
import gzip
import tempfile


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


DEFAULT_PLUGINS = ["clinvar", "spliceai", "cadd", "alphamissense", "dbnsfp"]
PLUGIN_SOURCE_ENV_VARS = {
    "clinvar": "VEPYR_DIFFLY_PLUGIN_CLINVAR_SOURCE",
    "spliceai": "VEPYR_DIFFLY_PLUGIN_SPLICEAI_SOURCE",
    "cadd_snv": "VEPYR_DIFFLY_PLUGIN_CADD_SNV_SOURCE",
    "cadd_indel": "VEPYR_DIFFLY_PLUGIN_CADD_INDEL_SOURCE",
    "alphamissense": "VEPYR_DIFFLY_PLUGIN_ALPHAMISSENSE_SOURCE",
    "dbnsfp": "VEPYR_DIFFLY_PLUGIN_DBNSFP_SOURCE",
}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    from vepyr_diffly.settings import env_int, env_int_or, env_path, load_repo_env

    load_repo_env()
    parser = argparse.ArgumentParser(
        description="Build a vepyr cache, optionally restricted to selected chromosomes."
    )
    parser.add_argument(
        "chromosomes",
        nargs="*",
        help="Optional chromosome list, e.g. `1 Y` or `chr1 chrY`. Omit for the full cache.",
    )
    parser.add_argument(
        "--chromosomes",
        dest="chromosomes_csv",
        help="Optional comma-separated chromosome list, e.g. `1,Y`.",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=env_path("VEPYR_DIFFLY_VEPYR_CACHE_OUTPUT_DIR"),
        help="Output cache root. Defaults to VEPYR_DIFFLY_VEPYR_CACHE_OUTPUT_DIR.",
    )
    parser.add_argument(
        "--vep-cache-dir",
        type=Path,
        default=env_path("VEPYR_DIFFLY_VEP_CACHE_DIR"),
        help="Ensembl VEP cache root. Defaults to VEPYR_DIFFLY_VEP_CACHE_DIR.",
    )
    parser.add_argument(
        "--local-cache",
        type=Path,
        help="Explicit local Ensembl cache directory. Overrides --vep-cache-dir derived path.",
    )
    parser.add_argument(
        "--vepyr-path",
        type=Path,
        default=env_path("VEPYR_DIFFLY_VEPYR_PATH"),
        help="Local vepyr checkout to install editable before building.",
    )
    parser.add_argument(
        "--release",
        type=int,
        default=env_int_or("VEPYR_DIFFLY_VEP_CACHE_VERSION", 115),
        help="Ensembl cache release. Defaults to VEPYR_DIFFLY_VEP_CACHE_VERSION or 115.",
    )
    parser.add_argument("--species", default="homo_sapiens")
    parser.add_argument("--assembly", default="GRCh38")
    parser.add_argument(
        "--cache-flavor",
        choices=("vep", "merged", "refseq"),
        default="vep",
        help="Suffix used for the Ensembl cache source and feature root.",
    )
    parser.add_argument("--partitions", type=int, default=8)
    parser.add_argument("--memory-limit-gb", type=int, default=32)
    parser.add_argument(
        "--plugins",
        default=",".join(DEFAULT_PLUGINS),
        help="Comma-separated plugin list. Use --no-plugins to skip plugin cache generation.",
    )
    parser.add_argument(
        "--no-plugins",
        action="store_true",
        help="Skip plugin cache build.",
    )
    parser.add_argument(
        "--force-plugin-source",
        action="store_true",
        help="Retained for CLI compatibility; local-source builds do not download plugin data.",
    )
    parser.add_argument(
        "--skip-install",
        action="store_true",
        help="Skip `pip install -e <vepyr_path> --no-build-isolation` before import.",
    )
    parser.add_argument(
        "--no-core-fjall",
        action="store_true",
        help="Skip building core `variation.fjall` and `translation_sift.fjall` from parquet.",
    )
    parser.add_argument(
        "--only-plugins",
        action="store_true",
        help="Skip core cache generation entirely and build only the requested plugins.",
    )
    parser.add_argument(
        "--clinvar-source",
        type=Path,
        default=env_path(PLUGIN_SOURCE_ENV_VARS["clinvar"]),
        help=f"Local ClinVar source file. Defaults to {PLUGIN_SOURCE_ENV_VARS['clinvar']}.",
    )
    parser.add_argument(
        "--spliceai-source",
        type=Path,
        default=env_path(PLUGIN_SOURCE_ENV_VARS["spliceai"]),
        help=f"Local SpliceAI source file. Defaults to {PLUGIN_SOURCE_ENV_VARS['spliceai']}.",
    )
    parser.add_argument(
        "--cadd-snv-source",
        type=Path,
        default=env_path(PLUGIN_SOURCE_ENV_VARS["cadd_snv"]),
        help=f"Local CADD SNV source file. Defaults to {PLUGIN_SOURCE_ENV_VARS['cadd_snv']}.",
    )
    parser.add_argument(
        "--cadd-indel-source",
        type=Path,
        default=env_path(PLUGIN_SOURCE_ENV_VARS["cadd_indel"]),
        help=(
            "Local CADD indel source file. "
            f"Defaults to {PLUGIN_SOURCE_ENV_VARS['cadd_indel']}."
        ),
    )
    parser.add_argument(
        "--alphamissense-source",
        type=Path,
        default=env_path(PLUGIN_SOURCE_ENV_VARS["alphamissense"]),
        help=(
            "Local AlphaMissense source file. "
            f"Defaults to {PLUGIN_SOURCE_ENV_VARS['alphamissense']}."
        ),
    )
    parser.add_argument(
        "--dbnsfp-source",
        type=Path,
        default=env_path(PLUGIN_SOURCE_ENV_VARS["dbnsfp"]),
        help=f"Local dbNSFP source file. Defaults to {PLUGIN_SOURCE_ENV_VARS['dbnsfp']}.",
    )
    preview_rows_default = env_int("VEPYR_DIFFLY_PLUGIN_PREVIEW_ROWS")
    parser.add_argument(
        "--preview-rows",
        type=int,
        default=preview_rows_default,
        help="Limit plugin conversion to the first N data rows from each source (headers still preserved).",
    )
    return parser.parse_args(argv)


def resolve_requested_chromosomes(positional: list[str], csv_value: str | None) -> list[str]:
    from vepyr_diffly.chromosomes import parse_chromosome_selection

    tokens = list(positional)
    if csv_value:
        tokens.extend(part.strip() for part in csv_value.split(","))
    if not tokens:
        return []
    normalized, _ = parse_chromosome_selection(",".join(tokens))
    return normalized


def resolve_local_cache_source(
    *,
    vep_cache_dir: Path,
    species: str,
    release: int,
    assembly: str,
    cache_flavor: str,
) -> Path:
    suffix = "" if cache_flavor == "vep" else f"_{cache_flavor}"
    return vep_cache_dir / species / f"{release}_{assembly}{suffix}"


def parse_plugins(raw_value: str) -> list[str]:
    plugins: list[str] = []
    for item in raw_value.split(","):
        value = item.strip().lower()
        if value and value not in plugins:
            plugins.append(value)
    return plugins


def plugin_source_argument_name(plugin_name: str) -> str:
    return plugin_name.lower().replace("-", "").replace("_", "")


def resolve_plugin_source(plugin_name: str, configured_path: Path | None) -> Path:
    normalized_name = plugin_name.lower()
    env_name = PLUGIN_SOURCE_ENV_VARS[normalized_name]
    if configured_path is None:
        raise SystemExit(
            f"missing local source for plugin '{normalized_name}': "
            f"use --{normalized_name}-source or set {env_name}"
        )
    resolved = configured_path.expanduser()
    if not resolved.is_file():
        raise SystemExit(
            f"plugin '{normalized_name}' source file not found: {resolved} "
            f"(from --{normalized_name}-source or {env_name})"
        )
    return resolved.resolve()


def resolve_plugin_sources(args: argparse.Namespace, plugins: list[str]) -> dict[str, object]:
    sources: dict[str, object] = {}
    for plugin in plugins:
        if plugin == "cadd":
            sources[plugin] = {
                "snv": resolve_plugin_source("cadd_snv", args.cadd_snv_source),
                "indel": resolve_plugin_source("cadd_indel", args.cadd_indel_source),
            }
            continue
        sources[plugin] = resolve_plugin_source(
            plugin, getattr(args, f"{plugin_source_argument_name(plugin)}_source")
        )
    return sources


def _open_text_reader(path: Path):
    if path.suffixes and path.suffixes[-1] == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def _open_text_writer(path: Path, suffixes: list[str]):
    if suffixes and suffixes[-1] == ".gz":
        return gzip.open(path, "wt", encoding="utf-8")
    return open(path, "wt", encoding="utf-8")


def _format_bytes(value: int) -> str:
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if value < 1024:
            return f"{value:.1f}{unit}"
        value /= 1024
    return f"{value:.1f}PB"


def _gzip_uncompressed_size(path: Path) -> int:
    if path.suffixes and path.suffixes[-1] == ".gz":
        result = subprocess.run(
            ["gzip", "-l", str(path)],
            capture_output=True,
            text=True,
            check=False,
        )
        if result.stdout:
            lines = result.stdout.strip().splitlines()
            if len(lines) >= 2:
                parts = lines[1].split()
                if len(parts) >= 2 and parts[1].isdigit():
                    return int(parts[1])
    return path.stat().st_size


def _rewrite_chrom_line(line: str, plugin_name: str) -> str:
    plugin_key = plugin_name.split('_')[0]
    if plugin_key in {"cadd", "dbnsfp"} and "\t" in line:
        chrom, rest = line.split("\t", 1)
        if not chrom.lower().startswith("chr"):
            chrom = f"chr{chrom}"
        return f"{chrom}\t{rest}"
    return line


def _create_preview_source(original: Path, plugin_name: str, preview_rows: int) -> tuple[Path, int]:
    suffix = "".join(original.suffixes) or original.suffix
    temp = tempfile.NamedTemporaryFile(prefix=f"{plugin_name}_preview_", suffix=suffix, delete=False)
    temp_path = Path(temp.name)
    temp.close()
    rows_written = 0
    raw_bytes = 0
    with _open_text_reader(original) as reader, _open_text_writer(temp_path, original.suffixes) as writer:
        for line in reader:
            formatted = line
            if not line.startswith("#") and not line.isspace():
                formatted = _rewrite_chrom_line(line, plugin_name)
            writer.write(formatted)
            if line.startswith("#"):
                continue
            rows_written += 1
            raw_bytes += len(formatted.encode("utf-8"))
            if rows_written >= preview_rows:
                break
    return temp_path, raw_bytes


def _prepare_preview_source(
    plugin_name: str,
    source_path: Path,
    preview_rows: int | None,
    temp_paths: list[Path],
) -> tuple[Path, dict[str, int]]:
    if preview_rows is None or preview_rows <= 0:
        return source_path, {
            "compressed_bytes": source_path.stat().st_size,
            "uncompressed_bytes": _gzip_uncompressed_size(source_path),
        }
    preview, raw_bytes = _create_preview_source(source_path, plugin_name, preview_rows)
    temp_paths.append(preview)
    return preview, {
        "compressed_bytes": preview.stat().st_size,
        "uncompressed_bytes": raw_bytes,
    }


def _combine_stats(*stats: dict[str, int]) -> dict[str, int]:
    return {
        "compressed_bytes": sum(stat.get("compressed_bytes", 0) for stat in stats),
        "uncompressed_bytes": sum(stat.get("uncompressed_bytes", 0) for stat in stats),
    }


def _parquet_size(plugin_dir: Path) -> int:
    if not plugin_dir.exists():
        return 0
    total = 0
    for file in plugin_dir.glob("*.parquet"):
        total += file.stat().st_size
    return total


def _log_plugin_stats(
    plugin_name: str,
    stats: dict[str, int],
    cache_dir: Path,
    preview_rows: int | None,
    source_label: str | None = None,
) -> None:
    label = source_label or plugin_name
    parquet_dir = cache_dir / plugin_name
    parquet_bytes = _parquet_size(parquet_dir)
    preview_note = f"preview={preview_rows}" if preview_rows else "full"
    print(
        f"{label} ({preview_note}): .gz={_format_bytes(stats['compressed_bytes'])}, "
        f"raw={_format_bytes(stats['uncompressed_bytes'])}, parquet={_format_bytes(parquet_bytes)}"
    )


def build_plugin_caches(
    vepyr_module,
    cache_dir: Path,
    plugins: list[str],
    plugin_sources: dict[str, object],
    partitions: int,
    chromosomes: list[str] | None,
    preview_rows: int | None,
) -> None:
    temp_preview_paths: list[Path] = []
    try:
        for plugin in plugins:
            if plugin == "cadd":
                cadd_sources = plugin_sources[plugin]
                assert isinstance(cadd_sources, dict)
                snv_source, snv_stats = _prepare_preview_source(
                    "cadd_snv",
                    cadd_sources["snv"],
                    preview_rows,
                    temp_preview_paths,
                )
                indel_source, indel_stats = _prepare_preview_source(
                    "cadd_indel",
                    cadd_sources["indel"],
                    preview_rows,
                    temp_preview_paths,
                )
                print(
                    "building plugin cache for cadd from "
                    f"{snv_source} and {indel_source}"
                )
                vepyr_module.build_plugin(
                    "cadd",
                    {"snv": str(snv_source), "indel": str(indel_source)},
                    str(cache_dir),
                    partitions=partitions,
                    chromosomes=chromosomes,
                )
                combined_stats = _combine_stats(snv_stats, indel_stats)
                _log_plugin_stats(
                    "cadd",
                    combined_stats,
                    cache_dir,
                    preview_rows,
                    source_label="cadd (snv+indel)",
                )
                _log_plugin_stats(
                    "cadd",
                    snv_stats,
                    cache_dir,
                    preview_rows,
                    source_label="cadd_snv",
                )
                _log_plugin_stats(
                    "cadd",
                    indel_stats,
                    cache_dir,
                    preview_rows,
                    source_label="cadd_indel",
                )
                continue

            source_path = plugin_sources[plugin]
            if not isinstance(source_path, Path):
                source_path = Path(str(source_path))
            limited_source, stats = _prepare_preview_source(
                plugin, source_path, preview_rows, temp_preview_paths
            )
            print(f"building plugin cache for {plugin} from {limited_source}")
            vepyr_module.build_plugin(
                plugin,
                str(limited_source),
                str(cache_dir),
                partitions=partitions,
                chromosomes=chromosomes,
            )
            _log_plugin_stats(plugin, stats, cache_dir, preview_rows)
    finally:
        for temp_path in temp_preview_paths:
            try:
                temp_path.unlink()
            except FileNotFoundError:
                pass


def ensure_vepyr_installed(vepyr_path: Path) -> None:
    subprocess.run(
        [
            sys.executable,
            "-m",
            "pip",
            "install",
            "--no-build-isolation",
            "-e",
            str(vepyr_path),
        ],
        check=True,
    )


def ensure_vepyr_import_path(vepyr_path: Path) -> None:
    source_root = vepyr_path / "src"
    if str(source_root) not in sys.path:
        sys.path.insert(0, str(source_root))


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.cache_dir is None:
        raise SystemExit("missing --cache-dir or VEPYR_DIFFLY_VEPYR_CACHE_OUTPUT_DIR")
    if not args.only_plugins and args.local_cache is None and args.vep_cache_dir is None:
        raise SystemExit("missing --local-cache or VEPYR_DIFFLY_VEP_CACHE_DIR")
    if args.vepyr_path is None and not args.skip_install:
        raise SystemExit("missing --vepyr-path or VEPYR_DIFFLY_VEPYR_PATH")

    requested_chromosomes = resolve_requested_chromosomes(args.chromosomes, args.chromosomes_csv)
    plugins = [] if args.no_plugins else parse_plugins(args.plugins)
    plugin_sources = resolve_plugin_sources(args, plugins)
    local_cache = None
    if not args.only_plugins:
        local_cache = args.local_cache or resolve_local_cache_source(
            vep_cache_dir=args.vep_cache_dir,
            species=args.species,
            release=args.release,
            assembly=args.assembly,
            cache_flavor=args.cache_flavor,
        )

    if not args.skip_install:
        ensure_vepyr_installed(args.vepyr_path)
    if args.vepyr_path is not None:
        ensure_vepyr_import_path(args.vepyr_path)

    import vepyr

    if not args.only_plugins:
        print(
            f"building core cache into {args.cache_dir} for "
            f"{requested_chromosomes or ['all chromosomes']}"
        )
        vepyr.build_cache(
            release=args.release,
            cache_dir=str(args.cache_dir),
            species=args.species,
            assembly=args.assembly,
            method=args.cache_flavor,
            local_cache=str(local_cache),
            partitions=args.partitions,
            memory_limit_gb=args.memory_limit_gb,
            chromosomes=requested_chromosomes or None,
        )

        if not args.no_core_fjall:
            print("building core fjall caches")
            vepyr.build_cache_fjall(
                str(local_cache),
                str(args.cache_dir),
                release=args.release,
                assembly=args.assembly,
                method=args.cache_flavor,
                partitions=args.partitions,
                chromosomes=requested_chromosomes or None,
            )
    else:
        print("skipping core cache; generating plugins only")

    if args.no_plugins:
        return 0

    if args.force_plugin_source:
        print("warning: --force-plugin-source is ignored for local-source plugin builds")

    build_plugin_caches(
        vepyr_module=vepyr,
        cache_dir=args.cache_dir,
        plugins=plugins,
        plugin_sources=plugin_sources,
        partitions=args.partitions,
        chromosomes=requested_chromosomes or None,
        preview_rows=args.preview_rows,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
