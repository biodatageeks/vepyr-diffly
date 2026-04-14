#!/usr/bin/env python3
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path
import gzip
import shutil
import tempfile
import re
import time


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
        "--chromosomes",
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
        help="Limit core/plugin conversion to preview-sized source slices (headers still preserved).",
    )
    parser.add_argument(
        "--remove-old-layout-cache",
        action="store_true",
        help=(
            "Remove recognized stale cache artifacts from the previous layout before building "
            "(old parquet/<version>, root-level *.fjall, and root-level plugin directories)."
        ),
    )
    parser.add_argument(
        "--assume-sorted-plugin-input",
        action="store_true",
        help=(
            "Opt-in: skip SQL ORDER BY for single-source plugin builds when the raw plugin input "
            "is already sorted by chrom,pos,ref,alt. This is not applied to CADD."
        ),
    )
    parser.add_argument(
        "--clean-plugin-output",
        action="store_true",
        help=(
            "Remove output directories for the requested plugins in the current cache layout "
            "before building them (e.g. <cache>/<version>/<plugin> and <plugin>.fjall)."
        ),
    )
    return parser.parse_args(argv)


def resolve_requested_chromosomes(csv_value: str | None) -> list[str]:
    from vepyr_diffly.chromosomes import parse_chromosome_selection

    tokens: list[str] = []
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


def resolve_partitioned_cache_dir(
    *,
    cache_dir: Path,
    release: int,
    assembly: str,
    cache_flavor: str,
) -> Path:
    return cache_dir / f"{release}_{assembly}_{cache_flavor}"


def resolve_legacy_partitioned_cache_dir(
    *,
    cache_dir: Path,
    release: int,
    assembly: str,
    cache_flavor: str,
) -> Path:
    return cache_dir / "parquet" / f"{release}_{assembly}_{cache_flavor}"


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


def _natural_file_key(path: Path) -> tuple[int, int, str]:
    name = path.name
    if name == "all_vars.gz":
        return (2, 0, name)
    if name.endswith(".csi"):
        return (3, 0, name)
    match = re.match(r"(\d+)-(\d+)(_reg)?\.gz$", name)
    if match:
        start = int(match.group(1))
        kind = 1 if match.group(3) else 0
        return (kind, start, name)
    return (4, 0, name)


def _copy_core_preview_dir(source_dir: Path, destination_dir: Path, preview_rows: int) -> None:
    destination_dir.mkdir(parents=True, exist_ok=True)
    shard_budget = max(1, (preview_rows + 1_000_000 - 1) // 1_000_000)
    copied_main = 0
    copied_reg = 0
    for path in sorted(source_dir.iterdir(), key=_natural_file_key):
        if path.is_dir():
            continue
        if path.suffix == ".csi":
            continue
        if path.name == "all_vars.gz":
            shutil.copy2(path, destination_dir / path.name)
            continue
        if path.name.endswith("_reg.gz"):
            if copied_reg < shard_budget:
                shutil.copy2(path, destination_dir / path.name)
                copied_reg += 1
            continue
        if path.name.endswith(".gz"):
            if copied_main < shard_budget:
                shutil.copy2(path, destination_dir / path.name)
                copied_main += 1


def create_core_preview_cache(
    local_cache: Path,
    preview_rows: int,
    chromosomes: list[str] | None,
) -> Path:
    preview_root = Path(tempfile.mkdtemp(prefix="vep-core-preview-"))
    for metadata_name in ("info.txt", "chr_synonyms.txt"):
        source = local_cache / metadata_name
        if source.exists():
            shutil.copy2(source, preview_root / metadata_name)

    requested = {chrom.removeprefix("chr") for chrom in chromosomes or []}
    source_dirs = sorted(
        [path for path in local_cache.iterdir() if path.is_dir()],
        key=lambda path: path.name,
    )
    for source_dir in source_dirs:
        if requested and source_dir.name not in requested:
            continue
        _copy_core_preview_dir(source_dir, preview_root / source_dir.name, preview_rows)
    return preview_root


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


def _normalize_preview_chrom(value: str) -> str:
    return value.strip().removeprefix("chr").upper()


def _requested_chromosome_set(chromosomes: list[str] | None) -> set[str]:
    return {_normalize_preview_chrom(chrom) for chrom in chromosomes or [] if chrom.strip()}


def _plugin_tabix_kind(plugin_name: str) -> str | None:
    plugin_key = plugin_name.split("_")[0]
    if plugin_key in {"clinvar", "spliceai"}:
        return "vcf"
    if plugin_key in {"cadd", "dbnsfp"}:
        return "tsv"
    return None


def _tabix_path_for(source_path: Path) -> Path:
    return source_path.with_suffix(source_path.suffix + ".tbi")


def _tabix_regions(chromosomes: list[str]) -> list[str]:
    regions: list[str] = []
    seen: set[str] = set()
    for chrom in chromosomes:
        normalized = chrom.removeprefix("chr")
        for value in (normalized, f"chr{normalized}"):
            if value not in seen:
                regions.append(value)
                seen.add(value)
    return regions


def _read_tsv_header_line(original: Path) -> str:
    with _open_text_reader(original) as reader:
        for line in reader:
            if line.strip() and "\t" in line and not line.startswith("##"):
                return line
    raise RuntimeError(f"failed to read TSV header from {original}")


def _extract_tabix_preview_source(
    original: Path,
    plugin_name: str,
    preview_rows: int | None,
    chromosomes: list[str],
) -> tuple[Path, int]:
    tabix_bin = shutil.which("tabix")
    if tabix_bin is None:
        raise RuntimeError("tabix not found on PATH")
    if not _tabix_path_for(original).exists():
        raise RuntimeError(f"tabix index not found for {original}")

    suffix = "".join(original.suffixes) or original.suffix
    temp = tempfile.NamedTemporaryFile(prefix=f"{plugin_name}_preview_", suffix=suffix, delete=False)
    temp_path = Path(temp.name)
    temp.close()

    raw_bytes = 0
    rows_written = 0
    requested = _requested_chromosome_set(chromosomes)
    command = [tabix_bin]
    tabix_kind = _plugin_tabix_kind(plugin_name)
    if tabix_kind == "vcf":
        command.append("-h")
    command.append(str(original))
    command.extend(_tabix_regions(chromosomes))

    with _open_text_writer(temp_path, original.suffixes) as writer:
        if tabix_kind == "tsv":
            writer.write(_read_tsv_header_line(original))
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
        assert process.stdout is not None
        try:
            for line in process.stdout:
                formatted = line
                if not line.startswith("#") and not line.isspace():
                    formatted = _rewrite_chrom_line(line, plugin_name)
                    chrom = formatted.split("\t", 1)[0]
                    if requested and _normalize_preview_chrom(chrom) not in requested:
                        continue
                writer.write(formatted)
                if line.startswith("#"):
                    continue
                rows_written += 1
                raw_bytes += len(formatted.encode("utf-8"))
                if preview_rows is not None and rows_written >= preview_rows:
                    process.terminate()
                    break
        finally:
            _, stderr = process.communicate()
        if process.returncode not in (0, -15):
            raise RuntimeError(f"tabix failed for {original}: {stderr.strip()}")
    return temp_path, raw_bytes


def _create_preview_source(
    original: Path,
    plugin_name: str,
    preview_rows: int | None,
    chromosomes: list[str] | None = None,
) -> tuple[Path, int]:
    suffix = "".join(original.suffixes) or original.suffix
    temp = tempfile.NamedTemporaryFile(prefix=f"{plugin_name}_preview_", suffix=suffix, delete=False)
    temp_path = Path(temp.name)
    temp.close()
    rows_written = 0
    raw_bytes = 0
    requested = _requested_chromosome_set(chromosomes)
    header_written = False
    is_tsv = _plugin_tabix_kind(plugin_name) == "tsv"
    with _open_text_reader(original) as reader, _open_text_writer(temp_path, original.suffixes) as writer:
        for line in reader:
            if is_tsv and not header_written:
                if line.strip() and "\t" in line and not line.startswith("##"):
                    writer.write(line)
                    header_written = True
                continue
            formatted = line
            if not line.startswith("#") and not line.isspace():
                formatted = _rewrite_chrom_line(line, plugin_name)
                if requested:
                    chrom = formatted.split("\t", 1)[0]
                    if _normalize_preview_chrom(chrom) not in requested:
                        continue
            writer.write(formatted)
            if line.startswith("#"):
                continue
            rows_written += 1
            raw_bytes += len(formatted.encode("utf-8"))
            if preview_rows is not None and rows_written >= preview_rows:
                break
    return temp_path, raw_bytes


def _prepare_preview_source(
    plugin_name: str,
    source_path: Path,
    preview_rows: int | None,
    chromosomes: list[str] | None,
    temp_paths: list[Path],
) -> tuple[Path, dict[str, int]]:
    apply_preview_limit = preview_rows is not None and preview_rows > 0
    if not apply_preview_limit and not chromosomes:
        return source_path, {
            "compressed_bytes": source_path.stat().st_size,
            "uncompressed_bytes": _gzip_uncompressed_size(source_path),
        }
    preview_error: RuntimeError | None = None
    tabix_used = False
    if chromosomes:
        try:
            preview, raw_bytes = _extract_tabix_preview_source(
                source_path,
                plugin_name,
                preview_rows if apply_preview_limit else None,
                chromosomes,
            )
            tabix_used = True
        except RuntimeError as exc:
            preview_error = exc
            preview, raw_bytes = _create_preview_source(
                source_path,
                plugin_name,
                preview_rows if apply_preview_limit else None,
                chromosomes,
            )
    else:
        preview, raw_bytes = _create_preview_source(
            source_path,
            plugin_name,
            preview_rows if apply_preview_limit else None,
            chromosomes,
        )
    temp_paths.append(preview)
    if tabix_used:
        mode = "preview" if apply_preview_limit else "full chromosome-sliced source"
        print(f"tabix source slice for {plugin_name}: using {mode} from {source_path.name}")
    if preview_error is not None:
        print(f"preview fallback for {plugin_name}: {preview_error}")
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


def remove_old_layout_cache(
    *,
    cache_dir: Path,
    release: int,
    assembly: str,
    cache_flavor: str,
    plugins: list[str],
) -> list[Path]:
    removed: list[Path] = []
    targets: list[Path] = [
        resolve_legacy_partitioned_cache_dir(
            cache_dir=cache_dir,
            release=release,
            assembly=assembly,
            cache_flavor=cache_flavor,
        ),
        cache_dir / "variation.fjall",
        cache_dir / "translation_sift.fjall",
    ]
    for plugin in sorted(set(DEFAULT_PLUGINS) | set(plugins)):
        targets.append(cache_dir / plugin)
        targets.append(cache_dir / f"{plugin}.fjall")

    for path in targets:
        if path.is_symlink() or path.is_file():
            path.unlink()
            removed.append(path)
            continue
        if path.is_dir():
            shutil.rmtree(path)
            removed.append(path)
    return removed


def remove_current_plugin_outputs(
    *,
    partitioned_cache_dir: Path,
    plugins: list[str],
) -> list[Path]:
    removed: list[Path] = []
    for plugin in plugins:
        for path in (
            partitioned_cache_dir / plugin,
            partitioned_cache_dir / f"{plugin}.fjall",
        ):
            if path.is_symlink() or path.is_file():
                path.unlink()
                removed.append(path)
                continue
            if path.is_dir():
                shutil.rmtree(path)
                removed.append(path)
    return removed


def _log_plugin_stats(
    plugin_name: str,
    stats: dict[str, int] | None,
    cache_dir: Path,
    preview_rows: int | None,
    preview_prep_seconds: float,
    convert_seconds: float,
    source_label: str | None = None,
) -> None:
    label = source_label or plugin_name
    parquet_dir = cache_dir / plugin_name
    parquet_bytes = _parquet_size(parquet_dir)
    preview_note = f"preview={preview_rows}" if preview_rows else "full"
    total_seconds = preview_prep_seconds + convert_seconds
    if stats is not None:
        print(
            f"{label} ({preview_note}): .gz={_format_bytes(stats['compressed_bytes'])}, "
            f"raw={_format_bytes(stats['uncompressed_bytes'])}, parquet={_format_bytes(parquet_bytes)}, "
            f"preview_prep={preview_prep_seconds:.1f}s, convert={convert_seconds:.1f}s, total={total_seconds:.1f}s"
        )
        return
    print(
        f"{label} ({preview_note}): parquet={_format_bytes(parquet_bytes)}, "
        f"preview_prep={preview_prep_seconds:.1f}s, convert={convert_seconds:.1f}s, total={total_seconds:.1f}s"
    )


def build_plugin_caches(
    vepyr_module,
    cache_dir: Path,
    plugins: list[str],
    plugin_sources: dict[str, object],
    partitions: int,
    chromosomes: list[str] | None,
    preview_rows: int | None,
    assume_sorted_input: bool,
) -> None:
    temp_preview_paths: list[Path] = []
    try:
        for plugin in plugins:
            if plugin == "cadd":
                cadd_sources = plugin_sources[plugin]
                assert isinstance(cadd_sources, dict)
                direct_builder_limit = None
                preview_started_at = time.perf_counter()
                snv_source, snv_stats = _prepare_preview_source(
                    "cadd_snv",
                    cadd_sources["snv"],
                    preview_rows,
                    chromosomes,
                    temp_preview_paths,
                )
                indel_source, indel_stats = _prepare_preview_source(
                    "cadd_indel",
                    cadd_sources["indel"],
                    preview_rows,
                    chromosomes,
                    temp_preview_paths,
                )
                preview_elapsed = time.perf_counter() - preview_started_at
                print(
                    "building plugin cache for cadd from "
                    f"{snv_source} and {indel_source}"
                )
                convert_started_at = time.perf_counter()
                vepyr_module.build_plugin(
                    "cadd",
                    {"snv": str(snv_source), "indel": str(indel_source)},
                    str(cache_dir),
                    partitions=partitions,
                    chromosomes=chromosomes,
                    assume_sorted_input=assume_sorted_input,
                    preview_rows=direct_builder_limit,
                )
                convert_elapsed = time.perf_counter() - convert_started_at
                combined_stats = (
                    _combine_stats(snv_stats, indel_stats)
                    if snv_stats is not None and indel_stats is not None
                    else None
                )
                _log_plugin_stats(
                    "cadd",
                    combined_stats,
                    cache_dir,
                    preview_rows,
                    preview_elapsed,
                    convert_elapsed,
                    source_label="cadd (snv+indel)",
                )
                _log_plugin_stats(
                    "cadd",
                    snv_stats,
                    cache_dir,
                    preview_rows,
                    preview_elapsed,
                    convert_elapsed,
                    source_label="cadd_snv",
                )
                _log_plugin_stats(
                    "cadd",
                    indel_stats,
                    cache_dir,
                    preview_rows,
                    preview_elapsed,
                    convert_elapsed,
                    source_label="cadd_indel",
                )
                continue

            source_path = plugin_sources[plugin]
            if not isinstance(source_path, Path):
                source_path = Path(str(source_path))
            direct_builder_limit = (
                preview_rows
                if preview_rows and not chromosomes and assume_sorted_input
                else None
            )
            if direct_builder_limit is not None:
                limited_source = source_path
                stats = None
                preview_elapsed = 0.0
                print(
                    f"builder row limit for {plugin}: preview_rows={preview_rows} (no temp preview file)"
                )
            else:
                preview_started_at = time.perf_counter()
                limited_source, stats = _prepare_preview_source(
                    plugin, source_path, preview_rows, chromosomes, temp_preview_paths
                )
                preview_elapsed = time.perf_counter() - preview_started_at
            print(f"building plugin cache for {plugin} from {limited_source}")
            convert_started_at = time.perf_counter()
            vepyr_module.build_plugin(
                plugin,
                str(limited_source),
                str(cache_dir),
                partitions=partitions,
                chromosomes=chromosomes,
                assume_sorted_input=assume_sorted_input,
                preview_rows=direct_builder_limit,
            )
            convert_elapsed = time.perf_counter() - convert_started_at
            _log_plugin_stats(
                plugin,
                stats,
                cache_dir,
                preview_rows,
                preview_elapsed,
                convert_elapsed,
            )
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

    requested_chromosomes = resolve_requested_chromosomes(args.chromosomes)
    plugins = [] if args.no_plugins else parse_plugins(args.plugins)
    plugin_sources = resolve_plugin_sources(args, plugins)
    if args.remove_old_layout_cache:
        removed = remove_old_layout_cache(
            cache_dir=args.cache_dir,
            release=args.release,
            assembly=args.assembly,
            cache_flavor=args.cache_flavor,
            plugins=plugins,
        )
        if removed:
            print("removed old-layout cache artifacts:")
            for path in removed:
                print(f"  - {path}")
        else:
            print("no old-layout cache artifacts found")
    partitioned_cache_dir = resolve_partitioned_cache_dir(
        cache_dir=args.cache_dir,
        release=args.release,
        assembly=args.assembly,
        cache_flavor=args.cache_flavor,
    )
    if args.clean_plugin_output and plugins:
        removed = remove_current_plugin_outputs(
            partitioned_cache_dir=partitioned_cache_dir,
            plugins=plugins,
        )
        if removed:
            print("removed current plugin outputs:")
            for path in removed:
                print(f"  - {path}")
        else:
            print("no current plugin outputs found")
    local_cache = None
    core_preview_cache: Path | None = None
    if not args.only_plugins:
        local_cache = args.local_cache or resolve_local_cache_source(
            vep_cache_dir=args.vep_cache_dir,
            species=args.species,
            release=args.release,
            assembly=args.assembly,
            cache_flavor=args.cache_flavor,
        )
        if args.preview_rows is not None and args.preview_rows > 0:
            core_preview_cache = create_core_preview_cache(
                local_cache,
                args.preview_rows,
                requested_chromosomes or None,
            )
            print(
                f"using core preview cache from {core_preview_cache} "
                f"(preview_rows={args.preview_rows}, copied full source shards)"
            )
            local_cache = core_preview_cache

    if not args.skip_install:
        ensure_vepyr_installed(args.vepyr_path)
    if args.vepyr_path is not None:
        ensure_vepyr_import_path(args.vepyr_path)

    import vepyr

    try:
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
            cache_dir=partitioned_cache_dir,
            plugins=plugins,
            plugin_sources=plugin_sources,
            partitions=args.partitions,
            chromosomes=requested_chromosomes or None,
            preview_rows=args.preview_rows,
            assume_sorted_input=args.assume_sorted_plugin_input,
        )
        return 0
    finally:
        if core_preview_cache is not None:
            shutil.rmtree(core_preview_cache, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
