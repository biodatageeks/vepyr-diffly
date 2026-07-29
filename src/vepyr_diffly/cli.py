from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
from time import perf_counter
import subprocess
import sys
import tempfile

from .presets import get_preset, load_presets
from .plugins import compare_plugin_field_aliases, compare_plugin_fields, parse_plugin_list
from .settings import env_bool, env_int, env_path, env_str, load_repo_env

try:
    from rich.console import Console
    from rich.table import Table
except ModuleNotFoundError:  # pragma: no cover

    class Console:  # type: ignore[override]
        def print(self, value: object = "") -> None:
            print(value)

    class Table:  # pragma: no cover
        def __init__(self, title: str | None = None) -> None:
            self.title = title
            self.columns: list[str] = []
            self.rows: list[list[str]] = []

        def add_column(self, name: str) -> None:
            self.columns.append(name)

        def add_row(self, *values: str) -> None:
            self.rows.append(list(values))

        def __str__(self) -> str:
            lines = []
            if self.title:
                lines.append(self.title)
            if self.columns:
                lines.append("\t".join(self.columns))
            lines.extend("\t".join(row) for row in self.rows)
            return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="vepyr-diff")
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("--preset", default=env_str("VEPYR_DIFFLY_PRESET"))
    run_parser.add_argument(
        "--execution-mode", default=env_str("VEPYR_DIFFLY_EXECUTION_MODE") or "docker"
    )
    run_parser.add_argument("--input-vcf", type=Path, default=env_path("VEPYR_DIFFLY_INPUT_VCF"))
    run_parser.add_argument("--output-dir", type=Path, default=env_path("VEPYR_DIFFLY_OUTPUT_DIR"))
    run_parser.add_argument(
        "--sample-first-n", type=int, default=env_int("VEPYR_DIFFLY_SAMPLE_FIRST_N")
    )
    run_parser.add_argument("--chromosomes", default=env_str("VEPYR_DIFFLY_CHROMOSOMES"))
    run_parser.add_argument("--plugins", default=env_str("VEPYR_DIFFLY_PLUGINS"))
    run_parser.add_argument(
        "--compare-only-plugins",
        action=argparse.BooleanOptionalAction,
        default=env_bool("VEPYR_DIFFLY_COMPARE_ONLY_PLUGINS") or False,
    )
    run_parser.add_argument("--vepyr-path", type=Path, default=env_path("VEPYR_DIFFLY_VEPYR_PATH"))
    run_parser.add_argument(
        "--vepyr-python", type=Path, default=env_path("VEPYR_DIFFLY_VEPYR_PYTHON")
    )
    run_parser.add_argument(
        "--vep-cache-dir", type=Path, default=env_path("VEPYR_DIFFLY_VEP_CACHE_DIR")
    )
    run_parser.add_argument(
        "--vepyr-cache-output-dir",
        type=Path,
        default=env_path("VEPYR_DIFFLY_VEPYR_CACHE_OUTPUT_DIR"),
    )
    run_parser.add_argument(
        "--reference-fasta", type=Path, default=env_path("VEPYR_DIFFLY_REFERENCE_FASTA")
    )
    run_parser.add_argument("--vep-bin", type=Path, default=env_path("VEPYR_DIFFLY_VEP_BIN"))
    run_parser.add_argument(
        "--vep-cache-version", default=env_str("VEPYR_DIFFLY_VEP_CACHE_VERSION")
    )
    run_parser.add_argument("--vep-perl5lib", default=env_str("VEPYR_DIFFLY_VEP_PERL5LIB"))
    run_parser.add_argument(
        "--use-fjall",
        action=argparse.BooleanOptionalAction,
        default=env_bool("VEPYR_DIFFLY_VEPYR_USE_FJALL") or False,
    )
    run_parser.add_argument(
        "--compare-mode",
        choices=("fast", "debug"),
        default=env_str("VEPYR_DIFFLY_COMPARE_MODE") or "fast",
    )
    run_parser.add_argument(
        "--bucket-count", type=int, default=env_int("VEPYR_DIFFLY_BUCKET_COUNT")
    )
    run_parser.add_argument(
        "--compare-workers", type=int, default=env_int("VEPYR_DIFFLY_COMPARE_WORKERS")
    )
    run_parser.add_argument(
        "--memory-budget-mb", type=int, default=env_int("VEPYR_DIFFLY_MEMORY_BUDGET_MB")
    )
    run_parser.add_argument(
        "--fingerprint-only",
        action=argparse.BooleanOptionalAction,
        default=env_bool("VEPYR_DIFFLY_FINGERPRINT_ONLY") or False,
    )

    vepyr_only_parser = subparsers.add_parser("annotate-vepyr")
    vepyr_only_parser.add_argument("--input-vcf", required=True, type=Path)
    vepyr_only_parser.add_argument("--output-vcf", required=True, type=Path)
    vepyr_only_parser.add_argument(
        "--cache-dir", type=Path, default=env_path("VEPYR_DIFFLY_VEPYR_CACHE_OUTPUT_DIR")
    )
    vepyr_only_parser.add_argument(
        "--reference-fasta", type=Path, default=env_path("VEPYR_DIFFLY_REFERENCE_FASTA")
    )
    vepyr_only_parser.add_argument(
        "--vepyr-python", type=Path, default=env_path("VEPYR_DIFFLY_VEPYR_PYTHON")
    )
    vepyr_only_parser.add_argument(
        "--use-fjall",
        action=argparse.BooleanOptionalAction,
        default=env_bool("VEPYR_DIFFLY_VEPYR_USE_FJALL") or False,
    )
    vepyr_only_parser.add_argument("--plugins", default=env_str("VEPYR_DIFFLY_PLUGINS"))
    vepyr_only_parser.add_argument("--log-path", type=Path)

    compare_parser = subparsers.add_parser("compare-existing")
    compare_parser.add_argument("--preset", default=env_str("VEPYR_DIFFLY_PRESET"))
    compare_parser.add_argument("--left-vcf", required=True, type=Path)
    compare_parser.add_argument("--right-vcf", required=True, type=Path)
    compare_parser.add_argument("--output-dir", required=True, type=Path)
    compare_parser.add_argument("--chromosomes", default=env_str("VEPYR_DIFFLY_CHROMOSOMES"))
    compare_parser.add_argument("--plugins", default=env_str("VEPYR_DIFFLY_PLUGINS"))
    compare_parser.add_argument(
        "--compare-only-plugins",
        action=argparse.BooleanOptionalAction,
        default=env_bool("VEPYR_DIFFLY_COMPARE_ONLY_PLUGINS") or False,
    )
    compare_parser.add_argument(
        "--compare-mode",
        choices=("fast", "debug"),
        default=env_str("VEPYR_DIFFLY_COMPARE_MODE") or "fast",
    )
    compare_parser.add_argument(
        "--bucket-count", type=int, default=env_int("VEPYR_DIFFLY_BUCKET_COUNT")
    )
    compare_parser.add_argument(
        "--compare-workers", type=int, default=env_int("VEPYR_DIFFLY_COMPARE_WORKERS")
    )
    compare_parser.add_argument(
        "--memory-budget-mb", type=int, default=env_int("VEPYR_DIFFLY_MEMORY_BUDGET_MB")
    )
    compare_parser.add_argument(
        "--fingerprint-only",
        action=argparse.BooleanOptionalAction,
        default=env_bool("VEPYR_DIFFLY_FINGERPRINT_ONLY") or False,
    )

    benchmark_parser = subparsers.add_parser("benchmark-compare")
    benchmark_parser.add_argument("--left-vcf", required=True, type=Path)
    benchmark_parser.add_argument("--right-vcf", required=True, type=Path)
    benchmark_parser.add_argument("--output-json", required=True, type=Path)
    benchmark_parser.add_argument("--modes", default="fast,debug")
    benchmark_parser.add_argument(
        "--bucket-count", type=int, default=env_int("VEPYR_DIFFLY_BUCKET_COUNT")
    )
    benchmark_parser.add_argument(
        "--compare-workers", type=int, default=env_int("VEPYR_DIFFLY_COMPARE_WORKERS")
    )
    benchmark_parser.add_argument(
        "--memory-budget-mb", type=int, default=env_int("VEPYR_DIFFLY_MEMORY_BUDGET_MB")
    )

    benchmark_plugins_parser = subparsers.add_parser("benchmark-plugins")
    benchmark_plugins_parser.add_argument("--preset", default=env_str("VEPYR_DIFFLY_PRESET"))
    benchmark_plugins_parser.add_argument("--input-vcf", required=True, type=Path)
    benchmark_plugins_parser.add_argument(
        "--reference-vcf",
        required=True,
        type=Path,
        help="Already-annotated VEP output (everything+hgvs+all plugins) to compare every vepyr config against.",
    )
    benchmark_plugins_parser.add_argument("--output-dir", required=True, type=Path)
    benchmark_plugins_parser.add_argument(
        "--cache-dir", required=True, type=Path, help="Core Ensembl parquet cache for vepyr."
    )
    benchmark_plugins_parser.add_argument(
        "--plugin-cache-root",
        required=True,
        type=Path,
        help="Root already containing plugin/<name>/ for every plugin in --plugins.",
    )
    benchmark_plugins_parser.add_argument(
        "--plugins",
        default=env_str("VEPYR_DIFFLY_PLUGINS") or "clinvar,spliceai,cadd,alphamissense,dbnsfp",
        help="Comma-separated plugin names to benchmark one at a time against a 'base' (no-plugin) run.",
    )
    benchmark_plugins_parser.add_argument(
        "--reference-fasta", default=env_str("VEPYR_DIFFLY_REFERENCE_FASTA")
    )
    benchmark_plugins_parser.add_argument(
        "--hgvs", action=argparse.BooleanOptionalAction, default=True
    )
    benchmark_plugins_parser.add_argument(
        "--vepyr-python", default=env_str("VEPYR_DIFFLY_VEPYR_PYTHON")
    )
    benchmark_plugins_parser.add_argument(
        "--rss-poll-interval-s",
        type=float,
        default=0.2,
        help="How often to sample the annotation subprocess's RSS while it runs.",
    )
    benchmark_plugins_parser.add_argument(
        "--compare-mode",
        choices=("fast", "debug"),
        default=env_str("VEPYR_DIFFLY_COMPARE_MODE") or "fast",
    )
    benchmark_plugins_parser.add_argument(
        "--bucket-count", type=int, default=env_int("VEPYR_DIFFLY_BUCKET_COUNT")
    )
    benchmark_plugins_parser.add_argument(
        "--compare-workers", type=int, default=env_int("VEPYR_DIFFLY_COMPARE_WORKERS")
    )
    benchmark_plugins_parser.add_argument(
        "--memory-budget-mb", type=int, default=env_int("VEPYR_DIFFLY_MEMORY_BUDGET_MB")
    )
    benchmark_plugins_parser.add_argument(
        "--gdrive-upload",
        action=argparse.BooleanOptionalAction,
        default=env_bool("VEPYR_DIFFLY_GDRIVE_UPLOAD") or False,
        help="Upload each config's annotated VCF to gdrive then delete the local copy "
        "once it has stopped being needed (disk-constrained full-genome runs).",
    )
    benchmark_plugins_parser.add_argument(
        "--gdrive-folder-id",
        default=env_str("VEPYR_DIFFLY_GDRIVE_FOLDER_ID") or "1Nkm66F6jStCLrZRZ0SvGQYgdqEr6_CGU",
    )
    benchmark_plugins_parser.add_argument(
        "--gdrive-dest-prefix",
        default=env_str("VEPYR_DIFFLY_GDRIVE_DEST_PREFIX") or "vepyr_results/116/benchmark-plugins",
    )
    benchmark_plugins_parser.add_argument(
        "--include-base",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Prepend the no-plugin 'base' config. Pass --no-include-base to run just "
        "the requested --plugins (e.g. one config at a time across separate invocations).",
    )

    subparsers.add_parser("list-presets")

    inspect_parser = subparsers.add_parser("inspect-run")
    inspect_parser.add_argument("--run-dir", required=True, type=Path)

    mismatch_parser = subparsers.add_parser("analyze-consequence-mismatches")
    mismatch_parser.add_argument("--run-dir", type=Path)
    mismatch_parser.add_argument("--mismatches-tsv", type=Path)
    mismatch_parser.add_argument("--output-json", type=Path)

    raw_csq_parser = subparsers.add_parser("extract-mismatch-csq-examples")
    raw_csq_parser.add_argument("--run-dir", type=Path)
    raw_csq_parser.add_argument("--analysis-json", type=Path)
    raw_csq_parser.add_argument("--left-vcf", type=Path)
    raw_csq_parser.add_argument("--right-vcf", type=Path)
    raw_csq_parser.add_argument("--output-json", required=True, type=Path)
    raw_csq_parser.add_argument("--per-category-limit", type=int, default=3)
    return parser


def _chromosome_equal(payload: dict[str, int]) -> bool:
    return (
        int(payload.get("left_only_rows", 0)) == 0
        and int(payload.get("right_only_rows", 0)) == 0
        and int(payload.get("unequal_rows", 0)) == 0
    )


def _resolve_compare_csq_fields(
    *,
    left_csq_fields: list[str],
    right_csq_fields: list[str],
    plugins: list[str],
    compare_only_plugins: bool,
) -> list[str]:
    if not compare_only_plugins:
        if left_csq_fields != right_csq_fields:
            raise ValueError("left/right CSQ headers differ and cannot be compared semantically")
        return left_csq_fields

    if not plugins:
        raise ValueError("--compare-only-plugins requires --plugins with at least one plugin")

    requested_fields = compare_plugin_fields(plugins)
    field_aliases = compare_plugin_field_aliases(plugins)
    left_selected = [
        field
        for field in requested_fields
        if any(alias in left_csq_fields for alias in field_aliases[field])
    ]
    right_selected = [
        field
        for field in requested_fields
        if any(alias in right_csq_fields for alias in field_aliases[field])
    ]
    if left_selected != right_selected:
        raise ValueError(
            "left/right plugin CSQ fields differ and cannot be compared with --compare-only-plugins"
        )
    if not left_selected:
        requested = ", ".join(requested_fields)
        raise ValueError(
            "no requested plugin CSQ fields were found in either header for "
            f"--compare-only-plugins: {requested}"
        )
    return left_selected


def _resolve_csq_field_indexes(
    *,
    selected_fields: list[str],
    header_fields: list[str],
    plugins: list[str],
    compare_only_plugins: bool,
) -> dict[str, int]:
    header_indexes = {field: index for index, field in enumerate(header_fields)}
    if not compare_only_plugins:
        return header_indexes

    field_aliases = compare_plugin_field_aliases(plugins)
    selected_indexes: dict[str, int] = {}
    for field in selected_fields:
        for alias in field_aliases[field]:
            if alias in header_indexes:
                selected_indexes[field] = header_indexes[alias]
                break
    return selected_indexes


def _attribute_stage_timings(
    *,
    timings: dict[str, float],
    variant: object,
    consequence: object,
) -> dict[str, object]:
    from .models import TierResult

    assert isinstance(variant, TierResult)
    assert isinstance(consequence, TierResult)
    chroms = sorted(set(variant.per_chromosome) | set(consequence.per_chromosome))
    if not chroms:
        return {
            "requested": [],
            "effective": [],
            "per_chromosome": {},
        }

    variant_total = sum(
        int(payload.get("left_rows", 0)) + int(payload.get("right_rows", 0))
        for payload in variant.per_chromosome.values()
    )
    consequence_total = sum(
        int(payload.get("left_rows", 0)) + int(payload.get("right_rows", 0))
        for payload in consequence.per_chromosome.values()
    )
    diff_variant_total = sum(
        int(payload.get("joined_equal_rows", 0))
        + int(payload.get("left_only_rows", 0))
        + int(payload.get("right_only_rows", 0))
        + int(payload.get("unequal_rows", 0))
        for payload in variant.per_chromosome.values()
    )
    diff_consequence_total = sum(
        int(payload.get("joined_equal_rows", 0))
        + int(payload.get("left_only_rows", 0))
        + int(payload.get("right_only_rows", 0))
        + int(payload.get("unequal_rows", 0))
        for payload in consequence.per_chromosome.values()
    )

    per_chromosome: dict[str, object] = {}
    for chrom in chroms:
        variant_payload = variant.per_chromosome.get(
            chrom,
            {
                "left_rows": 0,
                "right_rows": 0,
                "joined_equal_rows": 0,
                "left_only_rows": 0,
                "right_only_rows": 0,
                "unequal_rows": 0,
            },
        )
        consequence_payload = consequence.per_chromosome.get(
            chrom,
            {
                "left_rows": 0,
                "right_rows": 0,
                "joined_equal_rows": 0,
                "left_only_rows": 0,
                "right_only_rows": 0,
                "unequal_rows": 0,
            },
        )
        variant_weight = (
            int(variant_payload["left_rows"]) + int(variant_payload["right_rows"])
        ) / variant_total if variant_total else 0.0
        consequence_weight = (
            int(consequence_payload["left_rows"]) + int(consequence_payload["right_rows"])
        ) / consequence_total if consequence_total else 0.0
        variant_diff_weight = (
            int(variant_payload["joined_equal_rows"])
            + int(variant_payload["left_only_rows"])
            + int(variant_payload["right_only_rows"])
            + int(variant_payload["unequal_rows"])
        ) / diff_variant_total if diff_variant_total else 0.0
        consequence_diff_weight = (
            int(consequence_payload["joined_equal_rows"])
            + int(consequence_payload["left_only_rows"])
            + int(consequence_payload["right_only_rows"])
            + int(consequence_payload["unequal_rows"])
        ) / diff_consequence_total if diff_consequence_total else 0.0
        per_chromosome[chrom] = {
            "variant": {
                **variant_payload,
                "equal": _chromosome_equal(variant_payload),
            },
            "consequence": {
                **consequence_payload,
                "equal": _chromosome_equal(consequence_payload),
            },
            "timings": {
                "variant_summary_seconds": round(
                    timings.get("variant_summary_seconds", 0.0) * variant_weight,
                    3,
                ),
                "consequence_bucketization_seconds": round(
                    timings.get("consequence_bucketization_seconds", 0.0) * consequence_weight,
                    3,
                ),
                "consequence_summary_seconds": round(
                    timings.get("consequence_summary_seconds", 0.0) * consequence_weight,
                    3,
                ),
                "variant_diff_seconds": round(
                    timings.get("variant_diff_seconds", 0.0) * variant_diff_weight,
                    3,
                ),
                "consequence_diff_seconds": round(
                    timings.get("consequence_diff_seconds", 0.0) * consequence_diff_weight,
                    3,
                ),
            },
        }
    return {
        "per_chromosome": per_chromosome,
    }


def _cmd_list_presets(console: Console) -> int:
    presets = load_presets()
    table = Table(title="Presets")
    table.add_column("Name")
    table.add_column("Enabled")
    table.add_column("Cache")
    table.add_column("Description")
    for preset in presets.values():
        table.add_row(
            preset.name,
            "yes" if preset.enabled else "no",
            preset.cache_flavor,
            preset.description,
        )
    console.print(table)
    return 0


def _cmd_annotate_vepyr(args: argparse.Namespace, console: Console) -> int:
    from .runtime import run_vepyr_annotation

    cache_dir = args.cache_dir
    if cache_dir is None:
        raise ValueError("--cache-dir is required")
    # Allow passing either the root cache dir or the parquet feature root already.
    feature_root_name = f"{env_str('VEPYR_DIFFLY_VEP_CACHE_VERSION') or '115'}_GRCh38_vep"
    if cache_dir.name != feature_root_name:
        current_feature_root = cache_dir / feature_root_name
        legacy_feature_root = cache_dir / "parquet" / feature_root_name
        if current_feature_root.exists():
            cache_dir = current_feature_root
        elif legacy_feature_root.exists():
            cache_dir = legacy_feature_root
    log_path = args.log_path or args.output_vcf.parent / "vepyr.log"
    console.print(f"vepyr: annotating {args.input_vcf} -> {args.output_vcf}")
    run_vepyr_annotation(
        input_vcf=args.input_vcf,
        output_vcf=args.output_vcf,
        cache_dir=cache_dir,
        log_path=log_path,
        reference_fasta=args.reference_fasta,
        vepyr_python=args.vepyr_python,
        use_fjall=args.use_fjall,
        plugins=parse_plugin_list(args.plugins),
    )
    console.print(f"vepyr: completed, log at {log_path}")
    return 0


def _cmd_benchmark_compare(args: argparse.Namespace, console: Console) -> int:
    modes = [mode.strip() for mode in args.modes.split(",") if mode.strip()]
    output_json = args.output_json.resolve()
    output_json.parent.mkdir(parents=True, exist_ok=True)
    results: list[dict[str, object]] = []
    with tempfile.TemporaryDirectory(prefix="vepyr-diffly-bench-") as temp_root:
        for mode in modes:
            temp_dir = Path(temp_root) / mode
            command = [
                sys.executable,
                "-m",
                "vepyr_diffly.cli",
                "compare-existing",
                "--preset",
                env_str("VEPYR_DIFFLY_PRESET") or "ensembl_everything",
                "--left-vcf",
                str(args.left_vcf),
                "--right-vcf",
                str(args.right_vcf),
                "--output-dir",
                str(temp_dir),
                "--compare-mode",
                mode,
            ]
            if args.bucket_count is not None:
                command.extend(["--bucket-count", str(args.bucket_count)])
            if args.compare_workers is not None:
                command.extend(["--compare-workers", str(args.compare_workers)])
            if args.memory_budget_mb is not None:
                command.extend(["--memory-budget-mb", str(args.memory_budget_mb)])
            env = os.environ.copy()
            env["PYTHONPATH"] = str(Path(__file__).resolve().parents[1]) + (
                os.pathsep + env["PYTHONPATH"] if "PYTHONPATH" in env else ""
            )
            subprocess.run(
                command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, env=env
            )
            summary = json.loads((temp_dir / "summary.json").read_text())
            results.append(
                {
                    "mode": mode,
                    "variant_equal": summary["tiers"]["variant"]["equal"],
                    "consequence_equal": summary["tiers"]["consequence"]["equal"],
                    "timings": summary["timings"],
                    "consequence_details": summary["tiers"]["consequence"].get("details", {}),
                }
            )
    output_json.write_text(
        json.dumps(
            {"left_vcf": str(args.left_vcf), "right_vcf": str(args.right_vcf), "results": results},
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    console.print(output_json)
    return 0


def _dir_size(path: Path) -> int:
    if not path.exists():
        return 0
    total = 0
    for entry in os.walk(path, followlinks=True):
        dirpath, _dirnames, filenames = entry
        for filename in filenames:
            try:
                total += (Path(dirpath) / filename).stat().st_size
            except OSError:
                continue
    return total


def _snv_indel_breakdown(vcf_path: Path) -> dict[str, int]:
    total = 0
    snv = 0
    with vcf_path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            total += 1
            ref, alt = fields[3], fields[4]
            alts = alt.split(",")
            if len(ref) == 1 and all(len(a) == 1 for a in alts):
                snv += 1
    return {"total": total, "snv": snv, "indel_or_other": total - snv}


def _read_csq_format_fields(vcf_path: Path) -> list[str]:
    with vcf_path.open() as handle:
        for line in handle:
            if line.startswith("##INFO=<ID=CSQ"):
                marker = "Format: "
                marker_index = line.find(marker)
                if marker_index == -1:
                    return []
                tail = line[marker_index + len(marker) :].split('">')[0]
                return tail.split("|")
            if not line.startswith("#"):
                break
    return []


def _extract_csq_values(
    vcf_path: Path, field_index: int, *, first_entry_only: bool
) -> list[str]:
    values: list[str] = []
    with vcf_path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            csq_start = line.find("CSQ=")
            if csq_start == -1:
                continue
            value_start = csq_start + 4
            end_semi = line.find(";", value_start)
            end_tab = line.find("\t", value_start)
            ends = [e for e in (end_semi, end_tab) if e != -1]
            csq_value = line[value_start : min(ends)] if ends else line[value_start:]
            entries = csq_value.split(",")
            if first_entry_only:
                entries = entries[:1]
            for entry in entries:
                subfields = entry.split("|")
                if field_index < len(subfields) and subfields[field_index]:
                    values.append(subfields[field_index])
    return values


def _consequence_type_breakdown(vcf_path: Path) -> dict[str, int]:
    from collections import Counter

    header_fields = _read_csq_format_fields(vcf_path)
    if "Consequence" not in header_fields:
        return {}
    field_index = header_fields.index("Consequence")
    counts: Counter[str] = Counter()
    for value in _extract_csq_values(vcf_path, field_index, first_entry_only=False):
        for term in value.split("&"):
            if term:
                counts[term] += 1
    return dict(counts.most_common())


PRIMARY_PLUGIN_FIELD: dict[str, str] = {
    "clinvar": "ClinVar_CLNSIG",
    "spliceai": "ds_ag",
    "cadd": "cadd_phred",
    "alphamissense": "am_class",
    "dbnsfp": "revel_score",
}


def _primary_plugin_field_histogram(vcf_path: Path, plugin: str) -> dict[str, object] | None:
    from collections import Counter

    canonical_field = PRIMARY_PLUGIN_FIELD.get(plugin)
    if canonical_field is None:
        return None
    aliases = compare_plugin_field_aliases([plugin]).get(canonical_field, (canonical_field,))
    header_fields = _read_csq_format_fields(vcf_path)
    field_index = None
    matched_alias = None
    for alias in aliases:
        if alias in header_fields:
            field_index = header_fields.index(alias)
            matched_alias = alias
            break
    if field_index is None:
        return None

    values = _extract_csq_values(vcf_path, field_index, first_entry_only=True)
    if not values:
        return {"field": matched_alias, "kind": "empty", "non_null_count": 0}

    # dbNSFP fields carry one value per overlapping transcript, joined with "&"
    # ("." marks a missing per-transcript value) -- flatten those before
    # judging numeric-ness, otherwise every dbNSFP field falls back to an
    # unreadable categorical dump of raw joined strings.
    numeric_values: list[float] = []
    is_numeric = True
    for value in values:
        for token in value.split("&"):
            if token in ("", "."):
                continue
            try:
                numeric_values.append(float(token))
            except ValueError:
                is_numeric = False
                break
        if not is_numeric:
            break

    if is_numeric and numeric_values:
        lo, hi = min(numeric_values), max(numeric_values)
        bucket_count = 5
        buckets = [0] * bucket_count
        span = hi - lo
        for value in numeric_values:
            bucket_index = (
                0 if span == 0 else min(bucket_count - 1, int((value - lo) / span * bucket_count))
            )
            buckets[bucket_index] += 1
        edges = [round(lo + span * i / bucket_count, 4) for i in range(bucket_count + 1)]
        return {
            "field": matched_alias,
            "kind": "numeric",
            "non_null_count": len(numeric_values),
            "min": round(lo, 4),
            "max": round(hi, 4),
            "mean": round(sum(numeric_values) / len(numeric_values), 4),
            "bucket_edges": edges,
            "bucket_counts": buckets,
        }

    counts = Counter(values)
    return {
        "field": matched_alias,
        "kind": "categorical",
        "non_null_count": len(values),
        "top_values": counts.most_common(10),
    }


def _run_annotation_with_rss(
    *,
    vepyr_python: str,
    input_vcf: Path,
    output_vcf: Path,
    cache_dir: Path,
    reference_fasta: str | None,
    plugin_cache_root: Path | None,
    hgvs: bool,
    poll_interval_s: float,
) -> dict[str, object]:
    import psutil

    runner_path = Path(__file__).resolve().parent / "plugin_bench_runner.py"
    command = [
        vepyr_python,
        str(runner_path),
        "--input-vcf",
        str(input_vcf),
        "--output-vcf",
        str(output_vcf),
        "--cache-dir",
        str(cache_dir),
    ]
    if reference_fasta:
        command.extend(["--reference-fasta", reference_fasta])
    if plugin_cache_root is not None:
        command.extend(["--plugin-cache-root", str(plugin_cache_root)])
    if hgvs:
        command.append("--hgvs")
    env = os.environ.copy()
    src_dir = str(Path(__file__).resolve().parents[1])
    env["PYTHONPATH"] = src_dir + (os.pathsep + env["PYTHONPATH"] if "PYTHONPATH" in env else "")

    peak_rss_bytes = 0
    cpu_user_seconds = 0.0
    cpu_system_seconds = 0.0
    start = perf_counter()
    process = subprocess.Popen(command, env=env)
    proc_handle = psutil.Process(process.pid)
    while process.poll() is None:
        try:
            procs = [proc_handle, *proc_handle.children(recursive=True)]
            rss_total = 0
            cpu_user_total = 0.0
            cpu_system_total = 0.0
            for proc in procs:
                try:
                    rss_total += proc.memory_info().rss
                    cpu_times = proc.cpu_times()
                    cpu_user_total += cpu_times.user
                    cpu_system_total += cpu_times.system
                except psutil.NoSuchProcess:
                    continue
            peak_rss_bytes = max(peak_rss_bytes, rss_total)
            cpu_user_seconds = max(cpu_user_seconds, cpu_user_total)
            cpu_system_seconds = max(cpu_system_seconds, cpu_system_total)
        except psutil.NoSuchProcess:
            pass
        try:
            process.wait(timeout=poll_interval_s)
        except subprocess.TimeoutExpired:
            continue
    returncode = process.wait()
    wall_seconds = perf_counter() - start
    if returncode != 0:
        raise RuntimeError(f"annotation subprocess failed (rc={returncode}): {' '.join(command)}")
    return {
        "wall_seconds": round(wall_seconds, 3),
        "peak_rss_mb": round(peak_rss_bytes / (1024 * 1024), 1),
        "cpu_user_seconds": round(cpu_user_seconds, 3),
        "cpu_system_seconds": round(cpu_system_seconds, 3),
        "output_bytes": output_vcf.stat().st_size if output_vcf.exists() else 0,
    }


def _run_compare_subprocess(
    *,
    left_vcf: Path,
    right_vcf: Path,
    output_dir: Path,
    preset: str,
    plugins: list[str],
    compare_only_plugins: bool,
) -> dict[str, object]:
    command = [
        sys.executable,
        "-m",
        "vepyr_diffly.cli",
        "compare-existing",
        "--preset",
        preset,
        "--left-vcf",
        str(left_vcf),
        "--right-vcf",
        str(right_vcf),
        "--output-dir",
        str(output_dir),
    ]
    if plugins:
        command.extend(["--plugins", ",".join(plugins)])
    if compare_only_plugins:
        command.append("--compare-only-plugins")
    env = os.environ.copy()
    env["PYTHONPATH"] = str(Path(__file__).resolve().parents[1]) + (
        os.pathsep + env["PYTHONPATH"] if "PYTHONPATH" in env else ""
    )
    subprocess.run(
        command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, env=env
    )
    return json.loads((output_dir / "summary.json").read_text())


def _format_bytes_mb(num_bytes: object) -> str:
    return f"{float(num_bytes) / (1024 * 1024):.1f} MB"


def _format_delta_mb(num_bytes: object) -> str:
    value = float(num_bytes) / (1024 * 1024)
    if value == 0:
        return "—"
    return f"{'+' if value > 0 else ''}{value:.1f} MB"


def _histogram_md_lines(histogram: dict[str, object] | None) -> list[str]:
    if histogram is None:
        return []
    field = histogram.get("field", "?")
    kind = histogram.get("kind")
    if kind == "empty":
        return [f"**{field}**: no non-empty values", ""]
    if kind == "numeric":
        edges = histogram["bucket_edges"]
        counts = histogram["bucket_counts"]
        lines = [
            f"**{field}** (numeric, n={histogram['non_null_count']}, "
            f"min={histogram['min']}, max={histogram['max']}, mean={histogram['mean']}):",
            "",
            "| bucket | count |",
            "|---|---|",
        ]
        for i, count in enumerate(counts):
            lines.append(f"| [{edges[i]}, {edges[i + 1]}] | {count} |")
        lines.append("")
        return lines
    lines = [f"**{field}** (categorical, n={histogram['non_null_count']}, top values):", "", "| value | count |", "|---|---|"]
    for value, count in histogram.get("top_values", []):
        lines.append(f"| {value} | {count} |")
    lines.append("")
    return lines


def _consequence_md_lines(breakdown: dict[str, int], top_n: int = 10) -> list[str]:
    if not breakdown:
        return []
    lines = ["**Consequence types (top {}):**".format(top_n), "", "| type | count |", "|---|---|"]
    for term, count in list(breakdown.items())[:top_n]:
        lines.append(f"| {term} | {count} |")
    lines.append("")
    return lines


def _upload_and_remove(
    local_path: Path, *, folder_id: str, dest_prefix: str, console: Console
) -> None:
    remote_dir = f"gdrive:{dest_prefix}"
    console.print(f"[dim]uploading {local_path.name} -> {remote_dir}[/dim]")
    subprocess.run(
        [
            "rclone",
            "copy",
            str(local_path),
            remote_dir,
            "--drive-root-folder-id",
            folder_id,
            "--low-level-retries",
            "50",
        ],
        check=True,
    )
    listing = subprocess.run(
        [
            "rclone",
            "lsjson",
            remote_dir,
            "--drive-root-folder-id",
            folder_id,
            "--files-only",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    entries = json.loads(listing.stdout)
    local_size = local_path.stat().st_size
    matched = next((e for e in entries if e["Name"] == local_path.name), None)
    if matched is None or matched["Size"] != local_size:
        raise RuntimeError(
            f"upload verification failed for {local_path.name}: "
            f"remote entry {matched!r}, local size {local_size}"
        )
    local_path.unlink()
    console.print(f"[dim]verified + removed local {local_path.name}[/dim]")


def _cmd_benchmark_plugins(args: argparse.Namespace, console: Console) -> int:
    if args.preset is None:
        raise ValueError("--preset is required or set VEPYR_DIFFLY_PRESET in .env")
    vepyr_python = args.vepyr_python or sys.executable
    plugin_names = parse_plugin_list(args.plugins)
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    input_breakdown = _snv_indel_breakdown(args.input_vcf)
    total_variants = input_breakdown["total"]

    configs: list[tuple[str, str | None]] = (
        [("base", None)] if args.include_base else []
    ) + [(name, name) for name in plugin_names]
    results: list[dict[str, object]] = []
    for config_name, plugin in configs:
        console.print(f"[bold]benchmark-plugins: running config '{config_name}'[/bold]")
        config_dir = output_dir / config_name
        config_dir.mkdir(parents=True, exist_ok=True)

        if plugin is None:
            plugin_cache_root_for_run: Path | None = None
            cache_bytes = 0
        else:
            farm_root = config_dir / "plugin_cache"
            farm_plugin_dir = farm_root / "plugin"
            farm_plugin_dir.mkdir(parents=True, exist_ok=True)
            link_path = farm_plugin_dir / plugin
            if not link_path.exists():
                link_path.symlink_to(
                    (args.plugin_cache_root / "plugin" / plugin).resolve(), target_is_directory=True
                )
            plugin_cache_root_for_run = farm_root
            cache_bytes = _dir_size(args.plugin_cache_root / "plugin" / plugin)

        annotated_vcf = config_dir / f"vepyr.{config_name}.annotated.vcf"
        run_stats = _run_annotation_with_rss(
            vepyr_python=vepyr_python,
            input_vcf=args.input_vcf,
            output_vcf=annotated_vcf,
            cache_dir=args.cache_dir,
            reference_fasta=args.reference_fasta,
            plugin_cache_root=plugin_cache_root_for_run,
            hgvs=args.hgvs,
            poll_interval_s=args.rss_poll_interval_s,
        )
        run_stats["cache_bytes"] = cache_bytes
        run_stats["variants_per_second"] = (
            round(total_variants / run_stats["wall_seconds"], 1)
            if run_stats["wall_seconds"] > 0
            else None
        )
        run_stats["consequence_breakdown"] = _consequence_type_breakdown(annotated_vcf)
        run_stats["plugin_value_histogram"] = (
            _primary_plugin_field_histogram(annotated_vcf, plugin) if plugin is not None else None
        )

        concordance: dict[str, object] | None = None
        if plugin is not None:
            compare_summary = _run_compare_subprocess(
                left_vcf=args.reference_vcf,
                right_vcf=annotated_vcf,
                output_dir=config_dir / "compare",
                preset=args.preset,
                plugins=[plugin],
                compare_only_plugins=True,
            )
            plugin_tier = compare_summary["tiers"]["consequence"]
            concordance = {
                "plugin_fields_equal": plugin_tier["equal"],
                "joined_equal_rows": plugin_tier["joined_equal_rows"],
                "left_only_rows": plugin_tier["left_only_rows"],
                "right_only_rows": plugin_tier["right_only_rows"],
                "unequal_rows": plugin_tier["unequal_rows"],
                "per_chromosome": plugin_tier.get("per_chromosome", {}),
            }
        else:
            console.print(
                "[yellow]benchmark-plugins: skipping concordance for 'base' "
                "(no no-plugin VEP reference given)[/yellow]"
            )

        if args.gdrive_upload:
            _upload_and_remove(
                annotated_vcf,
                folder_id=args.gdrive_folder_id,
                dest_prefix=args.gdrive_dest_prefix,
                console=console,
            )

        results.append(
            {
                "config": config_name,
                "plugin": plugin,
                **run_stats,
                "concordance": concordance,
            }
        )

    base_output_bytes = results[0]["output_bytes"]
    for row in results:
        row["output_bytes_delta_vs_base"] = row["output_bytes"] - base_output_bytes

    summary = {
        "input_vcf": str(args.input_vcf),
        "reference_vcf": str(args.reference_vcf),
        "plugins": plugin_names,
        "input_breakdown": input_breakdown,
        "results": results,
    }
    summary_json = output_dir / "benchmark_plugins_summary.json"
    summary_json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    table = Table(title="vepyr plugin benchmark")
    for column in (
        "config",
        "wall (s)",
        "cpu (s)",
        "RSS (MB)",
        "var/s",
        "cache",
        "output (Δ base)",
        "plugin fields equal",
        "mismatches (L/R/U)",
    ):
        table.add_column(column)
    for row in results:
        concordance = row["concordance"] or {}
        mismatches = (
            f"{concordance['left_only_rows']}/{concordance['right_only_rows']}/{concordance['unequal_rows']}"
            if concordance
            else "n/a"
        )
        table.add_row(
            str(row["config"]),
            f"{row['wall_seconds']:.1f}",
            f"{row['cpu_user_seconds'] + row['cpu_system_seconds']:.1f}",
            f"{row['peak_rss_mb']:.0f}",
            str(row["variants_per_second"]),
            _format_bytes_mb(row["cache_bytes"]),
            f"{_format_bytes_mb(row['output_bytes'])} ({_format_delta_mb(row['output_bytes_delta_vs_base'])})",
            str(concordance.get("plugin_fields_equal", "n/a")),
            mismatches,
        )
    console.print(table)

    summary_md = output_dir / "benchmark_plugins_summary.md"
    lines = [
        "# vepyr plugin benchmark",
        "",
        f"- input: `{args.input_vcf}`",
        f"- VEP reference: `{args.reference_vcf}`",
        f"- variants: {input_breakdown['total']} "
        f"(SNV: {input_breakdown['snv']}, indel/other: {input_breakdown['indel_or_other']})",
        "",
        "## Summary",
        "",
        "| config | wall (s) | cpu (s) | peak RSS (MB) | throughput (var/s) | cache size | "
        "output size (Δ base) | plugin fields equal | mismatches (L/R/U) |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for row in results:
        concordance = row["concordance"] or {}
        mismatches = (
            f"{concordance['left_only_rows']}/{concordance['right_only_rows']}/{concordance['unequal_rows']}"
            if concordance
            else "n/a"
        )
        lines.append(
            "| {config} | {wall:.1f} | {cpu:.1f} | {rss:.0f} | {vps} | {cache} | {out} | {eq} | {mm} |".format(
                config=row["config"],
                wall=row["wall_seconds"],
                cpu=row["cpu_user_seconds"] + row["cpu_system_seconds"],
                rss=row["peak_rss_mb"],
                vps=row["variants_per_second"],
                cache=_format_bytes_mb(row["cache_bytes"]),
                out=f"{_format_bytes_mb(row['output_bytes'])} ({_format_delta_mb(row['output_bytes_delta_vs_base'])})",
                eq=concordance.get("plugin_fields_equal", "n/a"),
                mm=mismatches,
            )
        )
    lines.extend(["", "## Per-config detail", ""])
    for row in results:
        lines.append(f"### {row['config']}")
        lines.append("")
        lines.extend(_consequence_md_lines(row["consequence_breakdown"]))
        if row["plugin_value_histogram"] is not None:
            lines.extend(_histogram_md_lines(row["plugin_value_histogram"]))
    summary_md.write_text("\n".join(lines) + "\n", encoding="utf-8")

    console.print(summary_json)
    console.print(summary_md)
    return 0


def _run_comparison_pipeline(
    *,
    console: Console,
    config: object,
    artifacts: object,
    left_vcf: Path,
    right_vcf: Path,
    left_name: str,
    right_name: str,
) -> int:
    from .compare import (
        compare_bucketed_consequence_tier,
        compare_bucketed_variant_tier,
        compare_tier,
    )
    from .models import RuntimeConfig, RunArtifacts
    from .normalize import (
        VARIANT_KEY,
        _count_rows,
        materialize_consequence_summary,
        materialize_variant_summary,
    )
    from .progress import ProgressReporter
    from .report import print_run_summary, write_run_summary
    from .resources import plan_compare_resources
    from .schema import validate_variant_schema
    import polars as pl

    assert isinstance(config, RuntimeConfig)
    assert isinstance(artifacts, RunArtifacts)
    reporter = ProgressReporter(log_path=artifacts.progress_log_path, console=console)
    reporter.start()
    timings: dict[str, float] = {}
    try:
        resource_plan = plan_compare_resources(
            config=config,
            left_vcf=left_vcf,
            right_vcf=right_vcf,
        )
        reporter.log(
            "comparison: resource plan "
            + ", ".join(f"{key}={value}" for key, value in resource_plan.to_dict().items())
        )
        reporter.log(
            "comparison: chromosome filter "
            + (
                ",".join(config.selected_chromosomes)
                if config.selected_chromosomes
                else "all"
            )
        )
        reporter.stage("comparison: starting normalization")
        stage_start = perf_counter()
        variant_kwargs = [
            dict(
                vcf_path=left_vcf,
                variant_path=artifacts.left_variant_path,
                bucket_root=artifacts.left_variant_bucket_dir
                if resource_plan.use_bucketed_variant
                else None,
                reporter=reporter,
                side_label=left_name,
                bucket_count=resource_plan.bucket_count,
                chunk_variants=resource_plan.variant_chunk_rows,
                chromosome_aliases=set(config.selected_chromosome_aliases),
            ),
            dict(
                vcf_path=right_vcf,
                variant_path=artifacts.right_variant_path,
                bucket_root=artifacts.right_variant_bucket_dir
                if resource_plan.use_bucketed_variant
                else None,
                reporter=reporter,
                side_label=right_name,
                bucket_count=resource_plan.bucket_count,
                chunk_variants=resource_plan.variant_chunk_rows,
                chromosome_aliases=set(config.selected_chromosome_aliases),
            ),
        ]
        if resource_plan.parallelize_sides:
            with ThreadPoolExecutor(max_workers=2) as executor:
                left_future = executor.submit(materialize_variant_summary, **variant_kwargs[0])
                right_future = executor.submit(materialize_variant_summary, **variant_kwargs[1])
                left_csq_fields = left_future.result()
                right_csq_fields = right_future.result()
        else:
            left_csq_fields = materialize_variant_summary(**variant_kwargs[0])
            right_csq_fields = materialize_variant_summary(**variant_kwargs[1])
        timings["variant_summary_seconds"] = round(perf_counter() - stage_start, 3)
        left_variant_rows = _count_rows(artifacts.left_variant_path)
        right_variant_rows = _count_rows(artifacts.right_variant_path)
        compare_csq_fields = _resolve_compare_csq_fields(
            left_csq_fields=left_csq_fields,
            right_csq_fields=right_csq_fields,
            plugins=config.plugins,
            compare_only_plugins=config.compare_only_plugins,
        )
        left_field_indexes = _resolve_csq_field_indexes(
            selected_fields=compare_csq_fields,
            header_fields=left_csq_fields,
            plugins=config.plugins,
            compare_only_plugins=config.compare_only_plugins,
        )
        right_field_indexes = _resolve_csq_field_indexes(
            selected_fields=compare_csq_fields,
            header_fields=right_csq_fields,
            plugins=config.plugins,
            compare_only_plugins=config.compare_only_plugins,
        )

        if resource_plan.use_bucketed_consequence:
            reporter.stage("comparison: bucketizing consequence rows")
            reporter.log(
                f"comparison: using chunk_variants={resource_plan.consequence_chunk_rows} "
                f"bucket_count={resource_plan.bucket_count} compare_mode={config.compare_mode} "
                f"compare_workers={resource_plan.compare_workers} "
                f"fingerprint_only={'yes' if config.fingerprint_only else 'no'}"
            )
            stage_start = perf_counter()
            consequence_calls = [
                dict(
                    vcf_path=left_vcf,
                    csq_fields=compare_csq_fields,
                    field_indexes=left_field_indexes,
                    drop_empty_csq_rows=config.compare_only_plugins,
                    bucket_key_fields=VARIANT_KEY if config.compare_only_plugins else None,
                    bucket_root=artifacts.left_consequence_bucket_dir,
                    reporter=reporter,
                    side_label=left_name,
                    bucket_count=resource_plan.bucket_count,
                    chunk_variants=resource_plan.consequence_chunk_rows,
                    total_variants=left_variant_rows,
                    chromosome_aliases=set(config.selected_chromosome_aliases),
                ),
                dict(
                    vcf_path=right_vcf,
                    csq_fields=compare_csq_fields,
                    field_indexes=right_field_indexes,
                    drop_empty_csq_rows=config.compare_only_plugins,
                    bucket_key_fields=VARIANT_KEY if config.compare_only_plugins else None,
                    bucket_root=artifacts.right_consequence_bucket_dir,
                    reporter=reporter,
                    side_label=right_name,
                    bucket_count=resource_plan.bucket_count,
                    chunk_variants=resource_plan.consequence_chunk_rows,
                    total_variants=right_variant_rows,
                    chromosome_aliases=set(config.selected_chromosome_aliases),
                ),
            ]
            from .normalize import materialize_consequence_buckets

            if resource_plan.parallelize_sides:
                with ThreadPoolExecutor(max_workers=2) as executor:
                    left_future = executor.submit(
                        materialize_consequence_buckets, **consequence_calls[0]
                    )
                    right_future = executor.submit(
                        materialize_consequence_buckets, **consequence_calls[1]
                    )
                    left_future.result()
                    right_future.result()
            else:
                materialize_consequence_buckets(**consequence_calls[0])
                materialize_consequence_buckets(**consequence_calls[1])
            timings["consequence_bucketization_seconds"] = round(perf_counter() - stage_start, 3)
        else:
            stage_start = perf_counter()
            consequence_calls = [
                dict(
                    vcf_path=left_vcf,
                    consequence_path=artifacts.left_consequence_path,
                    csq_fields=compare_csq_fields,
                    field_indexes=left_field_indexes,
                    reporter=reporter,
                    side_label=left_name,
                    chromosome_aliases=set(config.selected_chromosome_aliases),
                ),
                dict(
                    vcf_path=right_vcf,
                    consequence_path=artifacts.right_consequence_path,
                    csq_fields=compare_csq_fields,
                    field_indexes=right_field_indexes,
                    reporter=reporter,
                    side_label=right_name,
                    chromosome_aliases=set(config.selected_chromosome_aliases),
                ),
            ]
            if resource_plan.parallelize_sides:
                with ThreadPoolExecutor(max_workers=2) as executor:
                    left_future = executor.submit(
                        materialize_consequence_summary, **consequence_calls[0]
                    )
                    right_future = executor.submit(
                        materialize_consequence_summary, **consequence_calls[1]
                    )
                    left_future.result()
                    right_future.result()
            else:
                materialize_consequence_summary(**consequence_calls[0])
                materialize_consequence_summary(**consequence_calls[1])
            timings["consequence_summary_seconds"] = round(perf_counter() - stage_start, 3)

        reporter.stage("comparison: validating normalized variant schemas")
        stage_start = perf_counter()
        left_variant = pl.scan_parquet(artifacts.left_variant_path).head(1000).collect()
        right_variant = pl.scan_parquet(artifacts.right_variant_path).head(1000).collect()
        validate_variant_schema(left_variant)
        validate_variant_schema(right_variant)
        timings["schema_validation_seconds"] = round(perf_counter() - stage_start, 3)

        reporter.stage("comparison: running variant tier diff")
        stage_start = perf_counter()
        if resource_plan.use_bucketed_variant:
            variant = compare_bucketed_variant_tier(
                left_bucket_dir=artifacts.left_variant_bucket_dir,
                right_bucket_dir=artifacts.right_variant_bucket_dir,
                left_name=left_name,
                right_name=right_name,
                diff_frame_path=artifacts.variant_diff_path,
                mismatches_tsv_path=artifacts.variant_mismatches_tsv_path,
                reporter=reporter,
                bucket_count=resource_plan.bucket_count,
                compare_mode=config.compare_mode,
                fingerprint_only=config.fingerprint_only,
            ).tier
        else:
            variant = compare_tier(
                name="variant",
                left=pl.scan_parquet(artifacts.left_variant_path),
                right=pl.scan_parquet(artifacts.right_variant_path),
                primary_key=VARIANT_KEY,
                left_name=left_name,
                right_name=right_name,
                diff_frame_path=artifacts.variant_diff_path,
                mismatches_tsv_path=artifacts.variant_mismatches_tsv_path,
                reporter=reporter,
            ).tier
        timings["variant_diff_seconds"] = round(perf_counter() - stage_start, 3)

        reporter.stage("comparison: running consequence tier diff")
        if resource_plan.use_bucketed_consequence:
            stage_start = perf_counter()
            consequence = compare_bucketed_consequence_tier(
                left_bucket_dir=artifacts.left_consequence_bucket_dir,
                right_bucket_dir=artifacts.right_consequence_bucket_dir,
                csq_fields=compare_csq_fields,
                left_name=left_name,
                right_name=right_name,
                diff_frame_path=artifacts.consequence_diff_path,
                mismatches_tsv_path=artifacts.consequence_mismatches_tsv_path,
                reporter=reporter,
                bucket_count=resource_plan.bucket_count,
                max_workers=resource_plan.compare_workers,
                compare_mode=config.compare_mode,
                fingerprint_only=config.fingerprint_only,
                compare_duplicate_count=not config.compare_only_plugins,
            ).tier
            timings["consequence_diff_seconds"] = round(perf_counter() - stage_start, 3)
        else:
            consequence_key = VARIANT_KEY + compare_csq_fields
            stage_start = perf_counter()
            left_consequence = pl.scan_parquet(artifacts.left_consequence_path)
            right_consequence = pl.scan_parquet(artifacts.right_consequence_path)
            if config.compare_only_plugins:
                left_consequence = left_consequence.select(*consequence_key)
                right_consequence = right_consequence.select(*consequence_key)
            consequence = compare_tier(
                name="consequence",
                left=left_consequence,
                right=right_consequence,
                primary_key=consequence_key,
                left_name=left_name,
                right_name=right_name,
                diff_frame_path=artifacts.consequence_diff_path,
                mismatches_tsv_path=artifacts.consequence_mismatches_tsv_path,
                reporter=reporter,
            ).tier
            timings["consequence_diff_seconds"] = round(perf_counter() - stage_start, 3)

        reporter.stage("comparison: writing summaries")
        chromosome_summary = _attribute_stage_timings(
            timings=timings,
            variant=variant,
            consequence=consequence,
        )
        chromosome_summary["requested"] = (
            config.chromosome_filter_raw.split(",") if config.chromosome_filter_raw else []
        )
        chromosome_summary["effective"] = config.selected_chromosomes
        print_run_summary(
            console=console,
            config=config,
            variant=variant,
            consequence=consequence,
            left_vcf=left_vcf,
            right_vcf=right_vcf,
            progress_log_path=artifacts.progress_log_path,
            resource_plan=resource_plan.to_dict(),
            timings=timings,
            chromosome_summary=chromosome_summary,
        )
        write_run_summary(
            config=config,
            artifacts=artifacts,
            variant=variant,
            consequence=consequence,
            left_vcf=left_vcf,
            right_vcf=right_vcf,
            resource_plan=resource_plan.to_dict(),
            timings=timings,
            chromosome_summary=chromosome_summary,
        )
        for chrom, payload in chromosome_summary.get("per_chromosome", {}).items():
            reporter.log(
                f"comparison: chromosome {chrom} variant_equal={'yes' if payload['variant']['equal'] else 'no'} "
                f"consequence_equal={'yes' if payload['consequence']['equal'] else 'no'}"
            )
        reporter.log("comparison: completed successfully")
        return 0
    finally:
        reporter.stop()


def _cmd_run(args: argparse.Namespace, console: Console) -> int:
    from .runtime import (
        execute_engines,
        prepare_artifacts,
        resolve_runtime_config,
        write_effective_config,
    )

    if args.preset is None:
        raise ValueError("--preset is required or set VEPYR_DIFFLY_PRESET in .env")
    if args.input_vcf is None:
        raise ValueError("--input-vcf is required or set VEPYR_DIFFLY_INPUT_VCF in .env")
    if args.output_dir is None:
        raise ValueError("--output-dir is required or set VEPYR_DIFFLY_OUTPUT_DIR in .env")

    preset = get_preset(args.preset)
    if not preset.enabled:
        raise ValueError(f"preset {preset.name} is disabled")
    config = resolve_runtime_config(
        preset=preset,
        input_vcf=args.input_vcf,
        output_dir=args.output_dir,
        sample_first_n=args.sample_first_n,
        execution_mode=args.execution_mode,
        vepyr_path=args.vepyr_path,
        vepyr_python=args.vepyr_python,
        vep_cache_dir=args.vep_cache_dir,
        vepyr_cache_output_dir=args.vepyr_cache_output_dir,
        reference_fasta=args.reference_fasta,
        vep_bin=args.vep_bin,
        vep_cache_version=args.vep_cache_version,
        vep_perl5lib=args.vep_perl5lib,
        vepyr_use_fjall=args.use_fjall,
        compare_mode=args.compare_mode,
        compare_bucket_count=args.bucket_count,
        compare_workers=args.compare_workers,
        memory_budget_mb=args.memory_budget_mb,
        fingerprint_only=args.fingerprint_only,
        chromosome_filter_raw=args.chromosomes,
        plugins=parse_plugin_list(args.plugins),
        compare_only_plugins=args.compare_only_plugins,
    )
    artifacts = prepare_artifacts(config.output_dir)
    write_effective_config(config, artifacts)
    outputs = execute_engines(config, artifacts)
    return _run_comparison_pipeline(
        console=console,
        config=config,
        artifacts=artifacts,
        left_vcf=outputs.left_vcf,
        right_vcf=outputs.right_vcf,
        left_name=outputs.left_name,
        right_name=outputs.right_name,
    )


def _cmd_compare_existing(args: argparse.Namespace, console: Console) -> int:
    from .runtime import prepare_artifacts, resolve_runtime_config, write_effective_config

    if args.preset is None:
        raise ValueError("--preset is required or set VEPYR_DIFFLY_PRESET in .env")

    preset = get_preset(args.preset)
    if not preset.enabled:
        raise ValueError(f"preset {preset.name} is disabled")
    config = resolve_runtime_config(
        preset=preset,
        input_vcf=args.left_vcf,
        output_dir=args.output_dir,
        sample_first_n=None,
        execution_mode="compare-only",
        vepyr_path=None,
        vepyr_python=None,
        vep_cache_dir=None,
        vepyr_cache_output_dir=None,
        reference_fasta=None,
        vep_bin=None,
        vep_cache_version=None,
        vep_perl5lib=None,
        vepyr_use_fjall=False,
        compare_mode=args.compare_mode,
        compare_bucket_count=args.bucket_count,
        compare_workers=args.compare_workers,
        memory_budget_mb=args.memory_budget_mb,
        fingerprint_only=args.fingerprint_only,
        annotated_left_vcf=args.left_vcf,
        annotated_right_vcf=args.right_vcf,
        chromosome_filter_raw=args.chromosomes,
        plugins=parse_plugin_list(args.plugins),
        compare_only_plugins=args.compare_only_plugins,
    )
    artifacts = prepare_artifacts(config.output_dir)
    write_effective_config(config, artifacts)
    return _run_comparison_pipeline(
        console=console,
        config=config,
        artifacts=artifacts,
        left_vcf=args.left_vcf,
        right_vcf=args.right_vcf,
        left_name="VEP",
        right_name="vepyr",
    )


def _cmd_inspect_run(args: argparse.Namespace, console: Console) -> int:
    summary_path = args.run_dir / "summary.json"
    if not summary_path.exists():
        raise FileNotFoundError(f"summary file not found: {summary_path}")
    payload = json.loads(summary_path.read_text(encoding="utf-8"))
    console.print(json.dumps(payload, indent=2))
    return 0


def _cmd_analyze_consequence_mismatches(args: argparse.Namespace, console: Console) -> int:
    from .mismatch_analysis import analyze_consequence_mismatches, write_mismatch_analysis

    mismatches_tsv = args.mismatches_tsv
    if mismatches_tsv is None:
        if args.run_dir is None:
            raise ValueError("pass either --mismatches-tsv or --run-dir")
        mismatches_tsv = args.run_dir / "consequence_mismatches.tsv"
    if not mismatches_tsv.exists():
        raise FileNotFoundError(f"mismatch file not found: {mismatches_tsv}")

    payload = analyze_consequence_mismatches(mismatches_tsv)
    if args.output_json is not None:
        write_mismatch_analysis(args.output_json, payload)
    console.print(json.dumps(payload, indent=2))
    return 0


def _cmd_extract_mismatch_csq_examples(args: argparse.Namespace, console: Console) -> int:
    from .mismatch_analysis import (
        analyze_consequence_mismatches,
        extract_csq_examples_from_analysis,
        write_csq_examples,
    )

    analysis_json = args.analysis_json
    left_vcf = args.left_vcf
    right_vcf = args.right_vcf
    if args.run_dir is not None:
        if analysis_json is None:
            analysis_json = args.run_dir / "consequence_mismatch_analysis.json"
        if left_vcf is None:
            left_vcf = args.run_dir / "runtime" / "vep.annotated.vcf"
        if right_vcf is None:
            right_vcf = args.run_dir / "runtime" / "vepyr.annotated.vcf"
    if left_vcf is None or right_vcf is None:
        raise ValueError("left/right annotated VCF paths are required")

    if analysis_json is not None and analysis_json.exists():
        analysis = json.loads(analysis_json.read_text(encoding="utf-8"))
    elif args.run_dir is not None:
        analysis = analyze_consequence_mismatches(args.run_dir / "consequence_mismatches.tsv")
    else:
        raise ValueError("pass --run-dir or an existing --analysis-json")

    payload = extract_csq_examples_from_analysis(
        analysis,
        left_vcf=left_vcf,
        right_vcf=right_vcf,
        per_category_limit=args.per_category_limit,
    )
    write_csq_examples(args.output_json, payload)
    console.print(json.dumps(payload, indent=2))
    return 0


def main(argv: list[str] | None = None) -> int:
    load_repo_env()
    parser = build_parser()
    args = parser.parse_args(argv)
    console = Console()
    if args.command == "list-presets":
        return _cmd_list_presets(console)
    if args.command == "run":
        return _cmd_run(args, console)
    if args.command == "annotate-vepyr":
        return _cmd_annotate_vepyr(args, console)
    if args.command == "compare-existing":
        return _cmd_compare_existing(args, console)
    if args.command == "benchmark-compare":
        return _cmd_benchmark_compare(args, console)
    if args.command == "benchmark-plugins":
        return _cmd_benchmark_plugins(args, console)
    if args.command == "inspect-run":
        return _cmd_inspect_run(args, console)
    if args.command == "analyze-consequence-mismatches":
        return _cmd_analyze_consequence_mismatches(args, console)
    if args.command == "extract-mismatch-csq-examples":
        return _cmd_extract_mismatch_csq_examples(args, console)
    raise AssertionError(f"unhandled command {args.command}")


if __name__ == "__main__":
    raise SystemExit(main())
