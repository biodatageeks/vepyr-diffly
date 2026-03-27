from __future__ import annotations

import argparse
import json
from pathlib import Path

from .presets import get_preset, load_presets
from .settings import env_int, env_path, env_str, load_repo_env

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
    run_parser.add_argument("--execution-mode", default=env_str("VEPYR_DIFFLY_EXECUTION_MODE") or "docker")
    run_parser.add_argument("--input-vcf", type=Path, default=env_path("VEPYR_DIFFLY_INPUT_VCF"))
    run_parser.add_argument("--output-dir", type=Path, default=env_path("VEPYR_DIFFLY_OUTPUT_DIR"))
    run_parser.add_argument("--sample-first-n", type=int, default=env_int("VEPYR_DIFFLY_SAMPLE_FIRST_N"))
    run_parser.add_argument("--vepyr-path", type=Path, default=env_path("VEPYR_DIFFLY_VEPYR_PATH"))
    run_parser.add_argument("--vepyr-python", type=Path, default=env_path("VEPYR_DIFFLY_VEPYR_PYTHON"))
    run_parser.add_argument("--vep-cache-dir", type=Path, default=env_path("VEPYR_DIFFLY_VEP_CACHE_DIR"))
    run_parser.add_argument("--vepyr-cache-output-dir", type=Path, default=env_path("VEPYR_DIFFLY_VEPYR_CACHE_OUTPUT_DIR"))
    run_parser.add_argument("--reference-fasta", type=Path, default=env_path("VEPYR_DIFFLY_REFERENCE_FASTA"))
    run_parser.add_argument("--vep-bin", type=Path, default=env_path("VEPYR_DIFFLY_VEP_BIN"))
    run_parser.add_argument("--vep-cache-version", default=env_str("VEPYR_DIFFLY_VEP_CACHE_VERSION"))
    run_parser.add_argument("--vep-perl5lib", default=env_str("VEPYR_DIFFLY_VEP_PERL5LIB"))

    subparsers.add_parser("list-presets")

    inspect_parser = subparsers.add_parser("inspect-run")
    inspect_parser.add_argument("--run-dir", required=True, type=Path)
    return parser


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


def _cmd_run(args: argparse.Namespace, console: Console) -> int:
    from .compare import compare_tier
    from .normalize import VARIANT_KEY, normalize_annotated_vcf
    from .report import print_run_summary, write_run_summary
    from .runtime import (
        execute_engines,
        prepare_artifacts,
        resolve_runtime_config,
        write_effective_config,
    )
    from .schema import validate_variant_schema

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
    )
    artifacts = prepare_artifacts(config.output_dir)
    write_effective_config(config, artifacts)
    outputs = execute_engines(config, artifacts)

    left = normalize_annotated_vcf(outputs.left_vcf)
    right = normalize_annotated_vcf(outputs.right_vcf)
    validate_variant_schema(left.variant)
    validate_variant_schema(right.variant)

    variant = compare_tier(
        name="variant",
        left=left.variant,
        right=right.variant,
        primary_key=VARIANT_KEY,
        left_name=outputs.left_name,
        right_name=outputs.right_name,
        diff_frame_path=artifacts.variant_diff_path,
        mismatches_tsv_path=artifacts.variant_mismatches_tsv_path,
    ).tier

    consequence_key = VARIANT_KEY + left.csq_fields
    consequence = compare_tier(
        name="consequence",
        left=left.consequence,
        right=right.consequence,
        primary_key=consequence_key,
        left_name=outputs.left_name,
        right_name=outputs.right_name,
        diff_frame_path=artifacts.consequence_diff_path,
        mismatches_tsv_path=artifacts.consequence_mismatches_tsv_path,
    ).tier

    print_run_summary(console=console, config=config, variant=variant, consequence=consequence)
    write_run_summary(config=config, artifacts=artifacts, variant=variant, consequence=consequence)
    return 0


def _cmd_inspect_run(args: argparse.Namespace, console: Console) -> int:
    summary_path = args.run_dir / "summary.json"
    if not summary_path.exists():
        raise FileNotFoundError(f"summary file not found: {summary_path}")
    payload = json.loads(summary_path.read_text(encoding="utf-8"))
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
    if args.command == "inspect-run":
        return _cmd_inspect_run(args, console)
    raise AssertionError(f"unhandled command {args.command}")


if __name__ == "__main__":
    raise SystemExit(main())
