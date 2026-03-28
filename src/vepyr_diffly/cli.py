from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import subprocess
import sys

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

    compare_parser = subparsers.add_parser("compare-existing")
    compare_parser.add_argument("--preset", default=env_str("VEPYR_DIFFLY_PRESET"))
    compare_parser.add_argument("--left-vcf", required=True, type=Path)
    compare_parser.add_argument("--right-vcf", required=True, type=Path)
    compare_parser.add_argument("--output-dir", required=True, type=Path)

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
    from .compare import compare_bucketed_consequence_tier, compare_tier
    from .normalize import (
        STREAMING_CONSEQUENCE_THRESHOLD_BYTES,
        VARIANT_KEY,
        materialize_consequence_summary,
        materialize_variant_summary,
        recommend_consequence_bucket_count,
    )
    from .progress import ProgressReporter
    from .report import print_run_summary, write_run_summary
    from .schema import validate_variant_schema
    import polars as pl

    reporter = ProgressReporter(log_path=artifacts.progress_log_path, console=console)
    reporter.start()
    try:
        reporter.stage("comparison: starting normalization")
        with ThreadPoolExecutor(max_workers=2) as executor:
            left_future = executor.submit(
                materialize_variant_summary,
                vcf_path=left_vcf,
                variant_path=artifacts.left_variant_path,
                reporter=reporter,
                side_label=left_name,
            )
            right_future = executor.submit(
                materialize_variant_summary,
                vcf_path=right_vcf,
                variant_path=artifacts.right_variant_path,
                reporter=reporter,
                side_label=right_name,
            )
            left_csq_fields = left_future.result()
            right_csq_fields = right_future.result()
        if left_csq_fields != right_csq_fields:
            raise ValueError("left/right CSQ headers differ and cannot be compared semantically")

        use_bucketed_consequence = (
            max(left_vcf.stat().st_size, right_vcf.stat().st_size)
            >= STREAMING_CONSEQUENCE_THRESHOLD_BYTES
        )
        if use_bucketed_consequence:
            bucket_count = recommend_consequence_bucket_count(left_vcf, right_vcf)
            reporter.stage("comparison: bucketizing consequence rows")
            reporter.log(f"comparison: using {bucket_count} consequence buckets")
            env = os.environ.copy()
            env["PYTHONPATH"] = str(Path(__file__).resolve().parents[1]) + (
                os.pathsep + env["PYTHONPATH"] if "PYTHONPATH" in env else ""
            )
            left_proc = subprocess.Popen(
                [
                    sys.executable,
                    "-m",
                    "vepyr_diffly.worker",
                    "bucketize-side",
                    "--vcf",
                    str(left_vcf),
                    "--bucket-root",
                    str(artifacts.left_consequence_bucket_dir),
                    "--side-label",
                    left_name,
                    "--csq-fields-json",
                    json.dumps(left_csq_fields),
                    "--bucket-count",
                    str(bucket_count),
                    "--progress-log",
                    str(artifacts.progress_log_path),
                ],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
                text=True,
                env=env,
            )
            right_proc = subprocess.Popen(
                [
                    sys.executable,
                    "-m",
                    "vepyr_diffly.worker",
                    "bucketize-side",
                    "--vcf",
                    str(right_vcf),
                    "--bucket-root",
                    str(artifacts.right_consequence_bucket_dir),
                    "--side-label",
                    right_name,
                    "--csq-fields-json",
                    json.dumps(right_csq_fields),
                    "--bucket-count",
                    str(bucket_count),
                    "--progress-log",
                    str(artifacts.progress_log_path),
                ],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
                text=True,
                env=env,
            )
            left_stderr = left_proc.communicate()[1]
            right_stderr = right_proc.communicate()[1]
            if left_proc.returncode != 0:
                raise RuntimeError(
                    f"left consequence bucketization failed with exit code {left_proc.returncode}: {left_stderr}"
                )
            if right_proc.returncode != 0:
                raise RuntimeError(
                    f"right consequence bucketization failed with exit code {right_proc.returncode}: {right_stderr}"
                )
        else:
            materialize_consequence_summary(
                vcf_path=left_vcf,
                consequence_path=artifacts.left_consequence_path,
                csq_fields=left_csq_fields,
                reporter=reporter,
                side_label=left_name,
            )
            materialize_consequence_summary(
                vcf_path=right_vcf,
                consequence_path=artifacts.right_consequence_path,
                csq_fields=right_csq_fields,
                reporter=reporter,
                side_label=right_name,
            )

        reporter.stage("comparison: validating normalized variant schemas")
        left_variant = pl.scan_parquet(artifacts.left_variant_path).head(1000).collect()
        right_variant = pl.scan_parquet(artifacts.right_variant_path).head(1000).collect()
        validate_variant_schema(left_variant)
        validate_variant_schema(right_variant)

        reporter.stage("comparison: running variant tier diff")
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

        reporter.stage("comparison: running consequence tier diff")
        if use_bucketed_consequence:
            consequence = compare_bucketed_consequence_tier(
                left_bucket_dir=artifacts.left_consequence_bucket_dir,
                right_bucket_dir=artifacts.right_consequence_bucket_dir,
                csq_fields=left_csq_fields,
                left_name=left_name,
                right_name=right_name,
                diff_frame_path=artifacts.consequence_diff_path,
                mismatches_tsv_path=artifacts.consequence_mismatches_tsv_path,
                reporter=reporter,
                bucket_count=bucket_count,
            ).tier
        else:
            consequence_key = VARIANT_KEY + left_csq_fields
            consequence = compare_tier(
                name="consequence",
                left=pl.scan_parquet(artifacts.left_consequence_path),
                right=pl.scan_parquet(artifacts.right_consequence_path),
                primary_key=consequence_key,
                left_name=left_name,
                right_name=right_name,
                diff_frame_path=artifacts.consequence_diff_path,
                mismatches_tsv_path=artifacts.consequence_mismatches_tsv_path,
                reporter=reporter,
            ).tier

        reporter.stage("comparison: writing summaries")
        print_run_summary(
            console=console,
            config=config,
            variant=variant,
            consequence=consequence,
            left_vcf=left_vcf,
            right_vcf=right_vcf,
            progress_log_path=artifacts.progress_log_path,
        )
        write_run_summary(
            config=config,
            artifacts=artifacts,
            variant=variant,
            consequence=consequence,
            left_vcf=left_vcf,
            right_vcf=right_vcf,
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
        annotated_left_vcf=args.left_vcf,
        annotated_right_vcf=args.right_vcf,
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


def main(argv: list[str] | None = None) -> int:
    load_repo_env()
    parser = build_parser()
    args = parser.parse_args(argv)
    console = Console()
    if args.command == "list-presets":
        return _cmd_list_presets(console)
    if args.command == "run":
        return _cmd_run(args, console)
    if args.command == "compare-existing":
        return _cmd_compare_existing(args, console)
    if args.command == "inspect-run":
        return _cmd_inspect_run(args, console)
    raise AssertionError(f"unhandled command {args.command}")


if __name__ == "__main__":
    raise SystemExit(main())
