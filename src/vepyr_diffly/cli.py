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
    vepyr_only_parser.add_argument("--log-path", type=Path)

    compare_parser = subparsers.add_parser("compare-existing")
    compare_parser.add_argument("--preset", default=env_str("VEPYR_DIFFLY_PRESET"))
    compare_parser.add_argument("--left-vcf", required=True, type=Path)
    compare_parser.add_argument("--right-vcf", required=True, type=Path)
    compare_parser.add_argument("--output-dir", required=True, type=Path)
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


def _cmd_annotate_vepyr(args: argparse.Namespace, console: Console) -> int:
    from .runtime import run_vepyr_annotation

    cache_dir = args.cache_dir
    if cache_dir is None:
        raise ValueError("--cache-dir is required")
    # Allow passing either the root cache dir or the parquet feature root already.
    if cache_dir.name != "115_GRCh38_vep" and (cache_dir / "parquet").exists():
        cache_dir = cache_dir / "parquet" / "115_GRCh38_vep"
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
        if left_csq_fields != right_csq_fields:
            raise ValueError("left/right CSQ headers differ and cannot be compared semantically")

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
                    csq_fields=left_csq_fields,
                    bucket_root=artifacts.left_consequence_bucket_dir,
                    reporter=reporter,
                    side_label=left_name,
                    bucket_count=resource_plan.bucket_count,
                    chunk_variants=resource_plan.consequence_chunk_rows,
                    total_variants=left_variant_rows,
                ),
                dict(
                    vcf_path=right_vcf,
                    csq_fields=right_csq_fields,
                    bucket_root=artifacts.right_consequence_bucket_dir,
                    reporter=reporter,
                    side_label=right_name,
                    bucket_count=resource_plan.bucket_count,
                    chunk_variants=resource_plan.consequence_chunk_rows,
                    total_variants=right_variant_rows,
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
                    csq_fields=left_csq_fields,
                    reporter=reporter,
                    side_label=left_name,
                ),
                dict(
                    vcf_path=right_vcf,
                    consequence_path=artifacts.right_consequence_path,
                    csq_fields=right_csq_fields,
                    reporter=reporter,
                    side_label=right_name,
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
                csq_fields=left_csq_fields,
                left_name=left_name,
                right_name=right_name,
                diff_frame_path=artifacts.consequence_diff_path,
                mismatches_tsv_path=artifacts.consequence_mismatches_tsv_path,
                reporter=reporter,
                bucket_count=resource_plan.bucket_count,
                max_workers=resource_plan.compare_workers,
                compare_mode=config.compare_mode,
                fingerprint_only=config.fingerprint_only,
            ).tier
            timings["consequence_diff_seconds"] = round(perf_counter() - stage_start, 3)
        else:
            consequence_key = VARIANT_KEY + left_csq_fields
            stage_start = perf_counter()
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
            timings["consequence_diff_seconds"] = round(perf_counter() - stage_start, 3)

        reporter.stage("comparison: writing summaries")
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
    if args.command == "annotate-vepyr":
        return _cmd_annotate_vepyr(args, console)
    if args.command == "compare-existing":
        return _cmd_compare_existing(args, console)
    if args.command == "benchmark-compare":
        return _cmd_benchmark_compare(args, console)
    if args.command == "inspect-run":
        return _cmd_inspect_run(args, console)
    raise AssertionError(f"unhandled command {args.command}")


if __name__ == "__main__":
    raise SystemExit(main())
