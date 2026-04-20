from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from .chromosomes import parse_chromosome_selection
from .models import AnnotatedOutputs, Preset, RunArtifacts, RuntimeConfig
from .plugins import vep_plugin_args
from .sampling import prepare_vcf_for_annotation


def _expand_path(value: Path | None) -> Path | None:
    return None if value is None else value.expanduser()


def prepare_artifacts(output_dir: Path) -> RunArtifacts:
    runtime_dir = output_dir / "runtime"
    normalized_dir = output_dir / "normalized"
    runtime_dir.mkdir(parents=True, exist_ok=True)
    normalized_dir.mkdir(parents=True, exist_ok=True)
    return RunArtifacts(
        runtime_dir=runtime_dir,
        normalized_dir=normalized_dir,
        summary_json_path=output_dir / "summary.json",
        summary_md_path=output_dir / "summary.md",
        variant_diff_path=output_dir / "variant_diff.parquet",
        consequence_diff_path=output_dir / "consequence_diff.parquet",
        variant_mismatches_tsv_path=output_dir / "variant_mismatches.tsv",
        consequence_mismatches_tsv_path=output_dir / "consequence_mismatches.tsv",
        left_variant_path=normalized_dir / "left.variant.parquet",
        right_variant_path=normalized_dir / "right.variant.parquet",
        left_variant_bucket_dir=normalized_dir / "left.variant_buckets",
        right_variant_bucket_dir=normalized_dir / "right.variant_buckets",
        left_consequence_path=normalized_dir / "left.consequence.parquet",
        right_consequence_path=normalized_dir / "right.consequence.parquet",
        left_consequence_bucket_dir=normalized_dir / "left.consequence_buckets",
        right_consequence_bucket_dir=normalized_dir / "right.consequence_buckets",
        progress_log_path=runtime_dir / "compare.progress.log",
        logs={
            "vep": runtime_dir / "vep.log",
            "vepyr": runtime_dir / "vepyr.log",
        },
    )


def resolve_runtime_config(
    *,
    preset: Preset,
    input_vcf: Path,
    output_dir: Path,
    sample_first_n: int | None,
    execution_mode: str,
    vepyr_path: Path | None,
    vepyr_python: Path | None,
    vep_cache_dir: Path | None,
    vepyr_cache_output_dir: Path | None,
    reference_fasta: Path | None,
    vep_bin: Path | None,
    vep_cache_version: str | None,
    vep_perl5lib: str | None,
    vepyr_use_fjall: bool = False,
    compare_mode: str = "fast",
    compare_bucket_count: int | None = None,
    compare_workers: int | None = None,
    memory_budget_mb: int | None = None,
    fingerprint_only: bool = False,
    compare_only_plugins: bool = False,
    annotated_left_vcf: Path | None = None,
    annotated_right_vcf: Path | None = None,
    chromosome_filter_raw: str | None = None,
    plugins: list[str] | None = None,
) -> RuntimeConfig:
    if sample_first_n is not None and not preset.supports_sampling:
        raise ValueError(f"preset {preset.name} does not support sampling")
    selected_chromosomes, chromosome_aliases = parse_chromosome_selection(chromosome_filter_raw)
    return RuntimeConfig(
        preset=preset,
        input_vcf=input_vcf.expanduser().resolve(),
        output_dir=output_dir.resolve(),
        sample_first_n=sample_first_n,
        chromosome_filter_raw=chromosome_filter_raw,
        selected_chromosomes=selected_chromosomes,
        selected_chromosome_aliases=sorted(chromosome_aliases),
        plugins=[] if plugins is None else plugins,
        compare_only_plugins=compare_only_plugins,
        annotated_left_vcf=_expand_path(annotated_left_vcf),
        annotated_right_vcf=_expand_path(annotated_right_vcf),
        execution_mode=execution_mode,
        vepyr_path=None if vepyr_path is None else vepyr_path.expanduser().resolve(),
        vepyr_python=_expand_path(vepyr_python),
        vep_cache_dir=None if vep_cache_dir is None else vep_cache_dir.expanduser().resolve(),
        vepyr_cache_output_dir=(
            None
            if vepyr_cache_output_dir is None
            else vepyr_cache_output_dir.expanduser().resolve()
        ),
        reference_fasta=(
            None if reference_fasta is None else reference_fasta.expanduser().resolve()
        ),
        vep_bin=None if vep_bin is None else vep_bin.expanduser().resolve(),
        vep_cache_version=vep_cache_version,
        vep_perl5lib=vep_perl5lib,
        vepyr_use_fjall=vepyr_use_fjall,
        compare_mode=compare_mode,
        compare_bucket_count=compare_bucket_count,
        compare_workers=compare_workers,
        memory_budget_mb=memory_budget_mb,
        fingerprint_only=fingerprint_only,
    )


def prepare_input(config: RuntimeConfig, artifacts: RunArtifacts) -> Path:
    prepared = artifacts.runtime_dir / "prepared_input.vcf"
    stats = prepare_vcf_for_annotation(
        config.input_vcf,
        prepared,
        first_n=config.sample_first_n,
        chromosome_aliases=set(config.selected_chromosome_aliases),
    )
    (artifacts.runtime_dir / "input_preparation.json").write_text(
        json.dumps(
            {
                "source_input_vcf": str(config.input_vcf),
                "requested_sample_first_n": config.sample_first_n,
                "requested_chromosomes": config.chromosome_filter_raw,
                "effective_chromosomes": config.selected_chromosomes,
                "sampled_records": stats.sampled_records,
                "filtered_records": stats.filtered_records,
                "prepared_output_records": stats.output_records,
                "split_source_records": stats.split_source_records,
                "source_records_per_chromosome": stats.source_records_per_chromosome,
                "sampled_records_per_chromosome": stats.sampled_records_per_chromosome,
                "filtered_records_per_chromosome": stats.filtered_records_per_chromosome,
                "prepared_output_records_per_chromosome": stats.output_records_per_chromosome,
                "prepared_input_vcf": str(prepared),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return prepared


def remove_stale_runtime_outputs(artifacts: RunArtifacts) -> None:
    stale_paths = [
        artifacts.summary_json_path,
        artifacts.summary_md_path,
        artifacts.variant_diff_path,
        artifacts.consequence_diff_path,
        artifacts.variant_mismatches_tsv_path,
        artifacts.consequence_mismatches_tsv_path,
        artifacts.progress_log_path,
        artifacts.logs["vep"],
        artifacts.logs["vepyr"],
        artifacts.runtime_dir / "prepared_input.vcf",
        artifacts.runtime_dir / "sampled_input.vcf",
        artifacts.runtime_dir / "filtered_input.vcf",
        artifacts.runtime_dir / "input_preparation.json",
        artifacts.runtime_dir / "vep.annotated.vcf",
        artifacts.runtime_dir / "vepyr.annotated.vcf",
    ]
    stale_dirs = [
        artifacts.normalized_dir,
        artifacts.left_variant_bucket_dir,
        artifacts.right_variant_bucket_dir,
        artifacts.left_consequence_bucket_dir,
        artifacts.right_consequence_bucket_dir,
    ]
    for path in stale_dirs:
        if path.exists():
            shutil.rmtree(path)
    artifacts.normalized_dir.mkdir(parents=True, exist_ok=True)
    for path in stale_paths:
        if path.exists():
            path.unlink()


def write_effective_config(config: RuntimeConfig, artifacts: RunArtifacts) -> None:
    target = artifacts.runtime_dir / "effective_config.json"
    target.write_text(json.dumps(config.to_dict(), indent=2) + "\n", encoding="utf-8")


def _run_command(command: list[str], log_path: Path) -> None:
    env = os.environ.copy()
    local_bin = Path(sys.executable).resolve().parent
    env["PATH"] = f"{local_bin}:{env.get('PATH', '')}"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as handle:
        started_at = datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")
        handle.write(f"$ {' '.join(command)}\n\n")
        handle.write(f"STARTED {started_at}\n")
        handle.write("STREAMED OUTPUT\n")
        handle.flush()
        completed = subprocess.run(
            command,
            stdout=handle,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
            env=env,
        )
        ended_at = datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")
        handle.write(f"\nENDED {ended_at}\n")
        handle.write(f"EXIT {completed.returncode}\n")
    if completed.returncode != 0:
        raise RuntimeError(
            f"command failed with exit code {completed.returncode}: {' '.join(command)}"
        )


def _run_command_env(command: list[str], log_path: Path, env: dict[str, str]) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as handle:
        started_at = datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")
        handle.write(f"$ {' '.join(command)}\n\n")
        handle.write(f"STARTED {started_at}\n")
        handle.write("STREAMED OUTPUT\n")
        handle.flush()
        completed = subprocess.run(
            command,
            stdout=handle,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
            env=env,
        )
        ended_at = datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")
        handle.write(f"\nENDED {ended_at}\n")
        handle.write(f"EXIT {completed.returncode}\n")
    if completed.returncode != 0:
        raise RuntimeError(
            f"command failed with exit code {completed.returncode}: {' '.join(command)}"
        )


def run_vepyr_annotation(
    *,
    input_vcf: Path,
    output_vcf: Path,
    cache_dir: Path,
    log_path: Path,
    reference_fasta: Path | None = None,
    vepyr_python: Path | None = None,
    use_fjall: bool = False,
    plugins: list[str] | None = None,
) -> None:
    env = os.environ.copy()
    if vepyr_python is not None:
        env["VEPYR_DIFFLY_VEPYR_PYTHON"] = str(vepyr_python)
    command = [
        str(vepyr_python or sys.executable),
        str(Path(__file__).with_name("vepyr_runner.py")),
        "--input-vcf",
        str(input_vcf),
        "--output-vcf",
        str(output_vcf),
        "--cache-dir",
        str(cache_dir),
        "--reference-fasta",
        "" if reference_fasta is None else str(reference_fasta),
    ]
    if use_fjall:
        command.append("--use-fjall")
    if plugins:
        command.extend(["--plugins", ",".join(plugins)])
    _run_command_env(command, log_path, env)


def _cache_method(preset: Preset) -> str:
    if preset.cache_flavor == "ensembl":
        return "vep"
    return preset.cache_flavor


def _resolve_local_cache_source(config: RuntimeConfig) -> Path:
    if config.vep_cache_dir is None or config.vep_cache_version is None:
        raise ValueError("vep cache dir and cache version are required")
    suffix = "" if _cache_method(config.preset) == "vep" else f"_{_cache_method(config.preset)}"
    return (
        config.vep_cache_dir
        / config.preset.species
        / f"{config.vep_cache_version}_{config.preset.assembly}{suffix}"
    )


def _resolve_vepyr_feature_root(config: RuntimeConfig) -> Path:
    if config.vepyr_cache_output_dir is None or config.vep_cache_version is None:
        raise ValueError("vepyr cache output dir and vep cache version are required")
    current = (
        config.vepyr_cache_output_dir
        / f"{config.vep_cache_version}_{config.preset.assembly}_{_cache_method(config.preset)}"
    )
    if current.exists():
        return current
    return (
        config.vepyr_cache_output_dir
        / "parquet"
        / f"{config.vep_cache_version}_{config.preset.assembly}_{_cache_method(config.preset)}"
    )


def _ensure_vepyr_local_ready(config: RuntimeConfig, artifacts: RunArtifacts) -> Path:
    if config.vepyr_python is not None:
        feature_root = _resolve_vepyr_feature_root(config)
        if not feature_root.exists():
            raise RuntimeError(
                f"configured vepyr cache feature root does not exist: {feature_root}"
            )
        return feature_root
    if config.vepyr_path is None:
        raise ValueError("--vepyr-path is required for local execution")
    if config.vepyr_cache_output_dir is None:
        raise ValueError("--vepyr-cache-output-dir is required for local execution")

    install_log = artifacts.runtime_dir / "vepyr_install.log"
    _run_command(
        [
            sys.executable,
            "-m",
            "pip",
            "install",
            "--no-build-isolation",
            "-e",
            str(config.vepyr_path),
        ],
        install_log,
    )
    feature_root = _resolve_vepyr_feature_root(config)
    if feature_root.exists():
        return feature_root

    build_log = artifacts.runtime_dir / "vepyr_cache_build.log"
    source_cache = _resolve_local_cache_source(config)
    command = [
        sys.executable,
        "-c",
        (
            "import vepyr;"
            f"vepyr.build_cache(release={int(config.vep_cache_version)},"
            f" cache_dir={str(config.vepyr_cache_output_dir)!r},"
            f" species={config.preset.species!r},"
            f" assembly={config.preset.assembly!r},"
            f" method={_cache_method(config.preset)!r},"
            f" local_cache={str(source_cache)!r})"
        ),
    ]
    _run_command(command, build_log)
    if not feature_root.exists():
        raise RuntimeError(f"vepyr feature root was not created: {feature_root}")
    return feature_root


def _execute_local(config: RuntimeConfig, artifacts: RunArtifacts) -> AnnotatedOutputs:
    remove_stale_runtime_outputs(artifacts)
    input_vcf = prepare_input(config, artifacts)
    left_vcf = artifacts.runtime_dir / "vep.annotated.vcf"
    right_vcf = artifacts.runtime_dir / "vepyr.annotated.vcf"
    feature_root = _ensure_vepyr_local_ready(config, artifacts)

    if config.vep_bin is None or config.vep_cache_dir is None:
        raise ValueError("--vep-bin and --vep-cache-dir are required for local execution")

    vep_env = os.environ.copy()
    if config.vep_perl5lib:
        vep_env["PERL5LIB"] = config.vep_perl5lib

    vep_command = [
        *_vep_command_prefix(config.vep_bin),
        "--input_file",
        str(input_vcf),
        "--output_file",
        str(left_vcf),
        "--dir_cache",
        str(config.vep_cache_dir),
        "--species",
        config.preset.species,
        "--assembly",
        config.preset.assembly,
        "--cache_version",
        str(config.vep_cache_version or ""),
        "--force_overwrite",
        "--no_stats",
        *config.preset.vep_args,
        *vep_plugin_args(config.plugins),
    ]
    if config.reference_fasta is not None:
        vep_command.extend(["--fasta", str(config.reference_fasta)])

    right_vcf.parent.mkdir(parents=True, exist_ok=True)
    _run_command_env(vep_command, artifacts.logs["vep"], vep_env)
    run_vepyr_annotation(
        input_vcf=input_vcf,
        output_vcf=right_vcf,
        cache_dir=feature_root,
        log_path=artifacts.logs["vepyr"],
        reference_fasta=config.reference_fasta,
        vepyr_python=config.vepyr_python,
        use_fjall=config.vepyr_use_fjall,
        plugins=config.plugins,
    )
    return AnnotatedOutputs(
        left_name="VEP",
        right_name="vepyr",
        left_vcf=left_vcf,
        right_vcf=right_vcf,
    )


def _vep_command_prefix(vep_bin: Path) -> list[str]:
    if vep_bin.name == "vep":
        return ["perl", "-MBio::DB::HTS::Tabix", str(vep_bin)]
    return [str(vep_bin)]


def execute_engines(config: RuntimeConfig, artifacts: RunArtifacts) -> AnnotatedOutputs:
    if config.execution_mode == "local":
        return _execute_local(config, artifacts)
    if config.plugins:
        raise ValueError("plugin execution is currently supported only with --execution-mode local")
    input_vcf = prepare_input(config, artifacts)
    left_vcf = artifacts.runtime_dir / "vep.annotated.vcf"
    right_vcf = artifacts.runtime_dir / "vepyr.annotated.vcf"

    if config.vep_cache_dir is None:
        raise ValueError("--vep-cache-dir is required for runtime execution")
    if config.vepyr_path is None:
        raise ValueError("--vepyr-path is required for runtime execution")

    vep_command = [
        "docker",
        "run",
        "--rm",
        "-v",
        f"{input_vcf.parent}:/input",
        "-v",
        f"{left_vcf.parent}:/output",
        "-v",
        f"{config.vep_cache_dir}:/cache",
        config.vep_container_image,
        "vep",
        "--input_file",
        f"/input/{input_vcf.name}",
        "--output_file",
        f"/output/{left_vcf.name}",
        "--dir_cache",
        "/cache",
        *config.preset.vep_args,
    ]
    vepyr_command = [
        "docker",
        "run",
        "--rm",
        "-v",
        f"{input_vcf.parent}:/input",
        "-v",
        f"{right_vcf.parent}:/output",
        "-v",
        f"{config.vepyr_path}:/workspace/vepyr",
        "-v",
        f"{config.vep_cache_dir}:/cache",
        config.vepyr_container_image,
        "--input-vcf",
        f"/input/{input_vcf.name}",
        "--output-vcf",
        f"/output/{right_vcf.name}",
        "--cache-dir",
        "/cache",
        *config.preset.vepyr_args,
    ]
    _run_command(vep_command, artifacts.logs["vep"])
    _run_command(vepyr_command, artifacts.logs["vepyr"])
    return AnnotatedOutputs(
        left_name="VEP",
        right_name="vepyr",
        left_vcf=left_vcf,
        right_vcf=right_vcf,
    )
