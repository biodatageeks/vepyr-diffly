from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

from .models import AnnotatedOutputs, Preset, RunArtifacts, RuntimeConfig
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
    annotated_left_vcf: Path | None = None,
    annotated_right_vcf: Path | None = None,
) -> RuntimeConfig:
    if sample_first_n is not None and not preset.supports_sampling:
        raise ValueError(f"preset {preset.name} does not support sampling")
    return RuntimeConfig(
        preset=preset,
        input_vcf=input_vcf.expanduser().resolve(),
        output_dir=output_dir.resolve(),
        sample_first_n=sample_first_n,
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
    )


def prepare_input(config: RuntimeConfig, artifacts: RunArtifacts) -> Path:
    prepared = artifacts.runtime_dir / "prepared_input.vcf"
    stats = prepare_vcf_for_annotation(
        config.input_vcf,
        prepared,
        first_n=config.sample_first_n,
    )
    (artifacts.runtime_dir / "input_preparation.json").write_text(
        json.dumps(
            {
                "source_input_vcf": str(config.input_vcf),
                "requested_sample_first_n": config.sample_first_n,
                "sampled_records": stats.sampled_records,
                "prepared_output_records": stats.output_records,
                "split_source_records": stats.split_source_records,
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
        artifacts.runtime_dir / "prepared_input.vcf",
        artifacts.runtime_dir / "vep.annotated.vcf",
        artifacts.runtime_dir / "vepyr.annotated.vcf",
    ]
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
    completed = subprocess.run(command, capture_output=True, text=True, check=False, env=env)
    log_path.write_text(
        f"$ {' '.join(command)}\n\nSTDOUT\n{completed.stdout}\n\nSTDERR\n{completed.stderr}\n",
        encoding="utf-8",
    )
    if completed.returncode != 0:
        raise RuntimeError(f"command failed with exit code {completed.returncode}: {' '.join(command)}")


def _run_command_env(command: list[str], log_path: Path, env: dict[str, str]) -> None:
    completed = subprocess.run(command, capture_output=True, text=True, check=False, env=env)
    log_path.write_text(
        f"$ {' '.join(command)}\n\nSTDOUT\n{completed.stdout}\n\nSTDERR\n{completed.stderr}\n",
        encoding="utf-8",
    )
    if completed.returncode != 0:
        raise RuntimeError(f"command failed with exit code {completed.returncode}: {' '.join(command)}")


def _cache_method(preset: Preset) -> str:
    if preset.cache_flavor == "ensembl":
        return "vep"
    return preset.cache_flavor


def _resolve_local_cache_source(config: RuntimeConfig) -> Path:
    if config.vep_cache_dir is None or config.vep_cache_version is None:
        raise ValueError("vep cache dir and cache version are required")
    suffix = "" if _cache_method(config.preset) == "vep" else f"_{_cache_method(config.preset)}"
    return config.vep_cache_dir / config.preset.species / f"{config.vep_cache_version}_{config.preset.assembly}{suffix}"


def _resolve_vepyr_feature_root(config: RuntimeConfig) -> Path:
    if config.vepyr_cache_output_dir is None or config.vep_cache_version is None:
        raise ValueError("vepyr cache output dir and vep cache version are required")
    return (
        config.vepyr_cache_output_dir
        / "parquet"
        / f"{config.vep_cache_version}_{config.preset.assembly}_{_cache_method(config.preset)}"
    )


def _ensure_vepyr_local_ready(config: RuntimeConfig, artifacts: RunArtifacts) -> Path:
    if config.vepyr_python is not None:
        feature_root = _resolve_vepyr_feature_root(config)
        if not feature_root.exists():
            raise RuntimeError(f"configured vepyr cache feature root does not exist: {feature_root}")
        return feature_root
    if config.vepyr_path is None:
        raise ValueError("--vepyr-path is required for local execution")
    if config.vepyr_cache_output_dir is None:
        raise ValueError("--vepyr-cache-output-dir is required for local execution")

    install_log = artifacts.runtime_dir / "vepyr_install.log"
    _run_command(
        [sys.executable, "-m", "pip", "install", "--no-build-isolation", "-e", str(config.vepyr_path)],
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
        str(config.vep_bin),
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
    ]
    if config.reference_fasta is not None:
        vep_command.extend(["--fasta", str(config.reference_fasta)])

    vepyr_command = [
        str(config.vepyr_python or sys.executable),
        str(Path(__file__).with_name("vepyr_runner.py")),
        "--input-vcf",
        str(input_vcf),
        "--output-vcf",
        str(right_vcf),
        "--cache-dir",
        str(feature_root),
        "--reference-fasta",
        "" if config.reference_fasta is None else str(config.reference_fasta),
    ]
    _run_command_env(vep_command, artifacts.logs["vep"], vep_env)
    _run_command(vepyr_command, artifacts.logs["vepyr"])
    return AnnotatedOutputs(
        left_name="VEP",
        right_name="vepyr",
        left_vcf=left_vcf,
        right_vcf=right_vcf,
    )


def execute_engines(config: RuntimeConfig, artifacts: RunArtifacts) -> AnnotatedOutputs:
    if config.execution_mode == "local":
        return _execute_local(config, artifacts)
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
