from __future__ import annotations

from pathlib import Path

from vepyr_diffly.presets import get_preset
from vepyr_diffly.runtime import (
    prepare_artifacts,
    remove_stale_runtime_outputs,
    resolve_runtime_config,
)


def test_resolve_runtime_config_preserves_vepyr_python_symlink(tmp_path: Path) -> None:
    real_python = tmp_path / "python3.12"
    real_python.write_text("", encoding="utf-8")
    venv_bin = tmp_path / ".vepyr" / "bin"
    venv_bin.mkdir(parents=True)
    symlink_python = venv_bin / "python"
    symlink_python.symlink_to(real_python)

    config = resolve_runtime_config(
        preset=get_preset("ensembl_everything"),
        input_vcf=tmp_path / "input.vcf",
        output_dir=tmp_path / "runs",
        sample_first_n=1000,
        execution_mode="local",
        vepyr_path=tmp_path / "vepyr",
        vepyr_python=symlink_python,
        vep_cache_dir=tmp_path / ".vep",
        vepyr_cache_output_dir=tmp_path / "vepyr_cache",
        reference_fasta=tmp_path / "reference.fa",
        vep_bin=tmp_path / "vep",
        vep_cache_version="115",
        vep_perl5lib=None,
        vepyr_use_fjall=True,
        memory_budget_mb=768,
    )

    assert config.vepyr_python == symlink_python
    assert config.vepyr_use_fjall is True
    assert config.compare_mode == "fast"
    assert config.memory_budget_mb == 768
    assert config.fingerprint_only is False


def test_remove_stale_runtime_outputs_removes_old_vcfs(tmp_path: Path) -> None:
    artifacts = prepare_artifacts(tmp_path / "run")
    prepared = artifacts.runtime_dir / "prepared_input.vcf"
    sampled = artifacts.runtime_dir / "sampled_input.vcf"
    left = artifacts.runtime_dir / "vep.annotated.vcf"
    right = artifacts.runtime_dir / "vepyr.annotated.vcf"
    summary = artifacts.summary_json_path
    diff = artifacts.variant_diff_path
    progress = artifacts.progress_log_path
    normalized = artifacts.left_variant_path

    for path in (prepared, sampled, left, right, summary, diff, progress, normalized):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("stale", encoding="utf-8")

    remove_stale_runtime_outputs(artifacts)

    assert not prepared.exists()
    assert not sampled.exists()
    assert not left.exists()
    assert not right.exists()
    assert not summary.exists()
    assert not diff.exists()
    assert not progress.exists()
    assert not normalized.exists()
    assert artifacts.normalized_dir.exists()
