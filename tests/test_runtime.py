from __future__ import annotations

from pathlib import Path

from vepyr_diffly.presets import get_preset
from vepyr_diffly.runtime import (
    _vep_command_prefix,
    prepare_artifacts,
    remove_stale_runtime_outputs,
    run_vepyr_annotation,
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
        chromosome_filter_raw="1,chr2",
        plugins=["clinvar"],
        compare_only_plugins=True,
    )

    assert config.vepyr_python == symlink_python
    assert config.vepyr_use_fjall is True
    assert config.compare_mode == "fast"
    assert config.memory_budget_mb == 768
    assert config.fingerprint_only is False
    assert config.plugins == ["clinvar"]
    assert config.compare_only_plugins is True
    assert config.selected_chromosomes == ["1", "2"]
    assert "chr1" in config.selected_chromosome_aliases


def test_run_vepyr_annotation_passes_plugins_to_runner(monkeypatch, tmp_path: Path) -> None:
    captured: dict[str, object] = {}

    def fake_run_command_env(command, log_path, env):
        captured["command"] = command
        captured["log_path"] = log_path
        captured["env"] = env

    monkeypatch.setattr("vepyr_diffly.runtime._run_command_env", fake_run_command_env)

    run_vepyr_annotation(
        input_vcf=tmp_path / "input.vcf",
        output_vcf=tmp_path / "out.vcf",
        cache_dir=tmp_path / "cache",
        log_path=tmp_path / "vepyr.log",
        plugins=["clinvar", "cadd"],
    )

    command = captured["command"]
    assert "--plugins" in command
    assert command[command.index("--plugins") + 1] == "clinvar,cadd"


def test_vep_command_prefix_preloads_tabix_for_vep_script() -> None:
    assert _vep_command_prefix(Path("/opt/ensembl-vep/vep")) == [
        "perl",
        "-MBio::DB::HTS::Tabix",
        "/opt/ensembl-vep/vep",
    ]
    assert _vep_command_prefix(Path("/usr/local/bin/custom-vep")) == [
        "/usr/local/bin/custom-vep"
    ]


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
