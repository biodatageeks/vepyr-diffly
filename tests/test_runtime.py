from __future__ import annotations

from pathlib import Path

from vepyr_diffly.presets import get_preset
from vepyr_diffly.runtime import resolve_runtime_config


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
    )

    assert config.vepyr_python == symlink_python
