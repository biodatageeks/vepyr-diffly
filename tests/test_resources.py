from pathlib import Path

from vepyr_diffly.models import Preset, RuntimeConfig
from vepyr_diffly.resources import plan_compare_resources


def _config(tmp_path: Path, *, memory_budget_mb: int) -> RuntimeConfig:
    preset = Preset(
        name="ensembl_everything",
        enabled=True,
        description="test",
        species="homo_sapiens",
        assembly="GRCh38",
        cache_flavor="ensembl",
        normalization_policy="csq_semantic_v1",
        supports_full_run=True,
        supports_sampling=True,
        vep_args=[],
        vepyr_args=[],
    )
    return RuntimeConfig(
        preset=preset,
        input_vcf=tmp_path / "left.vcf",
        output_dir=tmp_path / "run",
        sample_first_n=None,
        memory_budget_mb=memory_budget_mb,
    )


def test_plan_compare_resources_disables_parallel_sides_for_large_inputs(tmp_path: Path) -> None:
    left = tmp_path / "left.vcf"
    right = tmp_path / "right.vcf"
    left.write_bytes(b"x" * (16 * 1024 * 1024))
    right.write_bytes(b"x" * (16 * 1024 * 1024))

    plan = plan_compare_resources(
        config=_config(tmp_path, memory_budget_mb=512),
        left_vcf=left,
        right_vcf=right,
    )

    assert plan.parallelize_sides is False
    assert plan.use_bucketed_variant is True
    assert plan.use_bucketed_consequence is True
    assert plan.compare_workers >= 1
