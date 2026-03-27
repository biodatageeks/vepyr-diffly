from pathlib import Path

from vepyr_diffly.models import Preset, RunArtifacts, RuntimeConfig, TierResult
from vepyr_diffly.report import write_run_summary


def test_write_run_summary_creates_json_and_markdown(tmp_path: Path) -> None:
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
    config = RuntimeConfig(
        preset=preset,
        input_vcf=tmp_path / "input.vcf",
        output_dir=tmp_path,
        sample_first_n=1000,
        vepyr_path=tmp_path / "vepyr",
        vep_cache_dir=tmp_path / "cache",
        reference_fasta=None,
    )
    artifacts = RunArtifacts(
        runtime_dir=tmp_path / "runtime",
        summary_json_path=tmp_path / "summary.json",
        summary_md_path=tmp_path / "summary.md",
        variant_diff_path=tmp_path / "variant.parquet",
        consequence_diff_path=tmp_path / "consequence.parquet",
        variant_mismatches_tsv_path=tmp_path / "variant.tsv",
        consequence_mismatches_tsv_path=tmp_path / "consequence.tsv",
    )
    tier = TierResult(
        name="variant",
        summary="summary",
        equal=False,
        left_only_rows=1,
        right_only_rows=2,
        unequal_rows=3,
        joined_equal_rows=4,
        diff_frame_path=artifacts.variant_diff_path,
        mismatches_tsv_path=artifacts.variant_mismatches_tsv_path,
    )

    write_run_summary(config=config, artifacts=artifacts, variant=tier, consequence=tier)

    assert artifacts.summary_json_path.exists()
    assert artifacts.summary_md_path.exists()
    assert "ensembl_everything" in artifacts.summary_md_path.read_text(encoding="utf-8")

