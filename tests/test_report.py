from pathlib import Path
import json

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
        chromosome_filter_raw="1,2",
        selected_chromosomes=["1", "2"],
        selected_chromosome_aliases=["1", "2", "chr1", "chr2"],
        vepyr_path=tmp_path / "vepyr",
        vep_cache_dir=tmp_path / "cache",
        reference_fasta=None,
        memory_budget_mb=512,
    )
    artifacts = RunArtifacts(
        runtime_dir=tmp_path / "runtime",
        normalized_dir=tmp_path / "normalized",
        summary_json_path=tmp_path / "summary.json",
        summary_md_path=tmp_path / "summary.md",
        variant_diff_path=tmp_path / "variant.parquet",
        consequence_diff_path=tmp_path / "consequence.parquet",
        variant_mismatches_tsv_path=tmp_path / "variant.tsv",
        consequence_mismatches_tsv_path=tmp_path / "consequence.tsv",
        left_variant_path=tmp_path / "normalized" / "left.variant.parquet",
        right_variant_path=tmp_path / "normalized" / "right.variant.parquet",
        left_variant_bucket_dir=tmp_path / "normalized" / "left.variant_buckets",
        right_variant_bucket_dir=tmp_path / "normalized" / "right.variant_buckets",
        left_consequence_path=tmp_path / "normalized" / "left.consequence.parquet",
        right_consequence_path=tmp_path / "normalized" / "right.consequence.parquet",
        left_consequence_bucket_dir=tmp_path / "normalized" / "left.consequence_buckets",
        right_consequence_bucket_dir=tmp_path / "normalized" / "right.consequence_buckets",
        progress_log_path=tmp_path / "runtime" / "compare.progress.log",
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
        per_chromosome={
            "1": {
                "left_rows": 1,
                "right_rows": 1,
                "joined_equal_rows": 1,
                "left_only_rows": 0,
                "right_only_rows": 0,
                "unequal_rows": 0,
            }
        },
    )

    write_run_summary(
        config=config,
        artifacts=artifacts,
        variant=tier,
        consequence=tier,
        left_vcf=tmp_path / "left.vcf",
        right_vcf=tmp_path / "right.vcf",
        resource_plan={"memory_budget_mb": 512, "bucket_count": 8},
        timings={"variant_diff_seconds": 1.25},
        chromosome_summary={
            "requested": ["1", "2"],
            "effective": ["1", "2"],
            "per_chromosome": {
                "1": {
                    "variant": {"equal": True, "left_only_rows": 0, "right_only_rows": 0, "unequal_rows": 0},
                    "consequence": {"equal": True, "left_only_rows": 0, "right_only_rows": 0, "unequal_rows": 0},
                    "timings": {"variant_diff_seconds": 0.5},
                }
            },
        },
    )

    assert artifacts.summary_json_path.exists()
    assert artifacts.summary_md_path.exists()
    assert "ensembl_everything" in artifacts.summary_md_path.read_text(encoding="utf-8")
    payload = json.loads(artifacts.summary_json_path.read_text(encoding="utf-8"))
    assert payload["compare_mode"] == "fast"
    assert payload["memory_budget_mb"] == 512
    assert payload["timings"]["variant_diff_seconds"] == 1.25
    assert payload["effective_chromosomes"] == ["1", "2"]
    assert payload["chromosomes"]["per_chromosome"]["1"]["timings"]["variant_diff_seconds"] == 0.5
