from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class Preset:
    name: str
    enabled: bool
    description: str
    species: str
    assembly: str
    cache_flavor: str
    normalization_policy: str
    supports_full_run: bool
    supports_sampling: bool
    vep_args: list[str]
    vepyr_args: list[str]


@dataclass(frozen=True)
class RuntimeConfig:
    preset: Preset
    input_vcf: Path
    output_dir: Path
    sample_first_n: int | None
    annotated_left_vcf: Path | None = None
    annotated_right_vcf: Path | None = None
    vepyr_path: Path | None = None
    vep_cache_dir: Path | None = None
    vepyr_cache_output_dir: Path | None = None
    reference_fasta: Path | None = None
    vepyr_python: Path | None = None
    vep_bin: Path | None = None
    vep_perl5lib: str | None = None
    vep_cache_version: str | None = None
    execution_mode: str = "docker"
    vep_container_image: str = "vepyr-diffly-vep:local"
    vepyr_container_image: str = "vepyr-diffly-vepyr:local"

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["input_vcf"] = str(self.input_vcf)
        payload["output_dir"] = str(self.output_dir)
        payload["vepyr_path"] = None if self.vepyr_path is None else str(self.vepyr_path)
        payload["annotated_left_vcf"] = (
            None if self.annotated_left_vcf is None else str(self.annotated_left_vcf)
        )
        payload["annotated_right_vcf"] = (
            None if self.annotated_right_vcf is None else str(self.annotated_right_vcf)
        )
        payload["vep_cache_dir"] = (
            None if self.vep_cache_dir is None else str(self.vep_cache_dir)
        )
        payload["vepyr_cache_output_dir"] = (
            None if self.vepyr_cache_output_dir is None else str(self.vepyr_cache_output_dir)
        )
        payload["reference_fasta"] = (
            None if self.reference_fasta is None else str(self.reference_fasta)
        )
        payload["vepyr_python"] = None if self.vepyr_python is None else str(self.vepyr_python)
        payload["vep_bin"] = None if self.vep_bin is None else str(self.vep_bin)
        return payload


@dataclass(frozen=True)
class AnnotatedOutputs:
    left_name: str
    right_name: str
    left_vcf: Path
    right_vcf: Path


@dataclass(frozen=True)
class TierResult:
    name: str
    summary: str
    equal: bool
    left_only_rows: int
    right_only_rows: int
    unequal_rows: int
    joined_equal_rows: int
    diff_frame_path: Path
    mismatches_tsv_path: Path


@dataclass
class RunArtifacts:
    runtime_dir: Path
    normalized_dir: Path
    summary_json_path: Path
    summary_md_path: Path
    variant_diff_path: Path
    consequence_diff_path: Path
    variant_mismatches_tsv_path: Path
    consequence_mismatches_tsv_path: Path
    left_variant_path: Path
    right_variant_path: Path
    left_consequence_path: Path
    right_consequence_path: Path
    left_consequence_bucket_dir: Path
    right_consequence_bucket_dir: Path
    progress_log_path: Path
    logs: dict[str, Path] = field(default_factory=dict)
