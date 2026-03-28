from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class PreparedInputStats:
    source_records: int
    sampled_records: int
    output_records: int
    split_source_records: int


def sample_vcf_first_n(input_vcf: Path, output_vcf: Path, first_n: int) -> int:
    if first_n < 0:
        raise ValueError("first_n must be >= 0")

    output_vcf.parent.mkdir(parents=True, exist_ok=True)
    kept_rows = 0
    with input_vcf.open(encoding="utf-8") as src, output_vcf.open("w", encoding="utf-8") as dst:
        for line in src:
            if line.startswith("#"):
                dst.write(line)
                continue
            if kept_rows >= first_n:
                break
            dst.write(line)
            kept_rows += 1
    return kept_rows


def decompose_multiallelic_vcf(input_vcf: Path, output_vcf: Path) -> PreparedInputStats:
    output_vcf.parent.mkdir(parents=True, exist_ok=True)
    source_records = 0
    output_records = 0
    split_source_records = 0
    with input_vcf.open(encoding="utf-8") as src, output_vcf.open("w", encoding="utf-8") as dst:
        for line in src:
            if line.startswith("#"):
                dst.write(line)
                continue
            source_records += 1
            row = line.rstrip("\n")
            if not row:
                continue
            columns = row.split("\t")
            if len(columns) < 8:
                dst.write(line)
                output_records += 1
                continue
            alts = columns[4].split(",")
            if len(alts) <= 1:
                dst.write(line)
                output_records += 1
                continue
            split_source_records += 1
            for alt in alts:
                cloned = list(columns)
                cloned[4] = alt
                dst.write("\t".join(cloned) + "\n")
                output_records += 1
    return PreparedInputStats(
        source_records=source_records,
        sampled_records=source_records,
        output_records=output_records,
        split_source_records=split_source_records,
    )


def prepare_vcf_for_annotation(
    input_vcf: Path,
    output_vcf: Path,
    *,
    first_n: int | None = None,
) -> PreparedInputStats:
    staged_input = input_vcf
    sampled_records = 0
    if first_n is not None:
        staged_input = output_vcf.parent / "sampled_input.vcf"
        sampled_records = sample_vcf_first_n(input_vcf, staged_input, first_n)
    stats = decompose_multiallelic_vcf(staged_input, output_vcf)
    return PreparedInputStats(
        source_records=stats.source_records if first_n is None else sampled_records,
        sampled_records=stats.sampled_records if first_n is None else sampled_records,
        output_records=stats.output_records,
        split_source_records=stats.split_source_records,
    )
