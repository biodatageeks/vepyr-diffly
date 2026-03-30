from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .chromosomes import normalize_chromosome_name


@dataclass(frozen=True)
class PreparedInputStats:
    source_records: int
    sampled_records: int
    filtered_records: int
    output_records: int
    split_source_records: int
    source_records_per_chromosome: dict[str, int]
    sampled_records_per_chromosome: dict[str, int]
    filtered_records_per_chromosome: dict[str, int]
    output_records_per_chromosome: dict[str, int]


def sample_vcf_first_n(input_vcf: Path, output_vcf: Path, first_n: int) -> int:
    kept_rows, _ = _sample_vcf_first_n_with_counts(input_vcf, output_vcf, first_n)
    return kept_rows


def _sample_vcf_first_n_with_counts(
    input_vcf: Path,
    output_vcf: Path,
    first_n: int,
) -> tuple[int, dict[str, int]]:
    if first_n < 0:
        raise ValueError("first_n must be >= 0")

    output_vcf.parent.mkdir(parents=True, exist_ok=True)
    kept_rows = 0
    per_chromosome: dict[str, int] = {}
    with input_vcf.open(encoding="utf-8") as src, output_vcf.open("w", encoding="utf-8") as dst:
        for line in src:
            if line.startswith("#"):
                dst.write(line)
                continue
            if kept_rows >= first_n:
                break
            dst.write(line)
            kept_rows += 1
            chrom = line.split("\t", 1)[0]
            _increment_counter(per_chromosome, chrom)
    return kept_rows, per_chromosome


def _increment_counter(counter: dict[str, int], chrom: str, increment: int = 1) -> None:
    canonical = normalize_chromosome_name(chrom)
    counter[canonical] = counter.get(canonical, 0) + increment


def _filter_vcf_by_chromosomes(
    input_vcf: Path,
    output_vcf: Path,
    *,
    chromosome_aliases: set[str],
) -> tuple[int, dict[str, int]]:
    output_vcf.parent.mkdir(parents=True, exist_ok=True)
    kept_rows = 0
    per_chromosome: dict[str, int] = {}
    with input_vcf.open(encoding="utf-8") as src, output_vcf.open("w", encoding="utf-8") as dst:
        for line in src:
            if line.startswith("#"):
                dst.write(line)
                continue
            columns = line.split("\t", 1)
            if not columns or columns[0] not in chromosome_aliases:
                continue
            dst.write(line)
            kept_rows += 1
            _increment_counter(per_chromosome, columns[0])
    return kept_rows, per_chromosome


def _summarize_vcf_records(input_vcf: Path) -> tuple[int, dict[str, int]]:
    total = 0
    per_chromosome: dict[str, int] = {}
    with input_vcf.open(encoding="utf-8") as src:
        for line in src:
            if line.startswith("#"):
                continue
            total += 1
            chrom = line.split("\t", 1)[0]
            _increment_counter(per_chromosome, chrom)
    return total, per_chromosome


def decompose_multiallelic_vcf(input_vcf: Path, output_vcf: Path) -> PreparedInputStats:
    output_vcf.parent.mkdir(parents=True, exist_ok=True)
    source_records = 0
    output_records = 0
    split_source_records = 0
    source_per_chromosome: dict[str, int] = {}
    output_per_chromosome: dict[str, int] = {}
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
                if columns:
                    _increment_counter(source_per_chromosome, columns[0])
                    _increment_counter(output_per_chromosome, columns[0])
                continue
            _increment_counter(source_per_chromosome, columns[0])
            alts = columns[4].split(",")
            if len(alts) <= 1:
                dst.write(line)
                output_records += 1
                _increment_counter(output_per_chromosome, columns[0])
                continue
            split_source_records += 1
            for alt in alts:
                cloned = list(columns)
                cloned[4] = alt
                dst.write("\t".join(cloned) + "\n")
                output_records += 1
                _increment_counter(output_per_chromosome, columns[0])
    return PreparedInputStats(
        source_records=source_records,
        sampled_records=source_records,
        filtered_records=source_records,
        output_records=output_records,
        split_source_records=split_source_records,
        source_records_per_chromosome=source_per_chromosome,
        sampled_records_per_chromosome=source_per_chromosome,
        filtered_records_per_chromosome=source_per_chromosome,
        output_records_per_chromosome=output_per_chromosome,
    )


def prepare_vcf_for_annotation(
    input_vcf: Path,
    output_vcf: Path,
    *,
    first_n: int | None = None,
    chromosome_aliases: set[str] | None = None,
) -> PreparedInputStats:
    source_records, source_records_per_chromosome = _summarize_vcf_records(input_vcf)
    staged_input = input_vcf
    sampled_records = 0
    sampled_records_per_chromosome: dict[str, int] = {}
    filtered_records = source_records
    filtered_records_per_chromosome = dict(source_records_per_chromosome)
    if chromosome_aliases:
        filtered_input = output_vcf.parent / "filtered_input.vcf"
        filtered_records, filtered_records_per_chromosome = _filter_vcf_by_chromosomes(
            input_vcf,
            filtered_input,
            chromosome_aliases=chromosome_aliases,
        )
        staged_input = filtered_input
    if first_n is not None:
        sampled_input = output_vcf.parent / "sampled_input.vcf"
        sampled_records, sampled_records_per_chromosome = _sample_vcf_first_n_with_counts(
            staged_input,
            sampled_input,
            first_n,
        )
        staged_input = sampled_input
    else:
        sampled_records = filtered_records
        sampled_records_per_chromosome = dict(filtered_records_per_chromosome)
    stats = decompose_multiallelic_vcf(staged_input, output_vcf)
    return PreparedInputStats(
        source_records=source_records,
        sampled_records=sampled_records,
        filtered_records=filtered_records,
        output_records=stats.output_records,
        split_source_records=stats.split_source_records,
        source_records_per_chromosome=source_records_per_chromosome,
        sampled_records_per_chromosome=sampled_records_per_chromosome,
        filtered_records_per_chromosome=filtered_records_per_chromosome,
        output_records_per_chromosome=stats.output_records_per_chromosome,
    )
