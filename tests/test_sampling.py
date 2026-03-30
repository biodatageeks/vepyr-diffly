from pathlib import Path

from vepyr_diffly.chromosomes import parse_chromosome_selection
from vepyr_diffly.sampling import (
    decompose_multiallelic_vcf,
    prepare_vcf_for_annotation,
    sample_vcf_first_n,
)


def test_sample_vcf_first_n_preserves_header_and_first_rows(tmp_path: Path) -> None:
    src = tmp_path / "input.vcf"
    src.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t1\t.\tA\tG\t.\tPASS\t.\n"
        "1\t2\t.\tC\tT\t.\tPASS\t.\n"
        "1\t3\t.\tG\tA\t.\tPASS\t.\n",
        encoding="utf-8",
    )
    dst = tmp_path / "sample.vcf"

    kept = sample_vcf_first_n(src, dst, 2)

    assert kept == 2
    assert dst.read_text(encoding="utf-8").splitlines() == [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        "1\t1\t.\tA\tG\t.\tPASS\t.",
        "1\t2\t.\tC\tT\t.\tPASS\t.",
    ]


def test_decompose_multiallelic_vcf_splits_alt_rows(tmp_path: Path) -> None:
    src = tmp_path / "input.vcf"
    src.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t1\t.\tA\tG,T\t.\tPASS\t.\n"
        "1\t2\t.\tC\tT\t.\tPASS\t.\n",
        encoding="utf-8",
    )
    dst = tmp_path / "prepared.vcf"

    stats = decompose_multiallelic_vcf(src, dst)

    assert stats.source_records == 2
    assert stats.output_records == 3
    assert stats.split_source_records == 1
    assert dst.read_text(encoding="utf-8").splitlines() == [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        "1\t1\t.\tA\tG\t.\tPASS\t.",
        "1\t1\t.\tA\tT\t.\tPASS\t.",
        "1\t2\t.\tC\tT\t.\tPASS\t.",
    ]


def test_prepare_vcf_for_annotation_samples_then_splits(tmp_path: Path) -> None:
    src = tmp_path / "input.vcf"
    src.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t1\t.\tA\tG,T\t.\tPASS\t.\n"
        "1\t2\t.\tC\tT\t.\tPASS\t.\n"
        "1\t3\t.\tG\tA,C\t.\tPASS\t.\n",
        encoding="utf-8",
    )
    dst = tmp_path / "prepared.vcf"

    stats = prepare_vcf_for_annotation(src, dst, first_n=2)

    assert stats.sampled_records == 2
    assert stats.output_records == 3
    assert stats.split_source_records == 1
    assert dst.read_text(encoding="utf-8").splitlines() == [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        "1\t1\t.\tA\tG\t.\tPASS\t.",
        "1\t1\t.\tA\tT\t.\tPASS\t.",
        "1\t2\t.\tC\tT\t.\tPASS\t.",
    ]


def test_prepare_vcf_for_annotation_filters_before_sampling(tmp_path: Path) -> None:
    src = tmp_path / "input.vcf"
    src.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "2\t1\t.\tA\tG\t.\tPASS\t.\n"
        "1\t2\t.\tC\tT\t.\tPASS\t.\n"
        "2\t3\t.\tG\tA\t.\tPASS\t.\n"
        "chr1\t4\t.\tT\tC\t.\tPASS\t.\n"
        "1\t5\t.\tA\tG\t.\tPASS\t.\n",
        encoding="utf-8",
    )
    dst = tmp_path / "prepared.vcf"
    _, chr1_aliases = parse_chromosome_selection("1")

    stats = prepare_vcf_for_annotation(
        src,
        dst,
        first_n=2,
        chromosome_aliases=chr1_aliases,
    )

    assert stats.source_records == 5
    assert stats.filtered_records == 3
    assert stats.sampled_records == 2
    assert stats.filtered_records_per_chromosome == {"1": 3}
    assert stats.sampled_records_per_chromosome == {"1": 2}
    assert dst.read_text(encoding="utf-8").splitlines() == [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        "1\t2\t.\tC\tT\t.\tPASS\t.",
        "chr1\t4\t.\tT\tC\t.\tPASS\t.",
    ]
