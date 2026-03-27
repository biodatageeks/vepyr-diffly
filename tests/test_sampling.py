from pathlib import Path

from vepyr_diffly.sampling import sample_vcf_first_n


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

