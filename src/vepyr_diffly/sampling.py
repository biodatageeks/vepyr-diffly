from __future__ import annotations

from pathlib import Path


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

