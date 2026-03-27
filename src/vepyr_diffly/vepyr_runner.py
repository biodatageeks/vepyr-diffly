from __future__ import annotations

import argparse


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-vcf", required=True)
    parser.add_argument("--output-vcf", required=True)
    parser.add_argument("--cache-dir", required=True)
    parser.add_argument("--reference-fasta", default="")
    args = parser.parse_args()

    import vepyr

    vepyr.annotate(
        vcf=args.input_vcf,
        cache_dir=args.cache_dir,
        everything=True,
        output_vcf=args.output_vcf,
        reference_fasta=args.reference_fasta or None,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
