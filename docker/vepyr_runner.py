from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-vcf", required=True)
    parser.add_argument("--output-vcf", required=True)
    parser.add_argument("--cache-dir", required=True)
    args, _extra = parser.parse_known_args()

    install = subprocess.run(
        [sys.executable, "-m", "pip", "install", "/workspace/vepyr"],
        check=False,
        capture_output=True,
        text=True,
    )
    if install.returncode != 0:
        sys.stderr.write(install.stdout)
        sys.stderr.write(install.stderr)
        return install.returncode

    import vepyr

    annotate = getattr(vepyr, "annotate", None) or getattr(vepyr, "annotate_vcf", None)
    if annotate is None:
        raise RuntimeError("vepyr does not expose annotate/annotate_vcf")

    annotate(
        input_vcf=Path(args.input_vcf),
        output_vcf=Path(args.output_vcf),
        cache_dir=Path(args.cache_dir),
        everything=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
