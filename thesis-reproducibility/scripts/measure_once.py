#!/usr/bin/env python3
"""Measure a single VEPYR annotation run in an isolated subprocess."""

import argparse
import json
import resource
import sys
import time
from pathlib import Path

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-vcf", required=True)
    parser.add_argument("--cache-dir", required=True)
    parser.add_argument("--reference-fasta", default="")
    parser.add_argument("--plugin-cache-root", default="")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--output-vcf", required=True)
    args = parser.parse_args()

    import vepyr

    out_path = Path(args.output_vcf)
    if out_path.exists():
        out_path.unlink()

    t0 = time.perf_counter()
    vepyr.annotate(
        vcf=args.input_vcf,
        cache_dir=args.cache_dir,
        everything=True,
        hgvs=True,
        reference_fasta=args.reference_fasta or None,
        plugin_cache_root=args.plugin_cache_root or None,
        workers=args.workers,
        output_vcf=str(out_path),
        show_progress=False,
    )
    annotation_seconds = time.perf_counter() - t0

    usage = resource.getrusage(resource.RUSAGE_SELF)
    # On macOS, ru_maxrss is in bytes; on Linux, in kilobytes.
    rss_bytes = usage.ru_maxrss if sys.platform == "darwin" else usage.ru_maxrss * 1024

    if out_path.exists():
        out_path.unlink()

    result = {
        "annotation_seconds": annotation_seconds,
        "max_rss_bytes": rss_bytes,
        "max_rss_gb": round(rss_bytes / (1024**3), 3),
    }
    print(json.dumps(result))

if __name__ == "__main__":
    main()
