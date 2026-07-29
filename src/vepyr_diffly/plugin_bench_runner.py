"""Subprocess entrypoint for `benchmark-plugins`: one `vepyr.annotate()` call.

A separate process per plugin config (rather than importing `vepyr` in-process
in `cli.py`) so `psutil.Process(pid)` can sample this process's own RSS in
isolation, unpolluted by the CLI's own polars/comparison-pipeline memory.

Unlike `vepyr_runner.py`, this calls the current `vepyr.annotate()` signature
directly (`plugin_cache_root`, `hgvs`) -- `vepyr_runner.py`'s `use_fjall`/
`plugins` kwargs are stale against the live API (the fjall backend and the
plugin allow-list were both removed; plugins are now selected purely by
which subdirectories exist under `plugin_cache_root`).
"""

from __future__ import annotations

import argparse

from vepyr_diffly.vepyr_runner import _import_vepyr


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-vcf", required=True)
    parser.add_argument("--output-vcf", required=True)
    parser.add_argument("--cache-dir", required=True)
    parser.add_argument("--reference-fasta", default="")
    parser.add_argument("--plugin-cache-root", default="")
    parser.add_argument("--hgvs", action="store_true")
    args = parser.parse_args()

    vepyr = _import_vepyr()
    vepyr.annotate(
        vcf=args.input_vcf,
        cache_dir=args.cache_dir,
        everything=True,
        hgvs=args.hgvs,
        output_vcf=args.output_vcf,
        reference_fasta=args.reference_fasta or None,
        plugin_cache_root=args.plugin_cache_root or None,
        show_progress=False,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
