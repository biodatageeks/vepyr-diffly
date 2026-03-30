from __future__ import annotations

import argparse
import importlib
import os
from pathlib import Path
import sys


def _repo_vepyr_src() -> Path:
    return Path(__file__).resolve().parents[2] / "vepyr" / "src"


def _prepare_import_path() -> None:
    repo_src = _repo_vepyr_src()
    if repo_src.exists():
        sys.path.insert(0, str(repo_src))


def _validate_vepyr(module):
    if hasattr(module, "annotate"):
        return module
    raise RuntimeError(
        "Imported 'vepyr' but it does not expose annotate(). "
        "Use the compiled vepyr interpreter from VEPYR_DIFFLY_VEPYR_PYTHON or "
        "the configured .vepyr/bin/python environment."
    )


def _import_vepyr():
    try:
        return _validate_vepyr(importlib.import_module("vepyr"))
    except Exception as first_error:
        repo_src = _repo_vepyr_src()
        if not repo_src.exists():
            raise first_error
        original_path = list(sys.path)
        sys.modules.pop("vepyr", None)
        sys.path.insert(0, str(repo_src))
        try:
            return _validate_vepyr(importlib.import_module("vepyr"))
        except Exception:
            raise first_error
        finally:
            sys.path[:] = original_path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-vcf", required=True)
    parser.add_argument("--output-vcf", required=True)
    parser.add_argument("--cache-dir", required=True)
    parser.add_argument("--reference-fasta", default="")
    parser.add_argument("--use-fjall", action="store_true")
    args = parser.parse_args()
    try:
        vepyr = _import_vepyr()
    except Exception as exc:
        configured_python = os.environ.get("VEPYR_DIFFLY_VEPYR_PYTHON", "")
        hint = (
            f" Try rerunning with VEPYR_DIFFLY_VEPYR_PYTHON={configured_python}"
            if configured_python
            else ""
        )
        raise RuntimeError(f"Unable to import a working vepyr module.{hint}") from exc

    vepyr.annotate(
        vcf=args.input_vcf,
        cache_dir=args.cache_dir,
        everything=True,
        output_vcf=args.output_vcf,
        reference_fasta=args.reference_fasta or None,
        use_fjall=args.use_fjall,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
