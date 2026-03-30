from __future__ import annotations

import argparse
import json
from pathlib import Path

from .compare import compare_bucket_shard, dump_bucket_shard_summary
from .normalize import materialize_consequence_buckets
from .progress import ProgressReporter


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="vepyr-diff-worker")
    subparsers = parser.add_subparsers(dest="command", required=True)

    compare_parser = subparsers.add_parser("compare-bucket-shard")
    compare_parser.add_argument("--left-bucket-dir", required=True, type=Path)
    compare_parser.add_argument("--right-bucket-dir", required=True, type=Path)
    compare_parser.add_argument("--bucket-ids", required=True)
    compare_parser.add_argument("--csq-fields-json", required=True)
    compare_parser.add_argument("--temp-diff-dir", required=True, type=Path)
    compare_parser.add_argument("--temp-tsv-dir", required=True, type=Path)
    compare_parser.add_argument("--summary-path", required=True, type=Path)
    compare_parser.add_argument("--progress-log", default="")
    compare_parser.add_argument("--compare-mode", choices=("fast", "debug"), default="fast")
    compare_parser.add_argument("--fingerprint-only", action="store_true")

    bucketize_parser = subparsers.add_parser("bucketize-side")
    bucketize_parser.add_argument("--vcf", required=True, type=Path)
    bucketize_parser.add_argument("--bucket-root", required=True, type=Path)
    bucketize_parser.add_argument("--side-label", required=True)
    bucketize_parser.add_argument("--csq-fields-json", required=True)
    bucketize_parser.add_argument("--bucket-count", required=True, type=int)
    bucketize_parser.add_argument("--chunk-variants", required=True, type=int)
    bucketize_parser.add_argument("--total-variants", required=True, type=int)
    bucketize_parser.add_argument("--progress-log", default="")
    return parser


def _build_reporter(progress_log: str) -> ProgressReporter | None:
    if not progress_log:
        return None
    reporter = ProgressReporter(log_path=Path(progress_log), console=None)
    reporter.start()
    return reporter


def _cmd_compare_bucket_shard(args: argparse.Namespace) -> int:
    reporter = _build_reporter(args.progress_log)
    try:
        bucket_ids = [int(item) for item in args.bucket_ids.split(",") if item]

        def _log_bucket_progress(
            completed_in_shard: int,
            total_in_shard: int,
            result,
        ) -> None:
            if reporter is None:
                return
            reporter.log(
                f"worker: completed bucket {result.bucket_id:04d} "
                f"({completed_in_shard}/{total_in_shard} in shard) "
                f"left={result.left_rows} right={result.right_rows} "
                f"unequal={result.unequal_rows} left_only={result.left_only_rows} "
                f"right_only={result.right_only_rows}"
            )

        results = compare_bucket_shard(
            left_bucket_dir=args.left_bucket_dir,
            right_bucket_dir=args.right_bucket_dir,
            bucket_ids=bucket_ids,
            csq_fields=json.loads(args.csq_fields_json),
            temp_diff_dir=args.temp_diff_dir,
            temp_tsv_dir=args.temp_tsv_dir,
            compare_mode=args.compare_mode,
            fingerprint_only=args.fingerprint_only,
            on_bucket_complete=_log_bucket_progress,
        )
        dump_bucket_shard_summary(results, args.summary_path)
        return 0
    finally:
        if reporter is not None:
            reporter.stop()


def _cmd_bucketize_side(args: argparse.Namespace) -> int:
    reporter = _build_reporter(args.progress_log)
    try:
        materialize_consequence_buckets(
            vcf_path=args.vcf,
            csq_fields=json.loads(args.csq_fields_json),
            bucket_root=args.bucket_root,
            reporter=reporter,
            side_label=args.side_label,
            bucket_count=args.bucket_count,
            chunk_variants=args.chunk_variants,
            total_variants=args.total_variants,
        )
        return 0
    finally:
        if reporter is not None:
            reporter.stop()


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    if args.command == "compare-bucket-shard":
        return _cmd_compare_bucket_shard(args)
    if args.command == "bucketize-side":
        return _cmd_bucketize_side(args)
    raise ValueError(f"unsupported command: {args.command}")


if __name__ == "__main__":
    raise SystemExit(main())
