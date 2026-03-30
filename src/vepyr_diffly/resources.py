from __future__ import annotations

import os
from pathlib import Path

from .models import CompareResourcePlan, RuntimeConfig
from .normalize import (
    CONSEQUENCE_BUCKET_CHUNK_VARIANTS,
    CONSEQUENCE_BUCKET_COUNT,
    STREAMING_CONSEQUENCE_THRESHOLD_BYTES,
)


DEFAULT_MEMORY_BUDGET_MB = 2048
STREAMING_VARIANT_THRESHOLD_BYTES = STREAMING_CONSEQUENCE_THRESHOLD_BYTES
MIN_BUCKET_COUNT = 8
MAX_BUCKET_COUNT = 512
MIN_VARIANT_CHUNK_ROWS = 2_000
MIN_CONSEQUENCE_CHUNK_ROWS = 1_000


def _clamp(value: int, *, lower: int, upper: int) -> int:
    return max(lower, min(value, upper))


def plan_compare_resources(
    *,
    config: RuntimeConfig,
    left_vcf: Path,
    right_vcf: Path,
) -> CompareResourcePlan:
    memory_budget_mb = config.memory_budget_mb or DEFAULT_MEMORY_BUDGET_MB
    budget_bytes = memory_budget_mb * 1024 * 1024
    max_size = max(left_vcf.stat().st_size, right_vcf.stat().st_size)
    cpu_count = max(1, os.cpu_count() or 1)

    bucket_count = config.compare_bucket_count or _clamp(
        max(MIN_BUCKET_COUNT, int(max_size / (256 * 1024 * 1024))),
        lower=MIN_BUCKET_COUNT,
        upper=MAX_BUCKET_COUNT,
    )
    if max_size >= 8 * 1024 * 1024 * 1024:
        bucket_count = max(bucket_count, CONSEQUENCE_BUCKET_COUNT)

    variant_chunk_rows = _clamp(
        budget_bytes // (32 * 1024),
        lower=MIN_VARIANT_CHUNK_ROWS,
        upper=50_000,
    )
    consequence_chunk_rows = _clamp(
        budget_bytes // (128 * 1024),
        lower=MIN_CONSEQUENCE_CHUNK_ROWS,
        upper=CONSEQUENCE_BUCKET_CHUNK_VARIANTS,
    )
    if max_size >= 8 * 1024 * 1024 * 1024:
        consequence_chunk_rows = min(consequence_chunk_rows, 20_000)
        variant_chunk_rows = min(variant_chunk_rows, 10_000)
    elif max_size >= 2 * 1024 * 1024 * 1024:
        consequence_chunk_rows = min(consequence_chunk_rows, 50_000)
        variant_chunk_rows = min(variant_chunk_rows, 20_000)

    max_safe_workers = max(1, budget_bytes // (512 * 1024 * 1024))
    compare_workers = config.compare_workers or min(cpu_count, max_safe_workers)
    compare_workers = max(1, compare_workers)

    use_bucketed_variant = max_size >= STREAMING_VARIANT_THRESHOLD_BYTES or memory_budget_mb <= 4096
    use_bucketed_consequence = (
        max_size >= STREAMING_CONSEQUENCE_THRESHOLD_BYTES or memory_budget_mb <= 4096
    )
    parallelize_sides = (
        max_size < STREAMING_VARIANT_THRESHOLD_BYTES
        and memory_budget_mb >= 4096
        and compare_workers >= 2
    )

    return CompareResourcePlan(
        memory_budget_mb=memory_budget_mb,
        bucket_count=bucket_count,
        variant_chunk_rows=variant_chunk_rows,
        consequence_chunk_rows=consequence_chunk_rows,
        compare_workers=compare_workers,
        parallelize_sides=parallelize_sides,
        use_bucketed_variant=use_bucketed_variant,
        use_bucketed_consequence=use_bucketed_consequence,
    )
