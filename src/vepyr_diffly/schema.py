from __future__ import annotations

import dataframely as dy
import polars as pl


class VariantSchema(dy.Schema):
    chrom = dy.String(nullable=False)
    pos = dy.Int64(nullable=False)
    ref = dy.String(nullable=False)
    alt = dy.String(nullable=False)
    record_count = dy.Int64(nullable=False)
    consequence_count = dy.Int64(nullable=False)
    ids = dy.String(nullable=False)
    filters = dy.String(nullable=False)

    @dy.rule()
    def non_negative_counts(cls) -> pl.Expr:
        return (pl.col("record_count") >= 0) & (pl.col("consequence_count") >= 0)


def validate_variant_schema(frame: pl.DataFrame) -> pl.DataFrame:
    return VariantSchema.validate(frame, cast=False)
