from pathlib import Path

import pytest

from vepyr_diffly.plugins import _indexed_path, compare_plugin_fields


def test_compare_plugin_fields_preserves_plugin_order_and_deduplicates() -> None:
    fields = compare_plugin_fields(["spliceai", "clinvar", "spliceai"])

    assert fields[:3] == ["symbol", "ds_ag", "ds_al"]
    assert fields[-2:] == ["ClinVar_CLNVC", "ClinVar_CLNVI"]
    assert fields.count("symbol") == 1


def test_indexed_path_accepts_bgz_sibling_for_plain_gzip_sources(tmp_path: Path) -> None:
    source = tmp_path / "AlphaMissense_hg38.tsv.gz"
    source.write_text("placeholder", encoding="utf-8")
    sibling_bgz = tmp_path / "AlphaMissense_hg38.tsv.bgz"
    sibling_bgz.write_text("placeholder", encoding="utf-8")
    Path(str(sibling_bgz) + ".tbi").write_text("index", encoding="utf-8")

    assert _indexed_path(source, "alphamissense") == sibling_bgz


def test_indexed_path_raises_when_no_indexed_variant_exists(tmp_path: Path) -> None:
    source = tmp_path / "AlphaMissense_hg38.tsv.gz"
    source.write_text("placeholder", encoding="utf-8")

    with pytest.raises(ValueError, match="not tabix-indexed"):
        _indexed_path(source, "alphamissense")
