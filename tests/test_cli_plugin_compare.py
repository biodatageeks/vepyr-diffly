from vepyr_diffly.cli import _resolve_compare_csq_fields, _resolve_csq_field_indexes


def test_resolve_compare_csq_fields_uses_plugin_subset() -> None:
    left = ["Allele", "Consequence", "ClinVar", "ClinVar_CLNSIG", "ClinVar_CLNVI"]
    right = ["Allele", "SYMBOL", "ClinVar", "ClinVar_CLNSIG", "ClinVar_CLNVI"]

    resolved = _resolve_compare_csq_fields(
        left_csq_fields=left,
        right_csq_fields=right,
        plugins=["clinvar"],
        compare_only_plugins=True,
    )

    assert resolved == ["ClinVar", "ClinVar_CLNSIG", "ClinVar_CLNVI"]


def test_resolve_compare_csq_fields_requires_plugins_when_compare_only_enabled() -> None:
    try:
        _resolve_compare_csq_fields(
            left_csq_fields=["Allele"],
            right_csq_fields=["Allele"],
            plugins=[],
            compare_only_plugins=True,
        )
    except ValueError as exc:
        assert "--compare-only-plugins requires --plugins" in str(exc)
    else:
        raise AssertionError("expected ValueError")


def test_resolve_compare_csq_fields_accepts_dbnsfp_vep_aliases() -> None:
    left = ["Allele", "SIFT4G_score", "SIFT4G_pred", "CADD_phred"]
    right = ["Allele", "sift4g_score", "sift4g_pred", "cadd_phred"]

    resolved = _resolve_compare_csq_fields(
        left_csq_fields=left,
        right_csq_fields=right,
        plugins=["dbnsfp"],
        compare_only_plugins=True,
    )

    assert resolved == ["sift4g_score", "sift4g_pred", "cadd_phred"]


def test_resolve_csq_field_indexes_maps_dbnsfp_aliases_to_canonical_fields() -> None:
    result = _resolve_csq_field_indexes(
        selected_fields=["sift4g_score", "cadd_phred"],
        header_fields=["Allele", "SIFT4G_score", "CADD_phred"],
        plugins=["dbnsfp"],
        compare_only_plugins=True,
    )

    assert result == {"sift4g_score": 1, "cadd_phred": 2}


def test_resolve_compare_csq_fields_accepts_spliceai_vep_aliases() -> None:
    left = [
        "Allele",
        "SpliceAI_pred_DP_AG",
        "SpliceAI_pred_DS_AG",
        "SpliceAI_pred_SYMBOL",
    ]
    right = ["Allele", "symbol", "ds_ag", "dp_ag"]

    resolved = _resolve_compare_csq_fields(
        left_csq_fields=left,
        right_csq_fields=right,
        plugins=["spliceai"],
        compare_only_plugins=True,
    )

    assert resolved == ["symbol", "ds_ag", "dp_ag"]


def test_resolve_csq_field_indexes_maps_spliceai_aliases_to_canonical_fields() -> None:
    result = _resolve_csq_field_indexes(
        selected_fields=["symbol", "ds_ag", "dp_ag"],
        header_fields=["Allele", "SpliceAI_pred_SYMBOL", "SpliceAI_pred_DS_AG", "SpliceAI_pred_DP_AG"],
        plugins=["spliceai"],
        compare_only_plugins=True,
    )

    assert result == {"symbol": 1, "ds_ag": 2, "dp_ag": 3}


def test_resolve_compare_csq_fields_accepts_cadd_aliases() -> None:
    left = ["Allele", "CADD_RAW", "CADD_PHRED"]
    right = ["Allele", "raw_score", "phred_score"]

    resolved = _resolve_compare_csq_fields(
        left_csq_fields=left,
        right_csq_fields=right,
        plugins=["cadd"],
        compare_only_plugins=True,
    )

    assert resolved == ["cadd_raw", "cadd_phred"]


def test_resolve_csq_field_indexes_maps_cadd_aliases_to_canonical_fields() -> None:
    result = _resolve_csq_field_indexes(
        selected_fields=["cadd_raw", "cadd_phred"],
        header_fields=["Allele", "CADD_RAW", "CADD_PHRED"],
        plugins=["cadd"],
        compare_only_plugins=True,
    )

    assert result == {"cadd_raw": 1, "cadd_phred": 2}
