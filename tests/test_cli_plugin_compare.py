from vepyr_diffly.cli import _resolve_compare_csq_fields


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
