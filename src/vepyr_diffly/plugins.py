from __future__ import annotations

import os
from pathlib import Path


SUPPORTED_PLUGINS = ("clinvar", "spliceai", "cadd", "alphamissense", "dbnsfp")

PLUGIN_COMPARE_FIELDS: dict[str, tuple[str, ...]] = {
    "clinvar": (
        "ClinVar",
        "ClinVar_CLNSIG",
        "ClinVar_CLNREVSTAT",
        "ClinVar_CLNDN",
        "ClinVar_CLNVC",
        "ClinVar_CLNVI",
    ),
    "spliceai": (
        "symbol",
        "ds_ag",
        "ds_al",
        "ds_dg",
        "ds_dl",
        "dp_ag",
        "dp_al",
        "dp_dg",
        "dp_dl",
    ),
    "cadd": ("raw_score", "phred_score"),
    "alphamissense": ("am_pathogenicity", "am_class"),
    "dbnsfp": (
        "sift4g_score",
        "sift4g_pred",
        "polyphen2_hdiv_score",
        "polyphen2_hvar_score",
        "lrt_score",
        "lrt_pred",
        "mutationtaster_score",
        "mutationtaster_pred",
        "fathmm_score",
        "fathmm_pred",
        "provean_score",
        "provean_pred",
        "vest4_score",
        "metasvm_score",
        "metasvm_pred",
        "metalr_score",
        "metalr_pred",
        "revel_score",
        "gerp_rs",
        "phylop100way",
        "phylop30way",
        "phastcons100way",
        "phastcons30way",
        "siphy_29way",
        "cadd_raw",
        "cadd_phred",
    ),
}


def parse_plugin_list(raw_value: str | None) -> list[str]:
    if raw_value is None:
        return []
    plugins: list[str] = []
    for item in raw_value.split(","):
        plugin = item.strip().lower()
        if not plugin:
            continue
        if plugin not in SUPPORTED_PLUGINS:
            supported = ", ".join(SUPPORTED_PLUGINS)
            raise ValueError(f"unsupported plugin '{plugin}'; supported plugins: {supported}")
        if plugin not in plugins:
            plugins.append(plugin)
    return plugins


def compare_plugin_fields(plugins: list[str]) -> list[str]:
    fields: list[str] = []
    for plugin in plugins:
        for field in PLUGIN_COMPARE_FIELDS[plugin]:
            if field not in fields:
                fields.append(field)
    return fields


def _env_path(name: str) -> Path | None:
    value = os.environ.get(name)
    if value is None or not value.strip():
        return None
    return Path(value).expanduser()


def _required_env_path(name: str, plugin: str) -> Path:
    path = _env_path(name)
    if path is None:
        raise ValueError(f"missing source path for plugin '{plugin}': set {name}")
    if not path.exists():
        raise ValueError(f"source path for plugin '{plugin}' does not exist: {path}")
    return path


def _indexed_path(path: Path, plugin: str) -> Path:
    if Path(str(path) + ".tbi").exists():
        return path
    if path.suffix == ".gz":
        bgz = path.with_name(path.name[:-3] + ".bgz")
        if Path(str(bgz) + ".tbi").exists():
            return bgz
    raise ValueError(
        f"plugin '{plugin}' source is not tabix-indexed for VEP: {path}. "
        "Create an index first, e.g. scripts/create_plugin_indexes.py --plugins "
        f"{plugin} --recompress-plain-gzip when the source is plain gzip."
    )


def vep_plugin_args(plugins: list[str]) -> list[str]:
    args: list[str] = []
    for plugin in plugins:
        if plugin == "clinvar":
            source = _indexed_path(
                _required_env_path("VEPYR_DIFFLY_PLUGIN_CLINVAR_SOURCE", plugin), plugin
            )
            args.extend(
                [
                    "--custom",
                    f"{source},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN,CLNVC,CLNVI",
                ]
            )
        elif plugin == "spliceai":
            snv = _indexed_path(
                _required_env_path("VEPYR_DIFFLY_PLUGIN_SPLICEAI_SOURCE", plugin), plugin
            )
            indel = _env_path("VEPYR_DIFFLY_PLUGIN_SPLICEAI_INDEL_SOURCE") or snv
            indel = _indexed_path(indel, plugin)
            args.extend(["--plugin", f"SpliceAI,snv={snv},indel={indel}"])
        elif plugin == "cadd":
            snv = _indexed_path(
                _required_env_path("VEPYR_DIFFLY_PLUGIN_CADD_SNV_SOURCE", plugin), plugin
            )
            indel = _indexed_path(
                _required_env_path("VEPYR_DIFFLY_PLUGIN_CADD_INDEL_SOURCE", plugin), plugin
            )
            args.extend(["--plugin", f"CADD,snv={snv},indels={indel}"])
        elif plugin == "alphamissense":
            source = _indexed_path(
                _required_env_path("VEPYR_DIFFLY_PLUGIN_ALPHAMISSENSE_SOURCE", plugin),
                plugin,
            )
            args.extend(["--plugin", f"AlphaMissense,file={source}"])
        elif plugin == "dbnsfp":
            source = _indexed_path(
                _required_env_path("VEPYR_DIFFLY_PLUGIN_DBNSFP_SOURCE", plugin), plugin
            )
            args.extend(
                [
                    "--plugin",
                    f"dbNSFP,{source},SIFT4G_score,Polyphen2_HDIV_score,REVEL_score,CADD_phred",
                ]
            )
    return args
