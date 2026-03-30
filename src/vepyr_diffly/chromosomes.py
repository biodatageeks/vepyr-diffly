from __future__ import annotations

from collections.abc import Iterable


STANDARD_CHROMOSOMES = [str(index) for index in range(1, 23)] + ["X", "Y", "MT"]


def normalize_chromosome_name(name: str) -> str:
    value = name.strip()
    lower = value.lower()
    if lower.startswith("chr"):
        value = value[3:]
        lower = value.lower()
    if lower in {"m", "mt"}:
        return "MT"
    if lower == "x":
        return "X"
    if lower == "y":
        return "Y"
    if value.isdigit():
        stripped = value.lstrip("0")
        return stripped or "0"
    return value


def chromosome_aliases(name: str) -> set[str]:
    canonical = normalize_chromosome_name(name)
    if canonical in STANDARD_CHROMOSOMES:
        aliases = {canonical, f"chr{canonical}"}
        if canonical == "MT":
            aliases.update({"M", "chrM"})
        return aliases
    return {name}


def parse_chromosome_selection(raw_value: str | None) -> tuple[list[str], set[str]]:
    if raw_value is None or not raw_value.strip():
        return [], set()
    canonical: list[str] = []
    seen: set[str] = set()
    aliases: set[str] = set()
    for item in raw_value.split(","):
        token = item.strip()
        if not token:
            continue
        normalized = normalize_chromosome_name(token)
        if normalized not in seen:
            canonical.append(normalized)
            seen.add(normalized)
        aliases.update(chromosome_aliases(token))
        aliases.update(chromosome_aliases(normalized))
    return canonical, aliases


def canonicalize_chromosome_iter(values: Iterable[str]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for value in values:
        canonical = normalize_chromosome_name(value)
        counts[canonical] = counts.get(canonical, 0) + 1
    return counts
