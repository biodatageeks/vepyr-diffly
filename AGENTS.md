# Contributor Guide

This repository compares Ensembl VEP output against local-checkout `vepyr` output using Python, Polars, and semantic CSQ normalization.

## Mandatory Working Rules

- Read [`.agent/CONTINUITY.md`](/Users/lukaszjezapkowicz/Desktop/magisterka/praca/vepyr_diffly/.agent/CONTINUITY.md) at the start of each turn.
- Keep the repo container-first. Do not install system packages on the host unless explicitly requested.
- Default runtime path is `docker compose` plus bind mounts for inputs, caches, and output directories.
- Do not mutate a `vepyr` checkout automatically. Treat it as an input mount.
- Use semantic comparison of parsed `INFO/CSQ` as the canonical equality rule. Raw `CSQ` strings are supporting evidence only.
- Preserve deterministic output structure under the selected run directory.

## Project Structure

- `src/vepyr_diffly/`: package code.
- `tests/`: unit and integration tests.
- `presets/`: staged TOML preset matrix.
- `docker/`: container definitions and helper entrypoints.

## Runtime Contract

- `vepyr-diff run` must:
  - resolve one preset,
  - optionally sample the first `N` data rows from a source VCF while preserving headers,
  - run Ensembl VEP and `vepyr` with equivalent effective flags,
  - normalize annotated VCF outputs into variant and consequence tables,
  - compare both tiers using `diffly`,
  - emit clear console output plus file artifacts.

- `vepyr-diff list-presets` prints the staged preset matrix.
- `vepyr-diff inspect-run` summarizes an existing artifact directory.

## Validation Order

Run these from repo root when code changes:

```bash
uv run ruff format .
uv run ruff check .
uv run pytest
```

If runtime code changes substantially, also run:

```bash
uv run python -m vepyr_diffly.cli list-presets
```

## Testing Expectations

- Prefer small deterministic fixture VCFs for unit tests.
- Keep one integration-style sampled workflow test that exercises artifact generation end to end without requiring external VEP infrastructure.
- Validate normalized table schemas explicitly.

## Reporting Expectations

- Console output must remain concise and operator-friendly.
- File outputs must include machine-readable summary plus inspectable mismatch tables.
- When introducing a new comparison rule, add fixture coverage that demonstrates the failure mode it prevents.

## Run Artifact Hygiene

- After every smoke test or benchmark-style non-golden run, clean its artifacts from `runs/` once the result has been captured.
- Keep long-lived artifacts in `runs/` only for golden tests or explicitly requested retained runs.
