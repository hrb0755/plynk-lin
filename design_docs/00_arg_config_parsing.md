# Arg/Config Parsing Module Design

## Overview
This module parses and validates command-line arguments for the supported PLINK-like subset, then emits a normalized run configuration for downstream modules.

## Responsibilities
- Parse supported flags: `--linear`, `--vcf`, `--pheno`, `--maf`, `--allow-no-sex`, and `--out`.
- Normalize user inputs into a typed configuration object with defaults.
- Validate flag presence, value types, and allowed combinations.
- Emit actionable, user-facing errors for invalid invocations.

Not owned by this module:
- Reading input files.
- Sample/variant alignment.
- QC filtering or regression.
- Output file rendering.

## Inputs
- CLI argument vector (`argv`) from process entrypoint.
- Optional environment defaults if configured by entrypoint (none required by current proposal).

## Outputs
- `RunConfig` object, for example:
  - `linear_enabled: bool` (required; must be true)
  - `vcf_path: str` (required)
  - `pheno_path: str` (required)
  - `maf_threshold: float | None` (optional)
  - `allow_no_sex: bool` (default `false`)
  - `out_prefix: str` (required)
- Structured parse/validation error object for CLI failure path.

## Interface Contract
- `parse_args(argv: list[str]) -> RunConfig`
- `validate_config(cfg: RunConfig) -> None` (raises validation error)
- Validation rules:
  - `--linear` is required.
  - `--vcf`, `--pheno`, and `--out` are required.
  - `--maf` must parse as numeric and satisfy `0.0 <= maf <= 0.5`.
  - Unknown flags fail fast.

Downstream contract:
- Returned config is complete and type-safe for IO/alignment/QC modules.
- No downstream module should re-interpret raw CLI tokens.

## Edge Cases
- Missing required flags (`--linear`, `--vcf`, `--pheno`, `--out`).
- `--maf` malformed (non-numeric) or outside bounds.
- Repeated flags with conflicting values.
- Empty `--out` prefix.

## Failure Modes
- Parse failure on unsupported syntax or unknown option.
- Validation failure on missing/invalid combinations.
- Error output must identify failing flag and reason, then exit non-zero.

## Module Dependencies
- Upstream: process entrypoint/main.
- Downstream:
  - IO consumes `vcf_path` and `pheno_path`.
  - QC consumes `maf_threshold`.
  - Association/output consume `linear_enabled` and `out_prefix`.

## Test Scenarios
- Valid happy-path invocation with all supported flags.
- Invocation with shuffled flag ordering still parses identically.
- `--maf` at `0.0`, `0.5`, below `0`, above `0.5`.
- Unknown flag fails with deterministic message.
