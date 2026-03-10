# Association Testing Module Design

## Overview
This module performs per-variant OLS-based linear association testing against phenotype, with optional covariate adjustment, producing PLINK-like association result rows.

Implementation note:
- Core regression machinery will be built on `numpy` linear algebra primitives for design-matrix assembly and least-squares fitting.
- `scipy` will be used for statistical inference utilities such as test-statistic to p-value calculations.
- During development or validation, `statsmodels` may be used as a reference implementation to cross-check OLS behavior before finalizing the project-specific regression path.

## Responsibilities
- Build regression design per tested variant using aligned phenotype/covariate/genotype data.
- Execute OLS for each QC-passed variant.
- Compute and emit effect estimate and inference statistics required for output compatibility.
- Preserve NMISS and df-compatible fields for PLINK-like interpretation.
- Respect `hide-covar` behavior for emitted rows.

Not owned by this module:
- Input parsing and sample alignment.
- QC threshold selection/computation.
- Final file formatting/writing.
- External validation against PLINK.

## Inputs
- `AlignedCohort` (`y`, `X_covar`, retained sample IDs).
- QC-passed per-variant genotype vectors and metadata.
- Config flags:
  - `linear_enabled`
  - `hide_covar`

## Outputs
- `AssocResultRow` stream with fields compatible with `.assoc.linear` expectations, including:
  - Variant identifiers/metadata fields.
  - Test label (additive genotype test and covariate tests when not hidden).
  - `BETA`, `SE`, `T_STAT`, `P`, `NMISS`, `DF`.
- Optional variant-level status metadata for skipped/failed fits.

## Interface Contract
- `run_linear_assoc(cohort: AlignedCohort, variants: FilteredVariantStream, cfg: RunConfig) -> AssocResultStream`
- `fit_variant_ols(y, X_base, g) -> FitStats`
- Model policy:
  - OLS with intercept, genotype term, and optional covariate terms.
  - Variant-specific usable row set based on non-missing genotype and required covariates/phenotype.
  - `NMISS` reflects usable rows for the specific variant test.
  - `DF` derived from usable rows minus number of fitted parameters.
- `hide-covar` policy:
  - If enabled, emit genotype test rows only.
  - If disabled, emit genotype and covariate test rows per output convention.

## Edge Cases
- Singular design matrix due to collinearity.
- Too few usable rows for stable estimate/valid df.
- Near-zero genotype variance after missingness filtering.
- Large or extreme phenotype/genotype values affecting numeric stability.

## Failure Modes
- Per-variant fit failure should yield deterministic skip/NA row policy (implementation chooses one policy and applies consistently).
- Run-level failure only when core model infrastructure is invalid (for example, empty cohort).
- Failures should carry reason codes for downstream logging.

## Module Dependencies
- Upstream:
  - Sample/variant alignment for cohort matrices/vectors.
  - QC/filters for variant eligibility.
- Downstream:
  - Output writing module consumes association row stream.
- External dependencies:
  - `numpy` (design matrices, solves, residual calculations, and per-variant matrix operations).
  - `scipy` (distribution functions and other statistical inference helpers).
  - `statsmodels` (optional reference implementation for validation and debugging of OLS behavior).

## Test Scenarios
- Happy-path per-variant OLS with and without covariates.
- `hide-covar` toggles row emission correctly.
- NMISS changes appropriately with variant-specific genotype missingness.
- Singular matrix case handled with deterministic skip/NA behavior.
- Output stats fields are present and type-consistent for all processed variants.
