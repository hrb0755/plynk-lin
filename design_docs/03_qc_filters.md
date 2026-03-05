# QC/Filters Module Design

## Overview
This module applies variant-level quality control filters to aligned variants, primarily minor allele frequency filtering (`--maf`), and emits pass/fail decisions with traceable metadata.

## Responsibilities
- Compute per-variant QC metrics required by current scope.
- Apply configured filters, at minimum MAF threshold.
- Handle genotype missingness consistently during metric computation.
- Emit variant pass/fail decisions and filter reasons.

Not owned by this module:
- CLI parsing and threshold validation.
- Sample intersection/listwise deletion logic.
- Regression statistic computation.
- Output file writing.

## Inputs
- `AlignedVariant` stream from alignment module.
- `RunConfig.maf_threshold` (optional).
- Cohort context needed for metric computation (retained samples, non-missing counts).

## Outputs
- `FilteredVariant` stream:
  - Variants that pass QC with retained genotype vectors.
- `QcDecision` metadata per variant:
  - `passed: bool`
  - `reasons: list[str]`
  - `maf_value: float | None`
  - `non_missing_genotype_count: int`

## Interface Contract
- `compute_maf(g: vector[float | missing]) -> float | None`
- `apply_variant_filters(variant, cfg: RunConfig) -> QcDecision`
- `filter_variants(stream, cfg: RunConfig) -> FilteredVariantStream`
- MAF computation policy:
  - Use non-missing diploid genotype calls only.
  - Derive allele frequency from called allele counts.
  - `MAF = min(p, 1 - p)` where `p` is alternate-allele frequency.
  - If no non-missing calls, return `None` and fail variant.
- Threshold policy:
  - If `maf_threshold` is unset, MAF filter is skipped.
  - Pass condition is `MAF >= threshold`.

## Edge Cases
- `maf == threshold` boundary behavior (must pass).
- All genotype calls missing (`MAF = None`, fail).
- Monomorphic variants (`MAF = 0.0`) under positive threshold.
- Extremely small non-missing call count.
- Variants previously marked malformed by IO/alignment path.

## Failure Modes
- Inconsistent genotype encoding that blocks allele-count computation.
- Missing metric prerequisites for a variant.
- Module should fail variant, not entire run, unless stream integrity is broken.

## Module Dependencies
- Upstream: Sample/variant alignment.
- Downstream: Association testing consumes QC-passed variants and QC metadata as needed.

## Test Scenarios
- MAF threshold off: all valid variants pass this module.
- `maf` boundary tests (`==`, below, above threshold).
- All-missing genotype variant is rejected with explicit reason.
- Monomorphic variant handling is deterministic.
- Pass/fail counts match expected fixture outcomes.
