# Sample/Variant Alignment Module Design

## Overview
This module aligns samples across VCF and phenotype inputs using sample IDs, applies deterministic ordering, and builds analysis-ready structures for QC and association testing.

Implementation note:
- Alignment outputs will be materialized with `numpy` arrays for efficient downstream masking, ordering, and matrix construction.
- If table-style joins or reshaping become more complex during implementation, a lightweight `pandas` layer may be used at the file-to-array boundary before converting to the canonical in-memory structures used by the pipeline.

## Responsibilities
- Resolve sample identity across parsed inputs by sample ID.
- Define a deterministic retained-sample order.
- Construct aligned phenotype vector.
- Construct per-variant genotype vectors aligned to retained samples.
- Apply listwise deletion policy for analysis rows requiring complete data.

Not owned by this module:
- File parsing details.
- MAF computation/filter pass-fail logic.
- OLS coefficient/statistic computation.
- Output table formatting.

## Inputs
- `VcfDataset` from IO.
- `PhenoTable` from IO.

## Outputs
- `AlignedCohort`:
  - `sample_ids: list[str]` retained for analysis.
  - `y: vector[float]` phenotype values in retained order.
- `AlignedVariantStream`:
  - For each variant, `g: vector[float | missing]` in retained sample order.
  - Variant metadata passthrough.
- Alignment audit metadata:
  - Counts dropped for missing/mismatch/duplicate causes.

## Interface Contract
- `align_samples(parsed: ParsedInputs) -> AlignedInputs`
- `build_variant_genotype_view(variant, aligned_sample_ids) -> AlignedVariant`
- Join and ordering policy:
  - Sample keying by canonical sample ID string.
  - Deterministic ordering inherited from VCF sample order after intersection.
  - Required phenotype presence for retained samples.
- Listwise deletion policy:
  - Retain rows only when phenotype is non-missing.
  - Genotype missingness is preserved per variant for NMISS handling downstream.

## Edge Cases
- IDs present in phenotype but absent from VCF and vice versa.
- Duplicate IDs across any source.
- All samples dropped after intersection or listwise deletion.
- Variant with all genotype calls missing among retained samples.

## Failure Modes
- Unresolvable ID schema mismatch (no overlapping samples).
- Duplicate-ID ambiguity without deterministic resolution policy.
- Empty aligned cohort after required deletions.
- Failures include concise dropped-count diagnostics.

## Module Dependencies
- Upstream: IO for parsed datasets.
- Downstream:
  - QC consumes aligned variant genotype vectors and cohort size context.
  - Association testing consumes aligned `y` and per-variant `g`.
- External dependencies:
  - `numpy` (primary array container for aligned phenotype and genotype data).
  - `pandas` (optional helper for ID-based joins and table reshaping if needed during implementation).

## Test Scenarios
- Shuffled phenotype row order still yields correct alignment by ID.
- Known mismatched-ID dataset drops non-intersecting samples deterministically.
- Listwise deletion removes rows with missing phenotype only.
- Variant genotype vectors reflect retained sample order exactly.
- Alignment audit counts match expected drop reasons.
