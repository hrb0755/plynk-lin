# Output Writing Module Design

## Overview
This module serializes association result rows into PLINK-like `.assoc.linear` output using the configured `--out` prefix, with deterministic ordering and formatting rules.

Implementation note:
- Output writing is expected to rely primarily on Python standard-library file I/O.
- If result assembly benefits from tabular staging before serialization, an off-the-shelf table library such as `pandas` may be used as a convenience layer, but it is not required for the core writer.

## Responsibilities
- Materialize output path from `--out` prefix.
- Write `.assoc.linear` header and rows in defined column order.
- Apply consistent numeric/text formatting, including NA conventions.
- Enforce deterministic row ordering from upstream stream order.
- Surface write-time errors with actionable diagnostics.

Not owned by this module:
- Statistical computation.
- Variant filtering logic.
- Upstream parsing/alignment validation.

## Inputs
- `AssocResultStream` from association module.
- `RunConfig.out_prefix`.
- Optional run metadata for header comments if adopted (not required by current proposal).

## Outputs
- File artifact: `<out_prefix>.assoc.linear`.
- Write summary metadata:
  - output path
  - row count
  - success/failure status

## Interface Contract
- `write_assoc_linear(results: AssocResultStream, out_prefix: str) -> WriteSummary`
- Output schema policy:
  - Include required PLINK-like columns in fixed order: `CHR SNP BP A1 TEST NMISS BETA STAT P`.
  - Serialize the public association row fields only; internal fit diagnostics such as SE/DF are not written to `.assoc.linear`.
  - Use deterministic whitespace-delimited text with one header row.
- Formatting policy:
  - Numeric formatting is consistent run-wide.
  - Missing/unavailable values use a single NA token policy.
  - Line termination and encoding are deterministic.

## Edge Cases
- Empty result stream should still write a valid header-only output.
- Existing output file path collision policy (overwrite or fail) must be explicit and deterministic.
- Rows containing NA statistics after fit failures.
- Output prefix containing nested/nonexistent directories.

## Failure Modes
- Directory creation or file-open failure.
- Mid-write I/O interruption.
- Schema mismatch between expected columns and provided rows.
- On failure, module returns/raises with output path context and row-progress context when available.

## Module Dependencies
- Upstream: Association testing result stream and run config.
- Downstream: none (terminal stage of pipeline).
- External dependencies:
  - Standard library I/O is expected to be sufficient for the core writer.
  - `pandas` is optional if tabular result staging simplifies formatting or export logic.

## Test Scenarios
- Happy-path write creates `<out_prefix>.assoc.linear` with expected header and row count.
- Header-only output on no passing variants.
- NA/stat-failure rows serialize without column shifts.
- Write failure path surfaces clear filesystem error context.
