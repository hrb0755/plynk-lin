from __future__ import annotations

from typing import Iterator

import numpy as np

from plynk_lin.config import (
    AlignedCohort,
    AlignedInputs,
    AlignedVariant,
    AlignmentAudit,
    InputParseError,
    ParsedInputs,
)


def align_samples(parsed: ParsedInputs) -> AlignedInputs:
    vcf = parsed.vcf
    pheno = parsed.pheno
    covar = parsed.covar

    retained_ids: list[str] = []
    y_values: list[float] = []
    covar_rows: list[list[float]] = []
    not_in_pheno = 0
    not_in_covar = 0
    missing_pheno = 0
    missing_covar = 0

    for sid in vcf.sample_ids:
        if sid not in pheno.values_by_sample:
            not_in_pheno += 1
            continue

        pheno_value = pheno.values_by_sample[sid]
        if pheno_value is None:
            missing_pheno += 1
            continue

        covar_values: list[float] = []
        if covar is not None:
            covar_row = covar.values_by_sample.get(sid)
            if covar_row is None:
                not_in_covar += 1
                continue

            has_missing = False
            for name in covar.covar_columns:
                value = covar_row[name]
                if value is None:
                    has_missing = True
                    break
                covar_values.append(float(value))
            if has_missing:
                missing_covar += 1
                continue

        retained_ids.append(sid)
        y_values.append(float(pheno_value))
        if covar is not None:
            covar_rows.append(covar_values)

    if not retained_ids:
        raise InputParseError("No overlapping samples remain after alignment and listwise deletion")

    sample_index = {sid: idx for idx, sid in enumerate(retained_ids)}
    audit = AlignmentAudit(
        not_in_pheno=not_in_pheno,
        not_in_covar=not_in_covar,
        missing_pheno=missing_pheno,
        missing_covar=missing_covar,
        retained=len(retained_ids),
    )
    cohort = AlignedCohort(
        sample_ids=retained_ids,
        y=np.asarray(y_values, dtype=float),
        X_covar=np.asarray(covar_rows, dtype=float) if covar is not None else None,
        covar_names=list(covar.covar_columns) if covar is not None else [],
        sample_index=sample_index,
        audit=audit,
    )

    def iter_aligned_variants() -> Iterator[AlignedVariant]:
        for variant in vcf.iter_variants():
            g = np.full(len(retained_ids), np.nan, dtype=float)
            for sid, idx in sample_index.items():
                dosage = variant.genotypes_by_sample.get(sid)
                if dosage is not None:
                    g[idx] = float(dosage)
            yield AlignedVariant(
                chrom=variant.chrom,
                pos=variant.pos,
                variant_id=variant.variant_id,
                ref=variant.ref,
                alt=variant.alt,
                a1=variant.alt,
                g=g,
            )

    return AlignedInputs(cohort=cohort, _variant_iter_factory=iter_aligned_variants)
