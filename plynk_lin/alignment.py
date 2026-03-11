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

    retained_ids: list[str] = []
    y_values: list[float] = []
    not_in_pheno = 0
    missing_pheno = 0

    for sid in vcf.sample_ids:
        if sid not in pheno.values_by_sample:
            not_in_pheno += 1
            continue

        pheno_value = pheno.values_by_sample[sid]
        if pheno_value is None:
            missing_pheno += 1
            continue

        retained_ids.append(sid)
        y_values.append(float(pheno_value))

    if not retained_ids:
        raise InputParseError("No overlapping samples remain after alignment and listwise deletion")

    sample_index = {sid: idx for idx, sid in enumerate(retained_ids)}
    audit = AlignmentAudit(
        not_in_pheno=not_in_pheno,
        missing_pheno=missing_pheno,
        retained=len(retained_ids),
    )
    cohort = AlignedCohort(
        sample_ids=retained_ids,
        y=np.asarray(y_values, dtype=float),
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
