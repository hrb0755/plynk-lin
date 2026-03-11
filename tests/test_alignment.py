from __future__ import annotations

import pytest

from plynk_lin.alignment import align_samples
from plynk_lin.config import (
    InputParseError,
    ParseReport,
    ParsedInputs,
    PhenoTable,
    VariantRecord,
    VcfDataset,
)


def _make_vcf(sample_ids: list[str], variants: list[VariantRecord]) -> VcfDataset:
    return VcfDataset(
        path="synthetic.vcf",
        sample_ids=sample_ids,
        report=ParseReport(source="synthetic.vcf"),
        _variant_iter_factory=lambda: iter(variants),
    )


def test_align_samples_uses_vcf_order_and_listwise_deletion() -> None:
    parsed = ParsedInputs(
        vcf=_make_vcf(
            ["S1", "S2", "S3", "S4"],
            [
                VariantRecord(
                    chrom="1",
                    pos=100,
                    variant_id="rs1",
                    ref="A",
                    alt="G",
                    genotypes_by_sample={"S1": 0, "S2": 1, "S3": None, "S4": 2},
                )
            ],
        ),
        pheno=PhenoTable(
            path="pheno.txt",
            sample_ids=["S2", "S1", "S3"],
            values_by_sample={"S2": 2.0, "S1": 1.0, "S3": None},
            phenotype_column="PHENO",
            report=ParseReport(source="pheno.txt"),
        ),
    )

    aligned = align_samples(parsed)

    assert aligned.cohort.sample_ids == ["S1", "S2"]
    assert aligned.cohort.audit.not_in_pheno == 1
    assert aligned.cohort.audit.missing_pheno == 1
    assert aligned.cohort.audit.retained == 2

    variant = next(aligned.iter_variants())
    assert variant.g.tolist() == [0.0, 1.0]


def test_align_samples_ignores_pheno_only_missing_rows() -> None:
    parsed = ParsedInputs(
        vcf=_make_vcf(
            ["S1", "S2"],
            [
                VariantRecord(
                    chrom="1",
                    pos=100,
                    variant_id="rs1",
                    ref="A",
                    alt="G",
                    genotypes_by_sample={"S1": 0, "S2": 1},
                )
            ],
        ),
        pheno=PhenoTable(
            path="pheno.txt",
            sample_ids=["S1", "S2"],
            values_by_sample={"S1": None, "S2": 2.0},
            phenotype_column="PHENO",
            report=ParseReport(source="pheno.txt"),
        ),
    )

    aligned = align_samples(parsed)

    assert aligned.cohort.sample_ids == ["S2"]
    assert aligned.cohort.audit.missing_pheno == 1


def test_align_samples_empty_overlap_fails() -> None:
    parsed = ParsedInputs(
        vcf=_make_vcf(
            ["S1"],
            [
                VariantRecord(
                    chrom="1",
                    pos=100,
                    variant_id="rs1",
                    ref="A",
                    alt="G",
                    genotypes_by_sample={"S1": 0},
                )
            ],
        ),
        pheno=PhenoTable(
            path="pheno.txt",
            sample_ids=["X1"],
            values_by_sample={"X1": 1.0},
            phenotype_column="PHENO",
            report=ParseReport(source="pheno.txt"),
        ),
    )

    with pytest.raises(InputParseError, match="No overlapping samples"):
        align_samples(parsed)
