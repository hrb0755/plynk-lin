from __future__ import annotations

import math

import numpy as np

from plynk_lin.association import fit_variant_ols, run_linear_assoc
from plynk_lin.config import (
    AlignedCohort,
    AlignedVariant,
    AssocSummary,
    FilteredVariant,
    QcDecision,
    RunConfig,
)


def _cohort_with_covar() -> AlignedCohort:
    return AlignedCohort(
        sample_ids=["S1", "S2", "S3", "S4", "S5"],
        y=np.asarray([1.0, 1.4, 2.8, 3.7, 5.2], dtype=float),
        X_covar=np.asarray(
            [
                [30.0, 0.10],
                [31.0, 0.20],
                [29.0, -0.10],
                [35.0, 0.30],
                [36.0, 0.25],
            ],
            dtype=float,
        ),
        covar_names=["AGE", "PC1"],
        sample_index={"S1": 0, "S2": 1, "S3": 2, "S4": 3, "S5": 4},
        audit=None,  # type: ignore[arg-type]
    )


def test_fit_variant_ols_no_covar_matches_manual_slope() -> None:
    y = np.asarray([1.0, 1.2, 3.1, 2.9, 5.2], dtype=float)
    g = np.asarray([0.0, 0.0, 1.0, 1.0, 2.0], dtype=float)
    fit = fit_variant_ols(y, np.ones((5, 1), dtype=float), g, ["ADD"])

    x_mean = g.mean()
    y_mean = y.mean()
    beta = float(np.dot(g - x_mean, y - y_mean) / np.dot(g - x_mean, g - x_mean))

    assert fit.status == "ok"
    assert math.isclose(float(fit.beta[0]), beta, rel_tol=1e-9)
    assert fit.nmiss == 5


def test_run_linear_assoc_emits_covar_rows_when_visible() -> None:
    cohort = _cohort_with_covar()
    filtered = [
        FilteredVariant(
            variant=AlignedVariant(
                chrom="1",
                pos=100,
                variant_id="rs1",
                ref="A",
                alt="G",
                a1="G",
                g=np.asarray([0.0, 0.0, 1.0, 1.0, 2.0], dtype=float),
            ),
            qc=QcDecision(True, [], 0.4, 5),
        )
    ]
    cfg = RunConfig(True, False, "x.vcf", "p.txt", "c.txt", ["AGE", "PC1"], None, False, "out")
    summary = AssocSummary()

    rows = list(run_linear_assoc(cohort, filtered, cfg, summary=summary))

    assert [row.test for row in rows] == ["ADD", "AGE", "PC1"]
    assert all(row.nmiss == 5 for row in rows)
    assert summary.emitted_rows == 3


def test_run_linear_assoc_hide_covar_emits_only_add() -> None:
    cohort = _cohort_with_covar()
    filtered = [
        FilteredVariant(
            variant=AlignedVariant(
                chrom="1",
                pos=100,
                variant_id="rs1",
                ref="A",
                alt="G",
                a1="G",
                g=np.asarray([0.0, 0.0, 1.0, 1.0, 2.0], dtype=float),
            ),
            qc=QcDecision(True, [], 0.4, 4),
        )
    ]
    cfg = RunConfig(True, True, "x.vcf", "p.txt", "c.txt", ["AGE", "PC1"], None, False, "out")

    rows = list(run_linear_assoc(cohort, filtered, cfg))

    assert [row.test for row in rows] == ["ADD"]
    assert rows[0].nmiss == 5
