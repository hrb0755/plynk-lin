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
)


def _cohort() -> AlignedCohort:
    return AlignedCohort(
        sample_ids=["S1", "S2", "S3", "S4", "S5"],
        y=np.asarray([1.0, 1.4, 2.8, 3.7, 5.2], dtype=float),
        sample_index={"S1": 0, "S2": 1, "S3": 2, "S4": 3, "S5": 4},
        audit=None,  # type: ignore[arg-type]
    )


def test_fit_variant_ols_matches_manual_slope() -> None:
    y = np.asarray([1.0, 1.2, 3.1, 2.9, 5.2], dtype=float)
    g = np.asarray([0.0, 0.0, 1.0, 1.0, 2.0], dtype=float)
    fit = fit_variant_ols(y, g)

    x_mean = g.mean()
    y_mean = y.mean()
    beta = float(np.dot(g - x_mean, y - y_mean) / np.dot(g - x_mean, g - x_mean))

    assert fit.status == "ok"
    assert math.isclose(float(fit.beta[0]), beta, rel_tol=1e-9)
    assert fit.nmiss == 5


def test_run_linear_assoc_emits_add_row() -> None:
    cohort = _cohort()
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
    summary = AssocSummary()

    rows = list(run_linear_assoc(cohort, filtered, summary=summary))

    assert [row.test for row in rows] == ["ADD"]
    assert all(row.nmiss == 5 for row in rows)
    assert summary.emitted_rows == 1


def test_run_linear_assoc_flips_a1_when_alt_is_major() -> None:
    cohort = _cohort()
    filtered = [
        FilteredVariant(
            variant=AlignedVariant(
                chrom="1",
                pos=100,
                variant_id="rs1",
                ref="A",
                alt="G",
                a1="G",
                g=np.asarray([1.0, 1.0, 2.0, 2.0, 2.0], dtype=float),
            ),
            qc=QcDecision(True, [], 0.4, 4),
        )
    ]

    rows = list(run_linear_assoc(cohort, filtered))

    assert [row.test for row in rows] == ["ADD"]
    assert rows[0].a1 == "A"
