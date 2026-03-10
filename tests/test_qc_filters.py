from __future__ import annotations

import numpy as np

from plynk_lin.config import AlignedVariant, RunConfig
from plynk_lin.qc_filters import apply_variant_filters, compute_maf


def _cfg(maf: float | None) -> RunConfig:
    return RunConfig(True, False, "x.vcf", "p.txt", None, None, maf, False, "out")


def test_compute_maf_boundary() -> None:
    g = np.asarray([0.0, 1.0, 2.0, np.nan], dtype=float)
    assert compute_maf(g) == 0.5


def test_apply_variant_filters_handles_missing_and_threshold() -> None:
    variant = AlignedVariant("1", 100, "rs1", "A", "G", "G", np.asarray([0.0, 0.0, 0.0], dtype=float))
    decision = apply_variant_filters(variant, _cfg(0.1))
    assert decision.passed is False
    assert decision.maf_value == 0.0
    assert decision.non_missing_genotype_count == 3


def test_apply_variant_filters_all_missing() -> None:
    variant = AlignedVariant("1", 100, "rs1", "A", "G", "G", np.asarray([np.nan, np.nan], dtype=float))
    decision = apply_variant_filters(variant, _cfg(None))
    assert decision.passed is False
    assert decision.reasons == ["all_genotypes_missing"]
    assert decision.maf_value is None
