from __future__ import annotations

from typing import Iterable, Iterator

import numpy as np

from plynk_lin.config import FilteredVariant, QcDecision, QcSummary, RunConfig


def compute_maf(g: np.ndarray) -> float | None:
    mask = ~np.isnan(g)
    non_missing = int(mask.sum())
    if non_missing == 0:
        return None
    alt_freq = float(np.sum(g[mask]) / (2.0 * non_missing))
    return min(alt_freq, 1.0 - alt_freq)


def apply_variant_filters(variant, cfg: RunConfig) -> QcDecision:
    mask = ~np.isnan(variant.g)
    non_missing = int(mask.sum())
    maf = compute_maf(variant.g)
    reasons: list[str] = []

    if maf is None:
        reasons.append("all_genotypes_missing")
    if cfg.maf_threshold is not None and maf is not None and maf < cfg.maf_threshold:
        reasons.append(f"maf_below_threshold:{cfg.maf_threshold}")

    return QcDecision(
        passed=not reasons,
        reasons=reasons,
        maf_value=maf,
        non_missing_genotype_count=non_missing,
    )


def filter_variants(
    stream: Iterable, cfg: RunConfig, summary: QcSummary | None = None
) -> Iterator[FilteredVariant]:
    stats = summary if summary is not None else QcSummary()
    for variant in stream:
        stats.total_variants += 1
        decision = apply_variant_filters(variant, cfg)
        if decision.passed:
            stats.passed_variants += 1
            yield FilteredVariant(variant=variant, qc=decision)
        else:
            stats.failed_variants += 1
