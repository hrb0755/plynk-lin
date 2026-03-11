from __future__ import annotations

from typing import Iterable, Iterator

import numpy as np
from scipy.stats import t as student_t

from plynk_lin.config import AlignedCohort, AssocResultRow, AssocSummary, FilteredVariant, FitStats


def _fit_additive_ols(y: np.ndarray, g: np.ndarray) -> FitStats:
    mask = ~np.isnan(g)
    y_use = y[mask]
    g_use = g[mask]
    nmiss = int(mask.sum())
    df = nmiss - 2
    if nmiss < 3:
        return FitStats(["ADD"], np.array([]), np.array([]), np.array([]), np.array([]), nmiss, df, "too_few_rows")

    g_mean = float(np.mean(g_use))
    y_mean = float(np.mean(y_use))
    g_centered = g_use - g_mean
    sxx = float(np.dot(g_centered, g_centered))
    if sxx <= 0.0:
        return FitStats(["ADD"], np.array([]), np.array([]), np.array([]), np.array([]), nmiss, df, "zero_variance")
    if df <= 0:
        return FitStats(["ADD"], np.array([]), np.array([]), np.array([]), np.array([]), nmiss, df, "nonpositive_df")

    beta_add = float(np.dot(g_centered, y_use - y_mean) / sxx)
    intercept = y_mean - beta_add * g_mean
    resid = y_use - (intercept + beta_add * g_use)
    rss = float(np.dot(resid, resid))
    sigma2 = rss / df
    se_add = float(np.sqrt(sigma2 / sxx))
    if se_add == 0.0:
        return FitStats(["ADD"], np.array([]), np.array([]), np.array([]), np.array([]), nmiss, df, "zero_se")
    stat_add = beta_add / se_add
    p_add = float(2.0 * student_t.sf(abs(stat_add), df))

    return FitStats(
        term_names=["ADD"],
        beta=np.asarray([beta_add], dtype=float),
        se=np.asarray([se_add], dtype=float),
        stat=np.asarray([stat_add], dtype=float),
        p=np.asarray([p_add], dtype=float),
        nmiss=nmiss,
        df=df,
        status="ok",
    )


def fit_variant_ols(y: np.ndarray, g: np.ndarray) -> FitStats:
    return _fit_additive_ols(y, g)


def run_linear_assoc(
    cohort: AlignedCohort,
    variants: Iterable[FilteredVariant],
    summary: AssocSummary | None = None,
) -> Iterator[AssocResultRow]:
    stats = summary if summary is not None else AssocSummary()
    for filtered in variants:
        stats.processed_variants += 1
        fit = fit_variant_ols(cohort.y, filtered.variant.g)
        if fit.status != "ok":
            stats.fit_failures += 1
            continue

        mask = ~np.isnan(filtered.variant.g)
        alt_freq = float(np.sum(filtered.variant.g[mask]) / (2.0 * fit.nmiss))
        sign = 1.0
        a1 = filtered.variant.alt
        if alt_freq > 0.5:
            sign = -1.0
            a1 = filtered.variant.ref

        row = AssocResultRow(
            chrom=filtered.variant.chrom,
            snp=filtered.variant.variant_id,
            bp=filtered.variant.pos,
            a1=a1,
            test="ADD",
            nmiss=fit.nmiss,
            beta=float(sign * fit.beta[0]),
            stat=float(sign * fit.stat[0]),
            p=float(fit.p[0]),
        )
        stats.emitted_rows += 1
        yield row
