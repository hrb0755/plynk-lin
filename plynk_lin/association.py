from __future__ import annotations

from typing import Iterable, Iterator

import numpy as np
from scipy.stats import t as student_t

from plynk_lin.config import (
    AlignedCohort,
    AssocResultRow,
    AssocSummary,
    FilteredVariant,
    FitStats,
    RunConfig,
)


def _fit_no_covar(y: np.ndarray, g: np.ndarray) -> FitStats:
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


def fit_variant_ols(y: np.ndarray, X_base: np.ndarray, g: np.ndarray, term_names: list[str]) -> FitStats:
    mask = ~np.isnan(g)
    y_use = y[mask]
    g_use = g[mask]
    X_base_use = X_base[mask]

    nmiss = int(mask.sum())
    n_params = X_base_use.shape[1] + 1
    df = nmiss - n_params
    if nmiss <= n_params or df <= 0:
        return FitStats(term_names, np.array([]), np.array([]), np.array([]), np.array([]), nmiss, df, "nonpositive_df")

    if X_base_use.shape[1] == 1:
        return _fit_no_covar(y, g)

    X = np.column_stack([X_base_use[:, :1], g_use, X_base_use[:, 1:]])

    beta, _, rank, _ = np.linalg.lstsq(X, y_use, rcond=None)
    if rank < X.shape[1]:
        return FitStats(term_names, np.array([]), np.array([]), np.array([]), np.array([]), nmiss, df, "rank_deficient")

    fitted = X @ beta
    resid = y_use - fitted
    rss = float(np.dot(resid, resid))
    sigma2 = rss / df

    try:
        xtx_inv = np.linalg.inv(X.T @ X)
    except np.linalg.LinAlgError:
        return FitStats(term_names, np.array([]), np.array([]), np.array([]), np.array([]), nmiss, df, "singular")

    se = np.sqrt(np.diag(sigma2 * xtx_inv))
    term_beta = beta[1:]
    term_se = se[1:]
    if np.any(term_se == 0.0):
        return FitStats(term_names, np.array([]), np.array([]), np.array([]), np.array([]), nmiss, df, "zero_se")
    stat = term_beta / term_se
    p = 2.0 * student_t.sf(np.abs(stat), df)

    return FitStats(
        term_names=term_names,
        beta=np.asarray(term_beta, dtype=float),
        se=np.asarray(term_se, dtype=float),
        stat=np.asarray(stat, dtype=float),
        p=np.asarray(p, dtype=float),
        nmiss=nmiss,
        df=df,
        status="ok",
    )


def run_linear_assoc(
    cohort: AlignedCohort,
    variants: Iterable[FilteredVariant],
    cfg: RunConfig,
    summary: AssocSummary | None = None,
) -> Iterator[AssocResultRow]:
    stats = summary if summary is not None else AssocSummary()
    X_base = np.ones((len(cohort.sample_ids), 1), dtype=float)
    if cohort.X_covar is not None:
        X_base = np.column_stack([X_base, cohort.X_covar])

    term_names = ["ADD", *cohort.covar_names]
    for filtered in variants:
        stats.processed_variants += 1
        fit = fit_variant_ols(cohort.y, X_base, filtered.variant.g, term_names)
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

        rows_to_emit = range(len(fit.term_names))
        if cfg.hide_covar or not cohort.covar_names:
            rows_to_emit = range(1)

        for idx in rows_to_emit:
            row = AssocResultRow(
                chrom=filtered.variant.chrom,
                snp=filtered.variant.variant_id,
                bp=filtered.variant.pos,
                a1=a1,
                test=fit.term_names[idx],
                nmiss=fit.nmiss,
                beta=float(sign * fit.beta[idx]),
                stat=float(sign * fit.stat[idx]),
                p=float(fit.p[idx]),
            )
            stats.emitted_rows += 1
            yield row
