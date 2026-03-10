from __future__ import annotations

from plynk_lin.config import (
    AlignedCohort,
    AssocSummary,
    ParsedInputs,
    QcSummary,
    RunConfig,
    WriteSummary,
)


def build_debug_summary(
    cfg: RunConfig,
    parsed: ParsedInputs,
    cohort: AlignedCohort,
    qc_summary: QcSummary,
    assoc_summary: AssocSummary,
    write_summary: WriteSummary,
) -> str:
    covar = parsed.covar
    audit = cohort.audit
    lines = [
        "== plynk-lin pipeline summary ==",
        f"VCF path: {parsed.vcf.path}",
        f"VCF samples: {len(parsed.vcf.sample_ids)}",
        f"Pheno path: {parsed.pheno.path}",
        f"Pheno rows: {len(parsed.pheno.sample_ids)}",
        f"Covar path: {covar.path if covar else '(none)'}",
        f"Covar rows: {len(covar.sample_ids) if covar else 0}",
        f"Covar columns: {', '.join(cohort.covar_names) if cohort.covar_names else '(none)'}",
        f"Aligned samples: {audit.retained}",
        f"Dropped (not in pheno): {audit.not_in_pheno}",
        f"Dropped (not in covar): {audit.not_in_covar}",
        f"Dropped (missing pheno): {audit.missing_pheno}",
        f"Dropped (missing covar): {audit.missing_covar}",
        f"QC total variants: {qc_summary.total_variants}",
        f"QC passed variants: {qc_summary.passed_variants}",
        f"QC failed variants: {qc_summary.failed_variants}",
        f"Association processed variants: {assoc_summary.processed_variants}",
        f"Association fit failures: {assoc_summary.fit_failures}",
        f"Rows written: {write_summary.row_count}",
        f"Output path: {write_summary.output_path}",
        f"hide-covar: {cfg.hide_covar}",
        f"maf threshold: {cfg.maf_threshold if cfg.maf_threshold is not None else '(none)'}",
    ]
    return "\n".join(lines)
