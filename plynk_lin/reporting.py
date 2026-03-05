from __future__ import annotations

from plynk_lin.config import ParsedInputs


def _format_variant_preview(parsed: ParsedInputs, max_rows: int) -> list[str]:
    lines: list[str] = []
    for idx, variant in enumerate(parsed.vcf.variants[:max_rows], start=1):
        gvals = list(variant.genotypes_by_sample.values())
        missing = sum(val is None for val in gvals)
        lines.append(
            f"{idx}. {variant.chrom}:{variant.pos} {variant.variant_id} {variant.ref}>{variant.alt} "
            f"missing={missing}/{len(gvals)}"
        )
    return lines


def build_debug_summary(parsed: ParsedInputs, *, sample_preview: int = 5, variant_preview: int = 5) -> str:
    vcf = parsed.vcf
    pheno = parsed.pheno
    covar = parsed.covar

    total_gt = len(vcf.sample_ids) * len(vcf.variants)
    missing_gt = 0
    for variant in vcf.variants:
        missing_gt += sum(value is None for value in variant.genotypes_by_sample.values())

    lines: list[str] = [
        "== plynk-lin parse summary ==",
        f"VCF path: {vcf.path}",
        f"VCF samples: {len(vcf.sample_ids)}",
        f"VCF variants: {len(vcf.variants)}",
        f"VCF missing GT: {missing_gt}/{total_gt}" if total_gt else "VCF missing GT: 0/0",
        f"Pheno path: {pheno.path}",
        f"Pheno rows: {len(pheno.sample_ids)}",
        f"Pheno column: {pheno.phenotype_column}",
        f"Covar path: {covar.path if covar else '(none)'}",
        f"Covar rows: {len(covar.sample_ids) if covar else 0}",
        f"Covar columns: {', '.join(covar.covar_columns) if covar else '(none)'}",
        f"Warnings: {len(vcf.report.warnings) + len(pheno.report.warnings) + (len(covar.report.warnings) if covar else 0)}",
        f"Errors: {len(vcf.report.errors) + len(pheno.report.errors) + (len(covar.report.errors) if covar else 0)}",
    ]

    lines.append("Sample preview:")
    for idx, sid in enumerate(vcf.sample_ids[:sample_preview], start=1):
        lines.append(f"{idx}. {sid}")

    lines.append("Variant preview:")
    lines.extend(_format_variant_preview(parsed, variant_preview))

    return "\n".join(lines)

