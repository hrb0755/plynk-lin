from __future__ import annotations

from cyvcf2 import VCF

from plynk_lin.config import (
    CovarTable,
    InputParseError,
    MISSING_TOKENS,
    ParseReport,
    ParsedInputs,
    PhenoTable,
    RunConfig,
    VariantRecord,
    VcfDataset,
)


def _is_autosomal(chrom: str) -> bool:
    norm = chrom.removeprefix("chr").removeprefix("CHR")
    if not norm.isdigit():
        return False
    val = int(norm)
    return 1 <= val <= 22


def _parse_float(token: str) -> float | None:
    if token in MISSING_TOKENS:
        return None
    try:
        return float(token)
    except ValueError as exc:
        raise ValueError(f"Expected numeric token, got '{token}'") from exc


def _format_locus(path: str, chrom: str, pos: int) -> str:
    return f"{path}:{chrom}:{pos}"


def _parse_diploid_gt(
    gt_entry: list[int | bool], *, path: str, chrom: str, pos: int
) -> int | None:
    if len(gt_entry) < 2:
        raise InputParseError(
            f"GT must be diploid; got malformed GT entry at {_format_locus(path, chrom, pos)}"
        )

    a1 = gt_entry[0]
    a2 = gt_entry[1]
    if not isinstance(a1, int) or not isinstance(a2, int):
        raise InputParseError(
            f"GT alleles must be integer-coded at {_format_locus(path, chrom, pos)}"
        )
    if a1 < 0 or a2 < 0:
        return None
    if a1 not in {0, 1} or a2 not in {0, 1}:
        raise InputParseError(
            f"Only biallelic 0/1 GT supported; got '{a1}/{a2}' at {_format_locus(path, chrom, pos)}"
        )

    return a1 + a2


def load_vcf(path: str) -> VcfDataset:
    report = ParseReport(source=path)
    try:
        reader = VCF(path)
    except Exception as exc:
        raise InputParseError(f"Failed to open VCF '{path}': {exc}") from exc

    sample_ids = list(reader.samples or [])
    variants: list[VariantRecord] = []
    if not sample_ids:
        raise InputParseError(f"No sample columns found in VCF header for {path}")
    if len(set(sample_ids)) != len(sample_ids):
        raise InputParseError(f"Duplicate sample IDs in VCF header for {path}")

    try:
        for variant in reader:
            report.records_read += 1
            chrom = str(variant.CHROM)
            pos = int(variant.POS)
            if not _is_autosomal(chrom):
                raise InputParseError(
                    f"Unsupported non-autosomal chromosome '{chrom}' at {_format_locus(path, chrom, pos)}"
                )

            alt_list = list(variant.ALT or [])
            if len(alt_list) != 1 or alt_list[0] in {None, "."}:
                raise InputParseError(
                    f"Unsupported non-biallelic ALT '{variant.ALT}' at {_format_locus(path, chrom, pos)}"
                )

            fmt_fields = list(variant.FORMAT or [])
            if "GT" not in fmt_fields:
                raise InputParseError(f"FORMAT missing GT at {_format_locus(path, chrom, pos)}")

            genotypes_raw = list(variant.genotypes or [])
            if len(genotypes_raw) != len(sample_ids):
                raise InputParseError(
                    f"Sample genotype count mismatch at {_format_locus(path, chrom, pos)}"
                )

            genotypes: dict[str, int | None] = {}
            for sid, gt_entry in zip(sample_ids, genotypes_raw):
                dosage = _parse_diploid_gt(gt_entry, path=path, chrom=chrom, pos=pos)
                genotypes[sid] = dosage

            variants.append(
                VariantRecord(
                    chrom=chrom,
                    pos=pos,
                    variant_id=str(variant.ID or "."),
                    ref=str(variant.REF),
                    alt=str(alt_list[0]),
                    genotypes_by_sample=genotypes,
                )
            )
    except InputParseError:
        raise
    except Exception as exc:
        raise InputParseError(f"Failed while parsing VCF '{path}': {exc}") from exc
    finally:
        reader.close()

    return VcfDataset(path=path, sample_ids=sample_ids, variants=variants, report=report)


def _load_table_rows(path: str) -> tuple[list[str], list[tuple[int, list[str]]]]:
    header: list[str] | None = None
    rows: list[tuple[int, list[str]]] = []

    with open(path, "r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue
            cols = line.split()
            if header is None:
                header = cols
                continue
            rows.append((line_number, cols))

    if header is None:
        raise InputParseError(f"Empty file: {path}")
    if len(header) < 3:
        raise InputParseError(f"Expected at least 3 columns (FID IID VALUE...) in {path}")
    if header[0] != "FID" or header[1] != "IID":
        raise InputParseError(f"Expected header to start with 'FID IID' in {path}")

    return header, rows


def load_pheno(path: str) -> PhenoTable:
    report = ParseReport(source=path)
    header, rows = _load_table_rows(path)
    value_columns = header[2:]
    phenotype_column = "PHENO" if "PHENO" in value_columns else value_columns[0]
    pheno_idx = header.index(phenotype_column)

    sample_ids: list[str] = []
    values_by_sample: dict[str, float | None] = {}

    for line_number, cols in rows:
        report.records_read += 1
        if len(cols) != len(header):
            raise InputParseError(
                f"Column count mismatch at {path}:{line_number} (expected {len(header)}, got {len(cols)})"
            )
        iid = cols[1]
        if iid in values_by_sample:
            raise InputParseError(f"Duplicate IID '{iid}' at {path}:{line_number}")
        try:
            values_by_sample[iid] = _parse_float(cols[pheno_idx])
        except ValueError as exc:
            raise InputParseError(f"Invalid phenotype value at {path}:{line_number}: {exc}") from exc
        sample_ids.append(iid)

    return PhenoTable(
        path=path,
        sample_ids=sample_ids,
        values_by_sample=values_by_sample,
        phenotype_column=phenotype_column,
        report=report,
    )


def load_covar(path: str, covar_names: list[str] | None) -> CovarTable:
    report = ParseReport(source=path)
    header, rows = _load_table_rows(path)
    value_columns = header[2:]

    if covar_names is None:
        selected = list(value_columns)
    else:
        missing = [name for name in covar_names if name not in value_columns]
        if missing:
            raise InputParseError(
                f"Unknown covariate column(s) in {path}: {', '.join(missing)}"
            )
        selected = list(covar_names)

    selected_idx = {name: header.index(name) for name in selected}
    sample_ids: list[str] = []
    values_by_sample: dict[str, dict[str, float | None]] = {}

    for line_number, cols in rows:
        report.records_read += 1
        if len(cols) != len(header):
            raise InputParseError(
                f"Column count mismatch at {path}:{line_number} (expected {len(header)}, got {len(cols)})"
            )
        iid = cols[1]
        if iid in values_by_sample:
            raise InputParseError(f"Duplicate IID '{iid}' at {path}:{line_number}")

        row_vals: dict[str, float | None] = {}
        for name, idx in selected_idx.items():
            token = cols[idx]
            try:
                row_vals[name] = _parse_float(token)
            except ValueError as exc:
                raise InputParseError(
                    f"Invalid covariate value for '{name}' at {path}:{line_number}: {exc}"
                ) from exc

        values_by_sample[iid] = row_vals
        sample_ids.append(iid)

    return CovarTable(
        path=path,
        sample_ids=sample_ids,
        covar_columns=selected,
        values_by_sample=values_by_sample,
        report=report,
    )


def load_inputs(cfg: RunConfig) -> ParsedInputs:
    vcf = load_vcf(cfg.vcf_path)
    pheno = load_pheno(cfg.pheno_path)
    covar = load_covar(cfg.covar_path, cfg.covar_names) if cfg.covar_path else None
    return ParsedInputs(vcf=vcf, pheno=pheno, covar=covar)
