from __future__ import annotations

from cyvcf2 import VCF

from plynk_lin.config import (
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
    gt_entry: list[int | bool], *, path: str, chrom: str, pos: int, alt_count: int, alt_index: int
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
    if alt_count > 1 and ((a1 not in {0, alt_index}) or (a2 not in {0, alt_index})):
        return None
    if a1 not in {0, alt_index} or a2 not in {0, alt_index}:
        raise InputParseError(
            f"Only biallelic 0/1 GT supported; got '{a1}/{a2}' at {_format_locus(path, chrom, pos)}"
        )

    return int(a1 == alt_index) + int(a2 == alt_index)


def _variant_from_cyvcf2(variant, *, path: str, sample_ids: list[str]) -> VariantRecord | None:
    chrom = str(variant.CHROM)
    pos = int(variant.POS)
    if not _is_autosomal(chrom):
        return None

    alt_list = list(variant.ALT or [])
    if not alt_list or alt_list[0] in {None, "."}:
        return None

    fmt_fields = list(variant.FORMAT or [])
    if "GT" not in fmt_fields:
        raise InputParseError(f"FORMAT missing GT at {_format_locus(path, chrom, pos)}")

    genotypes_raw = list(variant.genotypes or [])
    if len(genotypes_raw) != len(sample_ids):
        raise InputParseError(f"Sample genotype count mismatch at {_format_locus(path, chrom, pos)}")

    alt_index = 1
    if len(alt_list) > 1:
        allele_counts = {idx: 0 for idx in range(1, len(alt_list) + 1)}
        for gt_entry in genotypes_raw:
            if len(gt_entry) < 2:
                continue
            for allele in gt_entry[:2]:
                if isinstance(allele, int) and allele > 0 and allele in allele_counts:
                    allele_counts[allele] += 1
        alt_index = max(allele_counts, key=lambda idx: (allele_counts[idx], -idx))

    genotypes: dict[str, int | None] = {}
    for sid, gt_entry in zip(sample_ids, genotypes_raw):
        dosage = _parse_diploid_gt(
            gt_entry,
            path=path,
            chrom=chrom,
            pos=pos,
            alt_count=len(alt_list),
            alt_index=alt_index,
        )
        genotypes[sid] = dosage

    return VariantRecord(
        chrom=chrom,
        pos=pos,
        variant_id=str(variant.ID or "."),
        ref=str(variant.REF),
        alt=str(alt_list[alt_index - 1]),
        genotypes_by_sample=genotypes,
    )


def _iter_vcf_variants(path: str, sample_ids: list[str]):
    reader = VCF(path)
    try:
        for variant in reader:
            record = _variant_from_cyvcf2(variant, path=path, sample_ids=sample_ids)
            if record is not None:
                yield record
    except InputParseError:
        raise
    except Exception as exc:
        raise InputParseError(f"Failed while parsing VCF '{path}': {exc}") from exc
    finally:
        reader.close()


def load_vcf(path: str) -> VcfDataset:
    report = ParseReport(source=path)
    try:
        reader = VCF(path)
    except Exception as exc:
        raise InputParseError(f"Failed to open VCF '{path}': {exc}") from exc

    try:
        sample_ids = list(reader.samples or [])
        if not sample_ids:
            raise InputParseError(f"No sample columns found in VCF header for {path}")
        if len(set(sample_ids)) != len(sample_ids):
            raise InputParseError(f"Duplicate sample IDs in VCF header for {path}")
    finally:
        reader.close()

    return VcfDataset(
        path=path,
        sample_ids=sample_ids,
        report=report,
        _variant_iter_factory=lambda: _iter_vcf_variants(path, sample_ids),
    )


def _read_nonempty_rows(path: str) -> list[tuple[int, list[str]]]:
    rows: list[tuple[int, list[str]]] = []
    with open(path, "r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue
            rows.append((line_number, line.split()))
    if not rows:
        raise InputParseError(f"Empty file: {path}")
    return rows


def load_pheno(path: str) -> PhenoTable:
    report = ParseReport(source=path)
    rows = _read_nonempty_rows(path)

    first_line, first_cols = rows[0]
    has_header = len(first_cols) >= 3 and first_cols[0] == "FID" and first_cols[1] == "IID"
    if has_header:
        header = first_cols
        data_rows = rows[1:]
        if len(header) < 3:
            raise InputParseError(f"Expected at least 3 columns (FID IID VALUE...) in {path}")
        value_columns = header[2:]
        phenotype_column = "PHENO" if "PHENO" in value_columns else value_columns[0]
        pheno_idx = header.index(phenotype_column)
    else:
        header = ["FID", "IID", "PHENO"]
        data_rows = rows
        phenotype_column = "PHENO"
        pheno_idx = 2

    sample_ids: list[str] = []
    values_by_sample: dict[str, float | None] = {}

    for line_number, cols in data_rows:
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


def load_inputs(cfg: RunConfig) -> ParsedInputs:
    vcf = load_vcf(cfg.vcf_path)
    pheno = load_pheno(cfg.pheno_path)
    return ParsedInputs(vcf=vcf, pheno=pheno)
