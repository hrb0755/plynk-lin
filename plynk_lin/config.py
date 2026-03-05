from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional


MISSING_TOKENS = {"", ".", "NA", "nan", "-9"}


class PlynkError(Exception):
    """Base error type for parser + IO slice."""


class ConfigError(PlynkError):
    """CLI/config validation failure."""


class InputParseError(PlynkError):
    """Input file parse failure with context."""


@dataclass(frozen=True)
class RunConfig:
    linear_enabled: bool
    hide_covar: bool
    vcf_path: str
    pheno_path: str
    covar_path: Optional[str]
    covar_names: Optional[list[str]]
    maf_threshold: Optional[float]
    allow_no_sex: bool
    out_prefix: str
    debug: bool = False
    debug_preview: int = 5
    debug_variants: int = 5


@dataclass(frozen=True)
class ParseWarning:
    message: str
    path: Optional[str] = None
    line_number: Optional[int] = None


@dataclass(frozen=True)
class ParseErrorInfo:
    message: str
    path: Optional[str] = None
    line_number: Optional[int] = None


@dataclass
class ParseReport:
    source: str
    records_read: int = 0
    warnings: list[ParseWarning] = field(default_factory=list)
    errors: list[ParseErrorInfo] = field(default_factory=list)


@dataclass(frozen=True)
class VariantRecord:
    chrom: str
    pos: int
    variant_id: str
    ref: str
    alt: str
    genotypes_by_sample: dict[str, Optional[int]]


@dataclass
class VcfDataset:
    path: str
    sample_ids: list[str]
    variants: list[VariantRecord]
    report: ParseReport


@dataclass
class PhenoTable:
    path: str
    sample_ids: list[str]
    values_by_sample: dict[str, Optional[float]]
    phenotype_column: str
    report: ParseReport


@dataclass
class CovarTable:
    path: str
    sample_ids: list[str]
    covar_columns: list[str]
    values_by_sample: dict[str, dict[str, Optional[float]]]
    report: ParseReport


@dataclass
class ParsedInputs:
    vcf: VcfDataset
    pheno: PhenoTable
    covar: Optional[CovarTable]

