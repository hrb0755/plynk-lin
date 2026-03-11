from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Iterator, Optional

import numpy as np


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
    vcf_path: str
    pheno_path: str
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
    report: ParseReport
    _variant_iter_factory: Callable[[], Iterator[VariantRecord]] = field(repr=False)
    _variants_cache: list[VariantRecord] | None = field(default=None, init=False, repr=False)

    def iter_variants(self) -> Iterator[VariantRecord]:
        return self._variant_iter_factory()

    @property
    def variants(self) -> list[VariantRecord]:
        if self._variants_cache is None:
            self._variants_cache = list(self.iter_variants())
        return self._variants_cache


@dataclass
class PhenoTable:
    path: str
    sample_ids: list[str]
    values_by_sample: dict[str, Optional[float]]
    phenotype_column: str
    report: ParseReport


@dataclass
class ParsedInputs:
    vcf: VcfDataset
    pheno: PhenoTable


@dataclass(frozen=True)
class AlignmentAudit:
    not_in_pheno: int
    missing_pheno: int
    retained: int


@dataclass
class AlignedCohort:
    sample_ids: list[str]
    y: np.ndarray
    sample_index: dict[str, int]
    audit: AlignmentAudit


@dataclass(frozen=True)
class AlignedVariant:
    chrom: str
    pos: int
    variant_id: str
    ref: str
    alt: str
    a1: str
    g: np.ndarray


@dataclass
class AlignedInputs:
    cohort: AlignedCohort
    _variant_iter_factory: Callable[[], Iterator[AlignedVariant]] = field(repr=False)

    def iter_variants(self) -> Iterator[AlignedVariant]:
        return self._variant_iter_factory()


@dataclass(frozen=True)
class QcDecision:
    passed: bool
    reasons: list[str]
    maf_value: float | None
    non_missing_genotype_count: int


@dataclass(frozen=True)
class FilteredVariant:
    variant: AlignedVariant
    qc: QcDecision


@dataclass
class QcSummary:
    total_variants: int = 0
    passed_variants: int = 0
    failed_variants: int = 0


@dataclass(frozen=True)
class AssocResultRow:
    chrom: str
    snp: str
    bp: int
    a1: str
    test: str
    nmiss: int
    beta: float
    stat: float
    p: float


@dataclass(frozen=True)
class FitStats:
    term_names: list[str]
    beta: np.ndarray
    se: np.ndarray
    stat: np.ndarray
    p: np.ndarray
    nmiss: int
    df: int
    status: str


@dataclass
class AssocSummary:
    processed_variants: int = 0
    fit_failures: int = 0
    emitted_rows: int = 0


@dataclass(frozen=True)
class WriteSummary:
    output_path: str
    row_count: int
    success: bool
