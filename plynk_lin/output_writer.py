from __future__ import annotations

from pathlib import Path
from typing import Iterable

from plynk_lin.config import AssocResultRow, WriteSummary


HEADER = ["CHR", "SNP", "BP", "A1", "TEST", "NMISS", "BETA", "STAT", "P"]
FIELD_WIDTHS = (4, 47, 10, 4, 10, 8, 10, 12, 12)


def _format_float(value: float) -> str:
    return f"{value:.4g}"


def _format_fields(fields: list[str], *, trailing_space: bool = False) -> str:
    text = " ".join(
        f"{field:>{width}}" for field, width in zip(fields, FIELD_WIDTHS, strict=True)
    )
    return f"{text} " if trailing_space else text


def _format_row(row: AssocResultRow) -> str:
    return _format_fields(
        [
            row.chrom,
            row.snp,
            str(row.bp),
            row.a1,
            row.test,
            str(row.nmiss),
            _format_float(row.beta),
            _format_float(row.stat),
            _format_float(row.p),
        ]
    )


def write_assoc_linear(results: Iterable[AssocResultRow], out_prefix: str) -> WriteSummary:
    output_path = Path(f"{out_prefix}.assoc.linear")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    row_count = 0
    with output_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(_format_fields(HEADER, trailing_space=True))
        handle.write("\n")
        for row in results:
            handle.write(_format_row(row))
            handle.write("\n")
            row_count += 1

    return WriteSummary(output_path=str(output_path), row_count=row_count, success=True)
