from __future__ import annotations

from pathlib import Path
from typing import Iterable

from plynk_lin.config import AssocResultRow, WriteSummary


HEADER = ["CHR", "SNP", "BP", "A1", "TEST", "NMISS", "BETA", "STAT", "P"]


def _format_float(value: float) -> str:
    return f"{value:.4g}"


def _format_row(row: AssocResultRow) -> str:
    fields = [
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
    return "\t".join(fields)


def write_assoc_linear(results: Iterable[AssocResultRow], out_prefix: str) -> WriteSummary:
    output_path = Path(f"{out_prefix}.assoc.linear")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    row_count = 0
    with output_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(HEADER))
        handle.write("\n")
        for row in results:
            handle.write(_format_row(row))
            handle.write("\n")
            row_count += 1

    return WriteSummary(output_path=str(output_path), row_count=row_count, success=True)
