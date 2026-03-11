from __future__ import annotations

from pathlib import Path

from plynk_lin.config import AssocResultRow
from plynk_lin.output_writer import write_assoc_linear


def test_write_assoc_linear_header_only(tmp_path: Path) -> None:
    summary = write_assoc_linear([], str(tmp_path / "nested" / "out"))
    output_path = tmp_path / "nested" / "out.assoc.linear"
    text = output_path.read_text(encoding="utf-8").splitlines()

    assert summary.row_count == 0
    assert text == [" CHR                                             SNP         BP   A1       TEST    NMISS       BETA         STAT            P "]


def test_write_assoc_linear_writes_rows(tmp_path: Path) -> None:
    rows = [
        AssocResultRow("1", "rs1", 100, "G", "ADD", 5, 0.5, 2.0, 0.04),
    ]
    summary = write_assoc_linear(rows, str(tmp_path / "out"))
    text = (tmp_path / "out.assoc.linear").read_text(encoding="utf-8").strip().splitlines()

    assert summary.row_count == 1
    assert text[0].split() == ["CHR", "SNP", "BP", "A1", "TEST", "NMISS", "BETA", "STAT", "P"]
    assert text[1] == "   1                                             rs1        100    G        ADD        5        0.5            2         0.04"
    assert text[1].split()[:5] == ["1", "rs1", "100", "G", "ADD"]
