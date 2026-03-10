from __future__ import annotations

import os
from pathlib import Path

import pytest

from plynk_lin.__main__ import main


def _parse_assoc(path: Path) -> tuple[list[str], list[list[str]]]:
    lines = path.read_text(encoding="utf-8").strip().splitlines()
    header = lines[0].split()
    rows = [line.split() for line in lines[1:]]
    return header, rows


def test_synthetic_end_to_end_pipeline(tmp_path: Path) -> None:
    out_prefix = tmp_path / "synthetic"
    exit_code = main(
        [
            "--linear",
            "--vcf",
            "tests/data/pipeline.vcf",
            "--pheno",
            "tests/data/pheno_complete.txt",
            "--out",
            str(out_prefix),
        ]
    )

    assert exit_code == 0
    header, rows = _parse_assoc(tmp_path / "synthetic.assoc.linear")
    assert header == ["CHR", "SNP", "BP", "A1", "TEST", "NMISS", "BETA", "STAT", "P"]
    assert len(rows) == 3
    assert [row[1] for row in rows] == ["rs1", "rs2", "rs3"]
    assert [row[4] for row in rows] == ["ADD", "ADD", "ADD"]


@pytest.mark.skipif(
    os.environ.get("PLYNK_RUN_REFERENCE") != "1",
    reason="reference regression is expensive; set PLYNK_RUN_REFERENCE=1 to enable",
)
def test_reference_dataset_regression(tmp_path: Path) -> None:
    out_prefix = tmp_path / "ps3_gwas"
    exit_code = main(
        [
            "--linear",
            "--vcf",
            "ref_data/ps3_gwas.vcf.gz",
            "--pheno",
            "ref_data/ps3_gwas.phen",
            "--maf",
            "0.05",
            "--allow-no-sex",
            "--out",
            str(out_prefix),
        ]
    )
    assert exit_code == 0

    actual_header, actual_rows = _parse_assoc(tmp_path / "ps3_gwas.assoc.linear")
    expected_header, expected_rows = _parse_assoc(Path("ref_data/ref_out/ps3_gwas.assoc.linear"))

    assert actual_header == expected_header
    assert len(actual_rows) == len(expected_rows)

    for actual, expected in zip(actual_rows[:1000], expected_rows[:1000], strict=True):
        assert actual[:6] == expected[:6]
        assert pytest.approx(float(expected[6]), abs=1e-4) == float(actual[6])
        assert pytest.approx(float(expected[7]), abs=1e-4) == float(actual[7])
        assert pytest.approx(float(expected[8]), abs=1e-6) == float(actual[8])
