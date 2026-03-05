from __future__ import annotations

import subprocess
import sys


def test_cli_parse_success_with_debug() -> None:
    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "plynk_lin",
            "--linear",
            "--vcf",
            "tests/data/test.vcf",
            "--pheno",
            "tests/data/pheno.txt",
            "--covar",
            "tests/data/covar.txt",
            "--covar-name",
            "AGE,PC1",
            "--out",
            "tmp/out",
            "--debug",
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0
    assert "parse summary" in proc.stdout
    assert "VCF samples" in proc.stdout
    assert "Variant preview" in proc.stdout


def test_cli_parse_failure_bad_vcf() -> None:
    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "plynk_lin",
            "--linear",
            "--vcf",
            "tests/data/test_bad_gt.vcf",
            "--pheno",
            "tests/data/pheno.txt",
            "--out",
            "tmp/out",
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    assert proc.returncode != 0
    assert "Parse error" in proc.stderr

