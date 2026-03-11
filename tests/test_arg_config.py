from __future__ import annotations

import pytest

from plynk_lin.arg_config import parse_args
from plynk_lin.config import ConfigError


def test_parse_happy_path() -> None:
    cfg = parse_args(
        [
            "--linear",
            "--vcf",
            "x.vcf",
            "--pheno",
            "p.txt",
            "--maf",
            "0.1",
            "--allow-no-sex",
            "--out",
            "res/out",
            "--debug",
        ]
    )
    assert cfg.linear_enabled is True
    assert cfg.maf_threshold == 0.1
    assert cfg.debug is True


@pytest.mark.parametrize(
    "argv,expected",
    [
        (["--vcf", "x", "--pheno", "p", "--out", "o"], "--linear is required"),
        (["--linear", "--pheno", "p", "--out", "o"], "--vcf is required"),
        (["--linear", "--vcf", "x", "--out", "o"], "--pheno is required"),
        (["--linear", "--vcf", "x", "--pheno", "p"], "--out is required"),
        (["--linear", "--vcf", "x", "--pheno", "p", "--out", "o", "--maf", "0.8"], "--maf"),
    ],
)
def test_parse_validation_errors(argv: list[str], expected: str) -> None:
    with pytest.raises(ConfigError, match=expected):
        parse_args(argv)


def test_unknown_modifier_fails() -> None:
    with pytest.raises(ConfigError, match="Unknown modifier"):
        parse_args(["--linear", "--vcf", "x", "--pheno", "p", "--out", "o", "bad-mod"])
