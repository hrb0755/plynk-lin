from __future__ import annotations

import pytest

from plynk_lin.config import InputParseError
from plynk_lin.io import load_covar, load_pheno


def test_load_pheno_happy_path() -> None:
    pheno = load_pheno("tests/data/pheno.txt")
    assert pheno.sample_ids == ["S1", "S2", "S3"]
    assert pheno.phenotype_column == "PHENO"
    assert pheno.values_by_sample["S1"] == 1.2
    assert pheno.values_by_sample["S3"] is None


def test_load_pheno_duplicate_iid_fails() -> None:
    with pytest.raises(InputParseError, match="Duplicate IID"):
        load_pheno("tests/data/pheno_dup.txt")


def test_load_pheno_headerless_happy_path() -> None:
    pheno = load_pheno("tests/data/pheno_headerless.txt")
    assert pheno.sample_ids == ["S1", "S2", "S3"]
    assert pheno.phenotype_column == "PHENO"
    assert pheno.values_by_sample["S3"] is None


def test_load_covar_subset_happy_path() -> None:
    covar = load_covar("tests/data/covar.txt", ["AGE", "PC1"])
    assert covar.covar_columns == ["AGE", "PC1"]
    assert covar.values_by_sample["S1"]["AGE"] == 30.0
    assert covar.values_by_sample["S3"]["PC1"] is None


def test_load_covar_unknown_column_fails() -> None:
    with pytest.raises(InputParseError, match="Unknown covariate"):
        load_covar("tests/data/covar.txt", ["AGE", "MISSING_COL"])


def test_load_covar_duplicate_iid_fails() -> None:
    with pytest.raises(InputParseError, match="Duplicate IID"):
        load_covar("tests/data/covar_dup.txt", None)
