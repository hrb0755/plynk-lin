from __future__ import annotations

import pytest

from plynk_lin.config import InputParseError
from plynk_lin.io import load_vcf


def test_load_vcf_happy_path() -> None:
    ds = load_vcf("tests/data/test.vcf")
    assert ds.sample_ids == ["S1", "S2", "S3"]
    assert len(ds.variants) == 2
    assert ds.variants[0].genotypes_by_sample["S1"] == 0
    assert ds.variants[0].genotypes_by_sample["S2"] == 1
    assert ds.variants[0].genotypes_by_sample["S3"] == 2
    assert ds.variants[1].genotypes_by_sample["S2"] is None


def test_load_vcf_rejects_multiallelic() -> None:
    with pytest.raises(InputParseError, match="non-biallelic"):
        load_vcf("tests/data/test_bad_multiallelic.vcf")


def test_load_vcf_rejects_bad_gt() -> None:
    with pytest.raises(InputParseError, match="Only biallelic 0/1 GT supported"):
        load_vcf("tests/data/test_bad_gt.vcf")

