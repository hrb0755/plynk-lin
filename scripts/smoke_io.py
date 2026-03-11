#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Ensure script execution works without package installation.
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from plynk_lin.config import InputParseError
from plynk_lin.io import load_pheno, load_vcf


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Smoke-test parser + IO stack")
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--pheno", required=True)
    return parser.parse_args()


def main() -> int:
    ns = parse_args()
    try:
        vcf = load_vcf(ns.vcf)
        pheno = load_pheno(ns.pheno)
    except (InputParseError, FileNotFoundError) as exc:
        print(f"Smoke test failed: {exc}", file=sys.stderr)
        return 1

    print("smoke_io summary")
    print(f"vcf_samples={len(vcf.sample_ids)}")
    print(f"vcf_variants={len(vcf.variants)}")
    print(f"pheno_rows={len(pheno.sample_ids)}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
