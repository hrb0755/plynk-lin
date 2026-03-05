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
from plynk_lin.io import load_covar, load_pheno, load_vcf


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Smoke-test parser + IO stack")
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--pheno", required=True)
    parser.add_argument("--covar")
    parser.add_argument("--covar-name", dest="covar_name")
    return parser.parse_args()


def main() -> int:
    ns = parse_args()
    covar_names = None
    if ns.covar_name:
        covar_names = [p.strip() for p in ns.covar_name.split(",") if p.strip()]

    try:
        vcf = load_vcf(ns.vcf)
        pheno = load_pheno(ns.pheno)
        covar = load_covar(ns.covar, covar_names) if ns.covar else None
    except (InputParseError, FileNotFoundError) as exc:
        print(f"Smoke test failed: {exc}", file=sys.stderr)
        return 1

    print("smoke_io summary")
    print(f"vcf_samples={len(vcf.sample_ids)}")
    print(f"vcf_variants={len(vcf.variants)}")
    print(f"pheno_rows={len(pheno.sample_ids)}")
    print(f"covar_rows={len(covar.sample_ids) if covar else 0}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
