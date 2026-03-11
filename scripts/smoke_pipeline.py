#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from plynk_lin.__main__ import main


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Smoke-test the full plynk-lin pipeline")
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--pheno", required=True)
    parser.add_argument("--maf")
    parser.add_argument("--out", required=True)
    parser.add_argument("--debug", action="store_true", dest="debug")
    return parser.parse_args()


def main_cli() -> int:
    ns = parse_args()
    argv = [
        "--linear",
        "--vcf",
        ns.vcf,
        "--pheno",
        ns.pheno,
        "--out",
        ns.out,
    ]
    if ns.maf:
        argv.extend(["--maf", ns.maf])
    if ns.debug:
        argv.append("--debug")
    return main(argv)


if __name__ == "__main__":
    raise SystemExit(main_cli())
