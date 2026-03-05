from __future__ import annotations

import sys

from plynk_lin.arg_config import parse_args
from plynk_lin.config import ConfigError, InputParseError
from plynk_lin.io import load_inputs
from plynk_lin.reporting import build_debug_summary


def main(argv: list[str] | None = None) -> int:
    args = sys.argv[1:] if argv is None else argv
    try:
        cfg = parse_args(args)
        parsed = load_inputs(cfg)
    except ConfigError as exc:
        print(f"Config error: {exc}", file=sys.stderr)
        return 2
    except FileNotFoundError as exc:
        print(f"File error: {exc}", file=sys.stderr)
        return 2
    except InputParseError as exc:
        print(f"Parse error: {exc}", file=sys.stderr)
        return 2

    if cfg.debug:
        print(
            build_debug_summary(
                parsed,
                sample_preview=cfg.debug_preview,
                variant_preview=cfg.debug_variants,
            )
        )
    else:
        print(
            f"Parsed VCF={len(parsed.vcf.variants)} variants, "
            f"samples={len(parsed.vcf.sample_ids)}, "
            f"pheno_rows={len(parsed.pheno.sample_ids)}, "
            f"covar_rows={len(parsed.covar.sample_ids) if parsed.covar else 0}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
