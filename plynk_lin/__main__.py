from __future__ import annotations

import sys

from plynk_lin.alignment import align_samples
from plynk_lin.arg_config import parse_args
from plynk_lin.association import run_linear_assoc
from plynk_lin.config import AssocSummary, ConfigError, InputParseError, QcSummary
from plynk_lin.io import load_inputs
from plynk_lin.output_writer import write_assoc_linear
from plynk_lin.qc_filters import filter_variants
from plynk_lin.reporting import build_debug_summary


def main(argv: list[str] | None = None) -> int:
    args = sys.argv[1:] if argv is None else argv
    try:
        cfg = parse_args(args)
        parsed = load_inputs(cfg)
        aligned = align_samples(parsed)
        qc_summary = QcSummary()
        assoc_summary = AssocSummary()
        filtered = filter_variants(aligned.iter_variants(), cfg, summary=qc_summary)
        results = run_linear_assoc(aligned.cohort, filtered, cfg, summary=assoc_summary)
        write_summary = write_assoc_linear(results, cfg.out_prefix)
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
                cfg,
                parsed,
                aligned.cohort,
                qc_summary,
                assoc_summary,
                write_summary,
            )
        )
    else:
        print(
            f"Wrote {write_summary.row_count} rows to {write_summary.output_path}; "
            f"aligned_samples={aligned.cohort.audit.retained}, "
            f"qc_passed_variants={qc_summary.passed_variants}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
