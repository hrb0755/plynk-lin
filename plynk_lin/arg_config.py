from __future__ import annotations

import argparse
from typing import Sequence

from plynk_lin.config import ConfigError, RunConfig


def _split_covar_names(raw: str | None) -> list[str] | None:
    if raw is None:
        return None
    names = [part.strip() for part in raw.split(",") if part.strip()]
    if not names:
        raise ConfigError("--covar-name must include at least one non-empty name")
    return names


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="plynk_lin")
    parser.add_argument("--linear", action="store_true", dest="linear_enabled")
    parser.add_argument("--vcf", dest="vcf_path")
    parser.add_argument("--pheno", dest="pheno_path")
    parser.add_argument("--covar", dest="covar_path")
    parser.add_argument("--covar-name", dest="covar_names")
    parser.add_argument("--maf", dest="maf_threshold", type=float)
    parser.add_argument("--allow-no-sex", action="store_true", dest="allow_no_sex")
    parser.add_argument("--out", dest="out_prefix")
    parser.add_argument("--debug", action="store_true", dest="debug")
    parser.add_argument("--debug-preview", type=int, default=5, dest="debug_preview")
    parser.add_argument("--debug-variants", type=int, default=5, dest="debug_variants")
    parser.add_argument("linear_modifiers", nargs="*")
    return parser


def parse_args(argv: Sequence[str]) -> RunConfig:
    parser = _build_parser()
    try:
        namespace = parser.parse_args(list(argv))
    except SystemExit as exc:
        if exc.code == 0:
            raise
        raise ConfigError("Invalid arguments") from exc

    modifiers = namespace.linear_modifiers or []
    allowed_modifiers = {"hide-covar"}
    unknown = [m for m in modifiers if m not in allowed_modifiers]
    if unknown:
        raise ConfigError(f"Unknown modifier(s): {', '.join(unknown)}")

    cfg = RunConfig(
        linear_enabled=bool(namespace.linear_enabled),
        hide_covar="hide-covar" in modifiers,
        vcf_path=namespace.vcf_path or "",
        pheno_path=namespace.pheno_path or "",
        covar_path=namespace.covar_path,
        covar_names=_split_covar_names(namespace.covar_names),
        maf_threshold=namespace.maf_threshold,
        allow_no_sex=bool(namespace.allow_no_sex),
        out_prefix=namespace.out_prefix or "",
        debug=bool(namespace.debug),
        debug_preview=namespace.debug_preview,
        debug_variants=namespace.debug_variants,
    )

    validate_config(cfg)
    return cfg


def validate_config(cfg: RunConfig) -> None:
    if not cfg.linear_enabled:
        raise ConfigError("--linear is required")
    if not cfg.vcf_path:
        raise ConfigError("--vcf is required")
    if not cfg.pheno_path:
        raise ConfigError("--pheno is required")
    if not cfg.out_prefix:
        raise ConfigError("--out is required and cannot be empty")

    if cfg.hide_covar and not cfg.linear_enabled:
        raise ConfigError("hide-covar is only allowed with --linear")

    if cfg.covar_names is not None and not cfg.covar_path:
        raise ConfigError("--covar-name requires --covar")

    if cfg.maf_threshold is not None and not (0.0 <= cfg.maf_threshold <= 0.5):
        raise ConfigError("--maf must satisfy 0.0 <= maf <= 0.5")

    if cfg.debug_preview < 0:
        raise ConfigError("--debug-preview must be >= 0")
    if cfg.debug_variants < 0:
        raise ConfigError("--debug-variants must be >= 0")
