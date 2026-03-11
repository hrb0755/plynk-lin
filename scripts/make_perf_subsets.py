#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import json
import random
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, TextIO

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from plynk_lin.config import InputParseError
from plynk_lin.io import load_pheno


RNG_SEED = 284
SAMPLE_SCALING_COUNTS = (25, 50, 100)
VARIANT_SCALING_COUNTS = (1_000, 10_000, 100_000)


@dataclass(frozen=True)
class DatasetSpec:
    family: str
    dataset_name: str
    sample_count: int
    variant_count: int
    use_full_samples: bool
    use_full_variants: bool


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate deterministic benchmark subsets from PS3 VCF/pheno inputs")
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--pheno", required=True)
    parser.add_argument("--out", required=True)
    return parser.parse_args()


def _open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def _open_text_write(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "wt", encoding="utf-8", newline="\n")
    return path.open("w", encoding="utf-8", newline="\n")


def read_vcf_header(vcf_path: Path) -> tuple[list[str], list[str]]:
    header_lines: list[str] = []
    sample_ids: list[str] = []
    with _open_text(vcf_path) as handle:
        for raw in handle:
            if raw.startswith("##"):
                header_lines.append(raw)
                continue
            if raw.startswith("#CHROM"):
                header_lines.append(raw)
                sample_ids = raw.rstrip("\n").split("\t")[9:]
                break
    if not sample_ids:
        raise InputParseError(f"No #CHROM header with sample IDs found in {vcf_path}")
    return header_lines, sample_ids


def count_variants(vcf_path: Path) -> int:
    count = 0
    with _open_text(vcf_path) as handle:
        for raw in handle:
            if not raw.startswith("#"):
                count += 1
    return count


def build_nested_order(items: Iterable[str], seed: int) -> list[str]:
    ordered = list(items)
    random.Random(seed).shuffle(ordered)
    return ordered


def build_nested_variant_order(total_variants: int, seed: int) -> list[int]:
    ordered = list(range(total_variants))
    random.Random(seed).shuffle(ordered)
    return ordered


def build_dataset_specs(total_samples: int, total_variants: int) -> list[DatasetSpec]:
    sample_sizes = [count for count in SAMPLE_SCALING_COUNTS if count < total_samples]
    sample_sizes.append(total_samples)
    variant_sizes = [count for count in VARIANT_SCALING_COUNTS if count < total_variants]
    variant_sizes.append(total_variants)

    specs: list[DatasetSpec] = []
    for sample_count in sample_sizes:
        specs.append(
            DatasetSpec(
                family="sample_scaling",
                dataset_name=f"n{sample_count:03d}_mfull",
                sample_count=sample_count,
                variant_count=total_variants,
                use_full_samples=sample_count == total_samples,
                use_full_variants=True,
            )
        )
    for variant_count in variant_sizes:
        specs.append(
            DatasetSpec(
                family="variant_scaling",
                dataset_name=f"nfull_m{'full' if variant_count == total_variants else f'{variant_count:06d}'}",
                sample_count=total_samples,
                variant_count=variant_count,
                use_full_samples=True,
                use_full_variants=variant_count == total_variants,
            )
        )
    return specs


def validate_output_dir(out_dir: Path) -> None:
    if out_dir.exists():
        if not out_dir.is_dir():
            raise FileExistsError(f"Output path exists and is not a directory: {out_dir}")
        if any(out_dir.iterdir()):
            raise FileExistsError(f"Output directory must be empty: {out_dir}")


def load_matching_phenotypes(pheno_path: Path, sample_ids: list[str]) -> dict[str, float | None]:
    pheno = load_pheno(str(pheno_path))
    missing = [sid for sid in sample_ids if sid not in pheno.values_by_sample]
    if missing:
        preview = ", ".join(missing[:5])
        raise InputParseError(
            f"Phenotype file {pheno_path} is missing {len(missing)} VCF sample(s); first missing: {preview}"
        )
    return {sid: pheno.values_by_sample[sid] for sid in sample_ids}


def write_subset_pheno(path: Path, ordered_sample_ids: list[str], phenotype_map: dict[str, float | None]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for sid in ordered_sample_ids:
            value = phenotype_map[sid]
            token = "NA" if value is None else repr(float(value))
            handle.write(f"{sid}\t{sid}\t{token}\n")


def write_subset_vcf(
    source_vcf: Path,
    out_vcf: Path,
    retained_sample_ids: list[str],
    retained_variant_indices: set[int],
) -> None:
    retained_sample_set = set(retained_sample_ids)
    with _open_text(source_vcf) as src, _open_text_write(out_vcf) as dst:
        variant_idx = 0
        sample_col_idx: list[int] | None = None
        for raw in src:
            if raw.startswith("##"):
                dst.write(raw)
                continue
            if raw.startswith("#CHROM"):
                cols = raw.rstrip("\n").split("\t")
                base = cols[:9]
                samples = cols[9:]
                sample_col_idx = [idx for idx, sid in enumerate(samples) if sid in retained_sample_set]
                if len(sample_col_idx) != len(retained_sample_ids):
                    raise InputParseError("Failed to resolve all retained sample IDs against VCF header")
                reordered_samples = [samples[idx] for idx in sample_col_idx]
                if reordered_samples != retained_sample_ids:
                    index_map = {sid: idx for idx, sid in enumerate(samples)}
                    sample_col_idx = [index_map[sid] for sid in retained_sample_ids]
                    reordered_samples = [samples[idx] for idx in sample_col_idx]
                dst.write("\t".join(base + reordered_samples) + "\n")
                continue

            if variant_idx in retained_variant_indices:
                cols = raw.rstrip("\n").split("\t")
                if sample_col_idx is None:
                    raise InputParseError("Encountered variant row before #CHROM header")
                kept = [cols[9 + idx] for idx in sample_col_idx]
                dst.write("\t".join(cols[:9] + kept) + "\n")
            variant_idx += 1


def build_manifest(
    source_vcf: Path,
    source_pheno: Path,
    out_dir: Path,
    total_samples: int,
    total_variants: int,
    specs: list[DatasetSpec],
) -> dict[str, object]:
    datasets = [
        {
            **asdict(spec),
            "output_dir": str(out_dir / spec.family / spec.dataset_name),
        }
        for spec in specs
    ]
    return {
        "source_vcf": str(source_vcf),
        "source_pheno": str(source_pheno),
        "seed": RNG_SEED,
        "source_sample_count": total_samples,
        "source_variant_count": total_variants,
        "datasets": datasets,
    }


def generate_subsets(vcf_path: Path, pheno_path: Path, out_dir: Path) -> dict[str, object]:
    validate_output_dir(out_dir)
    header_lines, source_sample_ids = read_vcf_header(vcf_path)
    del header_lines
    total_variants = count_variants(vcf_path)
    phenotype_map = load_matching_phenotypes(pheno_path, source_sample_ids)
    specs = build_dataset_specs(len(source_sample_ids), total_variants)

    sample_order = build_nested_order(source_sample_ids, RNG_SEED)
    variant_order = build_nested_variant_order(total_variants, RNG_SEED + 1)

    out_dir.mkdir(parents=True, exist_ok=False)

    for spec in specs:
        dataset_dir = out_dir / spec.family / spec.dataset_name
        dataset_dir.mkdir(parents=True, exist_ok=False)

        retained_sample_ids = sample_order[: spec.sample_count]
        retained_variant_indices = set(variant_order[: spec.variant_count])

        write_subset_vcf(
            source_vcf=vcf_path,
            out_vcf=dataset_dir / "subset.vcf.gz",
            retained_sample_ids=retained_sample_ids,
            retained_variant_indices=retained_variant_indices,
        )
        write_subset_pheno(dataset_dir / "subset.phen", retained_sample_ids, phenotype_map)

    manifest = build_manifest(
        source_vcf=vcf_path,
        source_pheno=pheno_path,
        out_dir=out_dir,
        total_samples=len(source_sample_ids),
        total_variants=total_variants,
        specs=specs,
    )
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def main() -> int:
    ns = parse_args()
    vcf_path = Path(ns.vcf)
    pheno_path = Path(ns.pheno)
    out_dir = Path(ns.out)
    try:
        manifest = generate_subsets(vcf_path, pheno_path, out_dir)
    except (FileNotFoundError, FileExistsError, InputParseError) as exc:
        print(f"Subset generation failed: {exc}", file=sys.stderr)
        return 1

    print("generated benchmark subsets")
    print(f"source_samples={manifest['source_sample_count']}")
    print(f"source_variants={manifest['source_variant_count']}")
    for dataset in manifest["datasets"]:
        print(
            f"{dataset['family']} {dataset['dataset_name']} "
            f"samples={dataset['sample_count']} variants={dataset['variant_count']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
