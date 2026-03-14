# plynk-lin: A Minimal, Extensible Drop-in Subset of `plink --linear` GWAS

`plynk-lin` is a Python command-line tool for a narrow but useful subset of the
PLINK 1.9 linear-regression workflow. The goal is to reproduce the core behavior
of `plink --linear` for the PS3/CSE 284 evaluation setup while keeping the code
modular and easy to extend.

Supported CLI surface:
- `--linear`
- `--vcf`
- `--pheno`
- `--maf`
- `--allow-no-sex`
- `--out`

Current implementation includes:
- argument/config parsing
- VCF and phenotype parsing
- sample alignment with listwise deletion
- variant QC with optional MAF filtering
- per-variant additive linear regression
- PLINK-like `.assoc.linear` output writing
- smoke scripts and benchmark helpers

The presentation slides for the project can be found on [Google Docs](https://docs.google.com/presentation/d/1pEZaKuNE84LdVk36zvZgyPYlGT96Uvoj2YKKw-D7Tgc/edit?usp=sharing)

## Evaluation Goals

- correctness against the PS3 PLINK reference output; mostly on row number, index and numerical agreement
- deterministic behavior on sample/variant subset families
- runtime and peak RAM measurement on generated benchmark inputs

## Dataset Notes

- The main reference dataset is the PS3 subset derived from 1000 Genomes.
- The PS3 phenotype file is simulated.
- Performance evaluation is based on deterministic subsets of the PS3 VCF/pheno pair.
- If larger 1000 Genomes subsets are added later, the same benchmark scheme can be reused.

## Project Structure

Core modules:
- [Arg/config parsing](design_docs/00_arg_config_parsing.md)
- [IO](design_docs/01_io.md)
- [Sample/variant alignment](design_docs/02_sample_variant_alignment.md)
- [QC/filters](design_docs/03_qc_filters.md)
- [Association testing](design_docs/04_association_testing.md)
- [Output writing](design_docs/05_output_writing.md)

## Dependencies

Python 3.12\
`numpy`, `scipy` and `cyvcf2`

Additional requirement for testing:
- `pytest`


## Installation

```bash
pip install -e . 
```
add `--no-build-isolation` to the command if necessary.

## Drop-in PLINK Usage

This project is intended to mirror the PS3-style PLINK command:

```bash
plink --vcf path_to_ref/ps3_gwas.vcf.gz \
  --pheno path_to_ref/ps3_gwas.phen \
  --linear \
  --maf 0.05 \
  --allow-no-sex \
  --out ps3_gwas
```

The equivalent `plynk-lin` command is just replacing `plink` with `python -m plynk_lin`:

```bash
python -m plynk_lin \
  --linear \
  --vcf path_to_ref/ps3_gwas.vcf.gz \
  --pheno path_to_ref/ps3_gwas.phen \
  --maf 0.05 \
  --allow-no-sex \
  --out ps3_gwas
```

Example using the small test data fixtures (in `tests/data`) in this repo:

```bash
python -m plynk_lin \
  --linear \
  --vcf tests/data/pipeline.vcf \
  --pheno tests/data/pheno_complete.txt \
  --out tmp/example_run
```

Output:
- `tmp/example_run.assoc.linear`

Optional debugging:
- add `--debug` for a pipeline summary


## Testing

### Smoke Checks

Quick parsing check using included small test data:

```bash
python scripts/smoke_io.py \
  --vcf tests/data/pipeline.vcf \
  --pheno tests/data/pheno_complete.txt
```

Quick end-to-end pipeline check:

```bash
python scripts/smoke_pipeline.py \
  --vcf tests/data/pipeline.vcf \
  --pheno tests/data/pheno_complete.txt \
  --out tmp/smoke_pipeline \
  --debug
```

### Testing With Pytest

Run the full local test suite:

```bash
python -m pytest -q
```

Run an individual module test file:

```bash
python -m pytest -q tests/some_individual_test.py
```

### PS3 Reference Regression Test

For the full PS3 comparison against a PLINK reference output:

1. Put `ps3_gwas.vcf.gz` in `ref_data/` (directory is hard-coded in the testing script)
2. Put `ps3_gwas.phen` in `ref_data/`
3. Put the PLINK reference output at `ref_data/ref_out/ps3_gwas.assoc.linear`
4. Run:

```bash
PLYNK_RUN_REFERENCE=1 python -m pytest -q tests/test_end_to_end.py -k reference_dataset_regression
```

----

## Benchmarking 

### Benchmark Subset Generation

The benchmark workflow uses deterministic nested subsets generated from the PS3
VCF/pheno pair.

Default evaluation families:
- sample scaling: `N = 25, 50, 100, 207` with full `M = 917,845`
- variant scaling: `M = 1,000, 10,000, 100,000, 917,845` with full `N = 207`

Generate benchmark subsets:

```bash
python scripts/make_perf_subsets.py \
  --vcf ref_data/ps3_gwas.vcf.gz \
  --pheno ref_data/ps3_gwas.phen \
  --out benchmark
```

This writes:
- nested subset directories under `benchmark/sample_scaling/` and `benchmark/variant_scaling/`
- a manifest at `benchmark/manifest.json`

### Running Benchmarks

Measure `plynk-lin` on the generated subsets:

```bash
bash scripts/run_benchmarks.sh --benchmark-dir benchmark --plynk-lin
```

Measure PLINK on VCF input:

```bash
bash scripts/run_benchmarks.sh --benchmark-dir benchmark --plink
```

Optional: convert benchmark subsets to PLINK bed/bim/fam first:

```bash
bash scripts/convert_benchmark_to_bfile.sh --benchmark-dir benchmark
```

Then benchmark PLINK on binary input:

```bash
bash scripts/run_benchmarks.sh --benchmark-dir benchmark --plink-bfile
```

Benchmark outputs:
- `benchmark/results.csv`
- `benchmark/runs/(run-type)/stdout.log`
- `benchmark/runs/(run-type)/stderr.log`

### Benchmark Results

Benchmark results and their visualization used for presentation and report could be found in `benchmark` folder of the repo