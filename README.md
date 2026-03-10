# plynk-lin: A Minimal, Extensible Drop-in Subset of `plink --linear` GWAS
This is solo project for CSE 284. This repo is as of now still work in progress.

---

## Method to Implement
The goal is to implement a minimal but expandable Python command-line tool that reproduces the core functionality of **PLINK 1.9's linear-regression GWAS** (`plink --linear`). 

The project aims to create a drop-in compatible CLI that matches PLINK behavior for a restricted subset of flags:
* `--linear` (with `hide-covar` modifier)
* `--vcf` (matching PS3 workflow, diploid autosomal GT with multiallelic records reduced to the tested ALT allele for PLINK compatibility)
* `--pheno`, `--covar`, `--covar-name`, `--maf`, `--allow-no-sex`
* `--out` (output prefix for PLINK-like `.assoc.linear` table)

### Technical Focus
1.  **Correctness:** Accurately reproduce the OLS regression for GWAS.
2.  **Extensibility:** Implement the tool using separate modules for:
    * [Arg/config parsing](design_docs/00_arg_config_parsing.md)
    * [IO](design_docs/01_io.md)
    * [Sample/variant alignment](design_docs/02_sample_variant_alignment.md)
    * [QC/filters](design_docs/03_qc_filters.md)
    * [Association testing](design_docs/04_association_testing.md)
    * [Output writing](design_docs/05_output_writing.md)

---

## Install and Test
Currently implemented:
* Arg/config parsing
* IO parsing (VCF/pheno/covar)
* Sample alignment and cohort listwise deletion
* Variant QC with optional `--maf`
* Per-variant linear association testing
* PLINK-like `.assoc.linear` output writing
* End-to-end CLI and smoke scripts

### Dependencies
`python 3.12`\
`numpy, scipy, cyvcf2`\
and also `pytest` for testing

### Install
```bash
pip install -e . --no-build-isolation
```


### Run CLI
There are sample VCF/pheno/covar files in `tests/data`.
```bash
python -m plynk_lin --linear --vcf /path/to/input.vcf --pheno /path/to/pheno.txt --out /tmp/out --debug
```

Optional covariates:
```bash
python -m plynk_lin --linear --vcf /path/to/input.vcf --pheno /path/to/pheno.txt --covar /path/to/covar.txt --covar-name AGE,PC1 --out /tmp/out --debug
```

### Smoke check scripts
```bash
python scripts/smoke_io.py --vcf /path/to/input.vcf --pheno /path/to/pheno.txt --covar /path/to/covar.txt --covar-name AGE,PC1

python scripts/smoke_pipeline.py --vcf /path/to/input.vcf --pheno /path/to/pheno.txt --out /tmp/out --debug
```

### Tests
Run all tests:
```bash
python -m pytest -q
```

Run a specific test module:
```bash
python -m pytest -q tests/test_io_vcf.py
```

---

## Evaluation Strategy
1.  **Correctness and Behavior Alignment:** 
    * Compare results against PLINK using the same PS3-like dataset.
    * Validate per-SNP Beta/t-stat/P agreement (correlation and max absolute differences).
    * Compare top hits (top-K overlap and genome-wide significant hits).
2.  **Data Integrity & Edge-case Alignment:**
    * Verify ID-based alignment by shuffling phenotype/covariate file orders.
    * Confirm listwise deletion and matching NMISS/df conventions by injecting missing values.
3.  **Performance Benchmark:**
    * Measure runtime and peak memory usage on varying datasets.

---

## Dataset
* Primarily the subset from the **1000 Genome dataset** used in PS3.
* Additional subsets to be derived from the complete 1000 Genome dataset available on DataHub.
* Phenotype data to be simulated as 1000 Genome dataset does not contain them.
