#!/usr/bin/env bash
set -euo pipefail

BENCHMARK_DIR="benchmark"
RUN_PLINK=0
RUN_PLINK_BFILE=0
RUN_PLYNK_LIN=0
PLINK_CMD="plink"
PLYNK_LIN_CMD=""
PYTHON_BIN="${PYTHON:-conda run -n plynk python}"
RESULTS_CSV=""
RUNS_DIR=""

run_python() {
  # Deliberately allow word splitting so launchers like "conda run -n plynk python" work.
  ${PYTHON_BIN} "$@"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --benchmark-dir)
      BENCHMARK_DIR="$2"
      shift 2
      ;;
    --plink)
      RUN_PLINK=1
      shift
      ;;
    --plink-bfile)
      RUN_PLINK_BFILE=1
      shift
      ;;
    --plynk-lin)
      RUN_PLYNK_LIN=1
      shift
      ;;
    --plink-cmd)
      PLINK_CMD="$2"
      shift 2
      ;;
    --plynk-lin-cmd)
      PLYNK_LIN_CMD="$2"
      shift 2
      ;;
    --python)
      PYTHON_BIN="$2"
      shift 2
      ;;
    --results-csv)
      RESULTS_CSV="$2"
      shift 2
      ;;
    --runs-dir)
      RUNS_DIR="$2"
      shift 2
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 2
      ;;
  esac
done

if [[ "$RUN_PLINK" -eq 0 && "$RUN_PLINK_BFILE" -eq 0 && "$RUN_PLYNK_LIN" -eq 0 ]]; then
  echo "Provide at least one of --plink, --plink-bfile, or --plynk-lin" >&2
  exit 2
fi

if [[ -z "$PLYNK_LIN_CMD" ]]; then
  PLYNK_LIN_CMD="${PYTHON_BIN} -m plynk_lin"
fi

MANIFEST_PATH="${BENCHMARK_DIR}/manifest.json"
if [[ ! -f "$MANIFEST_PATH" ]]; then
  echo "Benchmark manifest not found: ${MANIFEST_PATH}" >&2
  exit 1
fi

if [[ -z "$RESULTS_CSV" ]]; then
  RESULTS_CSV="${BENCHMARK_DIR}/results.csv"
fi
if [[ -z "$RUNS_DIR" ]]; then
  RUNS_DIR="${BENCHMARK_DIR}/runs"
fi

mkdir -p "$RUNS_DIR"

while IFS=$'\t' read -r family dataset_name output_dir sample_count variant_count; do
  [[ -z "$family" ]] && continue
  vcf_path="${output_dir}/subset.vcf.gz"
  pheno_path="${output_dir}/subset.phen"
  if [[ ! -f "$vcf_path" || ! -f "$pheno_path" ]]; then
    echo "Skipping ${family}/${dataset_name}: missing subset files" >&2
    continue
  fi

  if [[ "$RUN_PLINK" -eq 1 ]]; then
    run_dir="${RUNS_DIR}/plink/${family}/${dataset_name}"
    out_prefix="${run_dir}/result"
    printf -v shell_cmd '%s --vcf %q --pheno %q --linear --maf 0.05 --allow-no-sex --out %q' \
      "$PLINK_CMD" "$vcf_path" "$pheno_path" "$out_prefix"
    run_python scripts/measure_command.py \
      --tool plink \
      --dataset-family "$family" \
      --dataset-name "$dataset_name" \
      --sample-count "$sample_count" \
      --variant-count "$variant_count" \
      --csv "$RESULTS_CSV" \
      --run-dir "$run_dir" \
      --shell-command "$shell_cmd" || true
  fi

  if [[ "$RUN_PLINK_BFILE" -eq 1 ]]; then
    bfile_prefix="${output_dir}/subset_bfile"
    if [[ ! -f "${bfile_prefix}.bed" || ! -f "${bfile_prefix}.bim" || ! -f "${bfile_prefix}.fam" ]]; then
      echo "Skipping ${family}/${dataset_name}: missing bfile files for prefix ${bfile_prefix}" >&2
    else
      run_dir="${RUNS_DIR}/plink-bfile/${family}/${dataset_name}"
      out_prefix="${run_dir}/result"
      printf -v shell_cmd '%s --bfile %q --pheno %q --linear --maf 0.05 --allow-no-sex --out %q' \
        "$PLINK_CMD" "$bfile_prefix" "$pheno_path" "$out_prefix"
      run_python scripts/measure_command.py \
        --tool plink-bfile \
        --dataset-family "$family" \
        --dataset-name "$dataset_name" \
        --sample-count "$sample_count" \
        --variant-count "$variant_count" \
        --csv "$RESULTS_CSV" \
        --run-dir "$run_dir" \
        --shell-command "$shell_cmd" || true
    fi
  fi

  if [[ "$RUN_PLYNK_LIN" -eq 1 ]]; then
    run_dir="${RUNS_DIR}/plynk-lin/${family}/${dataset_name}"
    out_prefix="${run_dir}/result"
    printf -v shell_cmd '%s --linear --vcf %q --pheno %q --maf 0.05 --allow-no-sex --out %q' \
      "$PLYNK_LIN_CMD" "$vcf_path" "$pheno_path" "$out_prefix"
    run_python scripts/measure_command.py \
      --tool plynk-lin \
      --dataset-family "$family" \
      --dataset-name "$dataset_name" \
      --sample-count "$sample_count" \
      --variant-count "$variant_count" \
      --csv "$RESULTS_CSV" \
      --run-dir "$run_dir" \
      --shell-command "$shell_cmd" || true
  fi
done < <(
  run_python -c '
import json, sys
with open(sys.argv[1], "r", encoding="utf-8") as handle:
    manifest = json.load(handle)
for entry in manifest["datasets"]:
    print(
        "\t".join(
            [
                entry["family"],
                entry["dataset_name"],
                entry["output_dir"],
                str(entry["sample_count"]),
                str(entry["variant_count"]),
            ]
        )
    )
' "$MANIFEST_PATH"
)
