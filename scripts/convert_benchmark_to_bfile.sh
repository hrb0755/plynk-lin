#!/usr/bin/env bash
set -euo pipefail

BENCHMARK_DIR="benchmark"
PLINK_CMD="plink"
OVERWRITE=0

run_plink() {
  # Deliberately allow word splitting so launchers like "conda run -n env plink" work.
  ${PLINK_CMD} "$@"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --benchmark-dir)
      BENCHMARK_DIR="$2"
      shift 2
      ;;
    --plink-cmd)
      PLINK_CMD="$2"
      shift 2
      ;;
    --overwrite)
      OVERWRITE=1
      shift
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 2
      ;;
  esac
done

MANIFEST_PATH="${BENCHMARK_DIR}/manifest.json"
if [[ ! -f "$MANIFEST_PATH" ]]; then
  echo "Benchmark manifest not found: ${MANIFEST_PATH}" >&2
  exit 1
fi

while IFS=$'\t' read -r family dataset_name output_dir; do
  [[ -z "$family" ]] && continue
  vcf_path="${output_dir}/subset.vcf.gz"
  out_prefix="${output_dir}/subset_bfile"

  if [[ ! -f "$vcf_path" ]]; then
    echo "Skipping ${family}/${dataset_name}: missing VCF ${vcf_path}" >&2
    continue
  fi

  if [[ "$OVERWRITE" -eq 0 && -f "${out_prefix}.bed" && -f "${out_prefix}.bim" && -f "${out_prefix}.fam" ]]; then
    echo "Skipping ${family}/${dataset_name}: bfile already exists at ${out_prefix}" >&2
    continue
  fi

  echo "Converting ${family}/${dataset_name} -> ${out_prefix}"
  run_plink --vcf "$vcf_path" --make-bed --out "$out_prefix"
done < <(
  python3 -c '
import json, sys
with open(sys.argv[1], "r", encoding="utf-8") as handle:
    manifest = json.load(handle)
for entry in manifest["datasets"]:
    print("\t".join([entry["family"], entry["dataset_name"], entry["output_dir"]]))
' "$MANIFEST_PATH"
)
