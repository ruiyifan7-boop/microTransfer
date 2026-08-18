#!/usr/bin/env bash
set -euo pipefail

: "${RICE_DATA_DIR:?Set RICE_DATA_DIR to the downloaded Dryad files}"
: "${RICE_OUT_DIR:?Set RICE_OUT_DIR to the results directory}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

for required in \
  lc_study_otu_table.tsv.gz \
  lc_study_mapping_file.tsv \
  gg_otus_tax.rds \
  organelle.rds
do
  if [[ ! -f "${RICE_DATA_DIR}/${required}" ]]; then
    echo "ERROR: missing ${RICE_DATA_DIR}/${required}" >&2
    exit 2
  fi
done

mkdir -p "${RICE_OUT_DIR}"
R CMD INSTALL "${PACKAGE_DIR}"

Rscript "${SCRIPT_DIR}/run_edwards_rice_generalization.R" \
  "${RICE_DATA_DIR}/lc_study_otu_table.tsv.gz" \
  "${RICE_DATA_DIR}/lc_study_mapping_file.tsv" \
  "${RICE_OUT_DIR}" \
  "${RICE_DATA_DIR}/gg_otus_tax.rds" \
  "${RICE_DATA_DIR}/organelle.rds" \
  2>&1 | tee "${RICE_OUT_DIR}/rice_generalization.log"
