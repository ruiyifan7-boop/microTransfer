#!/usr/bin/env bash
set -euo pipefail

cd "${IMETA_ROOT_DIR:-$HOME/projects/blueberry_meta}"

Rscript --vanilla make_main_figures_v3.R 2>&1 | tee nature_main_figures.log
Rscript --vanilla make_Figure3_pH_transfer_v4.R 2>&1 | tee nature_figure3.log
Rscript --vanilla make_management_transfer_figure_v2.R 2>&1 | tee nature_figure4.log

echo
echo "Nature-style vector outputs:"
find global_harmonized/figures_nature \
     global_harmonized/enhancement/figures_nature \
     -maxdepth 1 -type f \
     \( -iname '*.pdf' -o -iname '*.png' -o -iname '*.tiff' \) \
     -printf '%TY-%Tm-%Td %TH:%TM:%TS  %s  %p\n' | sort
