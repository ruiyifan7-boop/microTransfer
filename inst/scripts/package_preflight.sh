#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
BUILD_DIR="$(cd "${PACKAGE_DIR}/.." && pwd)"

cd "${BUILD_DIR}"

Rscript -e "files <- list.files('${PACKAGE_DIR}/R', full.names=TRUE); invisible(lapply(files, parse)); cat('R_PARSE_OK\\n')"
R CMD build "${PACKAGE_DIR}"

PACKAGE_TARBALL="$(ls -1t microTransfer_*.tar.gz | head -1)"
R CMD check "${PACKAGE_TARBALL}" --no-manual

echo "PACKAGE_TARBALL=${BUILD_DIR}/${PACKAGE_TARBALL}"
echo "CHECK_DIR=${BUILD_DIR}/microTransfer.Rcheck"
