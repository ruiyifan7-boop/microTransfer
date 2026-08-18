#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="$ROOT/data/rice"
mkdir -p "$OUT"

# Dryad's current download endpoint is /api/v2/files/<id>/download
# (the old downloads/file_stream/<id> path now returns HTTP 403).
base="https://datadryad.org/api/v2/files"

curl -fL --retry 3 --retry-delay 5 \
  "$base/17335/download" -o "$OUT/lc_study_otu_table.tsv.gz"
curl -fL --retry 3 --retry-delay 5 \
  "$base/17337/download" -o "$OUT/gg_otus_tax.rds"
curl -fL --retry 3 --retry-delay 5 \
  "$base/17338/download" -o "$OUT/organelle.rds"

cd "$OUT"
echo "85bde2c7a81140907d2a9cceae5f55c4  lc_study_otu_table.tsv.gz" | md5sum -c -
echo "cfb234a8be3f3f195b6b5b903d16782f  gg_otus_tax.rds" | md5sum -c -
echo "0965d49f14366ed053e37e4943c61e27  organelle.rds" | md5sum -c -

echo "Rice files downloaded and verified."

