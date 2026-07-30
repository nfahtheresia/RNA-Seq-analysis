#!/bin/bash
set -euo pipefail

REF="data/reference/Mus_musculus.GRCm39.dna.primary_assembly.fa"
INDEX_DIR="data/reference/hisat2_index"
INDEX_PREFIX="${INDEX_DIR}/Mus_musculus_GRCm39"

mkdir -p "$INDEX_DIR" logs

echo "Building HISAT2 index for Mus musculus GRCm39..."
hisat2-build -p 4 "$REF" "$INDEX_PREFIX" 2>&1 | tee logs/hisat2_build_index.log

echo "HISAT2 index completed."
