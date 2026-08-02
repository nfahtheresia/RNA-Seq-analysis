#!/bin/bash
set -euo pipefail

GTF="data/reference/genome/Mus_musculus.GRCm39.116.gtf"
BAM_DIR="results/alignment/bam"
OUT_DIR="results/counts"

mkdir -p "$OUT_DIR" "$OUT_DIR/logs"

echo "Running featureCounts..."

featureCounts \
  -T 4 \
  -p \
  --countReadPairs \
  -t exon \
  -g gene_id \
  -a "$GTF" \
  -o "$OUT_DIR/gene_counts.txt" \
  "$BAM_DIR"/*.sorted.bam \
  2>&1 | tee "$OUT_DIR/logs/featurecounts.log"

echo "featureCounts completed."
