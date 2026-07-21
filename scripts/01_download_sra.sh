#!/bin/bash
set -euo pipefail

mkdir -p data/raw logs

while read SRR; do
    echo "Downloading $SRR"

    prefetch "$SRR" -O data/raw

    echo "Converting $SRR to FASTQ"

    fasterq-dump "data/raw/${SRR}/${SRR}.sra" \
        --split-files \
        -O data/raw \
        -e 4 \
        -p

    gzip -f data/raw/${SRR}*.fastq

done < data/metadata/SRR_accessions.txt

echo "Download and FASTQ conversion completed."

