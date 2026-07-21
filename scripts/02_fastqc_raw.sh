#!/bin/bash
set -euo pipefail

mkdir -p results/qc/fastqc_raw
mkdir -p results/qc/multiqc
mkdir -p logs

fastqc data/raw/*.fastq.gz \
    -o results/qc/fastqc_raw \
    -t 4 2>&1 | tee logs/fastqc_raw.log

multiqc results/qc/fastqc_raw \
    -o results/qc/multiqc \
    -n multiqc_raw_report.html

echo "Raw FastQC and MultiQC completed."

