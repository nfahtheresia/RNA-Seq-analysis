#!/bin/bash
set -euo pipefail

mkdir -p data/trimmed
mkdir -p results/qc/fastp
mkdir -p results/qc/fastqc_trimmed
mkdir -p results/qc/multiqc
mkdir -p logs

while read SRR; do
    echo "Trimming $SRR with fastp"

    fastp \
        -i data/raw/${SRR}_1.fastq.gz \
        -I data/raw/${SRR}_2.fastq.gz \
        -o data/trimmed/${SRR}_1.trimmed.fastq.gz \
        -O data/trimmed/${SRR}_2.trimmed.fastq.gz \
        --detect_adapter_for_pe \
        --thread 4 \
        --html results/qc/fastp/${SRR}_fastp.html \
        --json results/qc/fastp/${SRR}_fastp.json

done < data/metadata/SRR_accessions.txt

fastqc data/trimmed/*.fastq.gz \
    -o results/qc/fastqc_trimmed \
    -t 4 2>&1 | tee logs/fastqc_trimmed.log

multiqc results/qc \
    -o results/qc/multiqc \
    -n multiqc_final_report.html

echo "fastp trimming, trimmed FastQC and final MultiQC completed."
