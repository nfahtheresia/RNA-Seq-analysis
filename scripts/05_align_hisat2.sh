#!/bin/bash
set -euo pipefail

ACCESSIONS="data/metadata/SRR_accessions.txt"
INDEX="data/reference/hisat2_index/Mus_musculus_GRCm39"

TRIMMED_DIR="data/trimmed"
RAW_DIR="data/raw"

BAM_DIR="results/alignment/bam"
LOG_DIR="results/alignment/logs"

mkdir -p "$BAM_DIR" "$LOG_DIR"

while read SRR; do
    echo "Aligning ${SRR}..."

    # Prefer trimmed reads if available; otherwise use raw reads
    if [[ -f "${TRIMMED_DIR}/${SRR}_1.trimmed.fastq.gz" && -f "${TRIMMED_DIR}/${SRR}_2.trimmed.fastq.gz" ]]; then
        R1="${TRIMMED_DIR}/${SRR}_1.trimmed.fastq.gz"
        R2="${TRIMMED_DIR}/${SRR}_2.trimmed.fastq.gz"
        echo "Using trimmed reads for ${SRR}"
    else
        R1="${RAW_DIR}/${SRR}_1.fastq.gz"
        R2="${RAW_DIR}/${SRR}_2.fastq.gz"
        echo "Using raw reads for ${SRR}"
    fi

    hisat2 \
        -p 4 \
        --dta \
        -x "$INDEX" \
        -1 "$R1" \
        -2 "$R2" \
        --summary-file "${LOG_DIR}/${SRR}_hisat2_summary.txt" \
    | samtools view -@ 4 -bS - \
    | samtools sort -@ 4 -o "${BAM_DIR}/${SRR}.sorted.bam"

    samtools index "${BAM_DIR}/${SRR}.sorted.bam"

    echo "${SRR} alignment completed."

done < "$ACCESSIONS"

echo "All alignments completed."
