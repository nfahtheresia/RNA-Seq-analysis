#!/bin/bash
set -euo pipefail

ACCESSIONS="data/metadata/SRR_accessions.txt"
INDEX="data/reference/hisat2_index/Mus_musculus_GRCm39"

RAW_DIR="data/raw"
TRIMMED_DIR="data/trimmed"

BAM_DIR="results/alignment/bam"
LOG_DIR="results/alignment/logs"

mkdir -p "$BAM_DIR" "$LOG_DIR"

while read SRR; do
    echo "Aligning ${SRR}..."

    rm -f "${BAM_DIR}/${SRR}.sorted.bam" "${BAM_DIR}/${SRR}.sorted.bam.bai"

    if [[ -f "${TRIMMED_DIR}/${SRR}_1.trimmed.fastq.gz" && -f "${TRIMMED_DIR}/${SRR}_2.trimmed.fastq.gz" ]]; then
        echo "Using trimmed paired-end reads for ${SRR}"

        hisat2 -p 1 --dta \
            -x "$INDEX" \
            -1 "${TRIMMED_DIR}/${SRR}_1.trimmed.fastq.gz" \
            -2 "${TRIMMED_DIR}/${SRR}_2.trimmed.fastq.gz" \
            --summary-file "${LOG_DIR}/${SRR}_hisat2_summary.txt" \
        | samtools view -@ 1 -bS - \
        | samtools sort -@ 1 -m 512M -o "${BAM_DIR}/${SRR}.sorted.bam"

    elif [[ -f "${RAW_DIR}/${SRR}_1.fastq.gz" && -f "${RAW_DIR}/${SRR}_2.fastq.gz" ]]; then
        echo "Using raw paired-end reads for ${SRR}"

        hisat2 -p 1 --dta \
            -x "$INDEX" \
            -1 "${RAW_DIR}/${SRR}_1.fastq.gz" \
            -2 "${RAW_DIR}/${SRR}_2.fastq.gz" \
            --summary-file "${LOG_DIR}/${SRR}_hisat2_summary.txt" \
        | samtools view -@ 1 -bS - \
        | samtools sort -@ 1 -m 512M -o "${BAM_DIR}/${SRR}.sorted.bam"

    elif [[ -f "${RAW_DIR}/${SRR}.fastq.gz" ]]; then
        echo "Using raw single-end reads for ${SRR}"

        hisat2 -p 1 --dta \
            -x "$INDEX" \
            -U "${RAW_DIR}/${SRR}.fastq.gz" \
            --summary-file "${LOG_DIR}/${SRR}_hisat2_summary.txt" \
        | samtools view -@ 1 -bS - \
        | samtools sort -@ 1 -m 512M -o "${BAM_DIR}/${SRR}.sorted.bam"

    else
        echo "ERROR: FASTQ files not found for ${SRR}"
        exit 1
    fi

    samtools index "${BAM_DIR}/${SRR}.sorted.bam"

    echo "${SRR} alignment completed."

done < "$ACCESSIONS"

echo "All alignments completed."
