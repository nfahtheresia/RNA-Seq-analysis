#!/bin/bash
set -euo pipefail

OUT="results/alignment/summaries/alignment_rates.tsv"

mkdir -p results/alignment/summaries

echo -e "SampleID\tOverall_Alignment_Rate\tFlag" > "$OUT"

for LOG in results/alignment/logs/*_hisat2_summary.txt; do
    SAMPLE=$(basename "$LOG" _hisat2_summary.txt)
    RATE=$(grep "overall alignment rate" "$LOG" | awk '{print $1}' | sed 's/%//')

    if awk "BEGIN {exit !($RATE < 75)}"; then
        FLAG="FLAG_LOW_ALIGNMENT"
    else
        FLAG="PASS"
    fi

    echo -e "${SAMPLE}\t${RATE}%\t${FLAG}" >> "$OUT"
done

cat "$OUT"
