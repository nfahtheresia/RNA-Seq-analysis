#!/usr/bin/env python3

import pandas as pd
from pathlib import Path
import re

input_dir = Path("results/counts/galaxy_runs")
output_file = Path("results/counts/combined_gene_count_matrix.tsv")


files = []
for pattern in ["**/*.txt", "**/*.tsv", "**/*.tabular"]:
    files.extend(input_dir.glob(pattern))


if not files:
    raise SystemExit("ERROR: No uncompressed count files found in results/counts/galaxy_runs")

merged = None

def clean_sample_name(name):
    name = str(name)
    name = name.split("/")[-1]
    name = name.replace(".sorted.bam", "")
    name = name.replace(".bam", "")
    name = name.replace(".tabular", "")
    name = name.replace(".txt", "")
    name = name.replace(".tsv", "")
    name = re.sub(r"Galaxy[0-9]+[-_]*", "", name)
    name = re.sub(r"[^A-Za-z0-9_]+", "_", name)
    name = name.strip("_")
    return name

for file in sorted(files):
    print(f"Reading: {file}")

    df = pd.read_csv(file, sep="\t", comment="#")

    if "Geneid" not in df.columns:
        raise SystemExit(f"ERROR: Geneid column not found in {file}")

    # featureCounts files usually have annotation columns:
    # Geneid, Chr, Start, End, Strand, Length, then count columns
    if "Length" in df.columns:
        length_index = list(df.columns).index("Length")
        count_columns = list(df.columns)[length_index + 1:]
    else:
        count_columns = [col for col in df.columns if col != "Geneid"]

    keep = ["Geneid"] + count_columns
    df = df[keep]

    new_columns = ["Geneid"]
    for col in count_columns:
        cleaned = clean_sample_name(col)

        # If Galaxy gives a poor column name, use the filename
        if cleaned == "" or cleaned.lower().startswith("unnamed"):
            cleaned = clean_sample_name(file.stem)

        new_columns.append(cleaned)

    df.columns = new_columns

    if merged is None:
        merged = df
    else:
        merged = pd.merge(merged, df, on="Geneid", how="outer")

merged = merged.fillna(0)

# Convert count columns to integers
for col in merged.columns:
    if col != "Geneid":
        merged[col] = merged[col].astype(int)

output_file.parent.mkdir(parents=True, exist_ok=True)
merged.to_csv(output_file, sep="\t", index=False)

print(f"Combined count matrix saved to: {output_file}")
print(f"Number of genes: {merged.shape[0]}")
print(f"Number of samples: {merged.shape[1] - 1}")
