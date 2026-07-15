#!/bin/bash

# Folder setup script for RNA-seq project

mkdir -p data/raw
mkdir -p data/processed
mkdir -p data/reference
mkdir -p scripts
mkdir -p results/fastqc
mkdir -p results/multiqc
mkdir -p results/trimmed_reads
mkdir -p results/alignment
mkdir -p results/counts
mkdir -p docs
mkdir -p logs

# Keep empty folders visible to Git
touch data/raw/.gitkeep
touch data/processed/.gitkeep
touch data/reference/.gitkeep
touch scripts/.gitkeep
touch results/fastqc/.gitkeep
touch results/multiqc/.gitkeep
touch results/trimmed_reads/.gitkeep
touch results/alignment/.gitkeep
touch results/counts/.gitkeep
touch docs/.gitkeep
touch logs/.gitkeep

echo "RNA-seq project folder structure created successfully."
