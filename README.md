
##Name
Theresia Ngonda Nfah

## Program Title
Bioinformatics Training: RNA-seq Data Analysis Project

## Project Description
RNA-seq, also called RNA sequencing, is a next-generation sequencing method used to study the complete set of RNA molecules in a biological sample. RNA-seq is commonly described as a method for sequencing RNA molecules, converting RNA to DNA for sequencing, and determining both transcription and gene expression levels.It helps researchers identify which genes are active, measure gene expression levels, and compare transcript patterns between samples or treatment groups.In this project, GitHub version control is used to organize scripts, data folders, results, documentation, and workflow outputs for reproducible RNA-seq analysis.

This project uses six publicly available RNA-seq samples from NCBI GEO/SRA series GSE96870. The dataset examines transcriptomic changes in the central nervous system of Mus musculus following upper-respiratory Influenza A infection. The selected samples include cerebellum RNA-seq data from non-infected Day 0 mice and Influenza A-infected mice at post-infection time points. The sequencing data were generated using Illumina HiSeq 2500 paired-end RNA-seq reads. Raw reads were downloaded from SRA, assessed using FastQC, summarized with MultiQC, and trimmed with fastp where adapter or low-quality sequence issues were detected.

## HISAT2 Alignment Summary

The trimmed paired-end RNA-seq reads were aligned to the Mus musculus GRCm39 reference genome using HISAT2. Alignment summaries were generated for each sample, and the overall alignment rates are shown below.

| SampleID | Overall alignment rate | Flag |
|---|---:|---|
| SRR5364316 | paste_rate_here | PASS |
| SRR5364317 | paste_rate_here | PASS |
| SRR5364318 | paste_rate_here | PASS |
| SRR5364321 | paste_rate_here | PASS |
| SRR5364323 | paste_rate_here | PASS |
| SRR5364330 | paste_rate_here | PASS |

Samples with overall alignment rate below 75% were flagged as low-alignment samples. A low alignment rate may result from poor read quality, adapter contamination, wrong reference genome, sample contamination, or using reads from a different organism.

SRRXXXXXXX had an alignment rate below 75%. This may be due to poor sequencing quality, adapter contamination, sample contamination, or mismatch between the reads and the selected reference genome.

