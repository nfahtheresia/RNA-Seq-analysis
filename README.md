
## Name
Theresia Ngonda Nfah

## Program Title
Bioinformatics Training: RNA-seq Data Analysis Project

## Project Description
RNA-seq, also called RNA sequencing, is a next-generation sequencing method used to study the complete set of RNA molecules in a biological sample. RNA-seq is commonly described as a method for sequencing RNA molecules, converting RNA to DNA for sequencing, and determining both transcription and gene expression levels.It helps researchers identify which genes are active, measure gene expression levels, and compare transcript patterns between samples or treatment groups.In this project, GitHub version control is used to organize scripts, data folders, results, documentation, and workflow outputs for reproducible RNA-seq analysis.

This project uses six publicly available RNA-seq samples from NCBI GEO/SRA series GSE96870. The dataset examines transcriptomic changes in the central nervous system of Mus musculus following upper-respiratory Influenza A infection. The selected samples include cerebellum RNA-seq data from non-infected Day 0 mice and Influenza A-infected mice at post-infection time points. The sequencing data were generated using Illumina HiSeq 2500 paired-end RNA-seq reads. Raw reads were downloaded from SRA, assessed using FastQC, summarized with MultiQC, and trimmed with fastp where adapter or low-quality sequence issues were detected.

## HISAT2 Alignment Summary

The trimmed paired-end RNA-seq reads were aligned to the Mus musculus GRCm39 reference genome using HISAT2. Alignment summaries were generated for each sample, and the overall alignment rates are shown below.

Samples with overall alignment rate below 75% were flagged as low-alignment samples. A low alignment rate may result from poor read quality, adapter contamination, wrong reference genome, sample contamination, or using reads from a different organism.

## Combined Count Matrix

featureCounts was run on UseGalaxy in separate batches for samples belonging to the same BioProject. The individual Galaxy count outputs were downloaded, copied into `results/counts/galaxy_runs/`, and merged by `Geneid` using the script `scripts/08_combine_galaxy_counts.py`.

Final combined count matrix:

- `results/counts/combined_gene_count_matrix.tsv`
- `results/counts/combined_gene_count_matrix.tsv.gz`

## Differential Gene Expression Analysis

Differential gene expression analysis was performed using DESeq2. The combined gene count matrix was imported into R, matched with the sample metadata, filtered for low-count genes, and analysed using the design formula `~ Condition`.

### DEG Summary

| Contrast | Total genes tested | Significant DEGs | Upregulated | Downregulated |
|---|---:|---:|---:|---:|
| InfluenzaA_Day4_vs_NonInfected_Day0 | 29722 | 21 | 19 | 2 |
| InfluenzaA_Day8_vs_NonInfected_Day0 | 29722 | 291 | 174 | 117 |

Significant DEGs were defined as genes with adjusted p-value < 0.05 and absolute log2 fold change ≥ 1.

### Output files

- `scripts/09_deseq2_analysis.R`
- `results/deseq2/DEG_summary.csv`
- `results/deseq2/vst_normalized_count_matrix.csv`
- `results/deseq2/DESeq2_significant_DEGs_InfluenzaA_Day4_vs_NonInfected_Day0.csv`
- `results/deseq2/DESeq2_significant_DEGs_InfluenzaA_Day8_vs_NonInfected_Day0.csv`

## Visualisation, Pathway Analysis and Final Report

This repository contains a complete RNA-seq analysis workflow for mouse influenza A infection data. The project includes raw data acquisition, quality control, trimming, genome alignment, read counting, differential gene expression analysis, visualization, GO Biological Process enrichment, KEGG pathway enrichment, and final report generation.

### Project Topic

Differential gene expression analysis of mouse cerebellum during influenza A infection

### Methods Summary

RNA-seq reads were processed using FastQC, MultiQC, fastp, HISAT2, samtools, and featureCounts. Gene-level count data were imported into R and analysed using DESeq2. Variance-stabilising transformation was used for exploratory visualization. PCA, volcano plots, top DEG heatmaps, and sample distance heatmaps were generated. Significant genes were analysed for GO Biological Process and KEGG pathway enrichment using clusterProfiler and org.Mm.eg.db.

### Output Files

- `scripts/10_visualization_pathway_analysis.R`
- `results/visualization_pathway/visualization_DEG_GO_KEGG_summary.csv`
- `results/visualization_pathway/vst_normalized_count_matrix_visualization.csv`
- `results/visualization_pathway/figures/`
- `results/visualization_pathway/enrichment/`
- `reports/final_report.Rmd`
- `docs/final_report.html`

### Final Report

The final HTML report is available at:

`docs/final_report.html`

### DEG, GO and KEGG Summary

Please see:

`results/visualization_pathway/visualization_DEG_GO_KEGG_summary.csv`

This file reports the total genes tested, significant DEGs, upregulated genes, downregulated genes, top 5 GO Biological Process terms, and top 3 KEGG pathways for each contrast.
