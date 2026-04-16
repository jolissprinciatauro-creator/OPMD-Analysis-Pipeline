# OPMD Transcriptomics Analysis Pipeline

## Project Overview
This repository contains a robust bioinformatics pipeline developed for the analysis of **Oral Premalignant Disorders (OPMD)**. Using R and Bioconductor, the pipeline processes transcriptomic data to identify key genetic drivers and biological pathways involved in oral lesion progression.

This project was developed during my tenure at the **Mazumdar Shaw Medical Foundation (MSMF)**.

## Pipeline Features
* **Automated Data Retrieval:** Integrated with `GEOquery` to fetch raw expression data directly from NCBI GEO.
* **Quality Control (QC):** Principal Component Analysis (PCA) to validate sample clustering and data integrity.
* **Differential Gene Expression (DGE):** Statistical modeling using `limma` to identify significantly upregulated and downregulated genes.
* **Expression Profiling:** Visualization of top DEGs through high-dimensional **ComplexHeatmaps**.
* **Functional Enrichment:** Gene Ontology (GO) and Pathway analysis using `clusterProfiler`.
* **Network Visualization:** Gene-Concept Networks (Cnetplots) to visualize the intersection of genes and biological pathways.

## Tools & Libraries
* **Language:** R
* **Bioinformatics:** Bioconductor, limma, GEOquery, clusterProfiler, org.Hs.eg.db
* **Visualization:** ggplot2, ComplexHeatmap, enrichplot

## How to Use
1. Clone the repository.
2. Open the `analysis_pipeline.R` script in RStudio.
3. Update the `gse_id` variable to analyze different datasets (Default is set to GSE10174).
4. Run the script to generate statistical results and plots.

## Author
**Joliss Princia Tauro** Computational Biologist intern | MSMF