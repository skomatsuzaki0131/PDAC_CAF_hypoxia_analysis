# PDAC_CAF_hypoxia_analysis

Code for analysis of hypoxia-associated CAF subpopulations in PDAC.

## Overview

This repository contains analysis scripts used for bulk RNA-seq, Xenium spatial transcriptomic analysis, CAF subclustering, gene set enrichment analysis, hypoxia mapping, spatial niche analysis, and figure generation.

## Repository structure

```text
code/
  01_Xenium_initial_clustering/
  02_CAF_subclustering/
  03_spatial_analysis/
  04_pseudotime_analysis/
  05_CAF_isolation_assays/
  06_bulk_RNAseq/
  07_figures/

data/
  README_data.md

metadata/
  gene_sets/

results/
  README_results.md
```

## Requirements

Analyses were performed in R.

Main R packages include:

- Seurat
- harmony
- edgeR
- clusterProfiler
- GSVA
- msigdbr
- ComplexHeatmap
- ggplot2
- dplyr

Detailed package versions will be added before publication.

## Data availability

Raw and processed data availability will be described in the manuscript Data Availability Statement.

Large raw data files, patient-derived spatial transcriptomic objects, and intermediate Seurat objects are not included in this repository.

## Code availability

Custom R scripts used for bulk RNA-seq analysis, Xenium spatial transcriptomic analysis, CAF subclustering, gene set enrichment analysis, hypoxia mapping, spatial niche analysis, and figure generation will be made available in this repository upon publication.
