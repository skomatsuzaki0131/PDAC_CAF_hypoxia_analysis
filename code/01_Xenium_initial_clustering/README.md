# Xenium initial clustering

This folder contains scripts for Xenium spatial transcriptomic data preprocessing, quality control, integration, dimensionality reduction, clustering, and cell type annotation.

## Scripts

- `TX5K_AnalyzeCAFbyIntegration_raw.R`: Exploratory Xenium integration and CAF analysis script before final parameter fixation.
- `TX5K_AnalyzeCAFbyIntegration_Final_raw.R`: Final unrefactored script using the parameter settings applied in the manuscript analyses.

## Notes

The current scripts are archived as raw analysis records. Absolute paths and local file names will be replaced with project-relative paths before public release.

`TX5K_AnalyzeCAFbyIntegration_raw.R` contains exploratory code before final parameter fixation.

`TX5K_AnalyzeCAFbyIntegration_Final_raw.R` contains the code corresponding to the final parameter settings used for the manuscript analyses. This script is currently unrefactored and will be organized before publication. This script includes Xenium initial clustering, CAF subclustering, spatial analyses, pseudotime analyses, bulk RNA-seq-related downstream analyses, supplemental table generation, and figure-generation code. It will be refactored into analysis-specific scripts before publication.

