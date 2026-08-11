# Reproducibility check of Xenium initial clustering

We compared the Seurat object reconstructed from raw Xenium output with the processed Seurat object used in the manuscript.

## Summary

- Cell set: identical
- Cell order: identical after matching sample merge order
- Normalized data layer: identical in checked entries
- Variable feature set: identical
- Scaled data: identical in checked entries
- PCA embeddings: identical up to arbitrary sign flips
- Harmony embeddings: identical up to arbitrary sign flips and negligible numerical differences
- Nearest-neighbor/SNN graphs: slightly different
- Final graph-based clustering: differed slightly depending on graph construction

## Interpretation

The raw preprocessing script reproduces the analysis object up to the Harmony embedding level. Minor differences in nearest-neighbor graph construction and graph-based clustering may arise from package versions, approximate nearest-neighbor search, and numerical tie handling.

Therefore, downstream analyses in this repository use the processed Seurat object used in the manuscript.