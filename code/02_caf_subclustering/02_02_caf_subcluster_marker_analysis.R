# CAF subcluster marker analysis
#
# This script identifies differentially expressed genes for each CAF
# subcluster using the annotated CAF subclustering object.


## Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(presto)
  library(fs)
})


## Load settings -------------------------------------------------------------

project_dir <- getwd()

if (basename(project_dir) != "PDAC_CAF_hypoxia_analysis") {
  stop(
    "Please open PDAC_CAF_hypoxia_analysis.Rproj before running this script."
  )
}

source(file.path(project_dir, "code", "00_config", "00_01_analysis_parameters.R"))
source(file.path(project_dir, "code", "00_config", "00_02_paths.R"))

caf_subclustering_table_dir <- file.path(
  xenium_table_dir,
  "caf_subclustering"
)


## Read Xenium gene identifier mapping --------------------------------------

xenium_gene_identifier_mapping <- read.csv(
  file.path(
    xenium_table_dir,
    "xenium_gene_symbol_entrez_mapping.csv"
  ),
  stringsAsFactors = FALSE,
  na.strings = ""
)

required_mapping_columns <- c(
  "gene_symbol",
  "entrez_id"
)

missing_mapping_columns <- setdiff(
  required_mapping_columns,
  colnames(xenium_gene_identifier_mapping)
)

if (length(missing_mapping_columns) > 0) {
  stop(
    paste0(
      "The Xenium gene identifier mapping table is missing the following ",
      "required column(s): ",
      paste(missing_mapping_columns, collapse = ", "),
      "."
    ),
    call. = FALSE
  )
}

if (nrow(xenium_gene_identifier_mapping) == 0) {
  stop(
    "The Xenium gene identifier mapping table contains no genes.",
    call. = FALSE
  )
}

if (anyNA(xenium_gene_identifier_mapping$gene_symbol)) {
  stop(
    "Missing gene symbols were found in the Xenium gene identifier mapping table.",
    call. = FALSE
  )
}

if (any(xenium_gene_identifier_mapping$gene_symbol == "")) {
  stop(
    "Empty gene symbols were found in the Xenium gene identifier mapping table.",
    call. = FALSE
  )
}

if (anyDuplicated(xenium_gene_identifier_mapping$gene_symbol)) {
  stop(
    "Duplicated gene symbols were found in the Xenium gene identifier mapping table.",
    call. = FALSE
  )
}

xenium_gene_identifier_mapping <- xenium_gene_identifier_mapping %>%
  dplyr::select(
    gene_symbol,
    entrez_id
  )


## Read annotated CAF subclustering object -----------------------------------

CAFObj <- readRDS(caf_object_file)

CAFObj <- JoinLayers(CAFObj)

stopifnot("caf_subcluster" %in% colnames(CAFObj@meta.data))

CAFObj$caf_subcluster <- factor(
  as.character(CAFObj$caf_subcluster),
  levels = caf_subcluster_display_order
)
stopifnot(!any(is.na(CAFObj$caf_subcluster)))

Idents(CAFObj) <- CAFObj$caf_subcluster


## Find differentially expressed genes across CAF subclusters ----------------

set.seed(123)

caf_subcluster_markers_all <- FindAllMarkers(
  CAFObj,
  slot = "data",
  only.pos = FALSE,
  min.pct = 0.00,
  min.cell.feature = 0,
  min.cells.group = 0,
  logfc.threshold = 0,
  return.thresh = 1.00,
  test.use = "wilcox"
)

marker_gene_symbols_missing_from_mapping <- setdiff(
  unique(caf_subcluster_markers_all$gene),
  xenium_gene_identifier_mapping$gene_symbol
)

if (length(marker_gene_symbols_missing_from_mapping) > 0) {
  stop(
    paste0(
      "The following marker gene symbol(s) were not found in the Xenium gene ",
      "identifier mapping table: ",
      paste(marker_gene_symbols_missing_from_mapping, collapse = ", "),
      "."
    ),
    call. = FALSE
  )
}

number_of_marker_rows_before_identifier_join <- nrow(
  caf_subcluster_markers_all
)

caf_subcluster_markers_all <- caf_subcluster_markers_all %>%
  dplyr::left_join(
    xenium_gene_identifier_mapping,
    by = c("gene" = "gene_symbol")
  ) %>%
  dplyr::mutate(
    cluster = as.character(cluster),
    cluster = factor(cluster, levels = caf_subcluster_display_order)
  ) %>%
  dplyr::arrange(
    cluster,
    dplyr::desc(avg_log2FC)
  )

if (nrow(caf_subcluster_markers_all) !=
    number_of_marker_rows_before_identifier_join) {
  stop(
    "The number of CAF marker rows changed after joining Entrez IDs.",
    call. = FALSE
  )
}

## Save full CAF subcluster DEG table ----------------------------------------

fs::dir_create(caf_subclustering_table_dir)

write.csv(
  caf_subcluster_markers_all,
  file = file.path(
    caf_subclustering_table_dir,
    "xenium_caf_subcluster_markers_all.csv"
  ),
  row.names = FALSE
)


## Define CAF-8 marker gene set ----------------------------------------------

caf8_marker_table <- caf_subcluster_markers_all %>%
  dplyr::mutate(
    cluster = as.character(cluster)
  ) %>%
  dplyr::filter(
    cluster == proliferating_caf_subcluster,
    avg_log2FC > 1,
    p_val_adj < 0.05,
    pct.1 > 0.2
  ) %>%
  dplyr::arrange(dplyr::desc(avg_log2FC))

write.csv(
  caf8_marker_table,
  file = file.path(
    caf_subclustering_table_dir,
    "xenium_caf8_marker_genes.csv"
  ),
  row.names = FALSE
)
