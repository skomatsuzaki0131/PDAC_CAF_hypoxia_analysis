# CAF subcluster gene set scoring
#
# This script calculates single-cell ssGSEA scores and cell-cycle scores
# for CAFs, and summarizes these scores by CAF subcluster.


## Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(GSVA)
  library(msigdbr)
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
source(file.path(project_dir, "code", "00_config", "00_03_genes_and_marker_definitions.R"))

caf_subclustering_table_dir <- file.path(
  xenium_table_dir,
  "caf_subclustering"
)

## Read CAF subclustering object and metadata --------------------------------

CAFObj <- readRDS(caf_object_file)

CAFObj <- JoinLayers(CAFObj)

caf_metadata_df <- read.csv(
  caf_subcluster_metadata_file,
  stringsAsFactors = FALSE
)

stopifnot(all(c(
  "cell_id",
  "xenium_sample_id",
  "xenium_sample_number",
  "tx_id",
  "x",
  "y",
  "initial_cluster",
  "initial_cell_type",
  "caf_subcluster",
  "caf_umap_1",
  "caf_umap_2"
) %in% colnames(caf_metadata_df)))

stopifnot(all(colnames(CAFObj) %in% caf_metadata_df$cell_id))

CAFObj$caf_subcluster <- factor(
  as.character(CAFObj$caf_subcluster),
  levels = caf_subcluster_display_order
)

stopifnot(!any(is.na(CAFObj$caf_subcluster)))
stopifnot(identical(all_genes_symbol, rownames(CAFObj)))


## Define gene sets ----------------------------------------------------------

hallmark_gene_sets <- msigdbr::msigdbr(
  species = "Homo sapiens",
  category = "H"
) %>%
  split(.$gs_name) %>%
  lapply(function(x) unique(x$gene_symbol))

c2_gene_sets <- msigdbr::msigdbr(
  species = "Homo sapiens",
  category = "C2"
)

additional_gene_sets <- list(
  BuffaOrig = buffa_hypoxia_genes,
  WinterOrig = winter_hypoxia_genes,
  myCAF = myCAF_xenium_panel,
  iCAF = iCAF_xenium_panel,
  KEGG_Cell_cycle = c2_gene_sets %>%
    dplyr::filter(gs_name == "KEGG_CELL_CYCLE") %>%
    dplyr::pull(gene_symbol) %>%
    unique()
)

gene_sets <- c(additional_gene_sets, hallmark_gene_sets)

gene_sets <- lapply(
  gene_sets,
  function(gs) intersect(gs, all_genes_symbol)
)

gene_sets <- gene_sets[lengths(gene_sets) > 0]


## Calculate single-cell ssGSEA scores --------------------------------------

expr <- GetAssayData(
  CAFObj,
  assay = "RNA",
  layer = "data"
) %>%
  as.matrix()

set.seed(1234)

ssgsea_param <- GSVA::ssgseaParam(
  exprData = expr,
  geneSets = gene_sets,
  minSize = 1,
  maxSize = Inf,
  normalize = TRUE
)

ssgsea_matrix <- GSVA::gsva(
  ssgsea_param,
  verbose = TRUE
)

ssgsea_df <- as.data.frame(t(ssgsea_matrix)) %>%
  tibble::rownames_to_column(var = "cell_id")

colnames(ssgsea_df)[-1] <- paste0(colnames(ssgsea_df)[-1], "_ssGSEA")


## Add ssGSEA scores to CAF metadata -----------------------------------------

ssgsea_metadata_df <- caf_metadata_df %>%
  dplyr::left_join(
    ssgsea_df,
    by = "cell_id"
  )

score_columns <- colnames(ssgsea_metadata_df) %>%
  stringr::str_subset("_ssGSEA$")

stopifnot(length(score_columns) > 0)
stopifnot(!any(is.na(ssgsea_metadata_df[, score_columns])))


## Calculate subcluster-level mean and Z-score -------------------------------

subcluster_mean_scores <- ssgsea_metadata_df %>%
  dplyr::select(caf_subcluster, dplyr::all_of(score_columns)) %>%
  dplyr::group_by(caf_subcluster) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(score_columns),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    caf_subcluster = factor(
      caf_subcluster,
      levels = caf_subcluster_display_order
    )
  ) %>%
  dplyr::arrange(caf_subcluster) %>%
  tibble::column_to_rownames(var = "caf_subcluster") %>%
  dplyr::rename(
    Buffa_hypoxia_ssGSEA = BuffaOrig_ssGSEA,
    Winter_hypoxia_ssGSEA = WinterOrig_ssGSEA,
    myCAF_signature_ssGSEA = myCAF_ssGSEA,
    iCAF_signature_ssGSEA = iCAF_ssGSEA
  )

subcluster_mean_long <- subcluster_mean_scores %>%
  tibble::rownames_to_column(var = "caf_subcluster") %>%
  tidyr::pivot_longer(
    cols = -caf_subcluster,
    names_to = "gene_set",
    values_to = "mean_ssGSEA_score"
  )

subcluster_zscore_long <- base::scale(subcluster_mean_scores) %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "caf_subcluster") %>%
  tidyr::pivot_longer(
    cols = -caf_subcluster,
    names_to = "gene_set",
    values_to = "z_score_within_gene_set"
  )

ssgsea_summary_table <- subcluster_mean_long %>%
  dplyr::inner_join(
    subcluster_zscore_long,
    by = c("caf_subcluster", "gene_set")
  ) %>%
  dplyr::mutate(
    caf_subcluster = factor(
      caf_subcluster,
      levels = caf_subcluster_display_order
    )
  ) %>%
  dplyr::arrange(caf_subcluster, gene_set) %>%
  dplyr::mutate(
    caf_subcluster = as.character(caf_subcluster)
  )


## Calculate cell-cycle scores -----------------------------------------------

s_features <- intersect(cc.genes$s.genes, all_genes_symbol)
g2m_features <- intersect(cc.genes$g2m.genes, all_genes_symbol)

CAFObj <- CellCycleScoring(
  CAFObj,
  s.features = s_features,
  g2m.features = g2m_features,
  set.ident = FALSE
)

cell_cycle_score_table <- CAFObj@meta.data %>%
  dplyr::select(
    cell_id,
    S.Score,
    G2M.Score
  ) %>%
  dplyr::left_join(
    caf_metadata_df,
    by = "cell_id"
  ) %>%
  dplyr::select(
    cell_id,
    xenium_sample_id,
    xenium_sample_number,
    tx_id,
    x,
    y,
    initial_cluster,
    initial_cell_type,
    caf_subcluster,
    caf_umap_1,
    caf_umap_2,
    S.Score,
    G2M.Score
  )

cell_cycle_summary_table <- cell_cycle_score_table %>%
  dplyr::group_by(caf_subcluster) %>%
  dplyr::summarise(
    mean_S_score = mean(S.Score, na.rm = TRUE),
    mean_G2M_score = mean(G2M.Score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    caf_subcluster = factor(
      caf_subcluster,
      levels = caf_subcluster_display_order
    )
  ) %>%
  dplyr::arrange(caf_subcluster) %>%
  dplyr::mutate(
    caf_subcluster = as.character(caf_subcluster)
  )


## Save score outputs --------------------------------------------------------

fs::dir_create(caf_subclustering_table_dir)

write.csv(
  ssgsea_metadata_df,
  file = file.path(
    caf_subclustering_table_dir,
    "xenium_caf_single_cell_ssgsea_scores.csv"
  ),
  row.names = FALSE
)

write.csv(
  ssgsea_summary_table,
  file = file.path(
    caf_subclustering_table_dir,
    "xenium_caf_subcluster_ssgsea_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  cell_cycle_score_table,
  file = file.path(
    caf_subclustering_table_dir,
    "xenium_caf_single_cell_cell_cycle_scores.csv"
  ),
  row.names = FALSE
)

write.csv(
  cell_cycle_summary_table,
  file = file.path(
    caf_subclustering_table_dir,
    "xenium_caf_subcluster_cell_cycle_summary.csv"
  ),
  row.names = FALSE
)
