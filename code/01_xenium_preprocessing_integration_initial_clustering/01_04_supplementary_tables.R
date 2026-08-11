# Initial clustering supplementary tables
#
# This script generates summary tables for Xenium initial clustering,
# including initial cluster counts and cell type composition across samples.


## Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
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
source(file.path(project_dir, "code", "00_config", "00_04_annotation_definitions.R"))

initial_clustering_table_dir <- file.path(
  xenium_table_dir,
  "initial_clustering"
)

fs::dir_create(initial_clustering_table_dir)


## Read initial clustering metadata -----------------------------------------

initial_metadata_df <- read.csv(
  initial_metadata_file,
  stringsAsFactors = FALSE
)

initial_celltype_order <- names(initial_cluster_annotation_list)


## Initial cluster summary ---------------------------------------------------

initial_cluster_summary <- initial_metadata_df %>%
  dplyr::count(
    initial_cluster,
    initial_cell_type,
    name = "n_cell"
  ) %>%
  dplyr::mutate(
    pct_of_total = 100 * n_cell / sum(n_cell),
    initial_cluster = paste0("C", initial_cluster)
  ) %>%
  dplyr::arrange(
    as.numeric(stringr::str_remove(initial_cluster, "^C"))
  ) %>%
  dplyr::select(
    initial_cluster,
    initial_cell_type,
    n_cell,
    pct_of_total
  )

write.csv(
  initial_cluster_summary,
  file = file.path(
    initial_clustering_table_dir,
    "xenium_initial_cluster_summary.csv"
  ),
  row.names = FALSE
)


## Cell type summary ---------------------------------------------------------

initial_celltype_summary <- initial_metadata_df %>%
  dplyr::count(
    initial_cell_type,
    name = "n_cell"
  ) %>%
  dplyr::mutate(
    pct_of_total = 100 * n_cell / sum(n_cell),
    initial_cell_type = factor(
      initial_cell_type,
      levels = initial_celltype_order
    )
  ) %>%
  dplyr::arrange(initial_cell_type) %>%
  dplyr::mutate(
    initial_cell_type = as.character(initial_cell_type)
  )

write.csv(
  initial_celltype_summary,
  file = file.path(
    initial_clustering_table_dir,
    "xenium_initial_celltype_summary.csv"
  ),
  row.names = FALSE
)


## Cell type composition across samples --------------------------------------

initial_celltype_composition_across_samples <- 
  initial_metadata_df %>%
  dplyr::count(
    initial_cell_type,
    xenium_sample_id,
    name = "n_cell"
  ) %>%
  dplyr::group_by(xenium_sample_id) %>%
  dplyr::mutate(
    pct_within_sample = 100 * n_cell / sum(n_cell),
    cell_number_and_percent = paste0(
      formatC(n_cell, big.mark = ","),
      " (",
      sprintf("%.1f", pct_within_sample),
      "%)"
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    initial_cell_type,
    xenium_sample_id,
    cell_number_and_percent
  ) %>%
  tidyr::pivot_wider(
    names_from = xenium_sample_id,
    values_from = cell_number_and_percent
  )

initial_celltype_composition_across_samples <- 
  initial_celltype_composition_across_samples %>%
  dplyr::mutate(
    initial_cell_type = factor(
      initial_cell_type,
      levels = initial_celltype_order
    )
  ) %>%
  dplyr::arrange(initial_cell_type) %>%
  dplyr::mutate(
    initial_cell_type = as.character(initial_cell_type)
  )

sample_totals <- initial_metadata_df %>%
  dplyr::count(
    xenium_sample_id,
    name = "n_cell"
  ) %>%
  dplyr::mutate(
    total_label = paste0("n = ", formatC(n_cell, big.mark = ","))
  ) %>%
  dplyr::select(
    xenium_sample_id,
    total_label
  ) %>%
  tidyr::pivot_wider(
    names_from = xenium_sample_id,
    values_from = total_label
  ) %>%
  dplyr::mutate(
    initial_cell_type = "Total cells"
  ) %>%
  dplyr::select(
    initial_cell_type,
    dplyr::everything()
  )

initial_celltype_composition_across_samples <- dplyr::bind_rows(
  sample_totals,
  initial_celltype_composition_across_samples
)

write.csv(
  initial_celltype_composition_across_samples,
  file = file.path(
    initial_clustering_table_dir,
    "xenium_initial_celltype_composition_across_samples.csv"
  ),
  row.names = FALSE
)
