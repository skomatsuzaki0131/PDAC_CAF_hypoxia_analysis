## Cell cycle pathway mapping of bulk RNA-seq DEGs --------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(fs)
  library(KEGGREST)
  library(pathview)
  library(readr)
  library(tidyr)
})

## Load settings -------------------------------------------------------------

project_dir <- getwd()

if (basename(project_dir) != "PDAC_CAF_hypoxia_analysis") {
  stop(
    "Please open PDAC_CAF_hypoxia_analysis.Rproj before running this script."
  )
}

source(file.path(
  project_dir,
  "code",
  "00_config",
  "00_01_analysis_parameters.R"
))

source(file.path(
  project_dir,
  "code",
  "00_config",
  "00_02_paths.R"
))


## Input and output directories ---------------------------------------------

differential_expression_table_dir <- file.path(
  bulk_rnaseq_table_dir,
  "differential_expression"
)

cell_cycle_pathview_dir <- file.path(
  bulk_rnaseq_figure_dir,
  "08_07_05_cell_cycle_pathview"
)

fs::dir_create(cell_cycle_pathview_dir)

significant_degs_for_ora_file <- file.path(
  differential_expression_table_dir,
  "bulk_rnaseq_significant_degs_for_ora.csv"
)

if (!file.exists(significant_degs_for_ora_file)) {
  stop(
    "Required DEG input file was not found: ",
    significant_degs_for_ora_file,
    call. = FALSE
  )
}


## Read significant DEGs ----------------------------------------------------

significant_degs_for_ora <- readr::read_csv(
  significant_degs_for_ora_file,
  col_types = readr::cols(
    comparison_name = readr::col_character(),
    gene_symbol = readr::col_character(),
    entrez_id = readr::col_character(),
    logFC = readr::col_double(),
    FDR = readr::col_double(),
    regulation = readr::col_character()
  )
)

required_columns <- c(
  "comparison_name",
  "entrez_id",
  "logFC",
  "FDR"
)

missing_columns <- setdiff(
  required_columns,
  colnames(significant_degs_for_ora)
)

if (length(missing_columns) > 0) {
  stop(
    "The significant DEG table is missing required columns: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}


## Define the three comparison settings -------------------------------------

pathview_comparison_settings <- tibble::tribble(
  ~comparison_name, ~comparison_label,
  "comparison_1_hypocaf_hypoxia_vs_normocaf_normoxia", "Original",
  "comparison_2_hypocaf_normoxia_vs_normocaf_normoxia", "Under21",
  "comparison_3_hypocaf_hypoxia_vs_normocaf_hypoxia", "Under1"
)

missing_comparisons <- setdiff(
  pathview_comparison_settings$comparison_name,
  unique(significant_degs_for_ora$comparison_name)
)

if (length(missing_comparisons) > 0) {
  stop(
    "The significant DEG table is missing comparisons: ",
    paste(missing_comparisons, collapse = ", "),
    call. = FALSE
  )
}

if (!exists("bulk_rnaseq_comparison_plot_labels")) {
  stop(
    "bulk_rnaseq_comparison_plot_labels was not found in 00_01_analysis_parameters.R.",
    call. = FALSE
  )
}

missing_plot_labels <- setdiff(
  pathview_comparison_settings$comparison_label,
  names(bulk_rnaseq_comparison_plot_labels)
)

if (length(missing_plot_labels) > 0) {
  stop(
    "bulk_rnaseq_comparison_plot_labels is missing labels for: ",
    paste(missing_plot_labels, collapse = ", "),
    call. = FALSE
  )
}


## Retrieve genes in the KEGG cell cycle pathway ----------------------------

cell_cycle_pathway_id <- "hsa04110"

cell_cycle_pathway <- KEGGREST::keggGet(cell_cycle_pathway_id)[[1]]
cell_cycle_gene_entries <- cell_cycle_pathway$GENE

cell_cycle_entrez_ids <- cell_cycle_gene_entries[
  seq(1, length(cell_cycle_gene_entries), by = 2)
] |>
  as.character() |>
  unique()


## Build the three-column log2 fold-change matrix ----------------------------

pathview_gene_table <- significant_degs_for_ora |>
  dplyr::filter(
    .data$comparison_name %in% pathview_comparison_settings$comparison_name,
    !is.na(.data$entrez_id),
    .data$entrez_id %in% cell_cycle_entrez_ids
  ) |>
  dplyr::select(
    .data$comparison_name,
    .data$entrez_id,
    .data$logFC
  ) |>
  dplyr::distinct(
    .data$comparison_name,
    .data$entrez_id,
    .keep_all = TRUE
  ) |>
  tidyr::pivot_wider(
    names_from = .data$comparison_name,
    values_from = .data$logFC
  )

pathview_logfc_matrix <- matrix(
  NA_real_,
  nrow = length(cell_cycle_entrez_ids),
  ncol = nrow(pathview_comparison_settings),
  dimnames = list(
    cell_cycle_entrez_ids,
    pathview_comparison_settings$comparison_name
  )
)

for (comparison_name in pathview_comparison_settings$comparison_name) {
  current_values <- pathview_gene_table |>
    dplyr::select(
      .data$entrez_id,
      dplyr::all_of(comparison_name)
    ) |>
    dplyr::filter(!is.na(.data[[comparison_name]]))

  pathview_logfc_matrix[
    current_values$entrez_id,
    comparison_name
  ] <- current_values[[comparison_name]]
}

colnames(pathview_logfc_matrix) <- unname(
  bulk_rnaseq_comparison_plot_labels[
    pathview_comparison_settings$comparison_label
  ]
)


## Save the mapped values used for Pathview ----------------------------------

pathview_input_table <- pathview_logfc_matrix |>
  as.data.frame() |>
  tibble::rownames_to_column("entrez_id")

readr::write_csv(
  pathview_input_table,
  file.path(
    cell_cycle_pathview_dir,
    "bulk_rnaseq_cell_cycle_pathview_input.csv"
  )
)


## Generate the KEGG-native Pathview image ----------------------------------

original_working_directory <- getwd()
on.exit(setwd(original_working_directory), add = TRUE)
setwd(cell_cycle_pathview_dir)

pathview::pathview(
  gene.data = pathview_logfc_matrix,
  pathway.id = cell_cycle_pathway_id,
  species = "hsa",
  limit = list(gene = 3, cpd = 1),
  low = list(gene = "red", cpd = "black"),
  mid = list(gene = "white", cpd = "black"),
  high = list(gene = "blue", cpd = "black"),
  same.layer = FALSE,
  kegg.native = TRUE,
  plot.col.key = FALSE,
  key.pos = "bottomright",
  out.suffix = "bulk_rnaseq_cell_cycle"
)

pathview_output_file <- file.path(
  cell_cycle_pathview_dir,
  paste0(
    cell_cycle_pathway_id,
    ".bulk_rnaseq_cell_cycle.multi.png"
  )
)

final_output_file <- file.path(
  cell_cycle_pathview_dir,
  "bulk_rnaseq_cell_cycle_pathview.png"
)

if (!file.exists(pathview_output_file)) {
  stop(
    "The expected Pathview output file was not generated: ",
    pathview_output_file,
    call. = FALSE
  )
}

if (file.exists(final_output_file)) {
  file.remove(final_output_file)
}

file.rename(
  from = pathview_output_file,
  to = final_output_file
)
