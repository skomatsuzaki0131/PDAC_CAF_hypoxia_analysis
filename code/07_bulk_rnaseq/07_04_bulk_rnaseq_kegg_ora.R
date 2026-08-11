## KEGG over-representation analysis of bulk RNA-seq DEGs --------------------

## Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(dplyr)
  library(fs)
  library(purrr)
  library(readr)
  library(tibble)
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


## Input and output directories ----------------------------------------------

differential_expression_table_dir <- file.path(
  bulk_rnaseq_table_dir,
  "differential_expression"
)

gene_set_analysis_table_dir <- file.path(
  bulk_rnaseq_table_dir,
  "gene_set_analysis"
)

kegg_ora_table_dir <- file.path(
  gene_set_analysis_table_dir,
  "kegg_ora"
)

fs::dir_create(kegg_ora_table_dir)

significant_degs_for_ora_file <- file.path(
  differential_expression_table_dir,
  "bulk_rnaseq_significant_degs_for_ora.csv"
)

ora_background_genes_file <- file.path(
  differential_expression_table_dir,
  "bulk_rnaseq_ora_background_genes.csv"
)

required_input_files <- c(
  significant_degs_for_ora_file,
  ora_background_genes_file
)

missing_input_files <- required_input_files[
  !file.exists(required_input_files)
]

if (length(missing_input_files) > 0) {
  stop(
    "Required KEGG ORA input files were not found: ",
    paste(missing_input_files, collapse = ", "),
    call. = FALSE
  )
}


## Read KEGG ORA input tables -------------------------------------------------

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

ora_background_genes <- readr::read_csv(
  ora_background_genes_file,
  col_types = readr::cols(
    comparison_name = readr::col_character(),
    gene_symbol = readr::col_character(),
    entrez_id = readr::col_character()
  )
)

required_significant_deg_columns <- c(
  "comparison_name",
  "gene_symbol",
  "entrez_id",
  "logFC",
  "FDR",
  "regulation"
)

required_background_columns <- c(
  "comparison_name",
  "gene_symbol",
  "entrez_id"
)

missing_significant_deg_columns <- setdiff(
  required_significant_deg_columns,
  colnames(significant_degs_for_ora)
)

missing_background_columns <- setdiff(
  required_background_columns,
  colnames(ora_background_genes)
)

if (length(missing_significant_deg_columns) > 0) {
  stop(
    "The significant DEG table is missing required columns: ",
    paste(missing_significant_deg_columns, collapse = ", "),
    call. = FALSE
  )
}

if (length(missing_background_columns) > 0) {
  stop(
    "The ORA background table is missing required columns: ",
    paste(missing_background_columns, collapse = ", "),
    call. = FALSE
  )
}


## Define comparisons used in the original KEGG analysis --------------------

kegg_comparison_settings <- tibble::tribble(
  ~comparison_name, ~comparison_label,
  "comparison_1_hypocaf_hypoxia_vs_normocaf_normoxia", "Original",
  "comparison_2_hypocaf_normoxia_vs_normocaf_normoxia", "Under21",
  "comparison_3_hypocaf_hypoxia_vs_normocaf_hypoxia", "Under1"
)

missing_significant_deg_comparisons <- setdiff(
  kegg_comparison_settings$comparison_name,
  unique(significant_degs_for_ora$comparison_name)
)

missing_background_comparisons <- setdiff(
  kegg_comparison_settings$comparison_name,
  unique(ora_background_genes$comparison_name)
)

if (length(missing_significant_deg_comparisons) > 0) {
  stop(
    "The significant DEG table is missing comparisons: ",
    paste(missing_significant_deg_comparisons, collapse = ", "),
    call. = FALSE
  )
}

if (length(missing_background_comparisons) > 0) {
  stop(
    "The ORA background table is missing comparisons: ",
    paste(missing_background_comparisons, collapse = ", "),
    call. = FALSE
  )
}


## Run KEGG over-representation analysis -------------------------------------

run_kegg_ora <- function(
    significant_gene_ids,
    background_gene_ids,
    current_comparison_name,
    current_comparison_label,
    current_enrichment_direction
) {
  significant_gene_ids <- significant_gene_ids |>
    as.character() |>
    unique() |>
    stats::na.omit()

  background_gene_ids <- background_gene_ids |>
    as.character() |>
    unique() |>
    stats::na.omit()

  significant_gene_ids <- intersect(
    significant_gene_ids,
    background_gene_ids
  )

  if (length(significant_gene_ids) == 0) {
    warning(
      "No significant Entrez IDs were available for ",
      current_comparison_name,
      " (",
      current_enrichment_direction,
      ")."
    )

    return(tibble::tibble())
  }

  set.seed(1234)

  kegg_enrichment_result <- clusterProfiler::enrichKEGG(
    gene = significant_gene_ids,
    universe = background_gene_ids,
    organism = "hsa",
    keyType = "ncbi-geneid",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1
  )

  kegg_result_table <- kegg_enrichment_result |>
    as.data.frame() |>
    tibble::as_tibble()

  if (nrow(kegg_result_table) == 0) {
    return(tibble::tibble())
  }

  kegg_result_table |>
    dplyr::mutate(
      comparison_name = current_comparison_name,
      comparison_label = current_comparison_label,
      enrichment_direction = current_enrichment_direction,
      .before = 1
    )
}

kegg_ora_results <- purrr::map_dfr(
  seq_len(nrow(kegg_comparison_settings)),
  function(comparison_index) {
    current_comparison_name <-
      kegg_comparison_settings$comparison_name[[comparison_index]]

    current_comparison_label <-
      kegg_comparison_settings$comparison_label[[comparison_index]]

    current_significant_degs <- significant_degs_for_ora |>
      dplyr::filter(
        .data$comparison_name == .env$current_comparison_name
      )

    current_background_gene_ids <- ora_background_genes |>
      dplyr::filter(
        .data$comparison_name == .env$current_comparison_name
      ) |>
      dplyr::pull(.data$entrez_id)

    upregulated_gene_ids <- current_significant_degs |>
      dplyr::filter(.data$regulation == "upregulated") |>
      dplyr::pull(.data$entrez_id)

    downregulated_gene_ids <- current_significant_degs |>
      dplyr::filter(.data$regulation == "downregulated") |>
      dplyr::pull(.data$entrez_id)

    dplyr::bind_rows(
      run_kegg_ora(
        significant_gene_ids = upregulated_gene_ids,
        background_gene_ids = current_background_gene_ids,
        current_comparison_name = current_comparison_name,
        current_comparison_label = current_comparison_label,
        current_enrichment_direction = "Upregulated in Hypo-CAF"
      ),
      run_kegg_ora(
        significant_gene_ids = downregulated_gene_ids,
        background_gene_ids = current_background_gene_ids,
        current_comparison_name = current_comparison_name,
        current_comparison_label = current_comparison_label,
        current_enrichment_direction = "Downregulated in Hypo-CAF"
      )
    )
  }
)

if (nrow(kegg_ora_results) == 0) {
  stop(
    "KEGG ORA returned no enriched pathways for any comparison.",
    call. = FALSE
  )
}


## Add derived columns used by downstream figures and tables -----------------

parse_ratio <- function(ratio_text) {
  ratio_parts <- strsplit(ratio_text, "/", fixed = TRUE)

  vapply(
    ratio_parts,
    function(current_ratio_parts) {
      as.numeric(current_ratio_parts[[1]]) /
        as.numeric(current_ratio_parts[[2]])
    },
    numeric(1)
  )
}

kegg_ora_results <- kegg_ora_results |>
  dplyr::mutate(
    gene_ratio_numeric = parse_ratio(.data$GeneRatio),
    background_ratio_numeric = parse_ratio(.data$BgRatio),
    signed_gene_ratio = dplyr::if_else(
      .data$enrichment_direction == "Upregulated in Hypo-CAF",
      .data$gene_ratio_numeric,
      -.data$gene_ratio_numeric
    ),
    significance_label = dplyr::case_when(
      .data$p.adjust < 0.001 ~ "***",
      .data$p.adjust < 0.01 ~ "**",
      .data$p.adjust < 0.05 ~ "*",
      TRUE ~ "n.s."
    )
  )

significant_kegg_ora_results <- kegg_ora_results |>
  dplyr::filter(.data$p.adjust < 0.05)


## Save KEGG ORA results ------------------------------------------------------

readr::write_csv(
  kegg_ora_results,
  file = file.path(
    kegg_ora_table_dir,
    "bulk_rnaseq_kegg_ora_all_results.csv"
  ),
  na = ""
)

readr::write_csv(
  significant_kegg_ora_results,
  file = file.path(
    kegg_ora_table_dir,
    "bulk_rnaseq_kegg_ora_significant_results.csv"
  ),
  na = ""
)
