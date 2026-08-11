## Hallmark gene set enrichment analysis of bulk RNA-seq results -------------


## Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(dplyr)
  library(fs)
  library(msigdbr)
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

source(file.path(
  project_dir,
  "code",
  "00_config",
  "00_03_genes_and_marker_definitions.R"
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

hallmark_gsea_table_dir <- file.path(
  gene_set_analysis_table_dir,
  "hallmark_gsea"
)


## Define the five differential-expression comparisons -----------------------

hallmark_gsea_comparison_settings <- tibble::tribble(
  ~comparison_name, ~comparison_label, ~comparison_category,
  "comparison_1_hypocaf_hypoxia_vs_normocaf_normoxia",
  "Original",
  "Isolation oxygen comparison",
  "comparison_2_hypocaf_normoxia_vs_normocaf_normoxia",
  "Under21",
  "Isolation oxygen comparison",
  "comparison_3_hypocaf_hypoxia_vs_normocaf_hypoxia",
  "Under1",
  "Isolation oxygen comparison",
  "comparison_4_hypocaf_hypoxia_vs_hypocaf_normoxia",
  "Switch_HypoCAF",
  "Assay oxygen comparison",
  "comparison_5_normocaf_hypoxia_vs_normocaf_normoxia",
  "Switch_NormoCAF",
  "Assay oxygen comparison"
) |>
  dplyr::mutate(
    differential_expression_file = file.path(
      differential_expression_table_dir,
      paste0(.data$comparison_name, ".csv")
    )
  )

missing_input_files <-
  hallmark_gsea_comparison_settings$differential_expression_file[
    !file.exists(
      hallmark_gsea_comparison_settings$differential_expression_file
    )
  ]

if (length(missing_input_files) > 0) {
  stop(
    "Required differential-expression files were not found: ",
    paste(missing_input_files, collapse = ", "),
    call. = FALSE
  )
}


## Obtain Hallmark gene sets from msigdbr ------------------------------------

# The original analysis read the MSigDB Hallmark collection from a local GMT
# file. For the public workflow, the same collection is obtained directly with
# msigdbr so that no separately downloaded GMT file is required.
hallmark_term2gene <- get_hallmark_term2gene()

if (nrow(hallmark_term2gene) == 0) {
  stop(
    "No Hallmark gene sets were returned by msigdbr.",
    call. = FALSE
  )
}


## Read and validate one differential-expression table -----------------------

read_differential_expression_table <- function(
    differential_expression_file,
    comparison_name
) {
  differential_expression_table <- readr::read_csv(
    differential_expression_file,
    col_types = readr::cols(
      gene_symbol = readr::col_character(),
      entrez_id = readr::col_character(),
      logFC = readr::col_double(),
      logCPM = readr::col_double(),
      F = readr::col_double(),
      PValue = readr::col_double(),
      FDR = readr::col_double()
    )
  )

  required_columns <- c(
    "gene_symbol",
    "entrez_id",
    "logFC",
    "FDR"
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(differential_expression_table)
  )

  if (length(missing_columns) > 0) {
    stop(
      "The differential-expression table for ",
      comparison_name,
      " is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  differential_expression_table
}


## Create the ranked Entrez-ID vector used by GSEA ---------------------------

create_ranked_gene_list <- function(
    differential_expression_table,
    comparison_name
) {
  ranked_gene_table <- differential_expression_table |>
    dplyr::filter(
      !is.na(.data$entrez_id),
      .data$entrez_id != "",
      !is.na(.data$logFC),
      is.finite(.data$logFC)
    )

  duplicated_entrez_ids <- ranked_gene_table |>
    dplyr::count(.data$entrez_id, name = "number_of_rows") |>
    dplyr::filter(.data$number_of_rows > 1)

  if (nrow(duplicated_entrez_ids) > 0) {
    stop(
      "Duplicated Entrez IDs were found in ",
      comparison_name,
      ". Resolve duplicated mappings before running Hallmark GSEA: ",
      paste(duplicated_entrez_ids$entrez_id, collapse = ", "),
      call. = FALSE
    )
  }

  # The differential-expression tables are ranked by decreasing logFC.
  # Sorting explicitly here makes the ranked-list requirement clear.
  ranked_gene_table <- ranked_gene_table |>
    dplyr::arrange(dplyr::desc(.data$logFC))

  ranked_gene_list <- stats::setNames(
    ranked_gene_table$logFC,
    ranked_gene_table$entrez_id
  )

  if (length(ranked_gene_list) == 0) {
    stop(
      "No ranked genes with Entrez IDs were available for ",
      comparison_name,
      ".",
      call. = FALSE
    )
  }

  if (any(diff(ranked_gene_list) > 0)) {
    stop(
      "The ranked gene list is not ordered by decreasing logFC for ",
      comparison_name,
      ".",
      call. = FALSE
    )
  }

  ranked_gene_list
}


## Run Hallmark GSEA with the original analysis settings ---------------------

run_hallmark_gsea <- function(ranked_gene_list) {
  set.seed(1234)

  clusterProfiler::GSEA(
    geneList = ranked_gene_list,
    minGSSize = 10,
    TERM2GENE = hallmark_term2gene,
    pvalueCutoff = 1,
    verbose = FALSE,
    eps = 0
  )
}


## Format one Hallmark GSEA result table -------------------------------------

format_hallmark_gsea_result_table <- function(
    hallmark_gsea_result,
    comparison_name,
    comparison_label,
    comparison_category
) {
  result_table <- hallmark_gsea_result@result |>
    tibble::as_tibble()

  if (nrow(result_table) == 0) {
    warning(
      "Hallmark GSEA returned no pathways for ",
      comparison_name,
      "."
    )

    return(tibble::tibble())
  }

  result_table |>
    dplyr::mutate(
      comparison_name = comparison_name,
      comparison_label = comparison_label,
      comparison_category = comparison_category,
      enrichment_direction = dplyr::case_when(
        comparison_category == "Isolation oxygen comparison" & .data$NES > 0 ~
          "Upregulated in Hypo-CAF",
        comparison_category == "Isolation oxygen comparison" & .data$NES < 0 ~
          "Upregulated in Normo-CAF",
        comparison_category == "Assay oxygen comparison" & .data$NES > 0 ~
          "Upregulated under 1% O2",
        comparison_category == "Assay oxygen comparison" & .data$NES < 0 ~
          "Upregulated under 21% O2",
        TRUE ~ NA_character_
      ),
      minus_log10_adjusted_p_value = -log10(.data$p.adjust),
      significance = dplyr::case_when(
        .data$p.adjust < 0.001 ~ "***",
        .data$p.adjust < 0.01 ~ "**",
        .data$p.adjust < 0.05 ~ "*",
        TRUE ~ "n.s."
      ),
      .before = 1
    ) |>
    dplyr::group_by(.data$enrichment_direction) |>
    dplyr::arrange(.data$p.adjust, .by_group = TRUE) |>
    dplyr::ungroup()
}


## Run Hallmark GSEA for all five comparisons --------------------------------

hallmark_gsea_analysis_results <- purrr::map(
  seq_len(nrow(hallmark_gsea_comparison_settings)),
  function(comparison_index) {
    current_comparison_name <-
      hallmark_gsea_comparison_settings$comparison_name[[comparison_index]]

    current_comparison_label <-
      hallmark_gsea_comparison_settings$comparison_label[[comparison_index]]

    current_comparison_category <-
      hallmark_gsea_comparison_settings$comparison_category[[comparison_index]]

    current_differential_expression_file <-
      hallmark_gsea_comparison_settings$differential_expression_file[[
        comparison_index
      ]]

    current_differential_expression_table <-
      read_differential_expression_table(
        differential_expression_file = current_differential_expression_file,
        comparison_name = current_comparison_name
      )

    current_ranked_gene_list <- create_ranked_gene_list(
      differential_expression_table = current_differential_expression_table,
      comparison_name = current_comparison_name
    )

    current_hallmark_gsea_result <- run_hallmark_gsea(
      ranked_gene_list = current_ranked_gene_list
    )

    current_hallmark_gsea_result_table <-
      format_hallmark_gsea_result_table(
        hallmark_gsea_result = current_hallmark_gsea_result,
        comparison_name = current_comparison_name,
        comparison_label = current_comparison_label,
        comparison_category = current_comparison_category
      )

    list(
      result_object = current_hallmark_gsea_result,
      result_table = current_hallmark_gsea_result_table
    )
  }
) |>
  stats::setNames(
    hallmark_gsea_comparison_settings$comparison_name
  )

hallmark_gsea_result_objects <- purrr::map(
  hallmark_gsea_analysis_results,
  "result_object"
)

hallmark_gsea_result_tables <- purrr::map(
  hallmark_gsea_analysis_results,
  "result_table"
)

hallmark_gsea_result_table <- dplyr::bind_rows(
  hallmark_gsea_result_tables
)

if (nrow(hallmark_gsea_result_table) == 0) {
  stop(
    "Hallmark GSEA returned no pathways for any comparison.",
    call. = FALSE
  )
}


## Prepare output tables ------------------------------------------------------

hallmark_gsea_output_table <- hallmark_gsea_result_table |>
  dplyr::select(
    comparison_name,
    comparison_label,
    comparison_category,
    enrichment_direction,
    ID,
    Description,
    setSize,
    enrichmentScore,
    NES,
    pvalue,
    p.adjust,
    qvalue,
    rank,
    leading_edge,
    core_enrichment,
    minus_log10_adjusted_p_value,
    significance
  )

isolation_oxygen_hallmark_gsea_table <- hallmark_gsea_output_table |>
  dplyr::filter(
    .data$comparison_category == "Isolation oxygen comparison"
  )

assay_oxygen_hallmark_gsea_table <- hallmark_gsea_output_table |>
  dplyr::filter(
    .data$comparison_category == "Assay oxygen comparison"
  )


## Save Hallmark GSEA results -------------------------------------------------

saveRDS(
  hallmark_gsea_result_objects,
  file = file.path(
    bulk_rnaseq_object_dir,
    "bulk_rnaseq_hallmark_gsea_result_objects.rds"
  )
)

readr::write_csv(
  hallmark_gsea_output_table,
  file = file.path(
    hallmark_gsea_table_dir,
    "bulk_rnaseq_hallmark_gsea_all_results.csv"
  ),
  na = ""
)

readr::write_csv(
  isolation_oxygen_hallmark_gsea_table,
  file = file.path(
    hallmark_gsea_table_dir,
    "bulk_rnaseq_hallmark_gsea_isolation_oxygen_comparisons.csv"
  ),
  na = ""
)

readr::write_csv(
  assay_oxygen_hallmark_gsea_table,
  file = file.path(
    hallmark_gsea_table_dir,
    "bulk_rnaseq_hallmark_gsea_assay_oxygen_comparisons.csv"
  ),
  na = ""
)
