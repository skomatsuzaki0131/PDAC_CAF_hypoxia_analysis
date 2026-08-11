## GO gene set enrichment analysis of bulk RNA-seq results -------------------


## Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(dplyr)
  library(fs)
  library(org.Hs.eg.db)
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

go_gsea_table_dir <- file.path(
  gene_set_analysis_table_dir,
  "go_gsea"
)

fs::dir_create(go_gsea_table_dir)


## Define comparisons used in the original GO GSEA --------------------------

go_gsea_comparison_settings <- tibble::tribble(
  ~comparison_name, ~comparison_label,
  "comparison_1_hypocaf_hypoxia_vs_normocaf_normoxia",  "Original",
  "comparison_2_hypocaf_normoxia_vs_normocaf_normoxia", "Under21",
  "comparison_3_hypocaf_hypoxia_vs_normocaf_hypoxia",   "Under1"
) |>
  dplyr::mutate(
    differential_expression_file = file.path(
      differential_expression_table_dir,
      paste0(.data$comparison_name, ".csv")
    )
  )

missing_input_files <- go_gsea_comparison_settings$differential_expression_file[
  !file.exists(go_gsea_comparison_settings$differential_expression_file)
]

if (length(missing_input_files) > 0) {
  stop(
    "Required differential-expression files were not found: ",
    paste(missing_input_files, collapse = ", "),
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


## Create the ranked Entrez-ID vector used by gseGO --------------------------

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
      ". Resolve duplicated mappings before running GO GSEA: ",
      paste(duplicated_entrez_ids$entrez_id, collapse = ", "),
      call. = FALSE
    )
  }

  # The original differential-expression tables are ordered by decreasing
  # logFC. Sorting explicitly here preserves the ranked-list requirement of
  # gseGO and makes that dependency clear in the public code.
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


## Run GO GSEA with the original analysis settings ---------------------------

run_go_gsea <- function(ranked_gene_list) {
  set.seed(1234)

  clusterProfiler::gseGO(
    geneList = ranked_gene_list,
    OrgDb = org.Hs.eg.db,
    ont = "ALL",
    minGSSize = 50,
    pvalueCutoff = 0.05,
    verbose = FALSE,
    eps = 0
  )
}


## Run GO GSEA for the three CAF-origin comparisons --------------------------

go_gsea_results <- purrr::map(
  seq_len(nrow(go_gsea_comparison_settings)),
  function(comparison_index) {
    current_comparison_name <-
      go_gsea_comparison_settings$comparison_name[[comparison_index]]

    current_differential_expression_file <-
      go_gsea_comparison_settings$differential_expression_file[[comparison_index]]

    current_differential_expression_table <-
      read_differential_expression_table(
        differential_expression_file = current_differential_expression_file,
        comparison_name = current_comparison_name
      )

    current_ranked_gene_list <- create_ranked_gene_list(
      differential_expression_table = current_differential_expression_table,
      comparison_name = current_comparison_name
    )

    current_go_gsea_result <- run_go_gsea(
      ranked_gene_list = current_ranked_gene_list
    )

    current_simplified_go_gsea_result <- clusterProfiler::simplify(
      current_go_gsea_result
    )

    list(
      ranked_gene_list = current_ranked_gene_list,
      go_gsea_result = current_go_gsea_result,
      simplified_go_gsea_result = current_simplified_go_gsea_result
    )
  }
)

names(go_gsea_results) <- go_gsea_comparison_settings$comparison_name


## Convert simplified GO GSEA results to one table ---------------------------

format_go_gsea_result_table <- function(
    simplified_go_gsea_result,
    comparison_name,
    comparison_label
) {
  result_table <- simplified_go_gsea_result@result |>
    tibble::as_tibble()

  if (nrow(result_table) == 0) {
    warning(
      "GO GSEA returned no significant terms for ",
      comparison_name,
      "."
    )

    return(tibble::tibble())
  }

  result_table |>
    dplyr::mutate(
      comparison_name = comparison_name,
      comparison_label = comparison_label,
      enrichment_direction = dplyr::if_else(
        .data$NES > 0,
        "Upregulated in Hypo-CAF",
        "Upregulated in Normo-CAF"
      ),
      minus_log10_adjusted_p_value = -log10(.data$p.adjust),
      .before = 1
    ) |>
    dplyr::group_by(.data$enrichment_direction) |>
    dplyr::arrange(.data$p.adjust, .by_group = TRUE) |>
    dplyr::ungroup()
}

go_gsea_result_table <- purrr::map_dfr(
  seq_len(nrow(go_gsea_comparison_settings)),
  function(comparison_index) {
    current_comparison_name <-
      go_gsea_comparison_settings$comparison_name[[comparison_index]]

    current_comparison_label <-
      go_gsea_comparison_settings$comparison_label[[comparison_index]]

    format_go_gsea_result_table(
      simplified_go_gsea_result = go_gsea_results[[
        current_comparison_name
      ]]$simplified_go_gsea_result,
      comparison_name = current_comparison_name,
      comparison_label = current_comparison_label
    )
  }
)

if (nrow(go_gsea_result_table) == 0) {
  stop(
    "GO GSEA returned no significant terms for any comparison.",
    call. = FALSE
  )
}


## Identify BP terms significant in the same direction in all comparisons ----

find_shared_biological_process_terms <- function(
    go_gsea_table,
    direction
) {
  terms_by_comparison <- go_gsea_table |>
    dplyr::filter(
      .data$ONTOLOGY == "BP",
      .data$enrichment_direction == direction,
      .data$p.adjust < 0.05
    ) |>
    dplyr::group_by(.data$comparison_name) |>
    dplyr::summarise(
      terms = list(unique(.data$Description)),
      .groups = "drop"
    ) |>
    dplyr::right_join(
      go_gsea_comparison_settings |>
        dplyr::select(.data$comparison_name),
      by = "comparison_name"
    ) |>
    dplyr::mutate(
      terms = purrr::map(
        .data$terms,
        function(current_terms) {
          if (is.null(current_terms) || length(current_terms) == 0) {
            character(0)
          } else {
            current_terms
          }
        }
      )
    )

  if (nrow(terms_by_comparison) != nrow(go_gsea_comparison_settings)) {
    stop(
      "Could not evaluate shared GO terms across all comparisons.",
      call. = FALSE
    )
  }

  purrr::reduce(
    terms_by_comparison$terms,
    intersect
  )
}

shared_go_bp_terms <- dplyr::bind_rows(
  tibble::tibble(
    enrichment_direction = "Upregulated in Normo-CAF",
    Description = find_shared_biological_process_terms(
      go_gsea_table = go_gsea_result_table,
      direction = "Upregulated in Normo-CAF"
    )
  ),
  tibble::tibble(
    enrichment_direction = "Upregulated in Hypo-CAF",
    Description = find_shared_biological_process_terms(
      go_gsea_table = go_gsea_result_table,
      direction = "Upregulated in Hypo-CAF"
    )
  )
)


## Save GO GSEA results -------------------------------------------------------

readr::write_csv(
  go_gsea_result_table |>
    dplyr::select(
      .data$comparison_name,
      .data$comparison_label,
      .data$enrichment_direction,
      .data$ONTOLOGY,
      .data$ID,
      .data$Description,
      .data$setSize,
      .data$enrichmentScore,
      .data$NES,
      .data$pvalue,
      .data$p.adjust,
      .data$qvalue,
      .data$rank,
      .data$leading_edge,
      .data$core_enrichment,
      .data$minus_log10_adjusted_p_value
    ),
  file = file.path(
    go_gsea_table_dir,
    "bulk_rnaseq_go_gsea_simplified_results.csv"
  ),
  na = ""
)

readr::write_csv(
  shared_go_bp_terms,
  file = file.path(
    go_gsea_table_dir,
    "bulk_rnaseq_go_gsea_shared_bp_terms.csv"
  ),
  na = ""
)

