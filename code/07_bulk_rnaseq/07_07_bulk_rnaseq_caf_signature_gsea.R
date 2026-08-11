## CAF signature gene set enrichment analysis of bulk RNA-seq results --------

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


## Input and output paths -----------------------------------------------------

differential_expression_table_dir <- file.path(
  bulk_rnaseq_table_dir,
  "differential_expression"
)

caf_signature_file <- file.path(
  gene_set_input_dir,
  "elyada_supplementary_table_s22_human_orthologs.csv"
)

gene_set_analysis_table_dir <- file.path(
  bulk_rnaseq_table_dir,
  "gene_set_analysis"
)

caf_signature_gsea_table_dir <- file.path(
  gene_set_analysis_table_dir,
  "caf_signature_gsea"
)

gene_set_analysis_object_dir <- file.path(
  bulk_rnaseq_object_dir,
  "gene_set_analysis"
)

caf_signature_gsea_object_dir <- file.path(
  gene_set_analysis_object_dir,
  "caf_signature_gsea"
)

fs::dir_create(caf_signature_gsea_table_dir)
fs::dir_create(caf_signature_gsea_object_dir)

required_input_files <- c(
  caf_signature_file,
  caf8_signature_genes_table_file
)

missing_input_files <- required_input_files[!file.exists(required_input_files)]

if (length(missing_input_files) > 0) {
  stop(
    "Required input files were not found: ",
    paste(missing_input_files, collapse = ", "),
    call. = FALSE
  )
}


## Define the three isolation-oxygen comparisons -----------------------------

caf_signature_gsea_comparison_settings <- tibble::tribble(
  ~comparison_name, ~comparison_label,
  "comparison_1_hypocaf_hypoxia_vs_normocaf_normoxia", "Original",
  "comparison_2_hypocaf_normoxia_vs_normocaf_normoxia", "Under21",
  "comparison_3_hypocaf_hypoxia_vs_normocaf_hypoxia", "Under1"
) |>
  dplyr::mutate(
    differential_expression_file = file.path(
      differential_expression_table_dir,
      paste0(.data$comparison_name, ".csv")
    )
  )

missing_differential_expression_files <-
  caf_signature_gsea_comparison_settings$differential_expression_file[
    !file.exists(
      caf_signature_gsea_comparison_settings$differential_expression_file
    )
  ]

if (length(missing_differential_expression_files) > 0) {
  stop(
    "Required differential-expression files were not found: ",
    paste(missing_differential_expression_files, collapse = ", "),
    call. = FALSE
  )
}


## Read the Elyada CAF subtype signatures ------------------------------------

elyada_signature_table <- readr::read_csv(
  caf_signature_file,
  show_col_types = FALSE
)

required_elyada_columns <- c(
  "Subtype",
  "Ortholog_Human",
  "ENTREZID"
)

missing_elyada_columns <- setdiff(
  required_elyada_columns,
  colnames(elyada_signature_table)
)

if (length(missing_elyada_columns) > 0) {
  stop(
    "The Elyada CAF signature table is missing required columns: ",
    paste(missing_elyada_columns, collapse = ", "),
    call. = FALSE
  )
}

elyada_term2gene <- elyada_signature_table |>
  dplyr::transmute(
    term = as.character(.data$Subtype),
    gene = as.character(.data$ENTREZID)
  ) |>
  dplyr::filter(
    !is.na(.data$term),
    .data$term != "",
    !is.na(.data$gene),
    .data$gene != ""
  )

expected_elyada_signatures <- c("myCAF", "iCAF", "apCAF")
missing_elyada_signatures <- setdiff(
  expected_elyada_signatures,
  unique(elyada_term2gene$term)
)

if (length(missing_elyada_signatures) > 0) {
  stop(
    "The following expected Elyada CAF signatures were not found: ",
    paste(missing_elyada_signatures, collapse = ", "),
    call. = FALSE
  )
}


## Define the CAF-8 marker signature -----------------------------------------

caf8_marker_table <- readr::read_csv(
  caf8_signature_genes_table_file,
  show_col_types = FALSE
) %>% 
  dplyr::mutate(
    entrez_id = as.character(.data$entrez_id)
  )

required_marker_columns <- c(
  "cluster",
  "gene",
  "entrez_id",
  "avg_log2FC",
  "pct.1",
  "p_val_adj"
)

missing_marker_columns <- setdiff(
  required_marker_columns,
  colnames(caf8_marker_table)
)

if (length(missing_marker_columns) > 0) {
  stop(
    "The CAF-8 marker table is missing required columns: ",
    paste(missing_marker_columns, collapse = ", "),
    call. = FALSE
  )
}

caf8_marker_table <- caf8_marker_table |>
  dplyr::filter(
    !is.na(.data$entrez_id),
    .data$entrez_id != ""
  )

if (nrow(caf8_marker_table) < 10) {
  stop(
    "The CAF-8 marker signature contains fewer than 10 mapped genes.",
    call. = FALSE
  )
}

if (anyDuplicated(caf8_marker_table$entrez_id)) {
  stop(
    "Duplicated Entrez IDs were found in the CAF-8 marker signature.",
    call. = FALSE
  )
}

caf8_term2gene <- caf8_marker_table |>
  dplyr::transmute(
    term = "CAF8",
    gene = .data$entrez_id
  )

caf_signature_term2gene <- dplyr::bind_rows(
  elyada_term2gene,
  caf8_term2gene
)

if (nrow(caf_signature_term2gene) == 0) {
  stop(
    "No CAF signature genes were available for GSEA.",
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
    show_col_types = FALSE
  )

  required_columns <- c(
    "gene_symbol",
    "entrez_id",
    "logFC"
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
    dplyr::mutate(
      entrez_id = as.character(.data$entrez_id)
    ) |>
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
      ": ",
      paste(duplicated_entrez_ids$entrez_id, collapse = ", "),
      call. = FALSE
    )
  }

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


## Run CAF signature GSEA with the original analysis settings ----------------

run_caf_signature_gsea <- function(ranked_gene_list) {
  # The original analysis reset the same seed before every GSEA run.
  set.seed(1234)

  clusterProfiler::GSEA(
    geneList = ranked_gene_list,
    minGSSize = 10,
    TERM2GENE = caf_signature_term2gene,
    pvalueCutoff = 1,
    verbose = FALSE,
    eps = 0
  )
}


## Format one CAF signature GSEA result table --------------------------------

format_caf_signature_gsea_result_table <- function(
    caf_signature_gsea_result,
    comparison_name,
    comparison_label
) {
  result_table <- caf_signature_gsea_result@result |>
    tibble::as_tibble()

  if (nrow(result_table) == 0) {
    warning(
      "CAF signature GSEA returned no gene sets for ",
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
      significance = dplyr::case_when(
        .data$p.adjust < 0.001 ~ "***",
        .data$p.adjust < 0.01 ~ "**",
        .data$p.adjust < 0.05 ~ "*",
        TRUE ~ "n.s."
      ),
      .before = 1
    )
}


## Run CAF signature GSEA for all three comparisons --------------------------

caf_signature_gsea_analysis_results <- purrr::map(
  seq_len(nrow(caf_signature_gsea_comparison_settings)),
  function(comparison_index) {
    current_comparison_name <-
      caf_signature_gsea_comparison_settings$comparison_name[[
        comparison_index
      ]]

    current_comparison_label <-
      caf_signature_gsea_comparison_settings$comparison_label[[
        comparison_index
      ]]

    current_differential_expression_file <-
      caf_signature_gsea_comparison_settings$differential_expression_file[[
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

    current_caf_signature_gsea_result <- run_caf_signature_gsea(
      ranked_gene_list = current_ranked_gene_list
    )

    current_caf_signature_gsea_result_table <-
      format_caf_signature_gsea_result_table(
        caf_signature_gsea_result = current_caf_signature_gsea_result,
        comparison_name = current_comparison_name,
        comparison_label = current_comparison_label
      )

    list(
      result_object = current_caf_signature_gsea_result,
      result_table = current_caf_signature_gsea_result_table
    )
  }
) |>
  stats::setNames(
    caf_signature_gsea_comparison_settings$comparison_name
  )

caf_signature_gsea_result_objects <- purrr::map(
  caf_signature_gsea_analysis_results,
  "result_object"
)

caf_signature_gsea_result_table <- caf_signature_gsea_analysis_results |>
  purrr::map("result_table") |>
  dplyr::bind_rows()

if (nrow(caf_signature_gsea_result_table) == 0) {
  stop(
    "CAF signature GSEA returned no gene sets for any comparison.",
    call. = FALSE
  )
}


## Prepare output tables ------------------------------------------------------

caf_signature_order <- c("CAF8", "myCAF", "iCAF", "apCAF")
comparison_order <- c("Original", "Under21", "Under1")

caf_signature_gsea_output_table <- caf_signature_gsea_result_table |>
  dplyr::mutate(
    Description = factor(
      .data$Description,
      levels = caf_signature_order
    ),
    comparison_label = factor(
      .data$comparison_label,
      levels = comparison_order
    )
  ) |>
  dplyr::arrange(
    .data$comparison_label,
    .data$Description
  ) |>
  dplyr::select(
    comparison_name,
    comparison_label,
    ID,
    Description,
    enrichment_direction,
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
  ) |>
  dplyr::mutate(
    comparison_label = as.character(.data$comparison_label),
    Description = as.character(.data$Description)
  )


## Save CAF signature GSEA results -------------------------------------------

saveRDS(
  caf_signature_gsea_result_objects,
  file = file.path(
    caf_signature_gsea_object_dir,
    "bulk_rnaseq_caf_signature_gsea_result_objects.rds"
  )
)

readr::write_csv(
  caf_signature_gsea_output_table,
  file = file.path(
    caf_signature_gsea_table_dir,
    "bulk_rnaseq_caf_signature_gsea_results.csv"
  ),
  na = ""
)
