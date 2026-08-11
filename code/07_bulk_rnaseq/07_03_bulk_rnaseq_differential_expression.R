## Differential expression analysis of bulk RNA-seq data ---------------------

## Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(dplyr)
  library(edgeR)
  library(org.Hs.eg.db)
  library(readr)
  library(tibble)
  library(fs)
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


## Output directory -----------------------------------------------------------

fs::dir_create(
  bulk_rnaseq_object_dir
)

differential_expression_table_dir <- file.path(
  bulk_rnaseq_table_dir,
  "differential_expression"
)

fs::dir_create(
  differential_expression_table_dir
)


## Read raw count data --------------------------------------------------------

raw_count_table <- readr::read_csv(
  file = file.path(
    bulk_rnaseq_input_dir,
    "expression_data",
    "raw_counts.csv"
  ),
  show_col_types = FALSE
)

gene_symbol_column <- "gene_symbol"
ensembl_id_column <- "ensembl_id"

required_count_columns <- c(
  gene_symbol_column,
  ensembl_id_column
)

missing_count_columns <- setdiff(
  required_count_columns,
  colnames(raw_count_table)
)

if (length(missing_count_columns) > 0) {
  stop(
    "Bulk RNA-seq raw count table is missing required columns: ",
    paste(missing_count_columns, collapse = ", "),
    call. = FALSE
  )
}

gene_symbols <- raw_count_table[[gene_symbol_column]]

stopifnot(
  !anyNA(gene_symbols),
  !anyDuplicated(gene_symbols)
)

raw_count_matrix <- raw_count_table |>
  dplyr::select(-dplyr::all_of(c(gene_symbol_column, ensembl_id_column))) |>
  as.matrix()

if (
  anyNA(raw_count_matrix) ||
  !all(is.finite(raw_count_matrix)) ||
  any(raw_count_matrix < 0) ||
  any(raw_count_matrix != floor(raw_count_matrix))
) {
  stop(
    "The raw count matrix must contain finite, non-negative integers.",
    call. = FALSE
  )
}

storage.mode(raw_count_matrix) <- "integer"
rownames(raw_count_matrix) <- gene_symbols

stopifnot(
  all(is.finite(raw_count_matrix)),
  all(raw_count_matrix >= 0)
)


## Parse sample information from public sample names -------------------------

sample_information <- tibble::tibble(
  sample = colnames(raw_count_matrix)
) |>
  tidyr::extract(
    col = sample,
    into = c(
      "pair",
      "isolation_oxygen_condition",
      "assay_oxygen_condition",
      "biological_replicate"
    ),
    regex = paste0(
      "^pair_([0-9]+)_",
      "(hypocaf|normocaf)_",
      "(hypoxia|normoxia)_",
      "rep([0-9]+)$"
    ),
    remove = FALSE,
    convert = TRUE
  )

if (anyNA(sample_information[, c(
  "pair",
  "isolation_oxygen_condition",
  "assay_oxygen_condition",
  "biological_replicate"
)])) {
  stop(
    "One or more sample names do not follow the expected format: ",
    "pair_XX_{hypocaf|normocaf}_{hypoxia|normoxia}_repN"
  )
}

pairs_used_for_bulk_rnaseq <- c(2, 8, 9)

sample_information <- sample_information |>
  dplyr::filter(pair %in% pairs_used_for_bulk_rnaseq) |>
  dplyr::mutate(
    pair = factor(pair, levels = pairs_used_for_bulk_rnaseq),
    isolation_oxygen_condition = factor(
      isolation_oxygen_condition,
      levels = c("normocaf", "hypocaf")
    ),
    assay_oxygen_condition = factor(
      assay_oxygen_condition,
      levels = c("normoxia", "hypoxia")
    )
  )

expected_sample_information <- tidyr::expand_grid(
  pair = pairs_used_for_bulk_rnaseq,
  isolation_oxygen_condition = c("normocaf", "hypocaf"),
  assay_oxygen_condition = c("normoxia", "hypoxia"),
  biological_replicate = 1:3
)

observed_sample_information <- sample_information |>
  dplyr::transmute(
    pair = as.integer(as.character(pair)),
    isolation_oxygen_condition = as.character(isolation_oxygen_condition),
    assay_oxygen_condition = as.character(assay_oxygen_condition),
    biological_replicate = biological_replicate
  )

stopifnot(
  nrow(sample_information) == 36,
  nrow(dplyr::anti_join(
    expected_sample_information,
    observed_sample_information,
    by = c("pair", "isolation_oxygen_condition", "assay_oxygen_condition", "biological_replicate")
  )) == 0,
  nrow(dplyr::anti_join(
    observed_sample_information,
    expected_sample_information,
    by = c("pair", "isolation_oxygen_condition", "assay_oxygen_condition", "biological_replicate")
  )) == 0
)


## Define the five comparisons used in the original analysis -----------------

comparison_settings <- list(
  comparison_1_hypocaf_hypoxia_vs_normocaf_normoxia = list(
    positive = c(isolation_oxygen_condition = "hypocaf", assay_oxygen_condition = "hypoxia"),
    negative = c(isolation_oxygen_condition = "normocaf", assay_oxygen_condition = "normoxia")
  ),
  comparison_2_hypocaf_normoxia_vs_normocaf_normoxia = list(
    positive = c(isolation_oxygen_condition = "hypocaf", assay_oxygen_condition = "normoxia"),
    negative = c(isolation_oxygen_condition = "normocaf", assay_oxygen_condition = "normoxia")
  ),
  comparison_3_hypocaf_hypoxia_vs_normocaf_hypoxia = list(
    positive = c(isolation_oxygen_condition = "hypocaf", assay_oxygen_condition = "hypoxia"),
    negative = c(isolation_oxygen_condition = "normocaf", assay_oxygen_condition = "hypoxia")
  ),
  comparison_4_hypocaf_hypoxia_vs_hypocaf_normoxia = list(
    positive = c(isolation_oxygen_condition = "hypocaf", assay_oxygen_condition = "hypoxia"),
    negative = c(isolation_oxygen_condition = "hypocaf", assay_oxygen_condition = "normoxia")
  ),
  comparison_5_normocaf_hypoxia_vs_normocaf_normoxia = list(
    positive = c(isolation_oxygen_condition = "normocaf", assay_oxygen_condition = "hypoxia"),
    negative = c(isolation_oxygen_condition = "normocaf", assay_oxygen_condition = "normoxia")
  )
)


## Run one paired edgeR quasi-likelihood analysis -----------------------------

run_differential_expression <- function(
    count_matrix,
    sample_table,
    comparison_setting
) {
  positive_samples <- sample_table |>
    dplyr::filter(
      isolation_oxygen_condition == comparison_setting$positive[["isolation_oxygen_condition"]],
      assay_oxygen_condition == comparison_setting$positive[["assay_oxygen_condition"]]
    ) |>
    dplyr::mutate(comparison_group = "PositiveSide")

  negative_samples <- sample_table |>
    dplyr::filter(
      isolation_oxygen_condition == comparison_setting$negative[["isolation_oxygen_condition"]],
      assay_oxygen_condition == comparison_setting$negative[["assay_oxygen_condition"]]
    ) |>
    dplyr::mutate(comparison_group = "NegativeSide")

  comparison_sample_table <- dplyr::bind_rows(
    negative_samples,
    positive_samples
  ) |>
    dplyr::arrange(pair, comparison_group, biological_replicate) |>
    dplyr::mutate(
      comparison_group = factor(
        comparison_group,
        levels = c("NegativeSide", "PositiveSide")
      )
    )

  stopifnot(
    nrow(comparison_sample_table) == 18,
    all(table(comparison_sample_table$pair) == 6),
    all(table(comparison_sample_table$comparison_group) == 9)
  )

  comparison_count_matrix <- count_matrix[
    , comparison_sample_table$sample,
    drop = FALSE
  ]

  stopifnot(
    identical(
      colnames(comparison_count_matrix),
      comparison_sample_table$sample
    )
  )

  edge_object <- edgeR::DGEList(
    counts = comparison_count_matrix,
    samples = comparison_sample_table
  )

  design_matrix <- stats::model.matrix(
    ~ pair + comparison_group,
    data = comparison_sample_table
  )

  genes_to_keep <- edgeR::filterByExpr(
    edge_object,
    design = design_matrix
  )

  edge_object <- edge_object[
    genes_to_keep,
    ,
    keep.lib.sizes = FALSE
  ]

  edge_object <- edgeR::calcNormFactors(
    edge_object,
    method = "TMM"
  )

  edge_object <- edgeR::estimateDisp(
    edge_object,
    design = design_matrix,
    robust = TRUE
  )

  quasi_likelihood_fit <- edgeR::glmQLFit(
    edge_object,
    design = design_matrix,
    robust = TRUE
  )

  quasi_likelihood_test <- edgeR::glmQLFTest(
    quasi_likelihood_fit,
    coef = "comparison_groupPositiveSide"
  )

  differential_expression_table <- edgeR::topTags(
    quasi_likelihood_test,
    n = Inf
  )$table |>
    tibble::rownames_to_column(var = "gene_symbol") |>
    dplyr::arrange(dplyr::desc(logFC))

  list(
    differential_expression_table = differential_expression_table,
    sample_information = comparison_sample_table,
    design_matrix = design_matrix,
    edge_object = edge_object,
    quasi_likelihood_fit = quasi_likelihood_fit,
    quasi_likelihood_test = quasi_likelihood_test
  )
}


differential_expression_results <- lapply(
  comparison_settings,
  function(comparison_setting) {
    run_differential_expression(
      count_matrix = raw_count_matrix,
      sample_table = sample_information,
      comparison_setting = comparison_setting
    )
  }
)


## Add Entrez identifiers -----------------------------------------------------

# org.Hs.eg.db maps TEC to two Entrez identifiers in the annotation version
# used for the original analysis. ENTREZID 100124696 is removed explicitly so
# that each gene symbol remains represented by a single row.

add_entrez_identifiers <- function(differential_expression_table) {
  differential_expression_table <- differential_expression_table |>
    dplyr::select(
      -dplyr::any_of(c("EntrezID", "entrez_id"))
    )
  
  gene_symbols <- differential_expression_table$gene_symbol

  symbol_to_entrez <- AnnotationDbi::select(
    x = org.Hs.eg.db,
    keys = gene_symbols,
    columns = "ENTREZID",
    keytype = "SYMBOL"
  ) |>
    tibble::as_tibble() |>
    dplyr::filter(
      !(SYMBOL == "TEC" & ENTREZID == "100124696")
    )

  duplicated_symbols <- symbol_to_entrez |>
    dplyr::filter(!is.na(ENTREZID)) |>
    dplyr::count(SYMBOL, name = "number_of_entrez_ids") |>
    dplyr::filter(number_of_entrez_ids > 1)

  if (nrow(duplicated_symbols) > 0) {
    stop(
      "Multiple Entrez IDs remain for the following gene symbols after ",
      "removing TEC/100124696: ",
      paste(duplicated_symbols$SYMBOL, collapse = ", ")
    )
  }

  symbol_to_entrez <- symbol_to_entrez |>
    dplyr::group_by(SYMBOL) |>
    dplyr::summarise(
      ENTREZID = dplyr::first(ENTREZID[!is.na(ENTREZID)], default = NA_character_),
      .groups = "drop"
    ) |>
    dplyr::rename(
      gene_symbol = SYMBOL,
      entrez_id = ENTREZID
    )

  differential_expression_table |>
    dplyr::left_join(
      symbol_to_entrez,
      by = "gene_symbol"
    ) |>
    dplyr::relocate(entrez_id, .after = gene_symbol)
}

for (comparison_name in names(differential_expression_results)) {
  differential_expression_results[[comparison_name]]$differential_expression_table <-
    add_entrez_identifiers(
      differential_expression_results[[comparison_name]]$differential_expression_table
    )
}

## Prepare enrichment analysis input tables ----------------------------------

ora_logfc_threshold <- 1
ora_fdr_threshold <- 0.05

ranked_gene_lists <- purrr::map_dfr(
  names(differential_expression_results),
  function(current_comparison_name) {
    differential_expression_results[[
      current_comparison_name
    ]]$differential_expression_table |>
      dplyr::filter(
        !is.na(.data$entrez_id),
        .data$entrez_id != "",
        is.finite(.data$logFC)
      ) |>
      dplyr::transmute(
        comparison_name = .env$current_comparison_name,
        gene_symbol = .data$gene_symbol,
        entrez_id = .data$entrez_id,
        ranking_metric = .data$logFC
      ) |>
      dplyr::arrange(
        dplyr::desc(.data$ranking_metric)
      )
  }
)

significant_degs_for_ora <- purrr::map_dfr(
  names(differential_expression_results),
  function(current_comparison_name) {
    differential_expression_results[[
      current_comparison_name
    ]]$differential_expression_table |>
      dplyr::filter(
        !is.na(.data$entrez_id),
        .data$entrez_id != "",
        abs(.data$logFC) >= ora_logfc_threshold,
        .data$FDR < ora_fdr_threshold
      ) |>
      dplyr::mutate(
        regulation = dplyr::case_when(
          .data$logFC >= ora_logfc_threshold ~ "upregulated",
          .data$logFC <= -ora_logfc_threshold ~ "downregulated"
        )
      ) |>
      dplyr::transmute(
        comparison_name = .env$current_comparison_name,
        gene_symbol = .data$gene_symbol,
        entrez_id = .data$entrez_id,
        logFC = .data$logFC,
        FDR = .data$FDR,
        regulation = .data$regulation
      )
  }
)

ora_background_genes <- purrr::map_dfr(
  names(differential_expression_results),
  function(current_comparison_name) {
    differential_expression_results[[
      current_comparison_name
    ]]$differential_expression_table |>
      dplyr::filter(
        !is.na(.data$entrez_id),
        .data$entrez_id != ""
      ) |>
      dplyr::transmute(
        comparison_name = .env$current_comparison_name,
        gene_symbol = .data$gene_symbol,
        entrez_id = as.character(.data$entrez_id)
      ) |>
      dplyr::distinct(
        .data$comparison_name,
        .data$entrez_id,
        .keep_all = TRUE
      )
  }
)

## Save results ---------------------------------------------------------------

for (comparison_name in names(differential_expression_results)) {
  readr::write_csv(
    differential_expression_results[[comparison_name]]$differential_expression_table,
    file = file.path(
      differential_expression_table_dir,
      paste0(comparison_name, ".csv")
    ),
    na = ""
  )
}

saveRDS(
  object = differential_expression_results,
  file = file.path(
    bulk_rnaseq_object_dir,
    "bulk_rnaseq_differential_expression_results.rds"
  )
)

ranked_gene_lists_file <- file.path(
  bulk_rnaseq_table_dir,
  "differential_expression",
  "bulk_rnaseq_ranked_gene_lists.csv"
)

significant_degs_for_ora_file <- file.path(
  bulk_rnaseq_table_dir,
  "differential_expression",
  "bulk_rnaseq_significant_degs_for_ora.csv"
)

ora_background_genes_file <- file.path(
  bulk_rnaseq_table_dir,
  "differential_expression",
  "bulk_rnaseq_ora_background_genes.csv"
)

readr::write_csv(
  ranked_gene_lists,
  ranked_gene_lists_file
)

readr::write_csv(
  significant_degs_for_ora,
  significant_degs_for_ora_file
)

readr::write_csv(
  ora_background_genes,
  ora_background_genes_file
)
