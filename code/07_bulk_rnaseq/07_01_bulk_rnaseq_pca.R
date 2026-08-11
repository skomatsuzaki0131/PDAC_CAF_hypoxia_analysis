# Perform principal component analysis of bulk RNA-seq samples
#
# This script performs principal component analysis (PCA) from the TPM
# expression matrix and saves the resulting PCA object and sample scores for
# downstream visualization. Genes with TPM < 1 in any sample or zero variance
# across samples are excluded before PCA.


## Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
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
  "00_02_paths.R"
))


## Define input and output paths ---------------------------------------------

bulk_rnaseq_sample_metadata_file <- file.path(
  bulk_rnaseq_input_dir,
  "metadata",
  "sample_metadata.csv"
)

bulk_rnaseq_pca_results_file <- file.path(
  bulk_rnaseq_object_dir,
  "bulk_rnaseq_pca_results.rds"
)

fs::dir_create(bulk_rnaseq_object_dir)


## Define PCA settings -------------------------------------------------------

minimum_tpm <- 1

expected_pairs <- c("pair_02", "pair_08", "pair_09")

condition_display_order <- c(
  "hypocaf_normoxia",
  "hypocaf_hypoxia",
  "normocaf_normoxia",
  "normocaf_hypoxia"
)


## Helper functions ----------------------------------------------------------

read_required_csv <- function(file, description) {
  if (!file.exists(file)) {
    stop(description, " not found: ", file, call. = FALSE)
  }

  readr::read_csv(
    file,
    show_col_types = FALSE,
    name_repair = "minimal"
  )
}

validate_required_columns <- function(data, required_columns, description) {
  missing_columns <- setdiff(required_columns, colnames(data))

  if (length(missing_columns) > 0) {
    stop(
      description,
      " is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
}


## Load input data -----------------------------------------------------------

tpm_table <- read_required_csv(
  bulk_rnaseq_tpm_file,
  "Bulk RNA-seq TPM file"
)

sample_metadata <- read_required_csv(
  bulk_rnaseq_sample_metadata_file,
  "Bulk RNA-seq sample metadata file"
)


## Validate input data -------------------------------------------------------

validate_required_columns(
  tpm_table,
  c("gene_symbol", "ensembl_id"),
  "Bulk RNA-seq TPM table"
)

validate_required_columns(
  sample_metadata,
  c(
    "sample_name",
    "pair",
    "isolation_oxygen_condition",
    "assay_oxygen_condition",
    "biological_replicate"
  ),
  "Bulk RNA-seq sample metadata"
)

if (anyDuplicated(sample_metadata$sample_name) > 0) {
  stop(
    "Duplicated sample names were detected in sample metadata.",
    call. = FALSE
  )
}

expression_sample_names <- setdiff(
  colnames(tpm_table),
  c("gene_symbol", "ensembl_id")
)

missing_expression_samples <- setdiff(
  sample_metadata$sample_name,
  expression_sample_names
)

unexpected_expression_samples <- setdiff(
  expression_sample_names,
  sample_metadata$sample_name
)

if (length(missing_expression_samples) > 0) {
  stop(
    "Samples in metadata but not in TPM table: ",
    paste(missing_expression_samples, collapse = ", "),
    call. = FALSE
  )
}

if (length(unexpected_expression_samples) > 0) {
  stop(
    "Samples in TPM table but not in metadata: ",
    paste(unexpected_expression_samples, collapse = ", "),
    call. = FALSE
  )
}

if (anyDuplicated(tpm_table$ensembl_id) > 0) {
  stop(
    "Duplicated Ensembl IDs were detected in the TPM table.",
    call. = FALSE
  )
}

if (!setequal(unique(sample_metadata$pair), expected_pairs)) {
  stop(
    "Unexpected CAF pairs in sample metadata. Expected: ",
    paste(expected_pairs, collapse = ", "),
    call. = FALSE
  )
}

sample_metadata <- sample_metadata %>%
  dplyr::mutate(
    condition = paste(
      .data$isolation_oxygen_condition,
      .data$assay_oxygen_condition,
      sep = "_"
    ),
    condition = factor(
      .data$condition,
      levels = condition_display_order
    ),
    pair_number = stringr::str_remove(.data$pair, "^pair_0*")
  )

if (any(is.na(sample_metadata$condition))) {
  stop(
    paste0(
      "Unexpected isolation_oxygen_condition and ",
      "assay_oxygen_condition combination in sample metadata."
    ),
    call. = FALSE
  )
}

sample_group_counts <- sample_metadata %>%
  dplyr::count(
    .data$pair,
    .data$condition,
    name = "n_samples"
  )

if (
  nrow(sample_group_counts) !=
    length(expected_pairs) * length(condition_display_order) ||
    any(sample_group_counts$n_samples != 3)
) {
  stop(
    paste0(
      "Expected three biological replicates for each of three CAF pairs ",
      "and four experimental conditions."
    ),
    call. = FALSE
  )
}


## Prepare expression matrix -------------------------------------------------

sample_metadata <- sample_metadata %>%
  dplyr::arrange(
    factor(.data$pair, levels = expected_pairs),
    factor(.data$condition, levels = condition_display_order),
    .data$biological_replicate
  )

expression_matrix <- tpm_table %>%
  dplyr::select(dplyr::all_of(sample_metadata$sample_name)) %>%
  as.matrix()

storage.mode(expression_matrix) <- "double"
rownames(expression_matrix) <- tpm_table$ensembl_id

if (anyNA(expression_matrix)) {
  stop(
    "Missing or non-numeric TPM values were detected.",
    call. = FALSE
  )
}

number_of_input_genes <- nrow(expression_matrix)

expressed_in_all_samples <- apply(
  expression_matrix,
  1,
  function(x) all(x >= minimum_tpm)
)

expression_matrix_filtered <- expression_matrix[
  expressed_in_all_samples,
  ,
  drop = FALSE
]

nonzero_variance <- apply(
  expression_matrix_filtered,
  1,
  stats::var
) > 0

expression_matrix_filtered <- expression_matrix_filtered[
  nonzero_variance,
  ,
  drop = FALSE
]

if (nrow(expression_matrix_filtered) < 3) {
  stop(
    "Fewer than three genes remained after PCA filtering.",
    call. = FALSE
  )
}


## Perform principal component analysis -------------------------------------

pca_result <- stats::prcomp(
  t(expression_matrix_filtered),
  center = TRUE,
  scale. = TRUE
)

variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)

pca_scores <- pca_result$x %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sample_name") %>%
  dplyr::left_join(
    sample_metadata,
    by = "sample_name"
  )

if (anyNA(pca_scores$condition)) {
  stop(
    "Failed to match PCA coordinates to sample metadata.",
    call. = FALSE
  )
}


## Save PCA results ----------------------------------------------------------

bulk_rnaseq_pca_results <- list(
  pca = pca_result,
  scores = pca_scores,
  variance_explained = variance_explained,
  parameters = list(
    minimum_tpm = minimum_tpm,
    number_of_input_samples = ncol(expression_matrix),
    number_of_input_genes = number_of_input_genes,
    number_of_retained_genes = nrow(expression_matrix_filtered)
  )
)

saveRDS(
  bulk_rnaseq_pca_results,
  file = bulk_rnaseq_pca_results_file
)


## Print analysis summary ----------------------------------------------------

message("Bulk RNA-seq PCA completed.")
message("Input samples: ", ncol(expression_matrix))
message("Input genes: ", number_of_input_genes)
message(
  "Genes retained after TPM and variance filtering: ",
  nrow(expression_matrix_filtered)
)
message(
  "Variance explained by PC1, PC2, and PC3: ",
  paste0(
    round(variance_explained[1:3] * 100, 1),
    "%",
    collapse = ", "
  )
)
message("PCA results saved to: ", bulk_rnaseq_pca_results_file)
