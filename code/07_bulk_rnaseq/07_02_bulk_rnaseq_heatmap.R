# Prepare bulk RNA-seq data for heatmap visualization
#
# This script selects matched Normo-CAF and Hypo-CAF samples from CAF pairs
# 2, 8, and 9, identifies the 2,000 genes with the highest standard deviation
# in log2(TPM + 1), calculates gene-wise Z scores, and saves the processed
# results as an RDS object.


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

bulk_rnaseq_heatmap_results_file <- file.path(
  bulk_rnaseq_object_dir,
  "bulk_rnaseq_heatmap_results.rds"
)

fs::dir_create(bulk_rnaseq_object_dir)


## Define analysis settings --------------------------------------------------

heatmap_pair_order <- c(
  "pair_02",
  "pair_08",
  "pair_09"
)

heatmap_condition_order <- c(
  "normocaf_normoxia",
  "hypocaf_hypoxia"
)

number_of_high_sd_genes <- 2000


## Helper functions ----------------------------------------------------------

validate_required_columns <- function(
    data,
    required_columns,
    data_name
) {
  missing_columns <- setdiff(
    required_columns,
    colnames(data)
  )

  if (length(missing_columns) > 0) {
    stop(
      data_name,
      " is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
}

calculate_gene_wise_z_scores <- function(expression_matrix) {
  z_score_matrix <- scale(
    t(expression_matrix)
  ) %>%
    as.matrix() %>%
    t()

  if (anyNA(z_score_matrix)) {
    stop(
      paste0(
        "NA values were generated during Z-score calculation. ",
        "Check for zero-variance genes."
      ),
      call. = FALSE
    )
  }

  z_score_matrix
}


## Load input data -----------------------------------------------------------

if (!file.exists(bulk_rnaseq_tpm_file)) {
  stop(
    "Bulk RNA-seq TPM file not found: ",
    bulk_rnaseq_tpm_file,
    call. = FALSE
  )
}

if (!file.exists(bulk_rnaseq_sample_metadata_file)) {
  stop(
    "Bulk RNA-seq sample metadata file not found: ",
    bulk_rnaseq_sample_metadata_file,
    call. = FALSE
  )
}

tpm_table <- readr::read_csv(
  bulk_rnaseq_tpm_file,
  show_col_types = FALSE,
  name_repair = "minimal"
)

sample_metadata <- readr::read_csv(
  bulk_rnaseq_sample_metadata_file,
  show_col_types = FALSE
)

caf8_marker_table <- readr::read_csv(
  caf8_signature_genes_table_file,
  show_col_types = FALSE
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

if (anyDuplicated(tpm_table$gene_symbol) > 0) {
  stop(
    "Duplicated gene symbols were detected in the TPM table.",
    call. = FALSE
  )
}

if (anyDuplicated(sample_metadata$sample_name) > 0) {
  stop(
    "Duplicated sample names were detected in sample metadata.",
    call. = FALSE
  )
}


## Select and order heatmap samples ------------------------------------------

heatmap_sample_metadata <- sample_metadata %>%
  dplyr::mutate(
    condition = paste(
      .data$isolation_oxygen_condition,
      .data$assay_oxygen_condition,
      sep = "_"
    )
  ) %>%
  dplyr::filter(
    .data$pair %in% heatmap_pair_order,
    .data$condition %in% heatmap_condition_order
  ) %>%
  dplyr::mutate(
    pair = factor(
      .data$pair,
      levels = heatmap_pair_order
    ),
    condition = factor(
      .data$condition,
      levels = heatmap_condition_order
    )
  ) %>%
  dplyr::arrange(
    .data$pair,
    .data$condition,
    .data$biological_replicate
  )

expected_sample_count <- (
  length(heatmap_pair_order) *
    length(heatmap_condition_order) *
    3
)

if (nrow(heatmap_sample_metadata) != expected_sample_count) {
  stop(
    "Expected ",
    expected_sample_count,
    " heatmap samples, but found ",
    nrow(heatmap_sample_metadata),
    ".",
    call. = FALSE
  )
}

missing_samples <- setdiff(
  heatmap_sample_metadata$sample_name,
  colnames(tpm_table)
)

if (length(missing_samples) > 0) {
  stop(
    "Heatmap samples missing from the TPM table: ",
    paste(missing_samples, collapse = ", "),
    call. = FALSE
  )
}


## Select high-SD genes and calculate Z scores -------------------------------

tpm_matrix <- tpm_table %>%
  dplyr::select(
    .data$gene_symbol,
    dplyr::all_of(heatmap_sample_metadata$sample_name)
  ) %>%
  tibble::column_to_rownames("gene_symbol") %>%
  as.matrix()

storage.mode(tpm_matrix) <- "numeric"

log2_tpm_matrix <- log2(
  tpm_matrix + 1
)

gene_sd <- apply(
  log2_tpm_matrix,
  1,
  stats::sd,
  na.rm = TRUE
)

high_sd_gene_order <- order(
  gene_sd,
  decreasing = TRUE,
  na.last = NA
)

if (length(high_sd_gene_order) < number_of_high_sd_genes) {
  stop(
    "Fewer than ",
    number_of_high_sd_genes,
    " genes have valid standard deviations.",
    call. = FALSE
  )
}

selected_gene_indices <- high_sd_gene_order[
  seq_len(number_of_high_sd_genes)
]

selected_log2_tpm_matrix <- log2_tpm_matrix[
  selected_gene_indices,
  ,
  drop = FALSE
]

gene_z_score_matrix <- calculate_gene_wise_z_scores(
  selected_log2_tpm_matrix
)

selected_gene_table <- tibble::tibble(
  gene_symbol = rownames(selected_log2_tpm_matrix),
  sd_log2_tpm = gene_sd[selected_gene_indices],
  sd_rank = seq_len(number_of_high_sd_genes)
)

## Define genes to annotate ---------------------------------------------------

validate_required_columns(
  caf8_marker_table,
  "gene",
  "CAF-8 marker table"
)

all_caf8_markers <- unique(
  caf8_marker_table$gene
)

all_caf8_markers_in_high_sd <- intersect(
  all_caf8_markers,
  rownames(gene_z_score_matrix)
)

representative_caf8_markers <- c(
  "TOP2A",
  "RRM2",
  "MYBL2",
  "ANLN",
  "BIRC5",
  "TYMS",
  "LMNB1",
  "MCM5",
  "MCM7",
  "SMC4"
)

senescence_related_genes <- c(
  "CDKN1A",
  "CDKN2A"
)

# Two CAF-8 marker sets are prepared for downstream heatmap visualization.
# The representative-marker set is used for the manuscript figure.
# The all-marker set is retained as an exploratory output to show all CAF-8
# markers included among the high-SD genes.

all_caf8_and_senescence_labelled_genes <- c(
  all_caf8_markers_in_high_sd,
  senescence_related_genes
)

representative_caf8_and_senescence_labelled_genes <- c(
  representative_caf8_markers,
  senescence_related_genes
)

labelled_gene_sets <- list(
  all_caf8_and_senescence =
    all_caf8_and_senescence_labelled_genes,
  representative_caf8_and_senescence =
    representative_caf8_and_senescence_labelled_genes
)

for (label_set_name in names(labelled_gene_sets)) {
  genes_not_in_heatmap <- setdiff(
    labelled_gene_sets[[label_set_name]],
    rownames(gene_z_score_matrix)
  )

  if (length(genes_not_in_heatmap) > 0) {
    stop(
      "Genes in label set '",
      label_set_name,
      "' were not found among the top ",
      number_of_high_sd_genes,
      " high-SD genes: ",
      paste(genes_not_in_heatmap, collapse = ", "),
      call. = FALSE
    )
  }
}


## Save heatmap results ------------------------------------------------------

bulk_rnaseq_heatmap_results <- list(
  gene_z_score_matrix = gene_z_score_matrix,
  selected_log2_tpm_matrix = selected_log2_tpm_matrix,
  selected_gene_table = selected_gene_table,
  ordered_sample_metadata = heatmap_sample_metadata,
  labelled_gene_sets = labelled_gene_sets,
  marker_gene_sets = list(
    all_caf8_markers = all_caf8_markers,
    all_caf8_markers_in_high_sd = all_caf8_markers_in_high_sd,
    representative_caf8_markers = representative_caf8_markers,
    senescence_related_genes = senescence_related_genes
  ),
  analysis_settings = list(
    heatmap_pair_order = heatmap_pair_order,
    heatmap_condition_order = heatmap_condition_order,
    number_of_high_sd_genes = number_of_high_sd_genes
  )
)

saveRDS(
  bulk_rnaseq_heatmap_results,
  file = bulk_rnaseq_heatmap_results_file
)


## Print output summary ------------------------------------------------------

message("Bulk RNA-seq heatmap data preparation completed.")
message("Saved heatmap results: ", bulk_rnaseq_heatmap_results_file)
