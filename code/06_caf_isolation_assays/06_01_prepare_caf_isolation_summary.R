# Prepare CAF isolation summary table
#
# This script validates the CAF isolation metadata and generates a pair-level
# summary table for downstream analyses, figure generation, and supplementary
# tables.


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

source(file.path(project_dir, "code", "00_config", "00_02_paths.R"))


## Define output paths -------------------------------------------------------

caf_isolation_table_dir <- file.path(
  caf_isolation_output_dir,
  "tables"
)

caf_isolation_summary_table_file <- file.path(
  caf_isolation_table_dir,
  "caf_isolation_summary_table.csv"
)

fs::dir_create(caf_isolation_table_dir)


## Helper functions ----------------------------------------------------------

read_required_csv <- function(file, description) {
  if (!file.exists(file)) {
    stop(description, " not found: ", file, call. = FALSE)
  }

  readr::read_csv(
    file,
    show_col_types = FALSE
  )
}

check_required_columns <- function(data, required_columns, table_name) {
  missing_columns <- setdiff(required_columns, colnames(data))

  if (length(missing_columns) > 0) {
    stop(
      "Missing required columns in ",
      table_name,
      ": ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(NULL)
}

check_allowed_values <- function(x, allowed_values, column_name) {
  unexpected_values <- setdiff(unique(x), allowed_values)

  if (length(unexpected_values) > 0) {
    stop(
      "Unexpected values in ",
      column_name,
      ": ",
      paste(unexpected_values, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(NULL)
}

check_unique_within_pair <- function(data, column, column_name) {
  inconsistent_pairs <- data %>%
    dplyr::group_by(.data$caf_pair_no) %>%
    dplyr::summarise(
      n_unique = dplyr::n_distinct({{ column }}),
      .groups = "drop"
    ) %>%
    dplyr::filter(.data$n_unique != 1)

  if (nrow(inconsistent_pairs) > 0) {
    stop(
      column_name,
      " is not consistent within CAF pair(s): ",
      paste(inconsistent_pairs$caf_pair_no, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(NULL)
}


## Read CAF isolation metadata ----------------------------------------------

caf_isolation_metadata <- read_required_csv(
  caf_isolation_metadata_file,
  "CAF isolation metadata file"
)


## Validate CAF isolation metadata ------------------------------------------

required_columns <- c(
  "caf_pair_no",
  "neoadjuvant_therapy",
  "tumor_location",
  "isolation_oxygen_percent",
  "days_to_outgrowth",
  "used_for_rna_seq",
  "used_for_proliferation_assay",
  "used_for_morphology_assay"
)

check_required_columns(
  caf_isolation_metadata,
  required_columns,
  "CAF isolation metadata"
)

caf_isolation_metadata <- caf_isolation_metadata %>%
  dplyr::select(dplyr::all_of(required_columns))

if (anyNA(caf_isolation_metadata)) {
  stop("Missing values were detected in CAF isolation metadata.", call. = FALSE)
}

if (!is.numeric(caf_isolation_metadata$caf_pair_no)) {
  stop("caf_pair_no must be numeric.", call. = FALSE)
}

if (!is.numeric(caf_isolation_metadata$isolation_oxygen_percent)) {
  stop("isolation_oxygen_percent must be numeric.", call. = FALSE)
}

if (!is.numeric(caf_isolation_metadata$days_to_outgrowth)) {
  stop("days_to_outgrowth must be numeric.", call. = FALSE)
}

logical_columns <- c(
  "used_for_rna_seq",
  "used_for_proliferation_assay",
  "used_for_morphology_assay"
)

non_logical_columns <- logical_columns[
  !vapply(
    caf_isolation_metadata[logical_columns],
    is.logical,
    logical(1)
  )
]

if (length(non_logical_columns) > 0) {
  stop(
    "The following columns must contain TRUE/FALSE values: ",
    paste(non_logical_columns, collapse = ", "),
    call. = FALSE
  )
}

expected_pair_numbers <- 1:18
observed_pair_numbers <- sort(unique(caf_isolation_metadata$caf_pair_no))

if (!setequal(observed_pair_numbers, expected_pair_numbers)) {
  stop(
    "caf_pair_no must contain each integer from 1 to 18 exactly as a pair ID.",
    call. = FALSE
  )
}

check_allowed_values(
  caf_isolation_metadata$neoadjuvant_therapy,
  c("None", "GS", "GnP"),
  "neoadjuvant_therapy"
)

check_allowed_values(
  caf_isolation_metadata$tumor_location,
  c("Ph", "Pb", "Pt"),
  "tumor_location"
)

check_allowed_values(
  caf_isolation_metadata$isolation_oxygen_percent,
  c(1, 21),
  "isolation_oxygen_percent"
)

if (any(caf_isolation_metadata$days_to_outgrowth <= 0)) {
  stop("days_to_outgrowth must be greater than zero.", call. = FALSE)
}

if (any(duplicated(
  caf_isolation_metadata[c(
    "caf_pair_no",
    "isolation_oxygen_percent"
  )]
))) {
  stop(
    "Duplicated CAF pair and isolation oxygen combinations were detected.",
    call. = FALSE
  )
}

oxygen_count_by_pair <- caf_isolation_metadata %>%
  dplyr::count(
    .data$caf_pair_no,
    .data$isolation_oxygen_percent,
    name = "n"
  ) %>%
  tidyr::complete(
    caf_pair_no = expected_pair_numbers,
    isolation_oxygen_percent = c(1, 21),
    fill = list(n = 0)
  )

if (any(oxygen_count_by_pair$n != 1)) {
  invalid_pairs <- oxygen_count_by_pair %>%
    dplyr::filter(.data$n != 1) %>%
    dplyr::pull("caf_pair_no") %>%
    unique()

  stop(
    "Each CAF pair must contain one 1% O2 row and one 21% O2 row. ",
    "Invalid pair(s): ",
    paste(invalid_pairs, collapse = ", "),
    call. = FALSE
  )
}

check_unique_within_pair(
  caf_isolation_metadata,
  neoadjuvant_therapy,
  "neoadjuvant_therapy"
)

check_unique_within_pair(
  caf_isolation_metadata,
  tumor_location,
  "tumor_location"
)

check_unique_within_pair(
  caf_isolation_metadata,
  used_for_rna_seq,
  "used_for_rna_seq"
)

check_unique_within_pair(
  caf_isolation_metadata,
  used_for_proliferation_assay,
  "used_for_proliferation_assay"
)

check_unique_within_pair(
  caf_isolation_metadata,
  used_for_morphology_assay,
  "used_for_morphology_assay"
)

message(
  "Validated CAF isolation metadata: ",
  nrow(caf_isolation_metadata),
  " rows representing ",
  dplyr::n_distinct(caf_isolation_metadata$caf_pair_no),
  " paired CAF isolations."
)


## Create pair-level summary data -------------------------------------------

caf_isolation_summary_data <- caf_isolation_metadata %>%
  dplyr::arrange(
    .data$caf_pair_no,
    dplyr::desc(.data$isolation_oxygen_percent)
  ) %>%
  dplyr::mutate(
    outgrowth_column = paste0(
      "days_to_outgrowth_",
      .data$isolation_oxygen_percent
    )
  ) %>%
  tidyr::pivot_wider(
    id_cols = dplyr::all_of(c(
      "caf_pair_no",
      "neoadjuvant_therapy",
      "tumor_location",
      "used_for_rna_seq",
      "used_for_proliferation_assay",
      "used_for_morphology_assay"
    )),
    names_from = "outgrowth_column",
    values_from = "days_to_outgrowth"
  ) %>%
  dplyr::select(
    dplyr::all_of(c(
      "caf_pair_no",
      "neoadjuvant_therapy",
      "tumor_location",
      "days_to_outgrowth_21",
      "days_to_outgrowth_1",
      "used_for_rna_seq",
      "used_for_proliferation_assay",
      "used_for_morphology_assay"
    ))
  ) %>%
  dplyr::arrange(.data$caf_pair_no)

stopifnot(nrow(caf_isolation_summary_data) == 18)
stopifnot(!anyNA(caf_isolation_summary_data))


## Format supplementary table -----------------------------------------------

caf_isolation_summary_table <- caf_isolation_summary_data %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(c(
        "used_for_rna_seq",
        "used_for_proliferation_assay",
        "used_for_morphology_assay"
      )),
      ~ dplyr::if_else(.x, "Yes", "No")
    )
  ) %>%
  dplyr::rename(
    `CAF pair No.` = caf_pair_no,
    `Neoadjuvant therapy` = neoadjuvant_therapy,
    `Tumor location` = tumor_location,
    `Days to outgrowth (21% O2)` = days_to_outgrowth_21,
    `Days to outgrowth (1% O2)` = days_to_outgrowth_1,
    `RNA-seq` = used_for_rna_seq,
    `Proliferation assay` = used_for_proliferation_assay,
    `Morphology analysis` = used_for_morphology_assay
  )


## Save supplementary table --------------------------------------------------

readr::write_csv(
  caf_isolation_summary_table,
  caf_isolation_summary_table_file
)

message("Saved CAF isolation summary table: ", caf_isolation_summary_table_file)


## Print analysis summary ----------------------------------------------------

outgrowth_summary <- caf_isolation_metadata %>%
  dplyr::group_by(.data$isolation_oxygen_percent) %>%
  dplyr::summarise(
    median_days_to_outgrowth = median(.data$days_to_outgrowth),
    mean_days_to_outgrowth = mean(.data$days_to_outgrowth),
    sd_days_to_outgrowth = sd(.data$days_to_outgrowth),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(.data$isolation_oxygen_percent))

assay_summary <- caf_isolation_summary_data %>%
  dplyr::summarise(
    rna_seq_pairs = sum(.data$used_for_rna_seq),
    proliferation_assay_pairs = sum(.data$used_for_proliferation_assay),
    morphology_assay_pairs = sum(.data$used_for_morphology_assay)
  )

message("\nCAF isolation summary")
message("CAF pairs: ", nrow(caf_isolation_summary_data))
message("\nDays to outgrowth by isolation oxygen condition:")
print(outgrowth_summary)
message("\nPairs used for downstream assays:")
print(assay_summary)
