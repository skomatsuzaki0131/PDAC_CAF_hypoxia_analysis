# Analyze CAF outgrowth
#
# This script compares the time to CAF outgrowth between isolation oxygen
# conditions and evaluates associations with tumor location and neoadjuvant
# therapy. It generates analysis summary tables only; figure generation is
# handled separately in code/08_figures.


## Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(fs)
  library(exactRankTests)
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

caf_outgrowth_pair_summary_file <- file.path(
  caf_isolation_table_dir,
  "caf_outgrowth_pair_summary.csv"
)

caf_outgrowth_oxygen_comparison_file <- file.path(
  caf_isolation_table_dir,
  "caf_outgrowth_oxygen_comparison.csv"
)

caf_outgrowth_tumor_location_file <- file.path(
  caf_isolation_table_dir,
  "caf_outgrowth_by_tumor_location.csv"
)

caf_outgrowth_neoadjuvant_therapy_file <- file.path(
  caf_isolation_table_dir,
  "caf_outgrowth_by_neoadjuvant_therapy.csv"
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

format_number <- function(x) {
  format(
    x,
    trim = TRUE,
    scientific = FALSE
  )
}

format_median_iqr <- function(median_value, q1, q3) {
  paste0(
    format_number(median_value),
    " (",
    format_number(q1),
    "-",
    format_number(q3),
    ")"
  )
}


## Read and validate CAF isolation metadata ---------------------------------

caf_isolation_metadata <- read_required_csv(
  caf_isolation_metadata_file,
  "CAF isolation metadata file"
)

required_columns <- c(
  "caf_pair_no",
  "neoadjuvant_therapy",
  "tumor_location",
  "isolation_oxygen_percent",
  "days_to_outgrowth"
)

missing_columns <- setdiff(
  required_columns,
  colnames(caf_isolation_metadata)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns in CAF isolation metadata: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}


## Summarize outgrowth time by CAF pair --------------------------------------

# Tumor location and neoadjuvant therapy are pair-level characteristics.
# Therefore, the two oxygen-condition measurements are averaged within each
# CAF pair before evaluating their associations with these characteristics.

caf_outgrowth_pair_summary <-
  caf_isolation_metadata %>%
  dplyr::group_by(
    caf_pair_no,
    tumor_location,
    neoadjuvant_therapy
  ) %>%
  dplyr::summarise(
    mean_days_to_outgrowth = mean(
      days_to_outgrowth,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

if (nrow(caf_outgrowth_pair_summary) !=
    dplyr::n_distinct(caf_isolation_metadata$caf_pair_no)) {
  stop(
    "The pair-level outgrowth summary does not contain one row per CAF pair.",
    call. = FALSE
  )
}

caf_outgrowth_data <- caf_isolation_metadata %>%
  dplyr::select(dplyr::all_of(required_columns)) %>%
  dplyr::arrange(
    .data$caf_pair_no,
    dplyr::desc(.data$isolation_oxygen_percent)
  )

if (anyNA(caf_outgrowth_data)) {
  stop("Missing values were detected in CAF outgrowth data.", call. = FALSE)
}

if (any(duplicated(
  caf_outgrowth_data[c(
    "caf_pair_no",
    "isolation_oxygen_percent"
  )]
))) {
  stop(
    "Duplicated CAF pair and isolation oxygen combinations were detected.",
    call. = FALSE
  )
}

oxygen_values_by_pair <- caf_outgrowth_data %>%
  dplyr::group_by(.data$caf_pair_no) %>%
  dplyr::summarise(
    oxygen_conditions = list(sort(.data$isolation_oxygen_percent)),
    .groups = "drop"
  )

valid_oxygen_pairs <- vapply(
  oxygen_values_by_pair$oxygen_conditions,
  identical,
  logical(1),
  c(1, 21)
)

if (!all(valid_oxygen_pairs)) {
  stop(
    "Each CAF pair must contain one 1% O2 row and one 21% O2 row.",
    call. = FALSE
  )
}

message(
  "Validated CAF outgrowth data: ",
  nrow(caf_outgrowth_data),
  " cultures from ",
  dplyr::n_distinct(caf_outgrowth_data$caf_pair_no),
  " CAF pairs."
)


## Prepare paired outgrowth data ---------------------------------------------

caf_outgrowth_wide <- caf_outgrowth_data %>%
  dplyr::mutate(
    oxygen_column = paste0(
      "days_to_outgrowth_",
      .data$isolation_oxygen_percent
    )
  ) %>%
  tidyr::pivot_wider(
    id_cols = dplyr::all_of(c(
      "caf_pair_no",
      "neoadjuvant_therapy",
      "tumor_location"
    )),
    names_from = "oxygen_column",
    values_from = "days_to_outgrowth"
  ) %>%
  dplyr::arrange(.data$caf_pair_no)


## Assess distribution -------------------------------------------------------

shapiro_outgrowth <- stats::shapiro.test(
  caf_outgrowth_data$days_to_outgrowth
)

message(
  "Shapiro-Wilk test for all outgrowth measurements: W = ",
  signif(unname(shapiro_outgrowth$statistic), 4),
  ", P = ",
  signif(shapiro_outgrowth$p.value, 4)
)


## Compare isolation oxygen conditions --------------------------------------

oxygen_test <- exactRankTests::wilcox.exact(
  x = caf_outgrowth_wide$days_to_outgrowth_21,
  y = caf_outgrowth_wide$days_to_outgrowth_1,
  paired = TRUE
)

oxygen_descriptive_statistics <- caf_outgrowth_data %>%
  dplyr::group_by(.data$isolation_oxygen_percent) %>%
  dplyr::summarise(
    median_days_to_outgrowth = median(.data$days_to_outgrowth),
    first_quartile = unname(
      stats::quantile(.data$days_to_outgrowth, 0.25)
    ),
    third_quartile = unname(
      stats::quantile(.data$days_to_outgrowth, 0.75)
    ),
    .groups = "drop"
  )

oxygen_21_summary <- oxygen_descriptive_statistics %>%
  dplyr::filter(.data$isolation_oxygen_percent == 21)

oxygen_1_summary <- oxygen_descriptive_statistics %>%
  dplyr::filter(.data$isolation_oxygen_percent == 1)

caf_outgrowth_oxygen_comparison <- tibble::tibble(
  Comparison = "21% O2 vs 1% O2",
  `CAF pairs` = nrow(caf_outgrowth_wide),
  `Median days to outgrowth (21% O2)` =
    oxygen_21_summary$median_days_to_outgrowth,
  `IQR (21% O2)` = paste0(
    format_number(oxygen_21_summary$first_quartile),
    "-",
    format_number(oxygen_21_summary$third_quartile)
  ),
  `Median days to outgrowth (1% O2)` =
    oxygen_1_summary$median_days_to_outgrowth,
  `IQR (1% O2)` = paste0(
    format_number(oxygen_1_summary$first_quartile),
    "-",
    format_number(oxygen_1_summary$third_quartile)
  ),
  `P value` = oxygen_test$p.value,
  `Statistical test` = "Exact Wilcoxon signed-rank test"
)


## Analyze outgrowth time by tumor location ---------------------------------

caf_outgrowth_location_test <-
  stats::kruskal.test(
    mean_days_to_outgrowth ~ tumor_location,
    data = caf_outgrowth_pair_summary
  )

caf_outgrowth_by_tumor_location <-
  caf_outgrowth_pair_summary %>%
  dplyr::group_by(tumor_location) %>%
  dplyr::summarise(
    n_caf_pairs = dplyr::n(),
    median_days_to_outgrowth = stats::median(
      mean_days_to_outgrowth
    ),
    q1_days_to_outgrowth = stats::quantile(
      mean_days_to_outgrowth,
      0.25,
      names = FALSE
    ),
    q3_days_to_outgrowth = stats::quantile(
      mean_days_to_outgrowth,
      0.75,
      names = FALSE
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    "P value" = caf_outgrowth_location_test$p.value,
    statistical_test = "Kruskal-Wallis test"
  )

## Analyze outgrowth time by neoadjuvant therapy ----------------------------

caf_outgrowth_nat_test <-
  stats::kruskal.test(
    mean_days_to_outgrowth ~ neoadjuvant_therapy,
    data = caf_outgrowth_pair_summary
  )

caf_outgrowth_by_neoadjuvant_therapy <-
  caf_outgrowth_pair_summary %>%
  dplyr::group_by(neoadjuvant_therapy) %>%
  dplyr::summarise(
    n_caf_pairs = dplyr::n(),
    median_days_to_outgrowth = stats::median(
      mean_days_to_outgrowth
    ),
    q1_days_to_outgrowth = stats::quantile(
      mean_days_to_outgrowth,
      0.25,
      names = FALSE
    ),
    q3_days_to_outgrowth = stats::quantile(
      mean_days_to_outgrowth,
      0.75,
      names = FALSE
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    "P value" = caf_outgrowth_nat_test$p.value,
    statistical_test = "Kruskal-Wallis test"
  )


## Save analysis summary tables ----------------------------------------------

readr::write_csv(
  caf_outgrowth_pair_summary,
  caf_outgrowth_pair_summary_file
)

readr::write_csv(
  caf_outgrowth_oxygen_comparison,
  caf_outgrowth_oxygen_comparison_file
)

readr::write_csv(
  caf_outgrowth_by_tumor_location,
  caf_outgrowth_tumor_location_file
)

readr::write_csv(
  caf_outgrowth_by_neoadjuvant_therapy,
  caf_outgrowth_neoadjuvant_therapy_file
)


## Print analysis summary ----------------------------------------------------

message("\nCAF outgrowth analysis")

message(
  "Isolation oxygen comparison: P = ",
  signif(oxygen_test$p.value, 4),
  " (exact Wilcoxon signed-rank test)"
)

message(
  "Tumor location comparison: P = ",
  signif(caf_outgrowth_location_test$p.value, 4),
  " (Kruskal-Wallis test)"
)

message(
  "Neoadjuvant therapy comparison: P = ",
  signif(caf_outgrowth_nat_test$p.value, 4),
  " (Kruskal-Wallis test)"
)

message("\nSaved analysis tables:")
message("  ", caf_outgrowth_pair_summary_file)
message("  ", caf_outgrowth_oxygen_comparison_file)
message("  ", caf_outgrowth_tumor_location_file)
message("  ", caf_outgrowth_neoadjuvant_therapy_file)

