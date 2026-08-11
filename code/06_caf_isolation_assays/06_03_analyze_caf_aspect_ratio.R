# Analyze CAF aspect ratio
#
# This script analyzes single-cell aspect-ratio measurements from CAF pairs
# cultured under 21% and 1% O2 during isolation. It generates analysis summary
# tables only; figure generation is handled separately in code/08_figures.
#
# Statistical design:
#   - Thirty independently selected cells per CAF culture
#   - 21% O2 versus 1% O2 compared separately within each CAF pair
#   - Exact Wilcoxon rank-sum test (unpaired)


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

source(file.path(
  project_dir,
  "code",
  "00_config",
  "00_02_paths.R"
))


## Define input and output paths ---------------------------------------------

caf_aspect_ratio_measurement_file <- file.path(
  caf_isolation_input_dir,
  "measurements",
  "caf_aspect_ratio_measurement.csv"
)

caf_aspect_ratio_descriptive_statistics_file <- file.path(
  caf_isolation_table_dir,
  "caf_aspect_ratio_descriptive_statistics.csv"
)

caf_aspect_ratio_pairwise_comparison_file <- file.path(
  caf_isolation_table_dir,
  "caf_aspect_ratio_pairwise_comparison.csv"
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


## Read CAF aspect-ratio measurements ----------------------------------------

caf_aspect_ratio_data <- read_required_csv(
  caf_aspect_ratio_measurement_file,
  "CAF aspect-ratio measurement file"
)

required_measurement_columns <- c(
  "caf_pair_no",
  "isolation_oxygen_percent",
  "cell_number",
  "aspect_ratio"
)

missing_measurement_columns <- setdiff(
  required_measurement_columns,
  colnames(caf_aspect_ratio_data)
)

if (length(missing_measurement_columns) > 0) {
  stop(
    "Missing required columns in CAF aspect-ratio measurements: ",
    paste(missing_measurement_columns, collapse = ", "),
    call. = FALSE
  )
}

unexpected_measurement_columns <- setdiff(
  colnames(caf_aspect_ratio_data),
  required_measurement_columns
)

if (length(unexpected_measurement_columns) > 0) {
  stop(
    "Unexpected columns in CAF aspect-ratio measurements: ",
    paste(unexpected_measurement_columns, collapse = ", "),
    call. = FALSE
  )
}

caf_aspect_ratio_data <- caf_aspect_ratio_data %>%
  dplyr::transmute(
    caf_pair_no = as.character(.data$caf_pair_no),
    isolation_oxygen_percent = as.character(
      .data$isolation_oxygen_percent
    ),
    cell_number = as.integer(.data$cell_number),
    aspect_ratio = as.numeric(.data$aspect_ratio)
  ) %>%
  dplyr::arrange(
    as.integer(.data$caf_pair_no),
    factor(
      .data$isolation_oxygen_percent,
      levels = c("21", "1")
    ),
    .data$cell_number
  )


## Validate measurement data -------------------------------------------------

if (anyNA(caf_aspect_ratio_data)) {
  stop(
    "Missing or non-convertible values were detected in the aspect-ratio data.",
    call. = FALSE
  )
}

if (any(caf_aspect_ratio_data$aspect_ratio <= 0)) {
  stop("Aspect ratios must be greater than zero.", call. = FALSE)
}

if (any(duplicated(
  caf_aspect_ratio_data[c(
    "caf_pair_no",
    "isolation_oxygen_percent",
    "cell_number"
  )]
))) {
  stop(
    "Duplicated CAF pair, oxygen condition, and cell-number combinations ",
    "were detected.",
    call. = FALSE
  )
}

observed_oxygen_conditions <- sort(unique(
  caf_aspect_ratio_data$isolation_oxygen_percent
))

if (!identical(observed_oxygen_conditions, c("1", "21"))) {
  stop(
    "isolation_oxygen_percent must contain only '21' and '1'.",
    call. = FALSE
  )
}

culture_cell_counts <- caf_aspect_ratio_data %>%
  dplyr::count(
    .data$caf_pair_no,
    .data$isolation_oxygen_percent,
    name = "cell_count"
  )

if (!all(culture_cell_counts$cell_count == 30)) {
  stop(
    "Each CAF culture must contain exactly 30 measured cells.",
    call. = FALSE
  )
}

cell_number_validation <- caf_aspect_ratio_data %>%
  dplyr::group_by(
    .data$caf_pair_no,
    .data$isolation_oxygen_percent
  ) %>%
  dplyr::summarise(
    valid_cell_numbers = identical(
      sort(.data$cell_number),
      1:30
    ),
    .groups = "drop"
  )

if (!all(cell_number_validation$valid_cell_numbers)) {
  stop(
    "cell_number must run from 1 to 30 within each CAF culture.",
    call. = FALSE
  )
}


## Validate against CAF isolation metadata -----------------------------------

caf_isolation_metadata <- read_required_csv(
  caf_isolation_metadata_file,
  "CAF isolation metadata file"
)

required_metadata_columns <- c(
  "caf_pair_no",
  "isolation_oxygen_percent",
  "used_for_morphology_assay"
)

missing_metadata_columns <- setdiff(
  required_metadata_columns,
  colnames(caf_isolation_metadata)
)

if (length(missing_metadata_columns) > 0) {
  stop(
    "Missing required columns in CAF isolation metadata: ",
    paste(missing_metadata_columns, collapse = ", "),
    call. = FALSE
  )
}

morphology_metadata <- caf_isolation_metadata %>%
  dplyr::filter(.data$used_for_morphology_assay) %>%
  dplyr::transmute(
    caf_pair_no = as.character(.data$caf_pair_no),
    isolation_oxygen_percent = as.character(
      .data$isolation_oxygen_percent
    )
  ) %>%
  dplyr::distinct()

measurement_cultures <- caf_aspect_ratio_data %>%
  dplyr::distinct(
    .data$caf_pair_no,
    .data$isolation_oxygen_percent
  )

if (!setequal(
  morphology_metadata,
  measurement_cultures
)) {
  stop(
    "CAF cultures in the aspect-ratio file do not match cultures marked ",
    "used_for_morphology_assay in caf_isolation_metadata.csv.",
    call. = FALSE
  )
}

message(
  "Validated CAF aspect-ratio data: ",
  nrow(caf_aspect_ratio_data),
  " cells from ",
  nrow(culture_cell_counts),
  " cultures and ",
  dplyr::n_distinct(caf_aspect_ratio_data$caf_pair_no),
  " CAF pairs."
)


## Calculate descriptive statistics -----------------------------------------

# Aspect-ratio distributions were compared using exact Wilcoxon rank-sum
# tests because cell-level measurements were not assumed to be normally
# distributed.

caf_aspect_ratio_descriptive_statistics <- caf_aspect_ratio_data %>%
  dplyr::group_by(
    .data$caf_pair_no,
    .data$isolation_oxygen_percent
  ) %>%
  dplyr::summarise(
    cells = dplyr::n(),
    mean_aspect_ratio = mean(.data$aspect_ratio),
    standard_deviation = stats::sd(.data$aspect_ratio),
    median_aspect_ratio = median(.data$aspect_ratio),
    first_quartile = unname(
      stats::quantile(.data$aspect_ratio, 0.25)
    ),
    third_quartile = unname(
      stats::quantile(.data$aspect_ratio, 0.75)
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    median_iqr = format_median_iqr(
      .data$median_aspect_ratio,
      .data$first_quartile,
      .data$third_quartile
    ),
    oxygen_order = factor(
      .data$isolation_oxygen_percent,
      levels = c("21", "1")
    )
  ) %>%
  dplyr::arrange(
    as.integer(.data$caf_pair_no),
    .data$oxygen_order
  ) %>%
  dplyr::transmute(
    `CAF pair No.` = .data$caf_pair_no,
    `Isolation oxygen (%)` = .data$isolation_oxygen_percent,
    `Cells` = .data$cells,
    `Mean aspect ratio` = .data$mean_aspect_ratio,
    `Standard deviation` = .data$standard_deviation,
    `Median aspect ratio` = .data$median_aspect_ratio,
    `First quartile` = .data$first_quartile,
    `Third quartile` = .data$third_quartile,
    `Median (IQR)` = .data$median_iqr
  )


## Compare oxygen conditions within each CAF pair ----------------------------

analyze_pair <- function(pair_number) {
  pair_data <- caf_aspect_ratio_data %>%
    dplyr::filter(.data$caf_pair_no == pair_number)

  aspect_ratio_21 <- pair_data %>%
    dplyr::filter(.data$isolation_oxygen_percent == "21") %>%
    dplyr::pull(.data$aspect_ratio)

  aspect_ratio_1 <- pair_data %>%
    dplyr::filter(.data$isolation_oxygen_percent == "1") %>%
    dplyr::pull(.data$aspect_ratio)

  wilcox_result <- exactRankTests::wilcox.exact(
    x = aspect_ratio_21,
    y = aspect_ratio_1,
    paired = FALSE
  )

  tibble::tibble(
    caf_pair_no = pair_number,
    cells_21_percent_o2 = length(aspect_ratio_21),
    cells_1_percent_o2 = length(aspect_ratio_1),
    p_value = wilcox_result$p.value
  )
}

caf_aspect_ratio_pairwise_comparison <- sort(unique(
  caf_aspect_ratio_data$caf_pair_no
)) %>%
  purrr::map_dfr(analyze_pair) %>%
  dplyr::arrange(as.integer(.data$caf_pair_no)) %>%
  dplyr::transmute(
    `CAF pair No.` = .data$caf_pair_no,
    `Cells (21% O2)` = .data$cells_21_percent_o2,
    `Cells (1% O2)` = .data$cells_1_percent_o2,
    `P value` = .data$p_value,
    `Statistical test` = "Exact Wilcoxon rank-sum test"
  )


## Save analysis summary tables ----------------------------------------------

readr::write_csv(
  caf_aspect_ratio_descriptive_statistics,
  caf_aspect_ratio_descriptive_statistics_file
)

readr::write_csv(
  caf_aspect_ratio_pairwise_comparison,
  caf_aspect_ratio_pairwise_comparison_file
)


## Print analysis summary ----------------------------------------------------

message("\nCAF aspect-ratio analysis")

purrr::pwalk(
  caf_aspect_ratio_pairwise_comparison,
  function(
    `CAF pair No.`,
    `Cells (21% O2)`,
    `Cells (1% O2)`,
    `P value`,
    `Statistical test`
  ) {
    message(
      "CAF pair #",
      `CAF pair No.`,
      ": P = ",
      signif(`P value`, 4),
      " (",
      `Statistical test`,
      ")"
    )
  }
)

message("\nSaved analysis tables:")
message("  ", caf_aspect_ratio_descriptive_statistics_file)
message("  ", caf_aspect_ratio_pairwise_comparison_file)
