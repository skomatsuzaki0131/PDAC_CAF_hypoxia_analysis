# Analyze CAF proliferation
#
# This script analyzes longitudinal luminescence measurements from paired CAF
# cultures established under 21% and 1% O2 during isolation and subsequently
# assayed under 21% and 1% O2. Raw luminescence values are normalized to the
# mean Day 1 value across three replicate cultures within each CAF pair and
# oxygen-condition combination. It generates analysis summary tables only;
# figure generation is handled separately in code/08_figures.
#
# Statistical design:
#   - Three replicate cultures per CAF pair and oxygen-condition combination
#   - Normalized luminescence calculated relative to the corresponding mean
#     Day 1 raw luminescence
#   - Five prespecified two-sided Welch's t-tests at Day 7 within each CAF pair


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

caf_proliferation_luminescence_file <- file.path(
  caf_isolation_input_dir,
  "measurements",
  "caf_proliferation_luminescence.csv"
)

caf_proliferation_normalized_luminescence_file <- file.path(
  caf_isolation_table_dir,
  "caf_proliferation_normalized_luminescence.csv"
)

caf_proliferation_day7_descriptive_statistics_file <- file.path(
  caf_isolation_table_dir,
  "caf_proliferation_day7_descriptive_statistics.csv"
)

caf_proliferation_day7_pairwise_comparison_file <- file.path(
  caf_isolation_table_dir,
  "caf_proliferation_day7_pairwise_comparison.csv"
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

assign_proliferation_group <- function(
    isolation_oxygen_percent,
    assay_oxygen_percent
) {
  dplyr::case_when(
    isolation_oxygen_percent == "1" &
      assay_oxygen_percent == "21" ~ "Hypo-CAF in 21% O2",
    isolation_oxygen_percent == "1" &
      assay_oxygen_percent == "1" ~ "Hypo-CAF",
    isolation_oxygen_percent == "21" &
      assay_oxygen_percent == "21" ~ "Normo-CAF",
    isolation_oxygen_percent == "21" &
      assay_oxygen_percent == "1" ~ "Normo-CAF in 1% O2",
    TRUE ~ NA_character_
  )
}


## Read CAF proliferation measurements --------------------------------------

caf_proliferation_data <- read_required_csv(
  caf_proliferation_luminescence_file,
  "CAF proliferation luminescence file"
)

required_measurement_columns <- c(
  "caf_pair_no",
  "isolation_oxygen_percent",
  "assay_oxygen_percent",
  "day",
  "time_after_seeding_hours",
  "replicate_no",
  "raw_luminescence"
)

missing_measurement_columns <- setdiff(
  required_measurement_columns,
  colnames(caf_proliferation_data)
)

if (length(missing_measurement_columns) > 0) {
  stop(
    "Missing required columns in CAF proliferation measurements: ",
    paste(missing_measurement_columns, collapse = ", "),
    call. = FALSE
  )
}

unexpected_measurement_columns <- setdiff(
  colnames(caf_proliferation_data),
  required_measurement_columns
)

if (length(unexpected_measurement_columns) > 0) {
  stop(
    "Unexpected columns in CAF proliferation measurements: ",
    paste(unexpected_measurement_columns, collapse = ", "),
    call. = FALSE
  )
}

caf_proliferation_data <- caf_proliferation_data %>%
  dplyr::transmute(
    caf_pair_no = stringr::str_remove(
      as.character(.data$caf_pair_no),
      "^pair_"
    ),
    isolation_oxygen_percent = as.character(
      .data$isolation_oxygen_percent
    ),
    assay_oxygen_percent = as.character(
      .data$assay_oxygen_percent
    ),
    day = as.integer(.data$day),
    time_after_seeding_hours = as.numeric(
      .data$time_after_seeding_hours
    ),
    replicate_no = as.integer(.data$replicate_no),
    raw_luminescence = as.numeric(.data$raw_luminescence)
  ) %>%
  dplyr::mutate(
    group = assign_proliferation_group(
      .data$isolation_oxygen_percent,
      .data$assay_oxygen_percent
    )
  ) %>%
  dplyr::arrange(
    as.integer(.data$caf_pair_no),
    factor(
      .data$group,
      levels = c(
        "Hypo-CAF in 21% O2",
        "Hypo-CAF",
        "Normo-CAF",
        "Normo-CAF in 1% O2"
      )
    ),
    .data$day,
    .data$replicate_no
  )


## Validate measurement data -------------------------------------------------

if (anyNA(caf_proliferation_data)) {
  stop(
    "Missing, non-convertible, or unrecognized values were detected in the ",
    "CAF proliferation data.",
    call. = FALSE
  )
}

if (any(caf_proliferation_data$time_after_seeding_hours <= 0)) {
  stop(
    "Time after seeding must be greater than zero.",
    call. = FALSE
  )
}

if (any(caf_proliferation_data$raw_luminescence <= 0)) {
  stop("Raw luminescence values must be greater than zero.", call. = FALSE)
}

if (any(duplicated(
  caf_proliferation_data[c(
    "caf_pair_no",
    "isolation_oxygen_percent",
    "assay_oxygen_percent",
    "day",
    "replicate_no"
  )]
))) {
  stop(
    "Duplicated CAF pair, oxygen condition, day, and replicate combinations ",
    "were detected.",
    call. = FALSE
  )
}

observed_isolation_oxygen_conditions <- sort(unique(
  caf_proliferation_data$isolation_oxygen_percent
))

if (!identical(observed_isolation_oxygen_conditions, c("1", "21"))) {
  stop(
    "isolation_oxygen_percent must contain only '21' and '1'.",
    call. = FALSE
  )
}

observed_assay_oxygen_conditions <- sort(unique(
  caf_proliferation_data$assay_oxygen_percent
))

if (!identical(observed_assay_oxygen_conditions, c("1", "21"))) {
  stop(
    "assay_oxygen_percent must contain only '21' and '1'.",
    call. = FALSE
  )
}

if (!all(caf_proliferation_data$day %in% 1:7)) {
  stop("day must contain only integers from 1 to 7.", call. = FALSE)
}

if (!all(caf_proliferation_data$replicate_no %in% 1:3)) {
  stop("replicate_no must contain only integers from 1 to 3.", call. = FALSE)
}

measurement_counts <- caf_proliferation_data %>%
  dplyr::count(
    .data$caf_pair_no,
    .data$isolation_oxygen_percent,
    .data$assay_oxygen_percent,
    .data$day,
    name = "replicate_count"
  )

if (!all(measurement_counts$replicate_count == 3)) {
  stop(
    "Each CAF pair, oxygen-condition combination, and day must contain ",
    "exactly three replicate measurements.",
    call. = FALSE
  )
}

replicate_number_validation <- caf_proliferation_data %>%
  dplyr::group_by(
    .data$caf_pair_no,
    .data$isolation_oxygen_percent,
    .data$assay_oxygen_percent,
    .data$day
  ) %>%
  dplyr::summarise(
    valid_replicate_numbers = identical(
      sort(.data$replicate_no),
      1:3
    ),
    .groups = "drop"
  )

if (!all(replicate_number_validation$valid_replicate_numbers)) {
  stop(
    "replicate_no must run from 1 to 3 within each CAF pair, oxygen-condition ",
    "combination, and day.",
    call. = FALSE
  )
}

day_validation <- caf_proliferation_data %>%
  dplyr::distinct(
    .data$caf_pair_no,
    .data$isolation_oxygen_percent,
    .data$assay_oxygen_percent,
    .data$day
  ) %>%
  dplyr::count(
    .data$caf_pair_no,
    .data$isolation_oxygen_percent,
    .data$assay_oxygen_percent,
    name = "day_count"
  )

if (!all(day_validation$day_count == 7)) {
  stop(
    "Each CAF pair and oxygen-condition combination must contain exactly ",
    "seven measurement days.",
    call. = FALSE
  )
}

condition_validation <- caf_proliferation_data %>%
  dplyr::distinct(
    .data$caf_pair_no,
    .data$isolation_oxygen_percent,
    .data$assay_oxygen_percent
  ) %>%
  dplyr::count(
    .data$caf_pair_no,
    name = "condition_count"
  )

if (!all(condition_validation$condition_count == 4)) {
  stop(
    "Each CAF pair must contain all four isolation and assay oxygen-condition ",
    "combinations.",
    call. = FALSE
  )
}

time_validation <- caf_proliferation_data %>%
  dplyr::group_by(
    .data$caf_pair_no,
    .data$isolation_oxygen_percent,
    .data$assay_oxygen_percent,
    .data$day
  ) %>%
  dplyr::summarise(
    time_count = dplyr::n_distinct(.data$time_after_seeding_hours),
    .groups = "drop"
  )

if (!all(time_validation$time_count == 1)) {
  stop(
    "Replicates within the same CAF pair, oxygen-condition combination, and ",
    "day must have the same time after seeding.",
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
  "used_for_proliferation_assay"
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

proliferation_metadata <- caf_isolation_metadata %>%
  dplyr::filter(.data$used_for_proliferation_assay) %>%
  dplyr::transmute(
    caf_pair_no = as.character(.data$caf_pair_no),
    isolation_oxygen_percent = as.character(
      .data$isolation_oxygen_percent
    )
  ) %>%
  dplyr::distinct()

measurement_cultures <- caf_proliferation_data %>%
  dplyr::distinct(
    .data$caf_pair_no,
    .data$isolation_oxygen_percent
  )

if (!setequal(
  proliferation_metadata,
  measurement_cultures
)) {
  stop(
    "CAF cultures in the proliferation file do not match cultures marked ",
    "used_for_proliferation_assay in caf_isolation_metadata.csv.",
    call. = FALSE
  )
}

message(
  "Validated CAF proliferation data: ",
  nrow(caf_proliferation_data),
  " measurements from ",
  dplyr::n_distinct(caf_proliferation_data$caf_pair_no),
  " CAF pairs, four oxygen-condition combinations, seven days, and three ",
  "replicate cultures."
)


## Calculate normalized luminescence -----------------------------------------

caf_proliferation_normalized_luminescence <- caf_proliferation_data %>%
  dplyr::group_by(
    .data$caf_pair_no,
    .data$isolation_oxygen_percent,
    .data$assay_oxygen_percent
  ) %>%
  dplyr::mutate(
    mean_day1_luminescence = mean(
      .data$raw_luminescence[.data$day == 1]
    ),
    normalized_luminescence =
      .data$raw_luminescence / .data$mean_day1_luminescence,
    log2_normalized_luminescence = log2(
      .data$normalized_luminescence
    )
  ) %>%
  dplyr::ungroup()

if (any(!is.finite(
  caf_proliferation_normalized_luminescence$
    log2_normalized_luminescence
))) {
  stop(
    "Non-finite log2-transformed normalized luminescence values were ",
    "detected.",
    call. = FALSE
  )
}


## Calculate Day 7 descriptive statistics -----------------------------------

caf_proliferation_day7_data <-
  caf_proliferation_normalized_luminescence %>%
  dplyr::filter(.data$day == 7)

caf_proliferation_day7_descriptive_statistics <-
  caf_proliferation_day7_data %>%
  dplyr::group_by(
    .data$caf_pair_no,
    .data$isolation_oxygen_percent,
    .data$assay_oxygen_percent,
    .data$group
  ) %>%
  dplyr::summarise(
    replicates = dplyr::n(),
    mean_log2_normalized_luminescence = mean(
      .data$log2_normalized_luminescence
    ),
    standard_deviation = stats::sd(
      .data$log2_normalized_luminescence
    ),
    standard_error =
      .data$standard_deviation / sqrt(.data$replicates),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    group_order = factor(
      .data$group,
      levels = c(
        "Hypo-CAF in 21% O2",
        "Hypo-CAF",
        "Normo-CAF",
        "Normo-CAF in 1% O2"
      )
    )
  ) %>%
  dplyr::arrange(
    as.integer(.data$caf_pair_no),
    .data$group_order
  ) %>%
  dplyr::transmute(
    `CAF pair No.` = .data$caf_pair_no,
    `Isolation oxygen (%)` = .data$isolation_oxygen_percent,
    `Assay oxygen (%)` = .data$assay_oxygen_percent,
    `Group` = .data$group,
    `Replicates` = .data$replicates,
    `Mean log2(NL)` = .data$mean_log2_normalized_luminescence,
    `Standard deviation` = .data$standard_deviation,
    `Standard error` = .data$standard_error
  )


## Compare planned groups at Day 7 within each CAF pair ----------------------

planned_comparisons <- tibble::tribble(
  ~comparison_no, ~group_1,              ~group_2,
  1L, "Hypo-CAF in 21% O2",             "Hypo-CAF",
  2L, "Hypo-CAF in 21% O2",             "Normo-CAF",
  3L, "Hypo-CAF",                       "Normo-CAF",
  4L, "Hypo-CAF",                       "Normo-CAF in 1% O2",
  5L, "Normo-CAF",                      "Normo-CAF in 1% O2"
)

analyze_comparison <- function(
    pair_number,
    comparison_number,
    comparison_group_1,
    comparison_group_2
) {
  pair_data <- caf_proliferation_day7_data %>%
    dplyr::filter(.data$caf_pair_no == pair_number)

  values_group_1 <- pair_data %>%
    dplyr::filter(.data$group == comparison_group_1) %>%
    dplyr::pull(.data$log2_normalized_luminescence)

  values_group_2 <- pair_data %>%
    dplyr::filter(.data$group == comparison_group_2) %>%
    dplyr::pull(.data$log2_normalized_luminescence)

  if (length(values_group_1) != 3 || length(values_group_2) != 3) {
    stop(
      "Each group must contain exactly three Day 7 replicate values.",
      call. = FALSE
    )
  }

  t_test_result <- stats::t.test(
    x = values_group_1,
    y = values_group_2,
    alternative = "two.sided",
    var.equal = FALSE
  )

  tibble::tibble(
    caf_pair_no = pair_number,
    comparison_no = comparison_number,
    group_1 = comparison_group_1,
    group_2 = comparison_group_2,
    replicates_group_1 = length(values_group_1),
    replicates_group_2 = length(values_group_2),
    mean_group_1 = mean(values_group_1),
    mean_group_2 = mean(values_group_2),
    mean_difference = mean(values_group_1) - mean(values_group_2),
    t_statistic = unname(t_test_result$statistic),
    degrees_of_freedom = unname(t_test_result$parameter),
    p_value = t_test_result$p.value
  )
}

caf_proliferation_day7_pairwise_comparison <- tidyr::crossing(
  caf_pair_no = sort(unique(caf_proliferation_day7_data$caf_pair_no)),
  planned_comparisons
) %>%
  purrr::pmap_dfr(
    function(
        caf_pair_no,
        comparison_no,
        group_1,
        group_2
    ) {
      analyze_comparison(
        pair_number = caf_pair_no,
        comparison_number = comparison_no,
        comparison_group_1 = group_1,
        comparison_group_2 = group_2
      )
    }
  ) %>%
  dplyr::mutate(
    significance = dplyr::case_when(
      .data$p_value < 0.001 ~ "***",
      .data$p_value < 0.01 ~ "**",
      .data$p_value < 0.05 ~ "*",
      TRUE ~ "n.s."
    )
  ) %>%
  dplyr::arrange(
    as.integer(.data$caf_pair_no),
    .data$comparison_no
  ) %>%
  dplyr::transmute(
    `CAF pair No.` = .data$caf_pair_no,
    `Comparison No.` = .data$comparison_no,
    `Group 1` = .data$group_1,
    `Group 2` = .data$group_2,
    `Replicates (group 1)` = .data$replicates_group_1,
    `Replicates (group 2)` = .data$replicates_group_2,
    `Mean log2(NL) (group 1)` = .data$mean_group_1,
    `Mean log2(NL) (group 2)` = .data$mean_group_2,
    `Mean difference` = .data$mean_difference,
    `t statistic` = .data$t_statistic,
    `Degrees of freedom` = .data$degrees_of_freedom,
    `P value` = .data$p_value,
    `Significance` = .data$significance,
    `Statistical test` = "Two-sided Welch's t-test"
  )


## Format normalized luminescence table -------------------------------------

caf_proliferation_normalized_luminescence <-
  caf_proliferation_normalized_luminescence %>%
  dplyr::mutate(
    group_order = factor(
      .data$group,
      levels = c(
        "Hypo-CAF in 21% O2",
        "Hypo-CAF",
        "Normo-CAF",
        "Normo-CAF in 1% O2"
      )
    )
  ) %>%
  dplyr::arrange(
    as.integer(.data$caf_pair_no),
    .data$group_order,
    .data$day,
    .data$replicate_no
  ) %>%
  dplyr::transmute(
    `CAF pair No.` = .data$caf_pair_no,
    `Isolation oxygen (%)` = .data$isolation_oxygen_percent,
    `Assay oxygen (%)` = .data$assay_oxygen_percent,
    `Group` = .data$group,
    `Day` = .data$day,
    `Time after seeding (h)` = .data$time_after_seeding_hours,
    `Replicate No.` = .data$replicate_no,
    `Raw luminescence` = .data$raw_luminescence,
    `Mean Day 1 luminescence` = .data$mean_day1_luminescence,
    `Normalized luminescence` = .data$normalized_luminescence,
    `log2(NL)` = .data$log2_normalized_luminescence
  )


## Save analysis summary tables ----------------------------------------------

readr::write_csv(
  caf_proliferation_normalized_luminescence,
  caf_proliferation_normalized_luminescence_file
)

readr::write_csv(
  caf_proliferation_day7_descriptive_statistics,
  caf_proliferation_day7_descriptive_statistics_file
)

readr::write_csv(
  caf_proliferation_day7_pairwise_comparison,
  caf_proliferation_day7_pairwise_comparison_file
)


## Print analysis summary ----------------------------------------------------

message("\nCAF proliferation analysis")

purrr::pwalk(
  caf_proliferation_day7_pairwise_comparison,
  function(
    `CAF pair No.`,
    `Comparison No.`,
    `Group 1`,
    `Group 2`,
    `Replicates (group 1)`,
    `Replicates (group 2)`,
    `Mean log2(NL) (group 1)`,
    `Mean log2(NL) (group 2)`,
    `Mean difference`,
    `t statistic`,
    `Degrees of freedom`,
    `P value`,
    `Significance`,
    `Statistical test`
  ) {
    message(
      "CAF pair #",
      `CAF pair No.`,
      ", ",
      `Group 1`,
      " vs ",
      `Group 2`,
      ": P = ",
      signif(`P value`, 4),
      " (",
      `Significance`,
      "; ",
      `Statistical test`,
      ")"
    )
  }
)

message("\nSaved analysis tables:")
message("  ", caf_proliferation_normalized_luminescence_file)
message("  ", caf_proliferation_day7_descriptive_statistics_file)
message("  ", caf_proliferation_day7_pairwise_comparison_file)
