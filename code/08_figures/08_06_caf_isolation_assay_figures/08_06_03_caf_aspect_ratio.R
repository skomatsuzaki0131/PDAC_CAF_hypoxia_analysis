# CAF aspect-ratio figure
#
# This script generates the CAF morphology figure showing single-cell aspect
# ratios under 21% and 1% O2 for CAF pairs #13, #16, and #17.
#
# Statistical results are read from the output of
# 06_03_analyze_caf_aspect_ratio.R; statistical tests are not repeated here.
#
# Note: During manuscript preparation, project_dir is set to a local repository
# path. For public release, replace it with project_dir <- getwd(), and run this
# script from the repository root directory.


## Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggsignif)
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


## Define input and output paths ---------------------------------------------

caf_aspect_ratio_measurement_file <- file.path(
  caf_isolation_input_dir,
  "measurements",
  "caf_aspect_ratio_measurement.csv"
)

caf_aspect_ratio_pairwise_comparison_file <- file.path(
  caf_isolation_table_dir,
  "caf_aspect_ratio_pairwise_comparison.csv"
)

caf_aspect_ratio_figure_dir <- file.path(
  caf_isolation_figure_dir,
  "08_06_03_caf_aspect_ratio_figures"
)

caf_aspect_ratio_figure_file <- file.path(
  caf_aspect_ratio_figure_dir,
  "caf_aspect_ratio.pdf"
)

fs::dir_create(caf_aspect_ratio_figure_dir)


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

p_value_to_significance <- function(p_value) {
  dplyr::case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ "n.s."
  )
}


## Validate shared plotting parameters ---------------------------------------

required_caf_culture_colors <- c("N_CAF", "H_CAF")

if (!exists("caf_culture_colors")) {
  stop(
    "caf_culture_colors is not defined in 00_01_analysis_parameters.R.",
    call. = FALSE
  )
}

if (!all(required_caf_culture_colors %in% names(caf_culture_colors))) {
  stop(
    "caf_culture_colors must contain N_CAF and H_CAF.",
    call. = FALSE
  )
}

## Read data and statistical results -----------------------------------------

caf_aspect_ratio_data <- read_required_csv(
  caf_aspect_ratio_measurement_file,
  "CAF aspect-ratio measurement file"
)

caf_aspect_ratio_comparison <- read_required_csv(
  caf_aspect_ratio_pairwise_comparison_file,
  "CAF aspect-ratio pairwise-comparison table"
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

required_comparison_columns <- c(
  "CAF pair No.",
  "P value"
)

missing_comparison_columns <- setdiff(
  required_comparison_columns,
  colnames(caf_aspect_ratio_comparison)
)

if (length(missing_comparison_columns) > 0) {
  stop(
    "Missing required columns in CAF aspect-ratio comparison table: ",
    paste(missing_comparison_columns, collapse = ", "),
    call. = FALSE
  )
}


## Prepare plotting data -----------------------------------------------------

caf_pair_order <- c("13", "16", "17")

caf_aspect_ratio_plot_data <- caf_aspect_ratio_data %>%
  dplyr::transmute(
    caf_pair_no = factor(
      as.character(.data$caf_pair_no),
      levels = caf_pair_order
    ),
    isolation_oxygen_percent = factor(
      as.character(.data$isolation_oxygen_percent),
      levels = c("21", "1")
    ),
    caf_culture = dplyr::recode(
      as.character(.data$isolation_oxygen_percent),
      `21` = "N_CAF",
      `1` = "H_CAF"
    ),
    cell_number = as.integer(.data$cell_number),
    aspect_ratio = as.numeric(.data$aspect_ratio)
  )

if (anyNA(caf_aspect_ratio_plot_data)) {
  stop(
    "Missing or unrecognized values were detected in the plotting data.",
    call. = FALSE
  )
}

caf_aspect_ratio_significance <- caf_aspect_ratio_comparison %>%
  dplyr::transmute(
    caf_pair_no = factor(
      as.character(.data$`CAF pair No.`),
      levels = caf_pair_order
    ),
    group1 = "21",
    group2 = "1",
    y_position = c(20.2, 20.8, 18.5),
    annotation = p_value_to_significance(.data$`P value`)
  )

if (anyNA(caf_aspect_ratio_significance$caf_pair_no)) {
  stop(
    "Unexpected CAF pair numbers were detected in the comparison table.",
    call. = FALSE
  )
}


## Figure settings -----------------------------------------------------------

text_size <- 8

caf_culture_fill_colors <- c(
  "21" = unname(caf_culture_colors["N_CAF"]),
  "1" = unname(caf_culture_colors["H_CAF"])
)

caf_culture_outline_colors <- c(
  "21" = "#B2182B",
  "1" = "#2166AC"
)


## Generate CAF aspect-ratio figure ------------------------------------------

plot_caf_aspect_ratio <- ggplot(
  caf_aspect_ratio_plot_data,
  aes(
    x = .data$isolation_oxygen_percent,
    y = .data$aspect_ratio
  )
) +
  geom_boxplot(
    aes(
      fill = .data$isolation_oxygen_percent,
      color = .data$isolation_oxygen_percent
    ),
    width = 0.72,
    staplewidth = 0.9,
    linewidth = 1.1,
    median.linewidth = 0.5,
    outlier.shape = NA
  ) +
  geom_dotplot(
    binaxis = "y",
    stackdir = "center",
    binwidth = 0.4,
    dotsize = 0.85,
    color = "black",
    fill = "black",
    stroke = 1.3,
    alpha = 0.8
  ) +
  ggsignif::geom_signif(
    data = caf_aspect_ratio_significance,
    aes(
      xmin = .data$group1,
      xmax = .data$group2,
      annotations = .data$annotation,
      y_position = .data$y_position
    ),
    manual = TRUE,
    inherit.aes = FALSE,
    textsize = text_size,
    size = 1.0,
    tip_length = 0.03,
    vjust = 0.3
  ) +
  facet_grid(
    cols = vars(.data$caf_pair_no),
    labeller = labeller(
      caf_pair_no = function(x) paste0("CAF pair #", x)
    )
  ) +
  scale_x_discrete(
    labels = c(
      "21" = "Normo-CAFs",
      "1" = "Hypo-CAFs"
    )
  ) +
  scale_y_continuous(
    limits = c(0, 25),
    breaks = seq(0, 20, by = 5),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_fill_manual(
    values = caf_culture_fill_colors
  ) +
  scale_color_manual(
    values = caf_culture_outline_colors
  ) +
  labs(
    x = NULL,
    y = "Aspect ratio"
  ) +
  theme(
    aspect.ratio = 2.28,
    text = element_text(),
    legend.position = "none",
    axis.text = element_text(
      face = "plain",
      color = "black"
    ),
    axis.text.x = element_text(
      color = "black",
      size = text_size * 3,
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(
      face = "plain",
      color = "black",
      size = text_size * 3
    ),
    axis.title = element_text(
      face = "plain",
      color = "black",
      size = text_size * 3,
    ),
    strip.text = element_text(
      face = "plain",
      color = "black",
      size = text_size * 3,
    ),
    strip.background = element_blank(),
    plot.background = element_rect(
      fill = "transparent",
      color = NA
    ),
    panel.background = element_rect(
      fill = "transparent",
      color = NA
    ),
    panel.grid = element_blank(),
    panel.border = element_rect(
      fill = NA,
      color = "black"
    ),
    axis.line = element_blank(),
    panel.spacing.x = unit(0.8, "lines")
  )


## Save figure ---------------------------------------------------------------

ggsave(
  filename = caf_aspect_ratio_figure_file,
  plot = plot_caf_aspect_ratio,
  device = cairo_pdf,
  width = 9,
  height = 8.5,
  units = "in",
  bg = "transparent"
)

message("Saved CAF aspect-ratio figure:")
message("  ", caf_aspect_ratio_figure_file)
