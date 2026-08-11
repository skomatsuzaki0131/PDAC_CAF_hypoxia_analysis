## Hallmark GSEA before/after oxygen switching dot plot ----------------------


## Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(cowplot)
  library(dplyr)
  library(fs)
  library(ggplot2)
  library(readr)
  library(stringr)
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

source(file.path(
  project_dir,
  "code",
  "00_config",
  "00_05_helper_functions.R"
))


## Input and output files -----------------------------------------------------

gene_set_analysis_table_dir <- file.path(
  bulk_rnaseq_table_dir,
  "gene_set_analysis"
)

hallmark_gsea_table_dir <- file.path(
  gene_set_analysis_table_dir,
  "hallmark_gsea"
)

hallmark_gsea_table_dir <- file.path(
  gene_set_analysis_table_dir,
  "hallmark_gsea"
)

hallmark_gsea_assay_oxygen_result_file <- file.path(
  hallmark_gsea_table_dir,
  "bulk_rnaseq_hallmark_gsea_assay_oxygen_comparisons.csv"
)

hallmark_assay_oxygen_condition_figure_dir <- file.path(
  bulk_rnaseq_figure_dir,
  "08_07_10_bulk_rnaseq_hallmark_assay_oxygen_condition_dotplot_figures"
)

fs::dir_create(hallmark_assay_oxygen_condition_figure_dir)

if (!file.exists(hallmark_gsea_assay_oxygen_result_file)) {
  stop(
    "Required Hallmark GSEA result table was not found: ",
    hallmark_gsea_assay_oxygen_result_file,
    call. = FALSE
  )
}


## Read and validate Hallmark GSEA results -----------------------------------

hallmark_gsea_assay_oxygen_results <- readr::read_csv(
  hallmark_gsea_assay_oxygen_result_file,
  show_col_types = FALSE
)

required_columns <- c(
  "comparison_name",
  "Description",
  "NES",
  "p.adjust"
)

missing_columns <- setdiff(
  required_columns,
  names(hallmark_gsea_assay_oxygen_results)
)

if (length(missing_columns) > 0) {
  stop(
    "The Hallmark GSEA result table is missing required columns: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

required_comparison_names <- c(
  "comparison_4_hypocaf_hypoxia_vs_hypocaf_normoxia",
  "comparison_5_normocaf_hypoxia_vs_normocaf_normoxia"
)

missing_comparison_names <- setdiff(
  required_comparison_names,
  unique(hallmark_gsea_assay_oxygen_results$comparison_name)
)

if (length(missing_comparison_names) > 0) {
  stop(
    "The Hallmark GSEA result table is missing comparisons: ",
    paste(missing_comparison_names, collapse = ", "),
    call. = FALSE
  )
}


## Select the displayed Hallmark gene sets -----------------------------------

excluded_gene_sets <- c(
  "iCAF_signature",
  "Buffa_Hypoxia",
  "Winter_Hypoxia",
  "KEGG_Cell_cycle"
)

select_extreme_gene_sets <- function(
    gsea_result_table,
    comparison_name
) {
  comparison_result_table <-
    gsea_result_table |>
    dplyr::filter(
      .data$comparison_name == .env$comparison_name,
      !.data$Description %in% .env$excluded_gene_sets
    )

  dplyr::bind_rows(
    dplyr::slice_max(
      comparison_result_table,
      order_by = .data$NES,
      n = 5,
      with_ties = TRUE
    ),
    dplyr::slice_min(
      comparison_result_table,
      order_by = .data$NES,
      n = 5,
      with_ties = TRUE
    )
  ) |>
    dplyr::arrange(dplyr::desc(.data$NES)) |>
    dplyr::mutate(
      Description = factor(
        .data$Description,
        levels = unique(.data$Description)
      ),
      minus_log10_adjusted_p_value = -log10(.data$p.adjust)
    )
}

hypocaf_switching_results <- select_extreme_gene_sets(
  gsea_result_table = hallmark_gsea_assay_oxygen_results,
  comparison_name =
    "comparison_4_hypocaf_hypoxia_vs_hypocaf_normoxia"
)

normocaf_switching_results <- select_extreme_gene_sets(
  gsea_result_table = hallmark_gsea_assay_oxygen_results,
  comparison_name =
    "comparison_5_normocaf_hypoxia_vs_normocaf_normoxia"
)


## Define common scales and labels -------------------------------------------

nes_range <- range(
  c(
    hypocaf_switching_results$NES,
    normocaf_switching_results$NES
  ),
  na.rm = TRUE
)

symmetric_nes_limit <- max(abs(nes_range))

minus_log10_adjusted_p_value_range <- range(
  c(
    hypocaf_switching_results$minus_log10_adjusted_p_value,
    normocaf_switching_results$minus_log10_adjusted_p_value
  ),
  na.rm = TRUE
)

maximum_minus_log10_adjusted_p_value <-
  max(minus_log10_adjusted_p_value_range)

format_hallmark_label <- function(gene_set_name) {
  formatted_label <- format_gene_set_label(gene_set_name)

  stringr::str_replace(
    formatted_label,
    pattern = "Epithelial mesenchymal transition",
    replacement = "EMT"
  )
}

common_labs <- labs(x = "NES")

common_scale_x_continuous <- scale_x_continuous(
  breaks = c(-2, 0, 2),
  expand = expansion(
    mult = c(0.05, 0.05)
  )
)

common_coordinate_system <- coord_cartesian(
  xlim = c(
    -symmetric_nes_limit,
    symmetric_nes_limit
  ),
  clip = "off"
)

common_scale_y_discrete <- scale_y_discrete(
  expand = expansion(
    mult = c(0.1, 0.1)
  )
)

point_size_scale <- scale_size_continuous(
  name = expression(-log[10]("adjusted P value")),
  breaks = c(2, 10, 30),
  limits = c(0, maximum_minus_log10_adjusted_p_value),
  range = c(0, 5.8),
  guide = guide_legend(
    override.aes = list(
      fill = "white",
      color = "black"
    ),
    label.position = "bottom"
  )
)

common_vertical_line <- geom_vline(
  xintercept = 0
)

common_theme <- theme(
  aspect.ratio = 0.727,
  plot.margin = margin(
    t = 5.5,
    r = 55,
    b = 5.5,
    l = 55,
    unit = "pt"
  ),
  axis.title.y = element_blank(),
  axis.title.x = element_text(
    face = "plain",
    color = "black"
  ),
  axis.text.y = element_blank(),
  axis.text.x = element_text(
    face = "plain",
    color = "black",
    size = 9
  ),
  axis.line.x = element_line(
    color = "black"
  ),
  axis.line.y = element_blank(),
  axis.ticks.y = element_blank(),
  legend.key = element_blank(),
  legend.title = element_text(
    face = "plain",
    color = "black"
  ),
  legend.text = element_text(
    face = "plain",
    color = "black"
  ),
  plot.background = element_rect(
    fill = "transparent",
    color = NA
  ),
  legend.background = element_rect(
    fill = "transparent",
    color = NA
  ),
  panel.background = element_rect(
    fill = "transparent",
    color = NA
  ),
  panel.grid = element_blank(),
  panel.border = element_rect(
    fill = "transparent",
    color = NA
  ),
  legend.position = "bottom"
)

text_size_description <- 5

## Create the Normo-CAF and Hypo-CAF panels ----------------------------------

normocaf_switching_plot <- 
  ggplot(
  normocaf_switching_results,
  aes(
    x = .data$NES,
    y = .data$Description
  )
) +
  geom_segment(
    aes(
      x = 0,
      xend = .data$NES,
      yend = .data$Description
    )
  ) +
  geom_point(
    aes(
      size = .data$minus_log10_adjusted_p_value
    ),
    shape = 21,
    fill = "#F03B20"
  ) +
  geom_text(
    aes(
      x = ifelse(.data$NES > 0, -0.2, 0.2),
      label = format_hallmark_label(.data$Description),
      hjust = ifelse(.data$NES > 0, 1, 0)
    ),
    size = text_size_description
  ) +
  common_vertical_line +
  common_scale_x_continuous +
  common_scale_y_discrete +
  point_size_scale +
  common_labs +
  common_coordinate_system +
  common_theme

hypocaf_switching_plot <- ggplot(
  hypocaf_switching_results,
  aes(
    x = .data$NES,
    y = .data$Description
  )
) +
  geom_segment(
    aes(
      x = 0,
      xend = .data$NES,
      yend = .data$Description
    )
  ) +
  geom_point(
    aes(
      size = .data$minus_log10_adjusted_p_value
    ),
    shape = 21,
    fill = "#2C7FB8"
  ) +
  geom_text(
    aes(
      x = ifelse(.data$NES > 0, -0.2, 0.2),
      label = format_hallmark_label(.data$Description),
      hjust = ifelse(.data$NES > 0, 1, 0)
    ),
    size = text_size_description
  ) +
  common_vertical_line +
  common_scale_x_continuous +
  common_scale_y_discrete +
  point_size_scale +
  common_labs +
  common_coordinate_system +
  common_theme


## Save the figures -----------------------------------------------------------

ggsave(
  filename = file.path(
    hallmark_assay_oxygen_condition_figure_dir,
    "bulk_rnaseq_hallmark_switching_of_normocaf.pdf"
  ),
  plot = normocaf_switching_plot,
  device = grDevices::cairo_pdf,
  width = 4.5,
  height = 3.5,
  units = "in",
  bg = "transparent"
)

ggsave(
  filename = file.path(
    hallmark_assay_oxygen_condition_figure_dir,
    "bulk_rnaseq_hallmark_switching_of_hypocaf.pdf"
  ),
  plot = hypocaf_switching_plot,
  device = grDevices::cairo_pdf,
  width = 4.5,
  height = 3.5,
  units = "in",
  bg = "transparent"
)
