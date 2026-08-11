## Hallmark GSEA dot plot for bulk RNA-seq results ----------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(fs)
  library(ggplot2)
  library(readr)
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


## Input and output directories ----------------------------------------------

gene_set_analysis_table_dir <- file.path(
  bulk_rnaseq_table_dir,
  "gene_set_analysis"
)

hallmark_gsea_table_dir <- file.path(
  gene_set_analysis_table_dir,
  "hallmark_gsea"
)

hallmark_gsea_figure_dir <- file.path(
  bulk_rnaseq_figure_dir,
  "08_07_07_bulk_rnaseq_hallmark_gsea_main_dotplot_figures"
)

fs::dir_create(hallmark_gsea_figure_dir)

hallmark_gsea_input_file <- file.path(
  hallmark_gsea_table_dir,
  "bulk_rnaseq_hallmark_gsea_isolation_oxygen_comparisons.csv"
)

if (!file.exists(hallmark_gsea_input_file)) {
  stop(
    "Required Hallmark GSEA result file was not found: ",
    hallmark_gsea_input_file,
    call. = FALSE
  )
}


## Read and validate Hallmark GSEA results -----------------------------------

hallmark_gsea_table <- readr::read_csv(
  hallmark_gsea_input_file,
  show_col_types = FALSE
)

required_columns <- c(
  "comparison_label",
  "Description",
  "NES",
  "p.adjust"
)

missing_columns <- setdiff(
  required_columns,
  colnames(hallmark_gsea_table)
)

if (length(missing_columns) > 0) {
  stop(
    "The Hallmark GSEA table is missing required columns: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

comparison_display_order <- c(
  "Original",
  "Under21",
  "Under1"
)

hallmark_pathway_display_order <- c(
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_MITOTIC_SPINDLE",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  #"HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_COMPLEMENT"#,
  #"HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  #"HALLMARK_TGF_BETA_SIGNALING",
  #"HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  #"HALLMARK_MYOGENESIS"
)

missing_display_pathways <- setdiff(
  hallmark_pathway_display_order,
  unique(hallmark_gsea_table$Description)
)

if (length(missing_display_pathways) > 0) {
  stop(
    "Required Hallmark pathways were not found in the GSEA table: ",
    paste(missing_display_pathways, collapse = ", "),
    call. = FALSE
  )
}


## Prepare plotting table ----------------------------------------------------

hallmark_dotplot_table <- hallmark_gsea_table |>
  dplyr::filter(
    .data$Description %in% hallmark_pathway_display_order
  ) |>
  dplyr::mutate(
    comparison_label = factor(
      .data$comparison_label,
      levels = comparison_display_order
    ),
    Description = factor(
      .data$Description,
      levels = rev(hallmark_pathway_display_order)
    ),
    minus_log10_adjusted_p_value = -log10(.data$p.adjust)
  )

expected_number_of_rows <-
  length(hallmark_pathway_display_order) *
  length(comparison_display_order)

if (nrow(hallmark_dotplot_table) != expected_number_of_rows) {
  stop(
    "The Hallmark GSEA dot-plot table does not contain one result for every ",
    "selected pathway and comparison.",
    call. = FALSE
  )
}

if (any(!is.finite(hallmark_dotplot_table$minus_log10_adjusted_p_value))) {
  stop(
    "Non-finite -log10 adjusted P values were found.",
    call. = FALSE
  )
}

hallmark_nes_limit <- max(
  abs(hallmark_dotplot_table$NES),
  na.rm = TRUE
)

hallmark_pathway_labels <- stats::setNames(
  format_gene_set_label(hallmark_pathway_display_order),
  hallmark_pathway_display_order
)


## Create Hallmark GSEA dot plot ---------------------------------------------

text_size <- 10

hallmark_gsea_dotplot <- ggplot(
  hallmark_dotplot_table,
  aes(
    x = .data$comparison_label,
    y = .data$Description
  )
) +
  geom_point(
    aes(
      size = .data$minus_log10_adjusted_p_value,
      fill = .data$NES
    ),
    shape = 21
  ) +
  labs(
    x = NULL,
    y = NULL,
    size = "−log10\n(Adj.P)",
    fill = "NES"
  ) +
  scale_size_continuous(
    range = c(1.5, 8),
    breaks = c(2, 10, 40),
    guide = guide_legend(
      order = 1,
      direction = "vertical",
      title.position = "top",
      label.position = "left"
    )
  ) +
  scale_fill_distiller(
    palette = "RdBu",
    direction = 1,
    breaks = c(-2, 0, 2),
    limits = c(-hallmark_nes_limit, hallmark_nes_limit),
    guide = guide_colorbar(
      order = 2,
      direction = "vertical",
      title.position = "top",
      label.position = "left",
      reverse = FALSE,
      frame.colour = "black",
      ticks.colour = "black"
    )
  ) +
  scale_x_discrete(
    labels = bulk_rnaseq_comparison_plot_labels
  ) +
  scale_y_discrete(
    position = "right",
    labels = hallmark_pathway_labels
  ) +
  coord_fixed(ratio = 0.8) +
  theme(
    aspect.ratio = 3.880 * 8/13,
    plot.margin = margin(
      t = 0.83, r = 0.5 , b = 0.83, l = 0.5, unit = "in"
    ),
    axis.text.x = element_text(
      face = "plain",
      color = "black",
      angle = 45,
      hjust = 1.0,
      size = text_size * 1.2
    ),
    axis.text.y = element_text(
      face = "plain",
      color = "black",
      size = text_size * 1.3
    ),
    legend.box = "vertical",
    legend.box.just = "right",
    legend.justification = "right",
    legend.key = element_blank(),
    legend.text = element_text(
      face = "plain",
      color = "black",
      size = text_size * 1.2
    ),
    legend.title = element_text(
      face = "plain",
      color = "black",
      size = text_size * 1.2,
      hjust = 0.5
    ),
    legend.background = element_rect(
      fill = "transparent",
      color = NA
    ),
    plot.background = element_rect(
      fill = "transparent",
      color = NA
    ),
    panel.background = element_rect(
      fill = "transparent",
      color = NA
    ),
    panel.border = element_rect(
      fill = NA,
      color = "black"
    ),
    panel.grid = element_blank(),
    axis.line = element_blank()
  )


## Save figure ---------------------------------------------------------------

figure_width <- 6.2
figure_height <- 4.4

ggsave(
  filename = file.path(
    hallmark_gsea_figure_dir,
    "bulk_rnaseq_hallmark_gsea_main_dotplot2.pdf"
  ),
  plot = hallmark_gsea_dotplot,
  device = grDevices::cairo_pdf,
  width = figure_width,
  height = figure_height,
  units = "in",
  bg = "transparent"
)
