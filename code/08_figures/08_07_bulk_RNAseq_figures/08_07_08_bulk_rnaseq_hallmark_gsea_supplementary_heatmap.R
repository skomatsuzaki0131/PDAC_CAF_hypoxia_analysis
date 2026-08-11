## Hallmark GSEA heatmap for bulk RNA-seq results -----------------------------


## Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(fs)
  library(ggdendro)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(tidyr)
  library(tibble)
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

hallmark_gsea_supplementary_figure_dir <- file.path(
  bulk_rnaseq_figure_dir,
  "08_07_08_bulk_rnaseq_hallmark_gsea_supplementary_heatmap_figures"
)

fs::dir_create(hallmark_gsea_supplementary_figure_dir)

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

unexpected_comparison_labels <- setdiff(
  unique(hallmark_gsea_table$comparison_label),
  comparison_display_order
)

if (length(unexpected_comparison_labels) > 0) {
  stop(
    "Unexpected comparison labels were found: ",
    paste(unexpected_comparison_labels, collapse = ", "),
    call. = FALSE
  )
}

hallmark_heatmap_table <- hallmark_gsea_table |>
  dplyr::select(
    comparison_label,
    Description,
    NES,
    p.adjust
  ) |>
  dplyr::mutate(
    comparison_label = factor(
      .data$comparison_label,
      levels = comparison_display_order
    ),
    significance_label = dplyr::case_when(
      .data$p.adjust < 0.001 ~ "***",
      .data$p.adjust < 0.01 ~ "**",
      .data$p.adjust < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

expected_number_of_rows <-
  length(unique(hallmark_heatmap_table$Description)) *
  length(comparison_display_order)

if (nrow(hallmark_heatmap_table) != expected_number_of_rows) {
  stop(
    "The isolation-oxygen Hallmark GSEA table does not contain one result ",
    "for every pathway and comparison.",
    call. = FALSE
  )
}


## Determine pathway order by hierarchical clustering -----------------------

hallmark_nes_matrix <- hallmark_heatmap_table |>
  dplyr::select(
    comparison_label,
    Description,
    NES
  ) |>
  tidyr::pivot_wider(
    names_from = Description,
    values_from = NES
  ) |>
  tibble::column_to_rownames(var = "comparison_label") |>
  as.matrix()

if (anyNA(hallmark_nes_matrix)) {
  stop(
    "Missing NES values were found in the Hallmark GSEA matrix.",
    call. = FALSE
  )
}

hallmark_distance <- stats::dist(
  t(hallmark_nes_matrix),
  method = "euclidean"
)

hallmark_clustering <- stats::hclust(
  hallmark_distance,
  method = "ward.D2"
)

hallmark_pathway_order <- colnames(hallmark_nes_matrix)[
  hallmark_clustering$order
]

hallmark_heatmap_table <- hallmark_heatmap_table |>
  dplyr::mutate(
    Description = factor(
      .data$Description,
      levels = hallmark_pathway_order
    )
  )

hallmark_nes_limit <- max(
  abs(hallmark_heatmap_table$NES),
  na.rm = TRUE
)

hallmark_pathway_labels <- stats::setNames(
  format_gene_set_label(hallmark_pathway_order),
  hallmark_pathway_order
)


## Create dendrogram ----------------------------------------------------------

hallmark_dendrogram_data <- ggdendro::dendro_data(
  stats::as.dendrogram(hallmark_clustering)
)

hallmark_dendrogram_plot <- ggplot(
  hallmark_dendrogram_data$segments
) +
  geom_segment(
    aes(
      x = .data$y,
      y = .data$x,
      xend = .data$yend,
      yend = .data$xend
    )
  ) +
  scale_y_continuous(
    limits = c(0.5, length(hallmark_clustering$labels) + 0.5),
    breaks = seq_along(hallmark_clustering$labels),
    labels = hallmark_clustering$labels,
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_reverse(
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  theme_void()


## Create Hallmark GSEA heatmap ----------------------------------------------

text_size <- 10

hallmark_heatmap_plot <- ggplot(
  hallmark_heatmap_table,
  aes(
    x = .data$comparison_label,
    y = .data$Description,
    fill = .data$NES
  )
) +
  geom_tile(
    color = "gray50"
  ) +
  geom_text(
    aes(label = .data$significance_label),
    size = text_size * 0.5,
    vjust = 0.8
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = "NES"
  ) +
  scale_x_discrete(
    expand = expansion(mult = c(0, 0)),
    labels = bulk_rnaseq_comparison_plot_labels
  ) +
  scale_y_discrete(
    expand = expansion(mult = c(0, 0)),
    position = "right",
    labels = hallmark_pathway_labels
  ) +
  scale_fill_distiller(
    palette = "RdBu",
    direction = 1,
    limits = c(-hallmark_nes_limit, hallmark_nes_limit),
    breaks = c(-2, 0, 2),
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      frame.colour = "black",
      ticks.colour = "black",
      reverse = FALSE
    )
  ) +
  theme(
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0
    ),
    legend.position = "bottom",
    axis.text.x = element_text(
      face = "plain",
      color = "black",
      angle = 45,
      hjust = 1,
      size = text_size * 0.95
    ),
    axis.text.y = element_text(
      face = "plain",
      color = "black",
      size = text_size * 0.95
    ),
    legend.title = element_text(
      face = "plain",
      color = "black",
      size = text_size * 1.4
    ),
    legend.text = element_text(
      face = "plain",
      color = "black",
      size = text_size * 1.1
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
      fill = "transparent",
      color = "black"
    ),
    panel.grid = element_blank(),
    axis.line = element_blank()
  )

hallmark_gsea_heatmap <-
  hallmark_dendrogram_plot +
  hallmark_heatmap_plot +
  patchwork::plot_layout(widths = c(0.4, 0.2)) &
  theme(
    plot.background = element_rect(
      fill = "transparent",
      color = NA
    ),
    panel.background = element_rect(
      fill = "transparent",
      color = NA
    )
  )


## Save figure ---------------------------------------------------------------

figure_width <- 5.5
figure_height <- 9

ggsave(
  filename = file.path(
    hallmark_gsea_supplementary_figure_dir,
    "bulk_rnaseq_hallmark_gsea_supplementary_heatmap.pdf"
  ),
  plot = hallmark_gsea_heatmap,
  device = grDevices::cairo_pdf,
  width = figure_width,
  height = figure_height,
  units = "in",
  bg = "transparent"
)
