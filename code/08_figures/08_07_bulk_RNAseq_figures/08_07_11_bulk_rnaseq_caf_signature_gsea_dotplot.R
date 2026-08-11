## Visualization of CAF signature GSEA results from bulk RNA-seq -------------


## Load packages -------------------------------------------------------------

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


## Input and output paths -----------------------------------------------------

caf_signature_gsea_result_file <- file.path(
  bulk_rnaseq_table_dir,
  "gene_set_analysis",
  "caf_signature_gsea",
  "bulk_rnaseq_caf_signature_gsea_results.csv"
)

caf_signature_gsea_figure_dir <- file.path(
  bulk_rnaseq_figure_dir,
  "08_07_11_caf_signature_gsea"
)

fs::dir_create(caf_signature_gsea_figure_dir)

if (!file.exists(caf_signature_gsea_result_file)) {
  stop(
    "The CAF signature GSEA result file was not found: ",
    caf_signature_gsea_result_file,
    call. = FALSE
  )
}


## Read and validate GSEA results --------------------------------------------

caf_signature_gsea_result_table <- readr::read_csv(
  caf_signature_gsea_result_file,
  show_col_types = FALSE
)

required_columns <- c(
  "comparison_label",
  "Description",
  "NES",
  "p.adjust",
  "significance"
)

missing_columns <- setdiff(
  required_columns,
  colnames(caf_signature_gsea_result_table)
)

if (length(missing_columns) > 0) {
  stop(
    "The CAF signature GSEA result table is missing required columns: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

caf_signature_order <- c("CAF8", "myCAF", "iCAF", "apCAF")
comparison_order <- c("Original", "Under21", "Under1")

unexpected_signatures <- setdiff(
  unique(caf_signature_gsea_result_table$Description),
  caf_signature_order
)

if (length(unexpected_signatures) > 0) {
  stop(
    "Unexpected CAF signatures were found: ",
    paste(unexpected_signatures, collapse = ", "),
    call. = FALSE
  )
}

unexpected_comparisons <- setdiff(
  unique(caf_signature_gsea_result_table$comparison_label),
  comparison_order
)

if (length(unexpected_comparisons) > 0) {
  stop(
    "Unexpected comparison labels were found: ",
    paste(unexpected_comparisons, collapse = ", "),
    call. = FALSE
  )
}

if (anyNA(caf_signature_gsea_result_table$NES)) {
  stop(
    "Missing NES values were found in the CAF signature GSEA results.",
    call. = FALSE
  )
}

if (anyNA(caf_signature_gsea_result_table$p.adjust)) {
  stop(
    "Missing adjusted P values were found in the CAF signature GSEA results.",
    call. = FALSE
  )
}

expected_combinations <- tidyr::expand_grid(
  comparison_label = comparison_order,
  Description = caf_signature_order
)

observed_combinations <- caf_signature_gsea_result_table |>
  dplyr::distinct(
    .data$comparison_label,
    .data$Description
  )

missing_combinations <- dplyr::anti_join(
  expected_combinations,
  observed_combinations,
  by = c("comparison_label", "Description")
)

if (nrow(missing_combinations) > 0) {
  stop(
    "Some expected comparison-signature combinations are missing from the CAF signature GSEA results.",
    call. = FALSE
  )
}

caf_signature_gsea_plot_table <- caf_signature_gsea_result_table |>
  dplyr::transmute(
    Condition = factor(
      .data$comparison_label,
      levels = comparison_order
    ),
    Description = factor(
      .data$Description,
      levels = rev(caf_signature_order)
    ),
    NES = .data$NES,
    p.adjust = .data$p.adjust,
    Signif = .data$significance
  ) %>% 
  dplyr::mutate(
    minuslog10p = (-1)*log10(p.adjust)
  )

nes_range <- range(caf_signature_gsea_plot_table$NES, na.rm = TRUE)
maximum_absolute_nes <- max(abs(nes_range))


## Create the CAF signature GSEA tile plot -----------------------------------

caf_signature_gsea_plot <- 
  caf_signature_gsea_plot_table |>
  ggplot(
    aes(
      x = .data$Condition,
      y = .data$Description
    )
  ) +
  geom_point(
    aes(
      fill = .data$NES,
      size = .data$minuslog10p),
    color = "black",
    shape = 21
  ) +
  #geom_text(
  #  aes(
  #    label = .data$Signif,
  #    vjust = ifelse(.data$p.adjust < 0.05, 0.9, 0.4),
  #    size = ifelse(.data$p.adjust < 0.05, 8.0, 5.5)
  #  ),
  #  fontface = "plain"
  #) +
  labs(
    x = NULL,
    y = NULL,
    size = "−log10\n(adj.P)"
  ) +
  scale_x_discrete(
    #expand = expansion(mult = c(0, 0)),
    labels = function(x) {
      label_map <- c(
        Original = "plain('Isolation')",
        Under21 = "plain('21% ' * O[2])",
        Under1 = "plain('1% ' * O[2])"
      )
      parse(text = unname(label_map[x]))
    }
  ) +
  scale_y_discrete(
    #expand = expansion(mult = c(0, 0)),
    labels = c(CAF8 = "CAF-8"),
    position = "right"
  ) +
  scale_fill_distiller(
    palette = "RdBu",
    direction = 1,
    limits = c(-maximum_absolute_nes, maximum_absolute_nes),
    breaks = c(2, 0, -2),
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      frame.colour = "black",
      ticks.colour = "black",
      label.position = "left"
    )
  ) +
  scale_size(
    breaks = c(1, 4, 8),
    range = c(1, 10)
    ) +
  theme(
    aspect.ratio = 1.1936,
    plot.margin = margin(0.5, 0,  0.5, 0, unit = "in"),
    legend.position = "left",
    axis.text.x = element_text(
      face = "plain",
      color = "black",
      angle = 45,
      hjust = 1,
      #vjust = 0,
      size = 13
    ),
    axis.text.y = element_text(
      face = "plain",
      color = "black",
      size = 13
    ),
    legend.title = element_text(
      face = "plain",
      color = "black",
      size = 12
    ),
    legend.text = element_text(
      face = "plain",
      color = "black",
      size = 10
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
      color = "black"
    )
  )


## Save the figure ------------------------------------------------------------

ggsave(
  filename = file.path(
    caf_signature_gsea_figure_dir,
    "bulk_rnaseq_caf_signature_gsea_dotplot.pdf"
  ),
  plot = caf_signature_gsea_plot,
  device = grDevices::cairo_pdf,
  width = 4,
  height = 3,
  units = "in",
  bg = "transparent"
)
