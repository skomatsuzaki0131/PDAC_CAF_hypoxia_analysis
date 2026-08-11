## Visualization of CAF-8 marker gene expression in bulk RNA-seq -------------


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

differential_expression_table_dir <- file.path(
  bulk_rnaseq_table_dir,
  "differential_expression"
)

comparison_files <- c(
  Isolation = "comparison_1_hypocaf_hypoxia_vs_normocaf_normoxia.csv",
  Under21 = "comparison_2_hypocaf_normoxia_vs_normocaf_normoxia.csv",
  Under1 = "comparison_3_hypocaf_hypoxia_vs_normocaf_hypoxia.csv"
)

caf8_marker_expression_figure_dir <- file.path(
  bulk_rnaseq_figure_dir,
  "08_07_13_caf8_marker_gene_expression"
)

fs::dir_create(caf8_marker_expression_figure_dir)

required_input_files <- c(
  caf8_signature_genes_table_file,
  file.path(
    differential_expression_table_dir,
    unname(comparison_files)
  )
)

missing_input_files <- required_input_files[
  !file.exists(required_input_files)
]

if (length(missing_input_files) > 0) {
  stop(
    "The following required input files were not found: ",
    paste(missing_input_files, collapse = ", "),
    call. = FALSE
  )
}


## Read CAF-8 marker genes ----------------------------------------------------

caf8_marker_table <- readr::read_csv(
  caf8_signature_genes_table_file,
  show_col_types = FALSE
)

required_caf8_marker_columns <- "gene"
missing_caf8_marker_columns <- setdiff(
  required_caf8_marker_columns,
  colnames(caf8_marker_table)
)

if (length(missing_caf8_marker_columns) > 0) {
  stop(
    "The CAF-8 marker table is missing the required column: ",
    paste(missing_caf8_marker_columns, collapse = ", "),
    call. = FALSE
  )
}

caf8_marker_genes <- caf8_marker_table |>
  dplyr::pull(.data$gene) |>
  unique()

if (length(caf8_marker_genes) == 0) {
  stop(
    "No CAF-8 marker genes were found in the CAF-8 marker table.",
    call. = FALSE
  )
}

if (anyNA(caf8_marker_genes) || any(caf8_marker_genes == "")) {
  stop(
    "Missing or empty gene symbols were found in the CAF-8 marker table.",
    call. = FALSE
  )
}


## Read differential expression results -------------------------------------

read_comparison_table <- function(
    comparison_file_name,
    comparison_label
) {
  comparison_table <- readr::read_csv(
    file.path(
      differential_expression_table_dir,
      comparison_file_name
    ),
    show_col_types = FALSE
  )

  required_columns <- c(
    "gene_symbol",
    "logFC",
    "FDR"
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(comparison_table)
  )

  if (length(missing_columns) > 0) {
    stop(
      "The differential expression table for ",
      comparison_label,
      " is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  if (anyDuplicated(comparison_table$gene_symbol)) {
    stop(
      "Duplicated gene symbols were found in the differential expression table for ",
      comparison_label,
      ".",
      call. = FALSE
    )
  }

  comparison_table |>
    dplyr::filter(
      .data$gene_symbol %in% caf8_marker_genes
    ) |>
    dplyr::mutate(
      comparison = comparison_label
    )
}

comparison_tables <- Map(
  read_comparison_table,
  comparison_file_name = unname(comparison_files),
  comparison_label = names(comparison_files)
)

names(comparison_tables) <- names(comparison_files)

if (nrow(comparison_tables[["Isolation"]]) == 0) {
  stop(
    "None of the CAF-8 marker genes were retained in the Isolation comparison.",
    call. = FALSE
  )
}

# The original figure used the gene order from the Isolation comparison.
caf8_gene_display_order <- comparison_tables[["Isolation"]]$gene_symbol

caf8_marker_expression_plot_table <- dplyr::bind_rows(
  comparison_tables
) |>
  dplyr::filter(
    .data$gene_symbol %in% caf8_gene_display_order
  ) |>
  dplyr::mutate(
    gene_symbol = factor(
      .data$gene_symbol,
      levels = rev(caf8_gene_display_order)
    ),
    comparison = factor(
      .data$comparison,
      levels = names(comparison_files)
    ),
    significance = ifelse(
      .data$FDR < 0.05,
      "Signif",
      "ns"
    )
  )

if (anyNA(caf8_marker_expression_plot_table$logFC)) {
  stop(
    "Missing log2 fold-change values were found in the plotting table.",
    call. = FALSE
  )
}

if (anyNA(caf8_marker_expression_plot_table$FDR)) {
  stop(
    "Missing FDR values were found in the plotting table.",
    call. = FALSE
  )
}


## Create the CAF-8 marker gene lollipop plot --------------------------------

caf8_marker_expression_plot <- ggplot(
  caf8_marker_expression_plot_table,
  aes(
    x = .data$logFC,
    y = .data$gene_symbol
  )
) +
  geom_hline(
    yintercept = seq_along(caf8_gene_display_order),
    color = "gray90",
    linewidth = 0.1
  ) +
  geom_point(
    aes(color = .data$significance)
  ) +
  geom_segment(
    aes(
      x = 0,
      xend = .data$logFC,
      y = .data$gene_symbol,
      yend = .data$gene_symbol,
      color = .data$significance
    ),
    linewidth = 0.6
  ) +
  labs(
    x = expression(
      plain(log[2]~"fold change (Hypo-CAF / Normo-CAF)")
    ),
    y = NULL
  ) +
  geom_vline(
    xintercept = 0
  ) +
  scale_x_continuous(
    limits = c(-2.5, 2.5)
  ) +
  scale_color_manual(
    values = c(
      Signif = "black",
      ns = "gray80"
    ),
    labels = c(
      Signif = "FDR < 0.05",
      ns = "n.s."
    )
  ) +
  facet_grid(
    . ~ comparison,
    labeller = labeller(
      comparison = as_labeller(
        c(
          Isolation = "plain(Isolation)",
          Under21 = "plain(21*'%'~O[2])",
          Under1 = "plain(1*'%'~O[2])"
        ),
        default = label_parsed
      )
    )
  ) +
  theme(
    legend.position = "bottom",
    panel.spacing = grid::unit(5, "mm"),
    axis.text.x = element_text(
      face = "plain",
      color = "black",
      size = 11
    ),
    axis.text.y = element_text(
      face = "plain",
      color = "black",
      size = 11
    ),
    axis.title = element_text(
      face = "plain",
      color = "black",
      size = 15
    ),
    legend.title = element_blank(),
    legend.text = element_text(
      face = "plain",
      color = "black",
      size = 12
    ),
    plot.background = element_rect(
      fill = "transparent",
      color = NA
    ),
    strip.background = element_rect(
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
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.border = element_rect(
      fill = "transparent",
      color = NA
    ),
    axis.line = element_line(
      color = "black"
    )
  )


## Save the figure ------------------------------------------------------------

ggsave(
  filename = file.path(
    caf8_marker_expression_figure_dir,
    "bulk_rnaseq_caf8_marker_gene_expression_lollipop.pdf"
  ),
  plot = caf8_marker_expression_plot,
  device = grDevices::cairo_pdf,
  width = 8,
  height = 5,
  units = "in",
  bg = "transparent"
)
