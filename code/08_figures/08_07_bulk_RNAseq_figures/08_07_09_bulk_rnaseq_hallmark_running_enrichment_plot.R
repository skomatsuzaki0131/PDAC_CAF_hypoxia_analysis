## Hallmark GSEA running enrichment plots for bulk RNA-seq results ------------


## Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(enrichplot)
  library(fs)
  library(ggplot2)
  library(patchwork)
  library(stringr)
  library(tibble)
  library(ggrastr)
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

hallmark_gsea_result_object_file <- file.path(
  bulk_rnaseq_object_dir,
  "bulk_rnaseq_hallmark_gsea_result_objects.rds"
)

hallmark_running_enrichment_figure_dir <- file.path(
  bulk_rnaseq_figure_dir,
  "08_07_09_bulk_rnaseq_hallmark_running_enrichment_plot_figures"
)

hallmark_running_enrichment_panel_png_dir <- file.path(
  hallmark_running_enrichment_figure_dir,
  "panel_png"
)

fs::dir_create(hallmark_running_enrichment_figure_dir)
fs::dir_create(hallmark_running_enrichment_panel_png_dir)

if (!file.exists(hallmark_gsea_result_object_file)) {
  stop(
    "Required Hallmark GSEA result-object file was not found: ",
    hallmark_gsea_result_object_file,
    call. = FALSE
  )
}


## Read and validate Hallmark GSEA result objects -----------------------------

hallmark_gsea_result_objects <- readRDS(
  hallmark_gsea_result_object_file
)

if (!is.list(hallmark_gsea_result_objects)) {
  stop(
    "The Hallmark GSEA result-object file must contain a named list.",
    call. = FALSE
  )
}

required_comparison_names <- c(
  "comparison_1_hypocaf_hypoxia_vs_normocaf_normoxia",
  "comparison_2_hypocaf_normoxia_vs_normocaf_normoxia",
  "comparison_3_hypocaf_hypoxia_vs_normocaf_hypoxia",
  "comparison_4_hypocaf_hypoxia_vs_hypocaf_normoxia",
  "comparison_5_normocaf_hypoxia_vs_normocaf_normoxia"
)

missing_comparison_names <- setdiff(
  required_comparison_names,
  names(hallmark_gsea_result_objects)
)

if (length(missing_comparison_names) > 0) {
  stop(
    "The Hallmark GSEA result-object list is missing comparisons: ",
    paste(missing_comparison_names, collapse = ", "),
    call. = FALSE
  )
}


## Define pathway groups and plot settings -----------------------------------

proliferation_isolation_pathways <- c(
  "HALLMARK_E2F_TARGETS" = "#D55E00",
  "HALLMARK_G2M_CHECKPOINT" = "#4E79A7",
  "HALLMARK_MITOTIC_SPINDLE" = "#B3DE69"
)

inflammatory_and_stress_pathways <- c(
  "HALLMARK_INFLAMMATORY_RESPONSE" = "#2B83BA",
  "HALLMARK_P53_PATHWAY" = "#D7191C",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE" = "#FDAE61",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE" = "#006E00"
)

proliferation_assay_oxygen_pathways <- c(
  "HALLMARK_E2F_TARGETS" = "#D55E00",
  "HALLMARK_G2M_CHECKPOINT" = "#4E79A7",
  "HALLMARK_MYC_TARGETS_V2" = "#E69F00",
  "HALLMARK_MYC_TARGETS_V1" = "#7B3294"
)

running_enrichment_plot_settings <- tibble::tribble(
  ~comparison_name,
  ~pathway_group,
  ~output_file_stem,
  ~gene_rank_direction_label,

  "comparison_1_hypocaf_hypoxia_vs_normocaf_normoxia",
  "proliferation_isolation",
  "bulk_rnaseq_hallmark_running_enrichment_original_proliferation",
  "Upregulated \u2192 Downregulated",

  "comparison_1_hypocaf_hypoxia_vs_normocaf_normoxia",
  "inflammatory_and_stress",
  "bulk_rnaseq_hallmark_running_enrichment_original_inflammatory_and_stress",
  "Upregulated \u2192 Downregulated",

  "comparison_2_hypocaf_normoxia_vs_normocaf_normoxia",
  "proliferation_isolation",
  "bulk_rnaseq_hallmark_running_enrichment_under_21_percent_o2_proliferation",
  "Upregulated \u2192 Downregulated",

  "comparison_2_hypocaf_normoxia_vs_normocaf_normoxia",
  "inflammatory_and_stress",
  "bulk_rnaseq_hallmark_running_enrichment_under_21_percent_o2_inflammatory_and_stress",
  "Upregulated \u2192 Downregulated",

  "comparison_3_hypocaf_hypoxia_vs_normocaf_hypoxia",
  "proliferation_isolation",
  "bulk_rnaseq_hallmark_running_enrichment_under_1_percent_o2_proliferation",
  "Upregulated \u2192 Downregulated",

  "comparison_3_hypocaf_hypoxia_vs_normocaf_hypoxia",
  "inflammatory_and_stress",
  "bulk_rnaseq_hallmark_running_enrichment_under_1_percent_o2_inflammatory_and_stress",
  "Upregulated \u2192 Downregulated",

  "comparison_4_hypocaf_hypoxia_vs_hypocaf_normoxia",
  "proliferation_assay_oxygen",
  "bulk_rnaseq_hallmark_running_enrichment_hypocaf_assay_oxygen",
  "21% O2 \u2192 1% O2",

  "comparison_5_normocaf_hypoxia_vs_normocaf_normoxia",
  "proliferation_assay_oxygen",
  "bulk_rnaseq_hallmark_running_enrichment_normocaf_assay_oxygen",
  "21% O2 \u2192 1% O2"
)

hallmark_pathway_groups <- list(
  proliferation_isolation = proliferation_isolation_pathways,
  inflammatory_and_stress = inflammatory_and_stress_pathways,
  proliferation_assay_oxygen = proliferation_assay_oxygen_pathways
)


## Helper functions -----------------------------------------------------------

create_running_enrichment_plot <- function(
    hallmark_gsea_result,
    pathway_colors,
    gene_rank_direction_label
) {
  selected_pathways <- names(pathway_colors)

  available_pathways <- hallmark_gsea_result@result$ID

  missing_pathways <- setdiff(
    selected_pathways,
    available_pathways
  )

  if (length(missing_pathways) > 0) {
    stop(
      "Selected Hallmark pathways were not found in the GSEA result: ",
      paste(missing_pathways, collapse = ", "),
      call. = FALSE
    )
  }

  pathway_labels <- format_gene_set_label(selected_pathways) |>
    stringr::str_replace(
      pattern = "Interferon alpha",
      replacement = "IFN\u03b1"
    ) |>
    stringr::str_replace(
      pattern = "Interferon gamma",
      replacement = "IFN\u03b3"
    ) |>
    stats::setNames(selected_pathways)

  running_enrichment_plot_components <- enrichplot::gseaplot2(
    hallmark_gsea_result,
    geneSetID = selected_pathways,
    subplots = 1:2
  )

  enrichment_score_plot <-
    running_enrichment_plot_components[[1]] +
    labs(
      y = "Enrichment score",
      color = NULL
    ) +
    geom_hline(
      yintercept = 0,
      color = "gray70",
      alpha = 0.5
    ) +
    scale_color_manual(
      values = pathway_colors,
      labels = pathway_labels
    ) +
    theme(
      axis.title.y = element_text(
        face = "plain",
        hjust = 0.5,
        size = 18
      ),
      axis.text.y = element_text(
        face = "plain",
        size = 16
      ),
      legend.text = element_text(
        face = "plain",
        size = 15
      ),
      legend.background = element_rect(
        fill = "transparent",
        color = NA
      ),
      legend.key = element_blank(),
      plot.background = element_rect(
        fill = "transparent",
        color = NA
      ),
      panel.background = element_rect(
        fill = "transparent",
        color = NA
      ),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )

  # Preserve the running-line thickness used in the original figure code.
  if ("linewidth" %in% names(
    enrichment_score_plot$layers[[1]]$aes_params
  )) {
    enrichment_score_plot$layers[[1]]$aes_params$linewidth <- 1.3
  } else {
    enrichment_score_plot$layers[[1]]$aes_params$size <- 1.3
  }

  enrichment_score_plot$data$Description <- factor(
    enrichment_score_plot$data$Description,
    levels = rev(selected_pathways)
  )

  gene_rank_plot <-
    running_enrichment_plot_components[[2]] +
    labs(
      x = paste0(
        "Gene rank\n(",
        gene_rank_direction_label,
        ")"
      ),
      color = NULL
    ) +
    geom_hline(
      yintercept = seq_along(selected_pathways),
      color = "gray80"
    ) +
    scale_color_manual(
      values = pathway_colors,
      labels = pathway_labels
    ) +
    theme(
      legend.position = "none",
      axis.title.x = element_text(
        face = "plain",
        size = 18
      ),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.background = element_rect(
        fill = "transparent",
        color = NA
      ),
      panel.background = element_rect(
        fill = "transparent",
        color = NA
      )
    )
  
  # Panel-only version for a separate PNG. Because all axes, labels,
  # legends, and margins are removed, the saved device size equals the
  # physical panel size measured from the reference PDF.
  gene_rank_panel_plot <-
    gene_rank_plot +
    theme_void() +
    theme(
      legend.position = "none",
      plot.margin = margin(0, 0, 0, 0),
      plot.background = element_rect(
        fill = "transparent",
        color = NA
      ),
      panel.background = element_rect(
        fill = "transparent",
        color = NA
      )
    )

  # Rasterize only the hit-mark layer in the PDF version. Text, axes,
  # horizontal guides, and the running enrichment curves remain vector.
  gene_rank_plot_for_pdf <- ggrastr::rasterise(
    gene_rank_plot,
    layers = "Linerange",
    dpi = figure_dpi
  )

  running_enrichment_panels <-
    enrichment_score_plot /
    gene_rank_plot_for_pdf +
    patchwork::plot_layout(
      heights = c(1.8, 0.3)
    )

  running_enrichment_panels_with_legend <-
    running_enrichment_panels |
    patchwork::guide_area()

  combined_running_enrichment_plot <-
    running_enrichment_panels_with_legend +
    patchwork::plot_layout(
      guides = "collect",
      widths = c(1, 1)
    ) &
    theme(
      plot.background = element_rect(
        fill = "transparent",
        color = NA
      )
    )

  list(
    combined_plot = combined_running_enrichment_plot,
    gene_rank_panel_plot = gene_rank_panel_plot
  )
}


## Create and save running enrichment plots ----------------------------------

figure_width <- 8.5
figure_height <- 5

# Physical dimensions of the hit panel measured from
# bulk_rnaseq_hallmark_running_enrichment_hypocaf_assay_oxygen.pdf.
gene_rank_panel_width <- 3.68222613
gene_rank_panel_height <- 0.59869088

purrr::pwalk(
  running_enrichment_plot_settings,
  function(
      comparison_name,
      pathway_group,
      output_file_stem,
      gene_rank_direction_label
  ) {
    current_hallmark_gsea_result <-
      hallmark_gsea_result_objects[[comparison_name]]

    current_pathway_colors <-
      hallmark_pathway_groups[[pathway_group]]

    current_running_enrichment_plots <-
      create_running_enrichment_plot(
        hallmark_gsea_result = current_hallmark_gsea_result,
        pathway_colors = current_pathway_colors,
        gene_rank_direction_label = gene_rank_direction_label
      )

    ggsave(
      filename = file.path(
        hallmark_running_enrichment_figure_dir,
        paste0(output_file_stem, ".pdf")
      ),
      plot = current_running_enrichment_plots$combined_plot,
      device = grDevices::cairo_pdf,
      width = figure_width,
      height = figure_height,
      units = "in",
      bg = "transparent"
    )

    ggsave(
      filename = file.path(
        hallmark_running_enrichment_panel_png_dir,
        paste0(output_file_stem, "_hit_panel.png")
      ),
      plot = current_running_enrichment_plots$gene_rank_panel_plot,
      width = gene_rank_panel_width,
      height = gene_rank_panel_height,
      units = "in",
      dpi = figure_dpi,
      bg = "transparent"
    )
  }
)
