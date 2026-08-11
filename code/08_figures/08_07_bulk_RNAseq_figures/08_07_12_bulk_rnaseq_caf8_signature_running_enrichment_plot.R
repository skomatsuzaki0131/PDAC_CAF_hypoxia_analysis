## CAF-8 signature GSEA running enrichment plots for bulk RNA-seq results -----


## Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(enrichplot)
  library(fs)
  library(ggplot2)
  library(patchwork)
  library(purrr)
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


## Input and output files -----------------------------------------------------

caf_signature_gsea_result_object_file <- file.path(
  bulk_rnaseq_object_dir,
  "gene_set_analysis",
  "caf_signature_gsea",
  "bulk_rnaseq_caf_signature_gsea_result_objects.rds"
)

caf8_running_enrichment_figure_dir <- file.path(
  bulk_rnaseq_figure_dir,
  "08_07_12_bulk_rnaseq_caf8_signature_running_enrichment_plot_figures"
)

fs::dir_create(caf8_running_enrichment_figure_dir)

if (!file.exists(caf_signature_gsea_result_object_file)) {
  stop(
    "Required CAF signature GSEA result-object file was not found: ",
    caf_signature_gsea_result_object_file,
    call. = FALSE
  )
}


## Read and validate CAF signature GSEA result objects -----------------------

caf_signature_gsea_result_objects <- readRDS(
  caf_signature_gsea_result_object_file
)

if (!is.list(caf_signature_gsea_result_objects)) {
  stop(
    "The CAF signature GSEA result-object file must contain a named list.",
    call. = FALSE
  )
}

required_comparison_names <- c(
  "comparison_1_hypocaf_hypoxia_vs_normocaf_normoxia",
  "comparison_2_hypocaf_normoxia_vs_normocaf_normoxia",
  "comparison_3_hypocaf_hypoxia_vs_normocaf_hypoxia"
)

missing_comparison_names <- setdiff(
  required_comparison_names,
  names(caf_signature_gsea_result_objects)
)

if (length(missing_comparison_names) > 0) {
  stop(
    "The CAF signature GSEA result-object list is missing comparisons: ",
    paste(missing_comparison_names, collapse = ", "),
    call. = FALSE
  )
}


## Define comparison-specific plot settings ---------------------------------

running_enrichment_plot_settings <- tibble::tribble(
  ~comparison_name,
  ~output_file_stem,
  ~gene_rank_direction_label,

  "comparison_1_hypocaf_hypoxia_vs_normocaf_normoxia",
  "bulk_rnaseq_caf8_signature_running_enrichment_original",
  "Upregulated \u2192 Downregulated",

  "comparison_2_hypocaf_normoxia_vs_normocaf_normoxia",
  "bulk_rnaseq_caf8_signature_running_enrichment_under_21_percent_o2",
  "Upregulated \u2192 Downregulated",

  "comparison_3_hypocaf_hypoxia_vs_normocaf_hypoxia",
  "bulk_rnaseq_caf8_signature_running_enrichment_under_1_percent_o2",
  "Upregulated \u2192 Downregulated"
)

caf8_signature_id <- "CAF8"
caf8_signature_color <- "#2C7FB8"
axis_title_size_x <- 20
axis_text_size_y <- 20
axis_title_size_y <- 20
axis_linewidth <- 1.0

## Create one CAF-8 running enrichment plot ---------------------------------

create_caf8_running_enrichment_plot <- function(
    caf_signature_gsea_result,
    gene_rank_direction_label
) {
  available_gene_sets <- caf_signature_gsea_result@result$ID

  if (!caf8_signature_id %in% available_gene_sets) {
    stop(
      "The CAF-8 signature was not found in the GSEA result. ",
      "Expected gene-set ID: ",
      caf8_signature_id,
      call. = FALSE
    )
  }

  running_enrichment_plot_components <- enrichplot::gseaplot2(
    caf_signature_gsea_result,
    geneSetID = caf8_signature_id,
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
      values = stats::setNames(
        caf8_signature_color,
        caf8_signature_id
      )
    ) +
    theme(
      aspect.ratio = 0.829,
      legend.position = "none",
      axis.title.y = element_text(
        face = "plain",
        hjust = 0.5,
        size = axis_title_size_y
      ),
      axis.text.y = element_text(
        face = "plain",
        size = axis_text_size_y
      ),
      plot.background = element_rect(
        fill = "transparent",
        color = NA
      ),
      panel.background = element_rect(
        fill = "transparent",
        color = NA
      ),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.line.y = element_line(
        color = "black",
        linewidth = axis_linewidth
        ),
      axis.ticks.y = element_line(
        color = "black",
        linewidth = axis_linewidth
      ),
      axis.ticks.length.y = grid::unit(2, "mm"),
      axis.line.x = element_blank()
    )

  # Preserve the running-line thickness used in 08_07_09.
  if ("linewidth" %in% names(
    enrichment_score_plot$layers[[1]]$aes_params
  )) {
    enrichment_score_plot$layers[[1]]$aes_params$linewidth <- 1.3
  } else {
    enrichment_score_plot$layers[[1]]$aes_params$size <- 1.3
  }

  gene_rank_plot <- running_enrichment_plot_components[[2]]
  
  rect_layer_number <- which(
    vapply(
      gene_rank_plot$layers,
      function(layer) inherits(layer$geom, "GeomRect"),
      logical(1)
    )
  )
  
  gene_rank_plot$layers[[rect_layer_number]]$data$col <-
    rev(
      gene_rank_plot$layers[[rect_layer_number]]$data$col
    )
  
  gene_rank_plot <-
    gene_rank_plot +
    geom_hline(
      yintercept = 1, 
      color = "black"
      ) +
    labs(
      x = paste0("Gene rank"),
      color = NULL,
      fill = NULL
    ) +
    theme(
      aspect.ratio = 0.137,
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_text(
        size = axis_title_size_x,
        color = "black"
      ),
      legend.position = "none",
      axis.line = element_line(
        color = "black",
        linewidth = axis_linewidth
      )
    )

  enrichment_score_plot /
    gene_rank_plot &
    theme(
      plot.background = element_rect(
        fill = "transparent",
        color = NA
      )
    )
}


## Create and save the three running enrichment plots ------------------------

figure_width <- 5
figure_height <- 5

purrr::pwalk(
  running_enrichment_plot_settings,
  function(
      comparison_name,
      output_file_stem,
      gene_rank_direction_label
  ) {
    current_caf_signature_gsea_result <-
      caf_signature_gsea_result_objects[[comparison_name]]

    current_running_enrichment_plot <-
      create_caf8_running_enrichment_plot(
        caf_signature_gsea_result = current_caf_signature_gsea_result,
        gene_rank_direction_label = gene_rank_direction_label
      )

    ggsave(
      filename = file.path(
        caf8_running_enrichment_figure_dir,
        paste0(output_file_stem, ".pdf")
      ),
      plot = current_running_enrichment_plot,
      device = grDevices::pdf,
      width = figure_width,
      height = figure_height,
      units = "in",
      bg = "transparent"
    )
  }
)
