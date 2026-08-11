# Initial clustering marker feature plots
#
# This script generates UMAP feature plots for representative cell type marker
# genes using the integrated Xenium Seurat object.

## Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(fs)
  library(ggrastr)
  library(ragg)
})


## Load settings -------------------------------------------------------------

project_dir <- getwd()

if (basename(project_dir) != "PDAC_CAF_hypoxia_analysis") {
  stop(
    "Please open PDAC_CAF_hypoxia_analysis.Rproj before running this script."
  )
}

source(file.path(project_dir, "code", "00_config", "00_01_analysis_parameters.R"))
source(file.path(project_dir, "code", "00_config", "00_02_paths.R"))


## Read initial clustering data ----------------------------------------------

stopifnot(file.exists(initial_object_file))
stopifnot(file.exists(initial_metadata_file))

initial_obj <- readRDS(initial_object_file)
initial_metadata_df <- readr::read_csv(
  initial_metadata_file,
  show_col_types = FALSE
)

stopifnot(all(c(
  "cell_id",
  "umap_1",
  "umap_2",
  "initial_cell_type"
) %in% colnames(initial_metadata_df)))

stopifnot(all(initial_metadata_df$cell_id %in% colnames(initial_obj)))

DefaultAssay(initial_obj) <- "RNA"

if (inherits(initial_obj[["RNA"]], "Assay5")) {
  initial_obj <- JoinLayers(initial_obj)
}


## Define marker genes -------------------------------------------------------

marker_genes <- c(
  "EPCAM",
  "CD24",
  "AMY2A",
  "CFTR",
  "KRT19",
  "SCG2",
  "COL1A1", "COL3A1", "COL5A1",
  "NOTCH3",
  "CD34",
  "PTPRC",
  "CD3E",
  "CD19",
  "MZB1",
  "ITGAX",
  "KIT"
)

marker_gene_cell_type <- c(
  "EPCAM"   = "Pan-epithelial",
  "CD24"= "Pan-epithelial",
  "AMY2A"  = "Acinar",
  "CFTR"= "Ductal-like acinar",
  "KRT19"  = "Ductal",
  "SCG2"   = "Islet",
  "COL1A1" = "CAF",
  "COL3A1" = "CAF",
  "COL5A1" = "CAF",
  "NOTCH3" = "Mural",
  "CD34"   = "Endothelial",
  "PTPRC"  = "Pan-immune",
  "CD3E"   = "T lymphocyte",
  "CD19"   = "B lymphocyte",
  "MZB1"   = "Plasma",
  "ITGAX"  = "Myeloid",
  "KIT"    = "Mast"
)

marker_genes <- intersect(marker_genes, rownames(initial_obj))

if (length(marker_genes) == 0) {
  stop("No marker genes were found in the initial clustering object.")
}


## Prepare plotting metadata / parameters -------------------------------------

plot_metadata_df <- initial_metadata_df %>%
  dplyr::filter(initial_cell_type != "Unclassified") %>%
  dplyr::filter(cell_id %in% colnames(initial_obj)) %>%
  dplyr::select(
    cell_id,
    umap_1,
    umap_2,
    initial_cell_type
  )

stopifnot(nrow(plot_metadata_df) > 0)

# Use the same UMAP coordinate range and default ggplot2 expansion for both
# the PDF plots and the dot-only PNG panels.
umap_x_limits <- range(plot_metadata_df$umap_1, na.rm = TRUE)
umap_y_limits <- range(plot_metadata_df$umap_2, na.rm = TRUE)
umap_panel_aspect_ratio <-
  diff(umap_y_limits) / diff(umap_x_limits)

# Approximate width of the plotting panel measured from the existing
# 3 x 3 inch feature plot PDF.
feature_plot_panel_width <- 2.22

point_color_scale <- scale_color_viridis_c(
  option = "magma",
  trans = "sqrt",
  breaks = c(0, 2, 4, 6),
  guide = guide_colorbar(
    frame.colour = "black",
    ticks.colour = "black",
    barwidth = unit(4.5, "mm")
  )
)

pt_size_old <- 0.1
pt_shape_old <- 19
pt_alpha_old <- 0.2
pt_size <- 0.3
pt_shape <- 16
pt_alpha <- 0.2


## Extract expression matrix -------------------------------------------------

expression_matrix <- GetAssayData(
  initial_obj,
  assay = "RNA",
  layer = "data"
)


## Plot function -------------------------------------------------------------


plot_feature_gene <- function(target_gene, only_points = FALSE) {
  expression_df <- tibble::tibble(
    cell_id = colnames(initial_obj),
    expression = as.numeric(expression_matrix[target_gene, ])
  )
  
  plot_df <- plot_metadata_df %>%
    dplyr::inner_join(
      expression_df,
      by = "cell_id"
    )
  
  set.seed(123)
  
  plot_df <- plot_df %>%
    dplyr::slice_sample(prop = 1)

  if (only_points) {
    return(
      ggplot(
        plot_df,
        aes(x = umap_1, y = umap_2)
      ) +
        geom_point(
          aes(color = expression),
          size = pt_size,
          shape = pt_shape,
          alpha = pt_alpha
        ) +
        point_color_scale +
        scale_x_continuous(limits = umap_x_limits) +
        scale_y_continuous(limits = umap_y_limits) +
        coord_fixed(ratio = 1) +
        guides(color = "none") +
        theme_void() +
        theme(
          plot.margin = margin(0, 0, 0, 0, unit = "pt"),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent", color = NA)
        )
    )
  }

  ggplot(
    plot_df,
    aes(x = umap_1, y = umap_2)
  ) +
    ggrastr::rasterise(
      geom_point(
        aes(color = expression),
        size = pt_size,
        shape = pt_shape,
        alpha = pt_alpha
      ),
      dpi = figure_dpi,
      dev = "ragg"
    ) +
    point_color_scale +
    scale_x_continuous(limits = umap_x_limits) +
    scale_y_continuous(limits = umap_y_limits) +
    coord_fixed(ratio = 1) +
    labs(
      title = paste0(
        target_gene,
        " / ",
        marker_gene_cell_type[target_gene]
      ),
      color = "Exp.",
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    theme(
      legend.position = "right",
      plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "mm"),
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
      plot.title = element_text(
        face = ifelse(marker_gene_cell_type[target_gene] == "CAF", "bold", "plain"),
        color = "black",
        size = 15,
        margin = margin(b = 0)
      ),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_text(color = "black", size = 13),
      legend.text = element_text(color = "black", size = 11),
      plot.background = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = "transparent", color = "black")
    )

}


## Generate and save feature plots -------------------------------------------

featureplot_output_dir <- file.path(
  initial_clustering_figure_dir,
  "08_01_03_initial_clustering_marker_featureplots"
)

fs::dir_create(featureplot_output_dir)

for (target_gene in marker_genes) {
  feature_plot <- plot_feature_gene(target_gene)
  feature_plot_only_points <- plot_feature_gene(
    target_gene,
    only_points = TRUE
  )

  target_gene_file_label <- stringr::str_to_lower(target_gene)

  ggsave(
    plot = feature_plot,
    filename = file.path(
      featureplot_output_dir,
      paste0(
        "initial_clustering_featureplot_",
        target_gene_file_label,
        ".pdf"
      )
    ),
    width = 3,
    height = 3,
    bg = figure_bg
  )

  ggsave(
    plot = feature_plot_only_points,
    filename = file.path(
      featureplot_output_dir,
      paste0(
        "initial_clustering_featureplot_",
        target_gene_file_label,
        "_only_points.png"
      )
    ),
    width = feature_plot_panel_width,
    height = feature_plot_panel_width * umap_panel_aspect_ratio,
    units = "in",
    bg = figure_bg,
    dpi = figure_dpi
  )
}
