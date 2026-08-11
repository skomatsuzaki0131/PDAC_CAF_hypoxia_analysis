# Initial clustering UMAP and cell type composition plots
#
# This script generates UMAP plots for initial clustering and a stacked bar plot
# showing cell type composition within each Xenium sample.
# The main UMAP excludes unclassified clusters, whereas the supplementary UMAP
# includes all clusters and highlights the centers of unclassified clusters.


## Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(fs)
  library(ggrepel)
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
source(file.path(project_dir, "code", "00_config", "00_04_annotation_definitions.R"))
source(file.path(project_dir, "code", "00_config", "00_06_color_palettes.R"))


## Read initial clustering metadata -----------------------------------------

initial_metadata_df <- read.csv(
  initial_metadata_file,
  stringsAsFactors = FALSE
)


## Prepare plotting data -----------------------------------------------------

umap_df <- initial_metadata_df %>%
  dplyr::mutate(
    initial_cluster = as.character(initial_cluster)
  )

set.seed(123)

umap_df <- umap_df %>%
  dplyr::slice_sample(prop = 1)

umap_df_without_unclassified <- umap_df %>%
  dplyr::filter(!(initial_cell_type == "Unclassified"))

umap_center_df <- umap_df %>%
  dplyr::group_by(initial_cluster) %>%
  dplyr::summarise(
    umap_1 = median(umap_1, na.rm = TRUE),
    umap_2 = median(umap_2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    initial_cell_type = initial_cluster_annotation[as.character(initial_cluster)],
    cluster_label = paste0("C", initial_cluster)
  )


## Common plot settings ------------------------------------------------------

common_theme_initial_cluster <- theme(
  plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
  legend.position = "none",
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  axis.title = element_text(color = "black", size = 22),
  plot.background = element_rect(fill = "transparent", color = NA),
  panel.background = element_rect(fill = "transparent", color = NA),
  panel.grid = element_blank(),
  panel.border = element_rect(fill = "transparent", color = "black")
)

pt_size <- 0.51
pt_shape <- 16
pt_alpha <- 0.3

## UMAP coordinate range and panel aspect ratio --------------------------------

umap_x_limits <- range(umap_df$umap_1, na.rm = TRUE)
umap_y_limits <- range(umap_df$umap_2, na.rm = TRUE)

umap_panel_aspect_ratio <-
  diff(umap_y_limits) / diff(umap_x_limits)

# Approximate width of the UMAP panel area in the existing 6 x 6 inch PDFs.
initial_umap_panel_width <- 5.57

common_theme_initial_cluster_only_points <- theme(
  plot.margin = margin(0, 0, 0, 0, unit = "pt"),
  legend.position = "none",
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  axis.title = element_blank(),
  plot.background = element_rect(fill = "transparent", color = NA),
  panel.background = element_rect(fill = "transparent", color = NA),
  panel.grid = element_blank(),
  panel.border = element_blank()
)


## Main UMAP excluding unclassified clusters ---------------------------------

plot_initial_umap_without_unclassified_cells <- ggplot(
  umap_df_without_unclassified,
  aes(x = umap_1, y = umap_2)
) +
  ggrastr::rasterise(
    geom_point(
      aes(color = initial_cluster),
      size = pt_size,
      shape = pt_shape,
      alpha = pt_alpha
    ),
    dpi = figure_dpi,
    dev = "ragg"
  ) +
  labs(
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  scale_color_manual(values = initial_cluster_colors) +
  coord_fixed(ratio = 1) +
  common_theme_initial_cluster


## Supplementary UMAP excluding unclassified clusters with cluster labels -----

plot_initial_umap_without_unclassified_cells_with_cluster_labels <- ggplot(
  umap_df_without_unclassified,
  aes(x = umap_1, y = umap_2)
) +
  ggrastr::rasterise(
    geom_point(
      aes(color = initial_cluster),
      size = pt_size,
      shape = pt_shape,
      alpha = pt_alpha
    ),
    dpi = figure_dpi,
    dev = "ragg"
  ) +
  geom_point(
    data = umap_center_df %>%
      dplyr::filter(!(initial_cell_type == "Unclassified")),
    color = "black",
    size = 0.6
  ) +
  ggrepel::geom_text_repel(
    data = umap_center_df %>%
      dplyr::filter(!(initial_cell_type == "Unclassified")),
    aes(
      x = umap_1,
      y = umap_2,
      label = str_remove(cluster_label, pattern = "^C")
    ),
    inherit.aes = FALSE,
    color = "black",
    size = 7,
    box.padding = 0.25,
    point.padding = 0.1,
    min.segment.length = 0.5,
    segment.color = "black",
    segment.size = 0.3,
    max.overlaps = Inf,
    seed = 123
  ) +
  labs(
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  scale_color_manual(values = initial_cluster_colors) +
  coord_fixed(ratio = 1) +
  common_theme_initial_cluster


## Supplementary UMAP including unclassified clusters ------------------------

plot_initial_umap_all_qc_cells <- ggplot(
  umap_df,
  aes(x = umap_1, y = umap_2)
) +
  ggrastr::rasterise(
    geom_point(
      aes(color = initial_cluster),
      size = pt_size,
      shape = pt_shape,
      alpha = pt_alpha
    ),
    dpi = figure_dpi,
    dev = "ragg"
  ) +
  geom_point(
    data = umap_center_df %>%
      dplyr::filter(initial_cell_type == "Unclassified"),
    color = "black",
    size = 0.6
  ) +
  ggrepel::geom_text_repel(
    data = umap_center_df %>%
      dplyr::filter(initial_cell_type == "Unclassified"),
    aes(
      x = umap_1,
      y = umap_2,
      label = str_remove(cluster_label, pattern = "^C")
    ),
    inherit.aes = FALSE,
    color = "black",
    size = 7,
    box.padding = 0.25,
    point.padding = 0.1,
    min.segment.length = 0.5,
    segment.color = "black",
    segment.size = 0.3,
    max.overlaps = Inf,
    seed = 123
  ) +
  labs(
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  scale_color_manual(values = initial_cluster_colors) +
  coord_fixed(ratio = 1) +
  common_theme_initial_cluster


## Dot-only UMAP panels for figure assembly -----------------------------------

plot_initial_umap_without_unclassified_cells_only_points <- ggplot(
  umap_df_without_unclassified,
  aes(x = umap_1, y = umap_2)
) +
  ggrastr::rasterise(
    geom_point(
      aes(color = initial_cluster),
      size = pt_size,
      shape = pt_shape,
      alpha = pt_alpha
    ),
    dpi = figure_dpi,
    dev = "ragg"
  ) +
  scale_color_manual(values = initial_cluster_colors) +
  scale_x_continuous(limits = umap_x_limits) +
  scale_y_continuous(limits = umap_y_limits) +
  coord_fixed(ratio = 1) +
  common_theme_initial_cluster_only_points

plot_initial_umap_without_unclassified_cells_with_cluster_centers_only_points <- ggplot(
  umap_df_without_unclassified,
  aes(x = umap_1, y = umap_2)
) +
  ggrastr::rasterise(
    geom_point(
      aes(color = initial_cluster),
      size = pt_size,
      shape = pt_shape,
      alpha = pt_alpha
    ),
    dpi = figure_dpi,
    dev = "ragg"
  ) +
  geom_point(
    data = umap_center_df %>%
      dplyr::filter(!(
        initial_cell_type == "Unclassified")
        ),
    aes(x = umap_1, y = umap_2),
    inherit.aes = FALSE,
    color = "black",
    size = 0.6
  ) +
  scale_color_manual(values = initial_cluster_colors) +
  scale_x_continuous(limits = umap_x_limits) +
  scale_y_continuous(limits = umap_y_limits) +
  coord_fixed(ratio = 1) +
  common_theme_initial_cluster_only_points

plot_initial_umap_all_qc_cells_with_unclassified_centers_only_points <- ggplot(
  umap_df,
  aes(x = umap_1, y = umap_2)
) +
  ggrastr::rasterise(
    geom_point(
      aes(color = initial_cluster),
      size = pt_size,
      shape = pt_shape,
      alpha = pt_alpha
    ),
    dpi = figure_dpi,
    dev = "ragg"
  ) +
  geom_point(
    data = umap_center_df %>%
      dplyr::filter(initial_cell_type == "Unclassified"),
    aes(x = umap_1, y = umap_2),
    inherit.aes = FALSE,
    color = "black",
    size = 0.6
  ) +
  scale_color_manual(values = initial_cluster_colors) +
  scale_x_continuous(limits = umap_x_limits) +
  scale_y_continuous(limits = umap_y_limits) +
  coord_fixed(ratio = 1) +
  common_theme_initial_cluster_only_points


## Cell type composition within each sample ----------------------------------

celltype_group_labels <- c(
  "Acinar" = "Normal_acinar",
  "Ductal_like_acinar" = "Ductal",
  "Normal_ductal" = "Ductal",
  "PanIN" = "Ductal",
  "PDAC" = "Ductal",
  "Islet" = "Islet",
  "CAF" = "CAF",
  "Mural" = "Endothelial_mural",
  "Endothelial" = "Endothelial_mural",
  "Lymph_B" = "Lymphoid",
  "Lymph_T" = "Lymphoid",
  "Plasma" = "Lymphoid",
  "Myeloid" = "Myeloid",
  "Mast" = "Myeloid",
  "Nerve" = "Nerve",
  "Unclassified" = "Unclassified"
)

cell_count_by_sample_cluster <- table(
  initial_metadata_df$xenium_sample_id,
  initial_metadata_df$initial_cluster
) %>%
  as.data.frame() %>%
  magrittr::set_colnames(c("xenium_sample_id", "initial_cluster", "n_cell")) %>%
  dplyr::mutate(
    initial_cluster = as.character(initial_cluster),
    cell_type = initial_cluster_annotation[initial_cluster],
    celltype_group = celltype_group_labels[cell_type],
    celltype_group = factor(
      celltype_group,
      levels = unique(unname(celltype_group_labels))
    )
  )

cell_count_by_sample_celltype_group <- cell_count_by_sample_cluster %>%
  dplyr::group_by(xenium_sample_id, celltype_group) %>%
  dplyr::summarise(
    n_cell = sum(n_cell),
    .groups = "drop"
  )

plot_celltype_composition_within_sample <- ggplot(
  cell_count_by_sample_celltype_group,
  aes(x = xenium_sample_id, y = n_cell, fill = celltype_group)
) +
  geom_col(
    position = "fill",
    color = "gray30",
    linewidth = 0.3
  ) +
  labs(
    y = "Cell type proportion\nwithin sample (%)",
    x = NULL,
    fill = NULL
  ) +
  scale_x_discrete(
    expand = expansion(mult = c(0.12, 0.12)),
    labels = function(x) str_remove(x,pattern = "Xenium_")
  ) +
  scale_y_continuous(
    labels = function(x) x * 100,
    expand = expansion(mult = 0)
  ) +
  scale_fill_manual(
    values = c(
      "Normal_acinar" = "#ff00ff",
      "Ductal" = "#ddbcff",
      "Islet" = "#ffff00",
      "CAF" = "#00ff00",
      "Endothelial_mural" = "#990000",
      "Lymphoid" = "#00bfff",
      "Myeloid" = "#ff8800",
      "Nerve" = "#0000ff",
      "Unclassified" = "#D0D0D0"
    ),
    labels = c(
      "Normal_acinar" = "Normal acinar",
      "Ductal" = "Ductal",
      "Islet" = "Islet",
      "CAF" = "CAF",
      "Endothelial_mural" = "Endothelial and mural",
      "Lymphoid" = "Lymphoid",
      "Myeloid" = "Myeloid",
      "Nerve" = "Neuronal",
      "Unclassified" = "Unclassified"
    )
  ) +
  coord_fixed(ratio = 5.5) +
  theme(
    plot.margin = unit(c(0.5, 0.1, 0.1, 0.1), "cm"),
    legend.position = "right",
    legend.key = element_blank(),
    legend.text = element_text(color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 16),
    axis.text.x = element_text(color = "black", size = 16),
    axis.title = element_text(color = "black", size = 17),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
    panel.grid = element_blank()
  )

## Save figures --------------------------------------------------------------

initial_clustering_umap_and_cell_type_composition_figure_dir <- file.path(
  initial_clustering_figure_dir,
  "08_01_01_initial_clustering_umap_and_cell_type_composition"
)

fs::dir_create(initial_clustering_umap_and_cell_type_composition_figure_dir)

ggsave(
  plot = plot_initial_umap_without_unclassified_cells,
  filename = file.path(
    initial_clustering_umap_and_cell_type_composition_figure_dir,
    "xenium_initial_umap_without_unclassified_cells.pdf"
  ),
  width = 6,
  height = 6,
  bg = figure_bg
)

ggsave(
  plot = plot_initial_umap_without_unclassified_cells_with_cluster_labels,
  filename = file.path(
    initial_clustering_umap_and_cell_type_composition_figure_dir,
    "xenium_initial_umap_without_unclassified_cells_with_cluster_labels.pdf"
  ),
  width = 6,
  height = 6,
  bg = figure_bg
)

ggsave(
  plot = plot_initial_umap_all_qc_cells,
  filename = file.path(
    initial_clustering_umap_and_cell_type_composition_figure_dir,
    "xenium_initial_umap_all_qc_cells.pdf"
  ),
  width = 6,
  height = 6,
  bg = figure_bg
)

# Dot-only PNGs for replacement in Inkscape.
# The two labeled PDF variants retain their black cluster-center points here.
ggsave(
  plot = plot_initial_umap_without_unclassified_cells_only_points,
  filename = file.path(
    initial_clustering_umap_and_cell_type_composition_figure_dir,
    "xenium_initial_umap_without_unclassified_cells_only_points.png"
  ),
  width = initial_umap_panel_width,
  height = initial_umap_panel_width * umap_panel_aspect_ratio,
  units = "in",
  bg = figure_bg,
  dpi = figure_dpi
)

ggsave(
  plot = plot_initial_umap_without_unclassified_cells_with_cluster_centers_only_points,
  filename = file.path(
    initial_clustering_umap_and_cell_type_composition_figure_dir,
    paste0(
      "xenium_initial_umap_without_unclassified_cells_",
      "with_cluster_centers_only_points.png"
    )
  ),
  width = initial_umap_panel_width,
  height = initial_umap_panel_width * umap_panel_aspect_ratio,
  units = "in",
  bg = figure_bg,
  dpi = figure_dpi
)

ggsave(
  plot = plot_initial_umap_all_qc_cells_with_unclassified_centers_only_points,
  filename = file.path(
    initial_clustering_umap_and_cell_type_composition_figure_dir,
    paste0(
      "xenium_initial_umap_all_qc_cells_",
      "with_unclassified_centers_only_points.png"
    )
  ),
  width = initial_umap_panel_width,
  height = initial_umap_panel_width * umap_panel_aspect_ratio,
  units = "in",
  bg = figure_bg,
  dpi = figure_dpi
)


ggsave(
  plot = plot_celltype_composition_within_sample,
  filename = file.path(
    initial_clustering_umap_and_cell_type_composition_figure_dir,
    "xenium_initial_celltype_composition_within_sample.pdf"
  ),
  width = 6,
  height = 4,
  bg = figure_bg
)
