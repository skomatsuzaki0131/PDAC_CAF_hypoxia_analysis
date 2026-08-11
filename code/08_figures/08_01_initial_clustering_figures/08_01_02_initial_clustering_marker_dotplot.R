# Initial clustering marker dot plot (v3)
#
# This script generates a marker dot plot for initial cell type annotation
# using the annotated integrated Xenium Seurat object.
#
# Cell type labels, grouping segments, and cluster labels are drawn as three
# independent columns and combined with the dot plot using cowplot.


## Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(fs)
  library(cowplot)
  library(ggtext)
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


## Define input and output paths --------------------------------------------

initial_clustering_cluster_marker_dotplot_figure_dir <- file.path(
  initial_clustering_figure_dir,
  "08_01_02_initial_clustering_marker_dotplot"
)

fs::dir_create(initial_clustering_cluster_marker_dotplot_figure_dir)

stopifnot(file.exists(initial_object_file))


## Read initial clustering object -------------------------------------------

initial_obj <- readRDS(initial_object_file)

DefaultAssay(initial_obj) <- "RNA"

if (length(Layers(initial_obj[["RNA"]])) > 1) {
  initial_obj <- JoinLayers(initial_obj)
}

stopifnot("initial_cluster" %in% colnames(initial_obj@meta.data))
stopifnot("initial_cell_type" %in% colnames(initial_obj@meta.data))

initial_obj$initial_cluster <- as.character(initial_obj$initial_cluster)
initial_obj$initial_cluster <- stringr::str_remove(
  initial_obj$initial_cluster,
  pattern = "^C"
)


## Define marker genes -------------------------------------------------------

marker_genes <- c(
  "CD24",
  "AMY2A", "CPA1",
  "EPCAM", "CDH1", "SOX9",
  "CFTR","SERPINA5", "FGFR3",
  "KRT19", "MUC1", "TFF1", "MUC5B", 
  "CTSE", "CLDN18",  "MUC5AC",
  "LAMC2", "COL17A1", "FAM83A",
  "INS", "SCG2",
  "COL5A1", "COL3A1", "COL1A1",
  "MYH11", "NOTCH3", "RGS5",
  "PECAM1", "CD34", "FLT4",
  "PTPRC", "CD8A", "GZMA", "CD3E", "CD4",
  "CD19",
  "CD38", "TNFRSF17", "MZB1",
  "CSF1R", "CD68", "ITGAX",
  "LAMP3", "CCR7", "LY75",
  "KIT", "CMA1",
  "SCN7A", "MPZ", "NCAM1"
)

marker_genes <- marker_genes[marker_genes %in% rownames(initial_obj)]

if (length(marker_genes) == 0) {
  stop("No marker genes were found in the initial clustering object.")
}


## Define initial cluster display order and labels ---------------------------

initial_cluster_display_order <- c(
  0, 13, 24,
  1, 12, 33,
  30, 17, 23, 10, 27, 20,
  14,
  7, 6, 9,
  11,
  34, 3,
  2, 16,
  15,
  37, 18,
  5, 4, 31,
  26,
  29
) %>%
  as.character()

initial_cluster_display_order <- initial_cluster_display_order[
  initial_cluster_display_order %in% names(initial_cluster_annotation)
]

initial_cluster_display_order <- initial_cluster_display_order[
  !(initial_cluster_display_order %in%
      as.character(initial_cluster_annotation_list[["Unclassified"]]))
]

cell_type_label_map <- c(
  "Acinar" = "Acinar cell",
  "Ductal_like_acinar" = "Ductal-like acinar cell",
  "Normal_ductal" = "Ductal cell (normal)",
  "PanIN" = "Ductal cell (PanIN)",
  "PDAC" = "Ductal cell (PDAC)",
  "Islet" = "Islet cell",
  "CAF" = "CAF",
  "Mural" = "Mural cell",
  "Endothelial" = "Endothelial cell",
  "Lymph_T" = "T lymphocyte",
  "Lymph_B" = "B lymphocyte",
  "Plasma" = "Plasma cell",
  "Myeloid" = "Myeloid cell",
  "Mast" = "Mast cell",
  "Nerve" = "Neuronal cell"
)

cluster_label_df <- tibble::tibble(
  initial_cluster = initial_cluster_display_order,
  cluster_label = paste0("C", initial_cluster_display_order),
  cell_type = unname(
    initial_cluster_annotation[initial_cluster_display_order]
  ),
  y_position = match(
    initial_cluster_display_order,
    rev(initial_cluster_display_order)
  )
) %>%
  dplyr::mutate(
    cell_type_label = unname(cell_type_label_map[cell_type])
  )

stopifnot(!any(is.na(cluster_label_df$cell_type_label)))


## Define cell type grouping regions -----------------------------------------

cluster_region_half_height <- 0.3

cell_type_group_df <- cluster_label_df %>%
  dplyr::group_by(cell_type, cell_type_label) %>%
  dplyr::summarise(
    y_min = min(y_position),
    y_max = max(y_position),
    y_mid = mean(c(y_min, y_max)),
    line_y_min = y_min - cluster_region_half_height,
    line_y_max = y_max + cluster_region_half_height,
    .groups = "drop"
  )


## Calculate marker expression with Seurat DotPlot ---------------------------

initial_obj_annotated <- subset(
  initial_obj,
  subset = initial_cell_type != "Unclassified"
)

initial_obj_annotated$initial_cluster <- factor(
  initial_obj_annotated$initial_cluster,
  levels = initial_cluster_display_order
)

initial_obj_annotated <- subset(
  initial_obj_annotated,
  subset = !is.na(initial_cluster)
)

Idents(initial_obj_annotated) <- initial_obj_annotated$initial_cluster

marker_genes_cell_type_grouped <- setdiff(
  marker_genes,
  c("CD4", "LAMP3", "CCR7", "LY75")
)

plot_data_seurat <- DotPlot(
  initial_obj_annotated,
  features = marker_genes,
  scale = TRUE
)

marker_expression_df <- plot_data_seurat$data %>%
  dplyr::group_by(features.plot) %>%
  dplyr::mutate(
    scaled_expression = {
      expression_range <- max(avg.exp, na.rm = TRUE) -
        min(avg.exp, na.rm = TRUE)

      if (expression_range == 0) {
        rep(0, dplyr::n())
      } else {
        (avg.exp - min(avg.exp, na.rm = TRUE)) / expression_range
      }
    },
    id = as.character(id),
    y_position = match(
      id,
      rev(initial_cluster_display_order)
    ),
    features.plot = factor(
      features.plot,
      levels = marker_genes
    )
  ) %>%
  dplyr::ungroup()

stopifnot(!any(is.na(marker_expression_df$y_position)))


## Define figure settings ----------------------------------------------------

text_size <- 10
n_initial_clusters <- length(initial_cluster_display_order)

# Relative widths of the four columns.
cell_type_column_width <- 3.7
segment_column_width <- 0.4
cluster_label_column_width <- 0.5
dot_plot_column_width <- 15.5

# Horizontal spacing is controlled only by the segment column margins.
segment_margin_left_pt <- 4
segment_margin_right_pt <- 4

common_y_scale <- scale_y_continuous(
  limits = c(0.5, n_initial_clusters + 0.5),
  breaks = seq_len(n_initial_clusters),
  expand = expansion(mult = 0)
)


## Generate cell type label column -------------------------------------------

cell_type_label_plot <- ggplot(
  cell_type_group_df,
  aes(
    x = 1,
    y = y_mid,
    label = cell_type_label
  )
) +
  geom_text(
    hjust = 1,
    size = text_size * 2 / ggplot2::.pt
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    expand = expansion(mult = 0)
  ) +
  common_y_scale +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(
    plot.margin = margin(
      t = 5.5,
      r = 0,
      b = 5.5,
      l = 5.5
    ),
    plot.background = element_rect(
      fill = "transparent",
      color = NA
    )
  )


## Generate grouping segment column -----------------------------------------

segment_plot <- ggplot(
  cell_type_group_df,
  aes(
    x = 0.5,
    xend = 0.5,
    y = line_y_min,
    yend = line_y_max
  )
) +
  geom_segment(
    linewidth = 0.8,
    lineend = "butt"
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    expand = expansion(mult = 0)
  ) +
  common_y_scale +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(
    plot.margin = margin(
      t = 5.5,
      r = segment_margin_right_pt,
      b = 5.5,
      l = segment_margin_left_pt
    ),
    plot.background = element_rect(
      fill = "transparent",
      color = NA
    )
  )


## Generate cluster label column ---------------------------------------------

cluster_label_plot <- ggplot(
  cluster_label_df,
  aes(
    x = 1,
    y = y_position,
    label = cluster_label
  )
) +
  geom_text(
    hjust = 1,
    size = text_size * 2 / ggplot2::.pt
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    expand = expansion(mult = 0)
  ) +
  common_y_scale +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(
    plot.margin = margin(
      t = 5.5,
      r = 3,
      b = 5.5,
      l = 0
    ),
    plot.background = element_rect(
      fill = "transparent",
      color = NA
    )
  )


## Generate marker dot plot --------------------------------------------------

dot_plot <- marker_expression_df %>%
  ggplot(aes(x = features.plot, y = y_position)) +
  geom_point(
    shape = 21,
    aes(size = pct.exp, fill = scaled_expression)
  ) +
  labs(
    fill = "Scaled\nexpression\n(0-1 per gene)",
    size = "Percent\nexpressed"
  ) +
  common_y_scale +
  scale_fill_gradientn(
    breaks = c(0, 0.5, 1),
    colors = c("gray90", "white", "#fff7bc", "#fec44f", "#d95f0e"),
    guide = guide_colorbar(
      frame.colour = "black",
      ticks.colour = "black"
    )
  ) +
  scale_size(
    range = c(1, text_size),
    guide = guide_legend(
      title.position = "top",
      direction = "vertical"
    )
  ) +
  theme(
    plot.title = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_line(
      color = "black",
      linewidth = 0.4
    ),
    axis.ticks.length.y = grid::unit(3, "pt"),
    legend.position = "right",
    legend.key = element_blank(),
    axis.text.x = element_text(
      color = "black",
      size = text_size * 2,
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    legend.title = element_text(
      color = "black",
      face = "bold",
      size = text_size * 1.5,
      margin = margin(b = 10)
    ),
    legend.text = element_text(
      color = "black",
      size = text_size * 1.3
    ),
    plot.margin = margin(
      t = 5.5,
      r = 5.5,
      b = 5.5,
      l = 0
    ),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    panel.border = element_rect(fill = "transparent", color = "black"),
    panel.grid = element_blank()
  )


## Combine label columns and marker dot plot ---------------------------------

plot_initial_marker_dotplot <- cowplot::plot_grid(
  cell_type_label_plot,
  segment_plot,
  cluster_label_plot,
  dot_plot,
  nrow = 1,
  rel_widths = c(
    cell_type_column_width,
    segment_column_width,
    cluster_label_column_width,
    dot_plot_column_width
  ),
  align = "h",
  axis = "tb"
)


## Save figure ---------------------------------------------------------------

ggsave(
  plot = plot_initial_marker_dotplot,
  filename = file.path(
    initial_clustering_cluster_marker_dotplot_figure_dir,
    "initial_clustering_marker_dotplot_full.pdf"
  ),
  width = 21,
  height = 11,
  bg = figure_bg
)


## Generate cell type–grouped marker dot plot --------------------------------

# Assign each cell to the broader cell type defined by the initial cluster
# annotation. For example, clusters C0, C13, and C24 are merged as Acinar cell.
initial_obj_annotated$initial_cell_type_grouped <- unname(
  initial_cluster_annotation[
    as.character(initial_obj_annotated$initial_cluster)
  ]
)

stopifnot(!any(is.na(initial_obj_annotated$initial_cell_type_grouped)))

# Preserve the biological display order used in the cluster-level dot plot.
cell_type_display_order <- unique(cluster_label_df$cell_type)

initial_obj_annotated$initial_cell_type_grouped <- factor(
  initial_obj_annotated$initial_cell_type_grouped,
  levels = cell_type_display_order
)

Idents(initial_obj_annotated) <- initial_obj_annotated$initial_cell_type_grouped

cell_type_dotplot_seurat <- DotPlot(
  initial_obj_annotated,
  features = marker_genes_cell_type_grouped,
  scale = FALSE
)

cell_type_marker_expression_df <- cell_type_dotplot_seurat$data %>%
  dplyr::group_by(features.plot) %>%
  dplyr::mutate(
    expression_range = max(avg.exp, na.rm = TRUE) -
      min(avg.exp, na.rm = TRUE),
    scaled_expression = dplyr::case_when(
      expression_range == 0 ~ 0,
      TRUE ~ (avg.exp - min(avg.exp, na.rm = TRUE)) / expression_range
    ),
    id = factor(
      as.character(id),
      levels = rev(cell_type_display_order)
    ),
    features.plot = factor(
      as.character(features.plot),
      levels = marker_genes_cell_type_grouped
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(-expression_range)

cell_type_axis_labels <- unname(
  cell_type_label_map[cell_type_display_order]
)

names(cell_type_axis_labels) <- cell_type_display_order

cell_type_axis_labels[c(
  "Normal_ductal",
  "PanIN",
  "PDAC"
)] <- c(
  "Normal ductal cell",
  "PanIN",
  "PDAC"
)

cell_type_axis_labels["CAF"] <- "<b>CAF</b>"

plot_initial_marker_dotplot_by_cell_type <- ggplot(
  cell_type_marker_expression_df,
  aes(x = features.plot, y = id)
) +
  annotate(
    geom = "rect",
    xmin = 21.4,
    xmax = 24.6,
    ymin = 8.45,
    ymax = 9.55,
    fill = "#FFE699",
    alpha = 0.45,
    color = "#C99700",
    linewidth = 0.6
  ) +
  geom_point(
    shape = 21,
    aes(
      size = pct.exp,
      fill = scaled_expression
    )
  ) +
  labs(
    fill = "Scaled\nexpression\n(0-1 per gene)",
    size = "Percent\nexpressed",
    x = NULL,
    y = NULL
  ) +
  scale_y_discrete(
    labels = cell_type_axis_labels
  ) +
  scale_fill_gradientn(
    breaks = c(0, 0.5, 1),
    colors = c("gray90", "white", "#fff7bc", "#fec44f", "#d95f0e"),
    guide = guide_colorbar(
      frame.colour = "black",
      ticks.colour = "black"
    )
  ) +
  scale_size(
    range = c(1, text_size),
    guide = guide_legend(
      title.position = "top",
      direction = "vertical"
    )
  ) +
  theme(
    plot.title = element_blank(),
    axis.title = element_blank(),
    axis.text.y = ggtext::element_markdown(
      color = "black",
      size = text_size * 2
    ),
    axis.ticks.y = element_line(
      color = "black",
      linewidth = 0.4
    ),
    axis.ticks.length.y = grid::unit(3, "pt"),
    legend.position = "right",
    legend.key = element_blank(),
    axis.text.x = element_text(
      color = "black",
      size = text_size * 2,
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    legend.title = element_text(
      color = "black",
      face = "bold",
      size = text_size * 1.5,
      margin = margin(b = 10)
    ),
    legend.text = element_text(
      color = "black",
      size = text_size * 1.3
    ),
    plot.margin = margin(
      t = 5.5,
      r = 5.5,
      b = 5.5,
      l = 5.5
    ),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    panel.border = element_rect(fill = "transparent", color = "black"),
    panel.grid = element_blank()
  )

ggsave(
  plot = plot_initial_marker_dotplot_by_cell_type,
  filename = file.path(
    initial_clustering_cluster_marker_dotplot_figure_dir,
    "initial_clustering_marker_dotplot.pdf"
  ),
  width = 18,
  height = 7,
  bg = figure_bg
)

