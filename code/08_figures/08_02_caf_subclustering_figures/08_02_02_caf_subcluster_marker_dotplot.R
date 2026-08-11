# CAF subcluster marker dot plot
#
# This script generates a dot plot showing representative marker genes
# across CAF subclusters.


## Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(fs)
  library(cowplot)
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

caf_subclustering_table_dir <- file.path(
  xenium_table_dir,
  "caf_subclustering"
)

## Read CAF subclustering object and marker table ----------------------------

CAFObj <- readRDS(caf_object_file)

caf_subcluster_markers_all <- read.csv(
  caf_marker_table_file,
  stringsAsFactors = FALSE
)

stopifnot("caf_subcluster" %in% colnames(CAFObj@meta.data))
stopifnot(all(c("cluster", "gene", "avg_log2FC", "pct.1", "p_val_adj") %in%
                colnames(caf_subcluster_markers_all)))


## Prepare CAF subcluster identities ----------------------------------------

CAFObj$caf_subcluster <- factor(
  as.character(CAFObj$caf_subcluster),
  levels = caf_subcluster_display_order
)

stopifnot(!any(is.na(CAFObj$caf_subcluster)))

Idents(CAFObj) <- CAFObj$caf_subcluster


## Select CAF-8 marker genes -------------------------------------------------

caf8_marker_genes <- caf_subcluster_markers_all %>%
  dplyr::mutate(
    cluster = as.character(cluster)
  ) %>%
  dplyr::filter(
    cluster == proliferating_caf_subcluster,
    avg_log2FC > 1,
    pct.1 > 0.20,
    p_val_adj < 0.05
  ) %>%
  dplyr::arrange(dplyr::desc(avg_log2FC)) %>%
  dplyr::slice_head(n = 5) %>%
  dplyr::pull(gene)


## Define marker genes for dot plot ------------------------------------------

marker_gene_groups <- list(
  "CAF-4" = c("ACTA2", "COL11A1", "THBS2", "CDKN2B", "POSTN"),
  "CAF-8" = caf8_marker_genes,
  "CAF-0" = c("MMP11", "HOPX", "ANO1", "WNT5A", "NOTCH1"),
  "CAF-2" = c("CFD", "PI16", "PLA2G2A", "VEGFA", "CXCL12"),
  "CAF-5" = c("SOD3", "PTCH1", "TGM2", "CFH"),
  "CAF-3" = c("DPT"),
  "CAF-7" = c("CD74", "HLA-DRA", "IRF1", "STAT1"),
  "CAF-6" = c("NGFR", "APOD", "SCN7A", "SOX10")
)

marker_genes <- unname(unlist(marker_gene_groups))

stopifnot(!any(duplicated(marker_genes)))

missing_marker_genes <- setdiff(marker_genes, rownames(CAFObj))
stopifnot(length(missing_marker_genes) == 0)

marker_group_sizes <- lengths(marker_gene_groups)

marker_group_boundaries <- cumsum(marker_group_sizes)
marker_group_boundaries <- marker_group_boundaries[-length(marker_group_boundaries)] + 0.5


## Generate Seurat dot plot data --------------------------------------------

seurat_dotplot <- DotPlot(
  CAFObj,
  features = marker_genes,
  scale = FALSE
)

marker_expression_df <- seurat_dotplot$data %>%
  dplyr::mutate(
    id = factor(
      as.character(id),
      levels = rev(caf_subcluster_display_order)
    ),
    features.plot = factor(
      as.character(features.plot),
      levels = marker_genes
    )
  ) %>%
  dplyr::group_by(features.plot) %>%
  dplyr::mutate(
    expression_range = max(avg.exp, na.rm = TRUE) - min(avg.exp, na.rm = TRUE),
    scaled_expression = dplyr::case_when(
      expression_range == 0 ~ 0,
      TRUE ~ (avg.exp - min(avg.exp, na.rm = TRUE)) / expression_range
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(-expression_range)


## Define CAF subcluster labels and grouping regions -------------------------

caf_subcluster_axis_labels <- c(
  "CAF-4" = "CAF-4. Classical myCAF",
  "CAF-8" = "CAF-8. Proliferating CAF",
  "CAF-0" = "CAF-0. WNT/Notch-high CAF",
  "CAF-1" = "CAF-1. Attenuated myCAF",
  "CAF-2" = "CAF-2. Classical iCAF",
  "CAF-5" = "CAF-5. CFH-high CAF",
  "CAF-3" = "CAF-3. DPT-high CAF",
  "CAF-7" = "CAF-7. apCAF",
  "CAF-6" = "CAF-6. Peri-neural fibroblast"
)

caf_subcluster_label_df <- tibble::tibble(
  caf_subcluster = caf_subcluster_display_order,
  annotation_label = unname(
    caf_subcluster_axis_labels[caf_subcluster_display_order]
  ),
  y_position = match(
    caf_subcluster_display_order,
    rev(caf_subcluster_display_order)
  )
)

stopifnot(!any(is.na(caf_subcluster_label_df$annotation_label)))

cluster_region_half_height <- 0.3

caf_group_df <- tibble::tribble(
  ~caf_group,    ~first_subcluster, ~last_subcluster,
  "mCAF-like",   "CAF-4",           "CAF-1",
  "iCAF-like",   "CAF-2",           "CAF-3"
) %>%
  dplyr::mutate(
    y_min = pmin(
      caf_subcluster_label_df$y_position[
        match(first_subcluster, caf_subcluster_label_df$caf_subcluster)
      ],
      caf_subcluster_label_df$y_position[
        match(last_subcluster, caf_subcluster_label_df$caf_subcluster)
      ]
    ),
    y_max = pmax(
      caf_subcluster_label_df$y_position[
        match(first_subcluster, caf_subcluster_label_df$caf_subcluster)
      ],
      caf_subcluster_label_df$y_position[
        match(last_subcluster, caf_subcluster_label_df$caf_subcluster)
      ]
    ),
    y_mid = (y_min + y_max) / 2,
    line_y_min = y_min - cluster_region_half_height,
    line_y_max = y_max + cluster_region_half_height
  )


## Define figure settings ----------------------------------------------------

text_size <- 12
n_caf_subclusters <- length(caf_subcluster_display_order)

# Relative widths of the four columns.
caf_group_label_column_width <- 1.1
segment_column_width <- 0.2
annotation_label_column_width <- 3.0
dot_plot_column_width <- 12.5

# Horizontal spacing is controlled only by the segment column margins.
segment_margin_left_pt <- 4
segment_margin_right_pt <- 1

# The CAF annotation labels and y-axis ticks are drawn together in one plot.
# These x positions control the gap between the labels and the dot-plot panel.
annotation_text_x <- 0.955
annotation_tick_x_start <- 0.986
annotation_tick_x_end <- 1.00

common_y_scale <- scale_y_continuous(
  limits = c(0.5, n_caf_subclusters + 0.5),
  breaks = seq_len(n_caf_subclusters),
  expand = expansion(mult = 0)
)


## Generate CAF group label column -------------------------------------------

caf_group_label_plot <- ggplot(
  caf_group_df,
  aes(
    x = 1,
    y = y_mid,
    label = caf_group
  )
) +
  geom_text(
    hjust = 1,
    size = text_size * 1.5 / ggplot2::.pt
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
  caf_group_df,
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


## Generate CAF subcluster annotation and tick column -------------------------

caf_subcluster_annotation_plot <- ggplot(
  caf_subcluster_label_df,
  aes(y = y_position)
) +
  geom_text(
    aes(
      x = annotation_text_x,
      label = annotation_label
    ),
    hjust = 1,
    size = text_size * 1.5 / ggplot2::.pt
  ) +
  geom_segment(
    aes(
      x = annotation_tick_x_start,
      xend = annotation_tick_x_end,
      y = y_position,
      yend = y_position
    ),
    linewidth = 0.4,
    lineend = "butt",
    color = "black"
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
      l = 0
    ),
    plot.background = element_rect(
      fill = "transparent",
      color = NA
    )
  )


## Generate CAF subcluster marker dot plot -----------------------------------

dot_plot <- marker_expression_df %>%
  dplyr::mutate(
    y_position = match(
      as.character(id),
      rev(caf_subcluster_display_order)
    )
  ) %>%
  ggplot(
    aes(x = features.plot, y = y_position)
  ) +
  geom_vline(
    xintercept = marker_group_boundaries,
    linewidth = 0.5,
    color = "gray80"
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
  common_y_scale +
  scale_fill_gradientn(
    breaks = c(0, 0.5, 1.0),
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
    axis.ticks.y = element_blank(),
    axis.ticks.length.y = grid::unit(0, "pt"),
    legend.position = "right",
    legend.box.margin = margin(t = 70),
    legend.key = element_blank(),
    axis.text.x = element_text(
      color = "black",
      size = text_size * 1.5,
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    legend.title = element_text(
      color = "black",
      size = text_size * 1.5
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
    legend.background = element_rect(fill = "transparent", color = "transparent"),
    panel.background = element_rect(fill = "transparent", color = "transparent"),
    panel.border = element_rect(fill = "transparent", color = "black"),
    panel.grid = element_blank()
  )


## Combine label columns and marker dot plot ---------------------------------

plot_caf_subcluster_marker_dotplot <- cowplot::plot_grid(
  caf_group_label_plot,
  segment_plot,
  caf_subcluster_annotation_plot,
  dot_plot,
  nrow = 1,
  rel_widths = c(
    caf_group_label_column_width,
    segment_column_width,
    annotation_label_column_width,
    dot_plot_column_width
  ),
  align = "h",
  axis = "tb"
)


## Save figure ---------------------------------------------------------------

caf_subclustering_caf_subcluster_marker_dotplot_dir <- file.path(
  caf_subclustering_figure_dir,
  "08_02_02_caf_subcluster_marker_dotplot"
)

fs::dir_create(caf_subclustering_caf_subcluster_marker_dotplot_dir)

ggsave(
  plot = plot_caf_subcluster_marker_dotplot,
  filename = file.path(
    caf_subclustering_caf_subcluster_marker_dotplot_dir,
    "caf_subcluster_marker_dotplot.pdf"
  ),
  width = 20,
  height = 5.4,
  bg = figure_bg,
  device = grDevices::cairo_pdf
)

