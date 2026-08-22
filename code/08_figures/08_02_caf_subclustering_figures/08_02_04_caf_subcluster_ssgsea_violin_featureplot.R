# CAF subcluster ssGSEA violin and feature plots
#
# This script generates violin plots and UMAP feature plots for myCAF and iCAF
# ssGSEA scores calculated in 02_03_caf_subcluster_gene_set_scoring.R.


## Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
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
source(file.path(project_dir, "code", "00_config", "00_06_color_palettes.R"))

caf_subclustering_table_dir <- file.path(
  xenium_table_dir,
  "caf_subclustering"
)

## output directory -----------------------------------------------------------

caf_subclustering_mycaf_icaf_score_figure_dir <- file.path(
  caf_subclustering_figure_dir,
  "08_02_04_caf_subcluster_ssgsea_violin_featureplot"
)


## Read single-cell ssGSEA scores -------------------------------------------

ssgsea_score_file <- file.path(
  caf_subclustering_table_dir,
  "xenium_caf_single_cell_ssgsea_scores.csv"
)

stopifnot(file.exists(ssgsea_score_file))

ssgsea_score_df <- read.csv(
  ssgsea_score_file,
  stringsAsFactors = FALSE
)

required_columns <- c(
  "cell_id",
  "caf_subcluster",
  "caf_umap_1",
  "caf_umap_2",
  "myCAF_ssGSEA",
  "iCAF_ssGSEA"
)

stopifnot(all(required_columns %in% colnames(ssgsea_score_df)))

plot_df <- ssgsea_score_df %>%
  dplyr::transmute(
    cell_id,
    caf_subcluster = factor(
      caf_subcluster,
      levels = caf_subcluster_display_order
    ),
    umap_1 = caf_umap_1,
    umap_2 = caf_umap_2,
    myCAF_ssGSEA,
    iCAF_ssGSEA
  )

stopifnot(!any(is.na(plot_df$caf_subcluster)))


## Prepare colors ------------------------------------------------------------

caf_subcluster_plot_colors <- caf_subcluster_colors_label_named[
  caf_subcluster_display_order
]

stopifnot(identical(
  names(caf_subcluster_plot_colors),
  caf_subcluster_display_order
))
stopifnot(!any(is.na(caf_subcluster_plot_colors)))

feature_score_colors <- viridisLite::magma(14)[
  c(1, 2, 3, 4,  10, 11, 12, 13, 14)
]


feature_score_midpoints <- c(
  "myCAF_ssGSEA" = 0.85,
  "iCAF_ssGSEA" = 0.70
)


## Violin plots --------------------------------------------------------------

plot_violin_score <- function(data, score_column) {
  data %>%
    dplyr::select(caf_subcluster, dplyr::all_of(score_column)) %>%
    rlang::set_names(c("caf_subcluster", "score")) %>%
    ggplot(aes(x = caf_subcluster, y = score, fill = caf_subcluster)) +
    geom_violin(
      trim = TRUE,
      scale = "width",
      width = 0.6
    ) +
    geom_boxplot(
      width = 0.2,
      fill = "white",
      outlier.shape = NA
    ) +
    labs(
      x = "CAF subcluster",
      y = NULL,
      title = NULL
    ) +
    scale_fill_manual(values = caf_subcluster_plot_colors) +
    scale_x_discrete(labels = function(x) str_remove(x, pattern = "CAF-")) +
    coord_cartesian(
      ylim = c(NA, quantile(data[[score_column]], 0.999, na.rm = TRUE))
    ) +
    theme(
      aspect.ratio = 0.6,
      legend.position = "none",
      axis.text.x = element_text(
        color = "black",
        size = 12
      ),
      axis.text.y = element_text(color = "black", size = 10),
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black")
    )
}

plot_violin_mycaf <- plot_violin_score(
  data = plot_df,
  score_column = "myCAF_ssGSEA"
) +
  scale_y_continuous(
    breaks = c(0, 0.3, 0.6)
  )

plot_violin_icaf <- plot_violin_score(
  data = plot_df,
  score_column = "iCAF_ssGSEA"
) +
  scale_y_continuous(
    breaks = c(0.2, 0.4, 0.6)
  )

plot_violin_combined <- ggpubr::ggarrange(
  plot_violin_mycaf,
  plot_violin_icaf,
  ncol = 2,
  nrow = 1,
  align = "hv"
)


## UMAP feature plots --------------------------------------------------------

# Use identical UMAP limits and the default ggplot2 expansion for the PDF and
# dot-only PNG panels.

umap_x_limits <- range(plot_df$umap_1, na.rm = TRUE)
umap_y_limits <- range(plot_df$umap_2, na.rm = TRUE)

umap_scale_x <- scale_x_continuous(
  limits = umap_x_limits
)

umap_scale_y <- scale_y_continuous(
  limits = umap_y_limits
)

umap_panel_aspect_ratio <-
  diff(umap_y_limits) /
  diff(umap_x_limits)

# Approximate panel width measured from the existing 11 x 5 inch PDF.

feature_plot_panel_width <- 4.31

plot_feature_score <- function(data, score_column, midpoint) {
  feature_df <- data %>%
    dplyr::select(umap_1, umap_2, dplyr::all_of(score_column)) %>%
    rlang::set_names(c("umap_1", "umap_2", "score")) %>%
    dplyr::arrange(score)
  
  score_limits <- range(feature_df$score, na.rm = TRUE)
  
  color_scale_values <- c(
    seq(0, midpoint, length.out = 5)[1:4],
    seq(midpoint, 1, length.out = 5)
  )
  
  ggplot(feature_df, aes(x = umap_1, y = umap_2)) +
    ggrastr::rasterise(
      geom_point(
        aes(color = score),
        size = 0.20
      ),
      dpi = figure_dpi,
      dev = "ragg"
    ) +
    labs(
      x = "UMAP 1",
      y = "UMAP 2",
      color = "Score"
    ) +
    umap_scale_x +
    umap_scale_y +
    coord_fixed(ratio = 1) +
    scale_color_gradientn(
      colors = feature_score_colors,
      values = color_scale_values,
      limits = score_limits,
            breaks = dplyr::case_when(
              score_column == "myCAF_ssGSEA" ~ c(0.0, 0.3, 0.6),
              score_column == "iCAF_ssGSEA" ~ c(0.2, 0.4, 0.6)
            ),
            guide = guide_colorbar(
              title.position = "top",
              frame.colour = "black",
              ticks.colour = "white",
              barheight = 3.865,
              barwidth = 0.848
            )
    ) +
    theme(
      legend.position = "right",
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_text(color = "black", size = 20),
      legend.title = element_text(color = "black", size = 15, hjust = 0.5),
      legend.text = element_text(color = "black", size = 15),
      legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
      plot.background = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      plot.margin = margin(0, 1.8, 0, 1.8)
    )
}

plot_feature_score_only_points <- function(data, score_column, midpoint) {
  feature_df <- data %>%
    dplyr::select(umap_1, umap_2, dplyr::all_of(score_column)) %>%
    rlang::set_names(c("umap_1", "umap_2", "score")) %>%
    dplyr::arrange(score)

  score_limits <- range(feature_df$score, na.rm = TRUE)

  color_scale_values <- c(
    seq(0, midpoint, length.out = 5)[1:4],
    seq(midpoint, 1, length.out = 5)
  )

  ggplot(feature_df, aes(x = umap_1, y = umap_2)) +
    geom_point(
      aes(color = score),
      size = 0.20
    ) +
    umap_scale_x +
    umap_scale_y +
    coord_fixed(ratio = 1) +
    scale_color_gradientn(
      colors = feature_score_colors,
      values = color_scale_values,
      limits = score_limits
    ) +
    theme(
      plot.margin = margin(0, 0, 0, 0),
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA),
      panel.grid = element_blank(),
      legend.position = "none",
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.line = element_blank()
    )
}

plot_feature_mycaf <- 
  plot_feature_score(
  data = plot_df,
  score_column = "myCAF_ssGSEA",
  midpoint = feature_score_midpoints[["myCAF_ssGSEA"]]
)

plot_feature_icaf <- 
  plot_feature_score(
  data = plot_df,
  score_column = "iCAF_ssGSEA",
  midpoint = feature_score_midpoints[["iCAF_ssGSEA"]]
)

plot_feature_mycaf_only_points <- plot_feature_score_only_points(
  data = plot_df,
  score_column = "myCAF_ssGSEA",
  midpoint = feature_score_midpoints[["myCAF_ssGSEA"]]
)

plot_feature_icaf_only_points <- plot_feature_score_only_points(
  data = plot_df,
  score_column = "iCAF_ssGSEA",
  midpoint = feature_score_midpoints[["iCAF_ssGSEA"]]
)

plot_feature_combined <- ggpubr::ggarrange(
  plot_feature_mycaf,
  plot_feature_icaf,
  ncol = 2,
  nrow = 1,
  align = "hv"
)


## Save figures --------------------------------------------------------------

fs::dir_create(caf_subclustering_mycaf_icaf_score_figure_dir)

ggsave(
  plot = plot_violin_combined,
  filename = file.path(
    caf_subclustering_mycaf_icaf_score_figure_dir,
    "caf_subcluster_ssgsea_mycaf_icaf_violin.pdf"
  ),
  width = 8,
  height = 3,
  bg = figure_bg
)

ggsave(
  plot = plot_feature_combined,
  filename = file.path(
    caf_subclustering_mycaf_icaf_score_figure_dir,
    "caf_subcluster_ssgsea_mycaf_icaf_featureplot.pdf"
  ),
  width = 11,
  height = 5,
  bg = figure_bg
)

# Dot-only panels for replacement in Inkscape. The width is based on the
# measured panel area in the existing PDF; height preserves the UMAP aspect.

ggsave(
  plot = plot_feature_mycaf_only_points,
  filename = file.path(
    caf_subclustering_mycaf_icaf_score_figure_dir,
    "caf_subcluster_ssgsea_mycaf_featureplot_only_points.png"
  ),
  width = feature_plot_panel_width,
  height = feature_plot_panel_width * umap_panel_aspect_ratio,
  units = "in",
  bg = figure_bg,
  dpi = figure_dpi
)

ggsave(
  plot = plot_feature_icaf_only_points,
  filename = file.path(
    caf_subclustering_mycaf_icaf_score_figure_dir,
    "caf_subcluster_ssgsea_icaf_featureplot_only_points.png"
  ),
  width = feature_plot_panel_width,
  height = feature_plot_panel_width * umap_panel_aspect_ratio,
  units = "in",
  bg = figure_bg,
  dpi = figure_dpi
)


