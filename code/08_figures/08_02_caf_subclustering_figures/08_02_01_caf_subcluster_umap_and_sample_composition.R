# CAF subcluster UMAP and sample composition plots
#
# This script generates UMAP plots for CAF subclusters and stacked bar plots
# showing sample composition within each CAF subcluster.


## Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
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
source(file.path(project_dir, "code", "00_config", "00_06_color_palettes.R"))

caf_subclustering_table_dir <- file.path(
  xenium_table_dir,
  "caf_subclustering"
)

## Output directory ----------------------------------------------------------

caf_subcluster_umap_and_sample_composition_dir <- file.path(
  caf_subclustering_figure_dir, 
  "08_02_01_caf_subcluster_umap_and_sample_composition"
)

#caf_umap_feature_plots_dir <- 
#  file.path(caf_subclustering_figure_dir, "caf_subcluster_umap_feature_plots")


## Read CAF subclustering metadata -------------------------------------------

caf_metadata_df <- read.csv(
  caf_subcluster_metadata_file,
  stringsAsFactors = FALSE
)

stopifnot(all(c(
  "cell_id",
  "xenium_sample_id",
  "caf_subcluster",
  "caf_umap_1",
  "caf_umap_2"
) %in% colnames(caf_metadata_df)))


## Prepare CAF subcluster colors --------------------------------------------

caf_subcluster_plot_colors <- caf_subcluster_colors_label_named[
  caf_subcluster_display_order
]

stopifnot(identical(
  names(caf_subcluster_plot_colors),
  caf_subcluster_display_order
))
stopifnot(!any(is.na(caf_subcluster_plot_colors)))


## Prepare UMAP data ---------------------------------------------------------

umap_df <- caf_metadata_df %>%
  dplyr::transmute(
    cell_id,
    xenium_sample_id,
    caf_subcluster = factor(
      caf_subcluster,
      levels = caf_subcluster_display_order
    ),
    umap_1 = caf_umap_1,
    umap_2 = caf_umap_2
  ) %>%
  dplyr::filter(!is.na(caf_subcluster))

set.seed(123)

umap_df <- umap_df %>%
  dplyr::slice_sample(prop = 1)

umap_center_df <- umap_df %>%
  dplyr::group_by(caf_subcluster) %>%
  dplyr::summarise(
    umap_1 = median(umap_1, na.rm = TRUE),
    umap_2 = median(umap_2, na.rm = TRUE),
    .groups = "drop"
  )

umap_scale_x <- scale_x_continuous(
  limits = range(umap_df$umap_1, na.rm = TRUE)
)

umap_scale_y <- scale_y_continuous(
  limits = range(umap_df$umap_2, na.rm = TRUE)
)

umap_non_proliferating_df <- umap_df %>%
  dplyr::filter(
    caf_subcluster != proliferating_caf_subcluster
  )

umap_proliferating_df <- umap_df %>%
  dplyr::filter(
    caf_subcluster == proliferating_caf_subcluster
  )

common_theme_umap <- theme(
  plot.background = element_rect(fill = "transparent", color = NA),
  panel.background = element_rect(fill = "transparent", color = NA),
  panel.grid = element_blank(),
  legend.position = "none",
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  axis.title = element_text(color = "black"),
  axis.line = element_line(color = "black")
)

common_theme_only_points_umap <- theme(
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


## CAF subcluster UMAP -------------------------------------------------------

plot_caf_subcluster_umap_pdf <- ggplot(
  umap_df,
  aes(x = umap_1, y = umap_2)
) +
  ggrastr::rasterise(
    geom_point(
      data = umap_non_proliferating_df,
      aes(color = caf_subcluster),
      size = 0.6,
      alpha = 0.6,
      shape = 16
    ),
    dpi = figure_dpi,
    dev = "ragg"
  ) +
  ggrastr::rasterise(
    geom_point(
      data = umap_proliferating_df,
      color = caf_subcluster_plot_colors[[proliferating_caf_subcluster]],
      size = 0.9,
      alpha = 0.7,
      shape = 16
    ),
    dpi = figure_dpi,
    dev = "ragg"
  ) +
  geom_text(
    data = umap_center_df,
    aes(label = caf_subcluster),
    color = "black",
    size = 5
  ) +
  labs(
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  umap_scale_x +
  umap_scale_y +
  scale_color_manual(values = caf_subcluster_plot_colors) +
  coord_fixed(ratio = 1) +
  common_theme_umap +
  theme(
    axis.title = element_text(color = "black", size = 18)
  )


umap_panel_aspect_ratio <- 
  diff(range(umap_df$umap_2, na.rm = TRUE)) /  # range y
  diff(range(umap_df$umap_1, na.rm = TRUE))  # range x


plot_caf_subcluster_umap_only_points_png <- ggplot(
  umap_df,
  aes(x = umap_1, y = umap_2)
) +
  geom_point(
    data = umap_non_proliferating_df,
    aes(color = caf_subcluster),
    size = 0.6,
    alpha = 0.6,
    shape = 16
  ) +
  geom_point(
    data = umap_proliferating_df,
    color = caf_subcluster_plot_colors[[proliferating_caf_subcluster]],
    size = 0.9,
    alpha = 0.7,
    shape = 16
  ) +
  umap_scale_x +
  umap_scale_y +
  scale_color_manual(values = caf_subcluster_plot_colors) +
  coord_fixed(ratio = 1) +
  common_theme_only_points_umap


## CAF subcluster highlighted UMAPs -----------------------------------------

plot_highlighted_umap_list <- list()
plot_highlighted_umap_only_points_list <- list()

for (target_subcluster in caf_subcluster_display_order) {
  target_umap_df <- umap_df %>%
    dplyr::filter(caf_subcluster == target_subcluster)

  plot_highlighted_umap <- ggplot(
    umap_df,
    aes(x = umap_1, y = umap_2)
  ) +
    ggrastr::rasterise(
      geom_point(
        color = "gray90",
        size = 0.6,
        alpha = 0.6,
        shape = 16
      ),
      dpi = figure_dpi,
      dev = "ragg"
    ) +
    ggrastr::rasterise(
      geom_point(
        data = target_umap_df,
        color = caf_subcluster_plot_colors[[target_subcluster]],
        size = 0.6,
        alpha = 0.6,
        shape = 16
      ),
      dpi = figure_dpi,
      dev = "ragg"
    ) +
    geom_text(
      data = umap_center_df %>%
        dplyr::filter(caf_subcluster == target_subcluster),
      aes(label = caf_subcluster),
      color = "black",
      size = 4
    ) +
    labs(
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    umap_scale_x +
    umap_scale_y +
    coord_fixed(ratio = 1) +
    common_theme_umap

  plot_highlighted_umap_only_points <- ggplot(
    umap_df,
    aes(x = umap_1, y = umap_2)
  ) +
    geom_point(
      color = "gray90",
      size = 0.6,
      alpha = 0.6,
      shape = 16
    ) +
    geom_point(
      data = target_umap_df,
      color = caf_subcluster_plot_colors[[target_subcluster]],
      size = 0.6,
      alpha = 0.6,
      shape = 16
    ) +
    umap_scale_x +
    umap_scale_y +
    coord_fixed(ratio = 1) +
    common_theme_only_points_umap

  plot_highlighted_umap_list[[target_subcluster]] <-
    plot_highlighted_umap

  plot_highlighted_umap_only_points_list[[target_subcluster]] <-
    plot_highlighted_umap_only_points
}

plot_caf_subcluster_umap_highlighted <- patchwork::wrap_plots(
  plot_highlighted_umap_list,
  ncol = 5,
  nrow = 2
) &
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA)
  )


## Prepare sample composition data ------------------------------------------

cell_table <- umap_df %>%
  dplyr::count(
    caf_subcluster,
    xenium_sample_id,
    name = "cell"
  ) %>%
  tidyr::complete(
    caf_subcluster = caf_subcluster_display_order,
    xenium_sample_id = xenium_sample_id_display_order,
    fill = list(cell = 0)
  ) %>%
  dplyr::mutate(
    caf_subcluster = factor(
      caf_subcluster,
      levels = caf_subcluster_display_order
    ),
    xenium_sample_id = factor(
      xenium_sample_id,
      levels = xenium_sample_id_display_order
    )
  ) %>%
  dplyr::arrange(caf_subcluster, xenium_sample_id)

sample_composition_within_subcluster <- cell_table %>%
  dplyr::group_by(caf_subcluster) %>%
  dplyr::mutate(
    total = sum(cell),
    proportion = cell / total,
    percent_label = paste0(
      formatC(proportion * 100, format = "f", digits = 1),
      "%"
    ),
    cumulative_sum = cumsum(proportion),
    coord_y = 1 - cumulative_sum + proportion / 2
  ) %>%
  dplyr::ungroup()


## CAF subcluster sample composition plot -----------------------------------

text_size_composition <- 10

sample_colors <- c(
  "Xenium_01" = "#6E8FB8",
  "Xenium_02" = "#E5A857",
  "Xenium_03" = "#9FC9C1",
  "Xenium_04" = "#D98585",
  "Xenium_05" = "#8BBF7A",
  "Xenium_06" = "#E9D58A"
)

sample_colors <- sample_colors[xenium_sample_id_display_order]

stopifnot(identical(
  names(sample_colors),
  xenium_sample_id_display_order
))
stopifnot(!any(is.na(sample_colors)))

plot_sample_composition_within_subcluster <- ggplot(
  cell_table,
  aes(x = caf_subcluster, y = cell, fill = xenium_sample_id)
) +
  geom_col(
    position = "fill",
    color = "gray60",
    linewidth = 0.2
  ) +
  geom_text(
    data = sample_composition_within_subcluster,
    aes(
      y = coord_y,
      label = dplyr::case_when(
        proportion > 0.05 ~ percent_label,
        TRUE ~ ""
      )
    ),
    lineheight = 0.8,
    size = text_size_composition * 0.6
  ) +
  geom_text(
    data = sample_composition_within_subcluster %>%
      dplyr::group_by(caf_subcluster) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup(),
    aes(label = formatC(total, big.mark = ",")),
    y = 1.02,
    size = text_size_composition * 0.5
  ) +
  labs(
    y = "Sample proportion",
    x = NULL,
    fill = "Sample ID"
  ) +
  scale_y_continuous(
    labels = scales::percent_format(),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_fill_manual(
    values = sample_colors,
    breaks = xenium_sample_id_display_order
  ) +
  theme(
    axis.text.x = element_text(
      color = "black",
      size = text_size_composition * 1.3
    ),
    axis.text.y = element_text(
      color = "black",
      size = text_size_composition * 1.5
    ),
    axis.title = element_text(
      color = "black",
      size = text_size_composition
    ),
    legend.text = element_text(
      color = "black",
      size = text_size_composition * 1.3
    ),
    legend.title = element_text(
      color = "black",
      size = text_size_composition * 1.3
    ),
    legend.key = element_blank(),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black")
  )


## Save figures --------------------------------------------------------------

fs::dir_create(caf_subcluster_umap_and_sample_composition_dir)

ggsave(
  plot = plot_caf_subcluster_umap_pdf,
  filename = file.path(
    caf_subcluster_umap_and_sample_composition_dir,
    "caf_subcluster_umap.pdf"
  ),
  width = 6,
  height = 6,
  units = "in",
  bg = figure_bg,
  device = grDevices::cairo_pdf
)

ggsave(
  plot = plot_caf_subcluster_umap_only_points_png,
  filename = file.path(
    caf_subcluster_umap_and_sample_composition_dir,
    "caf_subcluster_umap_only_points.png"
  ),
  width = 5.57,
  height = 5.57 * umap_panel_aspect_ratio,
  units = "in",
  bg = figure_bg,
  dpi = figure_dpi
)

ggsave(
  plot = plot_caf_subcluster_umap_highlighted,
  filename = file.path(
    caf_subcluster_umap_and_sample_composition_dir,
    "caf_subcluster_umap_highlighted.pdf"
  ),
  width = 13.3,
  height = 5.3,
  units = "in",
  bg = figure_bg,
  device = grDevices::cairo_pdf
)

for (target_subcluster in caf_subcluster_display_order) {
  target_subcluster_file_label <- target_subcluster %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all("[^a-z0-9]+", "_") %>%
    stringr::str_remove("_$")

  ggsave(
    plot = plot_highlighted_umap_only_points_list[[target_subcluster]],
    filename = file.path(
      caf_subcluster_umap_and_sample_composition_dir,
      paste0(
        "caf_subcluster_umap_highlighted_",
        target_subcluster_file_label,
        "_only_points.png"
      )
    ),
    width = 2.28,
    height = 2.28 * umap_panel_aspect_ratio,
    units = "in",
    bg = figure_bg,
    dpi = figure_dpi
  )
}

ggsave(
  plot = plot_sample_composition_within_subcluster,
  filename = file.path(
    caf_subcluster_umap_and_sample_composition_dir,
    "caf_subcluster_sample_composition.pdf"
  ),
  width = 11,
  height = 5,
  bg = figure_bg
)
