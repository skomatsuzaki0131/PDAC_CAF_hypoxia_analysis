# Pseudotime trajectory figures
#
# This script generates the trajectory-coordinate plots included in the
# original pseudotime plotting code:
#   1. Pseudotime state
#   2. CAF subcluster
#   3. Pseudotime
#   4. Spatial niche
#
# The trajectory coordinates and minimum-spanning-tree edges are read from the
# tables exported by 05_02_pseudotime_extract_metadata_and_state_composition.R.
#
# Each plot is saved as:
#   - PDF with axes and legend
#   - panel-only PNG without axes or legend


## Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(viridis)
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
source(file.path(project_dir, "code", "00_config", "00_04_annotation_definitions.R"))
source(file.path(project_dir, "code", "00_config", "00_06_color_palettes.R"))


## Define input and output directories ---------------------------------------

pseudotime_table_dir <- file.path(
  xenium_table_dir,
  "pseudotime_analysis"
)

pseudotime_figure_dir <- file.path(
  xenium_figure_dir,
  "08_05_pseudotime_analysis"
)

pseudotime_trajectory_figure_dir <- file.path(
  pseudotime_figure_dir,
  "08_05_01_pseudotime_trajectory"
)

pseudotime_trajectory_panel_png_dir <- file.path(
  pseudotime_trajectory_figure_dir,
  "panel_png"
)


fs::dir_create(pseudotime_trajectory_figure_dir)
fs::dir_create(pseudotime_trajectory_panel_png_dir)


## Define input files --------------------------------------------------------

pseudotime_metadata_file <- file.path(
  pseudotime_table_dir,
  "pseudotime_analysis_cell_metadata.csv"
)

trajectory_edge_file <- file.path(
  pseudotime_table_dir,
  "pseudotime_analysis_trajectory_edges.csv"
)

pseudotime_niche_metadata_file <- file.path(
  pseudotime_table_dir,
  "pseudotime_analysis_cell_metadata_with_niche.csv"
)

required_input_files <- c(
  pseudotime_metadata_file,
  trajectory_edge_file,
  pseudotime_niche_metadata_file
)

missing_input_files <- required_input_files[
  !file.exists(required_input_files)
]

if (length(missing_input_files) > 0) {
  stop(
    "Required input files were not found:\n",
    paste(missing_input_files, collapse = "\n"),
    call. = FALSE
  )
}


## Figure settings -----------------------------------------------------------

# Match the physical size allocated to each trajectory plot in the original
# combined figure: left-column width = 5 inches and row height = 10 / 3 inches.
figure_width <- 5.00
figure_height <- 3.33

panel_width <- 3.36
panel_height <- 2.81

caf_subcluster_colors <- stats::setNames(
  unname(caf_subcluster_colors_number_named),
  paste0("CAF-", names(caf_subcluster_colors_number_named))
)


## Read and validate input tables -------------------------------------------

pseudotime_metadata_table <- readr::read_csv(
  pseudotime_metadata_file,
  show_col_types = FALSE
)

trajectory_edge_table <- readr::read_csv(
  trajectory_edge_file,
  show_col_types = FALSE
)

pseudotime_niche_metadata_table <- readr::read_csv(
  pseudotime_niche_metadata_file,
  show_col_types = FALSE
)

required_metadata_columns <- c(
  "full_id",
  "caf_subcluster",
  "pseudotime",
  "state",
  "state_number",
  "component_1",
  "component_2"
)

missing_metadata_columns <- setdiff(
  required_metadata_columns,
  colnames(pseudotime_metadata_table)
)

if (length(missing_metadata_columns) > 0) {
  stop(
    "Required pseudotime metadata columns are missing: ",
    paste(missing_metadata_columns, collapse = ", "),
    call. = FALSE
  )
}

required_edge_columns <- c(
  "x",
  "y",
  "xend",
  "yend"
)

missing_edge_columns <- setdiff(
  required_edge_columns,
  colnames(trajectory_edge_table)
)

if (length(missing_edge_columns) > 0) {
  stop(
    "Required trajectory edge columns are missing: ",
    paste(missing_edge_columns, collapse = ", "),
    call. = FALSE
  )
}

if (anyDuplicated(pseudotime_metadata_table$cell_id) > 0) {
  stop(
    "Duplicated cell_id values were found in the pseudotime metadata.",
    call. = FALSE
  )
}

if (
  any(!is.finite(pseudotime_metadata_table$component_1)) ||
    any(!is.finite(pseudotime_metadata_table$component_2))
) {
  stop(
    "Non-finite trajectory coordinates were found in the cell metadata.",
    call. = FALSE
  )
}

if (
  any(!is.finite(trajectory_edge_table$x)) ||
    any(!is.finite(trajectory_edge_table$y)) ||
    any(!is.finite(trajectory_edge_table$xend)) ||
    any(!is.finite(trajectory_edge_table$yend))
) {
  stop(
    "Non-finite coordinates were found in the trajectory edge table.",
    call. = FALSE
  )
}

state_levels <- paste0(
  "State ",
  sort(unique(pseudotime_metadata_table$state_number))
)

caf_subcluster_levels <- paste0(
  "CAF-",
  setdiff(
    caf_subcluster_number_order,
    pseudotime_excluded_caf_subcluster
  )
)

pseudotime_metadata_table <- pseudotime_metadata_table %>%
  dplyr::mutate(
    state = factor(state, levels = state_levels),
    caf_subcluster = factor(
      caf_subcluster,
      levels = caf_subcluster_levels
    )
  )

if (any(is.na(pseudotime_metadata_table$state))) {
  stop(
    "Some state labels could not be matched to state_number.",
    call. = FALSE
  )
}

if (any(is.na(pseudotime_metadata_table$caf_subcluster))) {
  stop(
    "Some CAF subcluster labels were not represented in the expected order.",
    call. = FALSE
  )
}

missing_caf_colors <- setdiff(
  caf_subcluster_levels,
  names(caf_subcluster_colors)
)

if (length(missing_caf_colors) > 0) {
  stop(
    "CAF subcluster colors were not defined for: ",
    paste(missing_caf_colors, collapse = ", "),
    call. = FALSE
  )
}


caf_subcluster_plot_table <- dplyr::bind_rows(
  dplyr::filter(
    pseudotime_metadata_table,
    caf_subcluster != proliferating_caf_subcluster
  ),
  dplyr::filter(
    pseudotime_metadata_table,
    caf_subcluster == proliferating_caf_subcluster
  )
)


required_niche_metadata_columns <- c(
  "cell_id",
  "component_1",
  "component_2",
  "niche_annotation"
)

missing_niche_metadata_columns <- setdiff(
  required_niche_metadata_columns,
  colnames(pseudotime_niche_metadata_table)
)

if (length(missing_niche_metadata_columns) > 0) {
  stop(
    "Required niche metadata columns are missing: ",
    paste(missing_niche_metadata_columns, collapse = ", "),
    call. = FALSE
  )
}

if (anyDuplicated(pseudotime_niche_metadata_table$cell_id) > 0) {
  stop(
    "Duplicated cell_id values were found in the niche metadata.",
    call. = FALSE
  )
}

niche_trajectory_plot_table <- pseudotime_niche_metadata_table %>%
  filter(!is.na(niche_annotation))

if (nrow(niche_trajectory_plot_table) == 0) {
  stop(
    "No cells with spatial niche annotations were found in ",
    basename(pseudotime_niche_metadata_file),
    ".",
    call. = FALSE
  )
}

niche_levels <- niche_trajectory_plot_table$niche_annotation %>%
  unique() %>%
  str_sort(numeric = TRUE)

normal_acinar_niche <- "N1. Normal acinar niche"

if (!normal_acinar_niche %in% niche_levels) {
  stop(
    "The normal acinar niche label was not found: ",
    normal_acinar_niche,
    call. = FALSE
  )
}

niche_trajectory_plot_table <- bind_rows(
  niche_trajectory_plot_table %>%
    filter(niche_annotation != normal_acinar_niche),
  niche_trajectory_plot_table %>%
    filter(niche_annotation == normal_acinar_niche)
) %>%
  mutate(
    niche_annotation = factor(
      niche_annotation,
      levels = niche_levels
    )
  )

niche_highlight_colors <- setNames(
  ifelse(
    niche_levels == normal_acinar_niche,
    "red",
    "gray"
  ),
  niche_levels
)


## Define common plot elements ----------------------------------------------

common_theme <- theme(
  axis.title = element_text(
    color = "black",
    size = 10
  ),
  legend.title = element_text(
    color = "black",
    size = 10
  ),
  axis.text.x = element_text(
    color = "black",
    size = 10
  ),
  axis.text.y = element_text(
    color = "black",
    size = 10
  ),
  legend.text = element_text(
    color = "black",
    size = 10
  ),
  plot.background = element_rect(
    fill = NA,
    color = NA
  ),
  panel.background = element_rect(
    fill = NA,
    color = NA
  ),
  legend.background = element_rect(
    fill = NA,
    color = NA
  ),
  legend.key = element_blank(),
  panel.border = element_rect(
    fill = NA,
    color = "black"
  ),
  panel.grid = element_blank(),
  aspect.ratio = 0.84
)

panel_only_theme <- theme(
  plot.margin = margin(0, 0, 0, 0),
  plot.background = element_rect(
    fill = "transparent",
    color = NA
  ),
  panel.background = element_rect(
    fill = "transparent",
    color = NA
  ),
  panel.grid = element_blank(),
  legend.position = "none",
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  axis.title = element_blank(),
  axis.line = element_blank(),
  panel.border = element_blank()
)

common_labels <- labs(
  x = "Component 1",
  y = "Component 2",
  color = NULL
)


## Define output function ----------------------------------------------------

save_trajectory_plot <- function(
    pdf_plot,
    panel_png_plot,
    output_stem
) {
  ggsave(
    filename = file.path(
      pseudotime_trajectory_figure_dir,
      paste0(output_stem, ".pdf")
    ),
    plot = pdf_plot,
    device = grDevices::cairo_pdf,
    width = figure_width,
    height = figure_height,
    units = "in",
    bg = "transparent"
  )

  ggsave(
    filename = file.path(
      pseudotime_trajectory_panel_png_dir,
      paste0(output_stem, "_panel.png")
    ),
    plot = panel_png_plot,
    device = ragg::agg_png,
    width = panel_width,
    height = panel_height,
    units = "in",
    dpi = figure_dpi,
    bg = "transparent"
  )
}


## State-colored trajectory --------------------------------------------------

state_plot_pdf <- ggplot(
  pseudotime_metadata_table,
  aes(
    x = component_1,
    y = component_2
  )
) +
  geom_segment(
    data = trajectory_edge_table,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend
    ),
    inherit.aes = FALSE
  ) +
  ggrastr::rasterise(
    geom_point(
      aes(color = state)
    ),
    dpi = figure_dpi,
    dev = "ragg"
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 5)
    )
  ) +
  common_labels +
  common_theme

state_plot_panel_png <- ggplot(
  pseudotime_metadata_table,
  aes(
    x = component_1,
    y = component_2
  )
) +
  geom_segment(
    data = trajectory_edge_table,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend
    ),
    inherit.aes = FALSE
  ) +
  geom_point(
    aes(color = state)
  ) +
  panel_only_theme

save_trajectory_plot(
  pdf_plot = state_plot_pdf,
  panel_png_plot = state_plot_panel_png,
  output_stem = "pseudotime_analysis_trajectory_state"
)


## CAF-subcluster-colored trajectory -----------------------------------------

caf_subcluster_plot_pdf <- ggplot(
  caf_subcluster_plot_table,
  aes(
    x = component_1,
    y = component_2
  )
) +
  ggrastr::rasterise(
    geom_segment(
      data = trajectory_edge_table,
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend
      ),
      inherit.aes = FALSE
    ),
    dpi = figure_dpi,
    dev = "ragg"
  ) +
  ggrastr::rasterise(
    geom_point(
      aes(color = caf_subcluster),
      size = 0.8
    ),
    dpi = figure_dpi,
    dev = "ragg"
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 5)
    )
  ) +
  scale_color_manual(
    values = caf_subcluster_colors[caf_subcluster_levels],
    drop = FALSE
  ) +
  common_labels +
  common_theme

caf_subcluster_plot_panel_png <- ggplot(
  caf_subcluster_plot_table,
  aes(
    x = component_1,
    y = component_2
  )
) +
  geom_segment(
    data = trajectory_edge_table,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend
    ),
    inherit.aes = FALSE
  ) +
  geom_point(
    aes(color = caf_subcluster),
    size = 0.8
  ) +
  scale_color_manual(
    values = caf_subcluster_colors[caf_subcluster_levels],
    drop = FALSE
  ) +
  panel_only_theme

save_trajectory_plot(
  pdf_plot = caf_subcluster_plot_pdf,
  panel_png_plot = caf_subcluster_plot_panel_png,
  output_stem = "pseudotime_analysis_trajectory_caf_subcluster"
)


## Pseudotime-colored trajectory ---------------------------------------------

pseudotime_plot_pdf <- ggplot(
  pseudotime_metadata_table,
  aes(
    x = component_1,
    y = component_2
  )
) +
  ggrastr::rasterise(
    geom_segment(
      data = trajectory_edge_table,
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend
      ),
      inherit.aes = FALSE
    ),
    dpi = figure_dpi,
    dev = "ragg"
  ) +
  ggrastr::rasterise(
    geom_point(
      aes(color = pseudotime)
    ),
    dpi = figure_dpi,
    dev = "ragg"
  ) +
  labs(
    x = "Component 1",
    y = "Component 2",
    color = "Pseudotime"
  ) +
  scale_color_viridis(
    option = "inferno",
    breaks = c(0, 4, 8),
    guide = guide_colorbar(
      direction = "vertical",
      title.position = "top",
      title.hjust = 0.5,
      frame.colour = "black",
      ticks.colour = "black"
    )
  ) +
  common_theme

pseudotime_plot_panel_png <- ggplot(
  pseudotime_metadata_table,
  aes(
    x = component_1,
    y = component_2
  )
) +
  geom_point(
    aes(color = pseudotime)
  ) +
  geom_segment(
    data = trajectory_edge_table,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend
    ),
    inherit.aes = FALSE
  ) +
  scale_color_viridis(
    option = "inferno"
  ) +
  panel_only_theme

save_trajectory_plot(
  pdf_plot = pseudotime_plot_pdf,
  panel_png_plot = pseudotime_plot_panel_png,
  output_stem = "pseudotime_analysis_trajectory_pseudotime"
)


## Spatial-niche-colored trajectory coordinates -----------------------------

niche_plot_pdf <- ggplot(
  niche_trajectory_plot_table,
  aes(
    x = component_1,
    y = component_2
  )
) +
  ggrastr::rasterise(
    geom_point(
      aes(color = niche_annotation),
      size = 0.45
    ),
    dpi = figure_dpi,
    dev = "ragg"
  ) +
  labs(
    x = "Component 1",
    y = "Component 2",
    color = NULL
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 5)
    )
  ) +
  scale_color_manual(
    values = niche_highlight_colors[niche_levels],
    drop = FALSE
  ) +
  common_theme

niche_plot_panel_png <- ggplot(
  niche_trajectory_plot_table,
  aes(
    x = component_1,
    y = component_2
  )
) +
  geom_point(
    aes(color = niche_annotation)
  ) +
  scale_color_manual(
    values = niche_highlight_colors[niche_levels],
    drop = FALSE
  ) +
  panel_only_theme

save_trajectory_plot(
  pdf_plot = niche_plot_pdf,
  panel_png_plot = niche_plot_panel_png,
  output_stem = "pseudotime_analysis_trajectory_normal_acinar_niche_highlight"
)


## Print output summary ------------------------------------------------------

message(
  "Pseudotime trajectory figures were saved to:\n",
  pseudotime_trajectory_figure_dir
)

