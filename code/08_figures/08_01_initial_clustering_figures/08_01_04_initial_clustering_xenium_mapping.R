# Initial clustering Xenium spatial mapping
#
# This script generates Xenium spatial maps for initial clustering labels.
# For each sample, it outputs a whole-tissue map with rectangles indicating
# two selected zoom-in regions and separate zoom-in maps for those regions.


## Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(arrow)
  library(fs)
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


## Define zoom-in coordinates ------------------------------------------------

coordinate_area1 <- list(
  "01" = c(2500, 2900),
  "02" = c(1300, 9500),
  "11" = c(2100, 3350),
  "16" = c(4600, 5200),
  "15" = c(1250, 10550),
  "19" = c(7770, 6670)
)

coordinate_area2 <- list(
  "01" = c(6200, 3100),
  "02" = c(3800, 2500),
  "11" = c(2650, 3000),
  "16" = c(4800, 3400),
  "15" = c(5400, 7450),
  "19" = c(8000, 750)
)

side_length <- 400


## Read initial clustering metadata -----------------------------------------

initial_metadata_df <- read.csv(
  initial_metadata_file,
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    initial_cluster = as.character(initial_cluster),
    cell_id_no_prefix = stringr::str_remove(
      cell_id,
      pattern = "^Xenium_[0-9]+_"
    )
  )


## Prepare colors ------------------------------------------------------------

unclassified_clusters <- as.character(initial_cluster_annotation_list[["Unclassified"]])

initial_cluster_colors_for_final <- c(
  initial_cluster_colors,
  "Excluded_byQC" = "#3A3A3A"
)

initial_cluster_colors_for_final[unclassified_clusters] <- "#3A3A3A"

stopifnot(!any(is.na(initial_cluster_colors_for_final)))


## Helper functions ----------------------------------------------------------

calculate_image_size <- function(range_x, range_y, base_size = 10) {
  length_x <- diff(range_x)
  length_y <- diff(range_y)
  
  if (length_y > length_x) {
    c(width = base_size, height = base_size * length_y / length_x)
  } else {
    c(width = base_size * length_x / length_y, height = base_size)
  }
}

calculate_scale_bar <- function(range_x, range_y, scale_bar_length = 1000) {
  length_x <- diff(range_x)
  length_y <- diff(range_y)
  
  if (length_y > length_x) {
    list(
      x = c(
        start = range_x[1] + length_x * 0.94,
        end = range_x[1] + length_x * 0.94 - scale_bar_length
      ),
      y = c(
        start = range_y[1] + length_y * 0.96,
        end = range_y[1] + length_y * 0.96
      )
    )
  } else {
    list(
      x = c(
        start = range_x[1] + length_x * 0.04,
        end = range_x[1] + length_x * 0.04
      ),
      y = c(
        start = range_y[1] + length_y * 0.94 - scale_bar_length,
        end = range_y[1] + length_y * 0.94
      )
    )
  }
}

plot_initial_xenium_map <- function(boundary_df, range_x, range_y) {
  ggplot(
    boundary_df,
    aes(x = vertex_x, y = vertex_y, group = label_id)
  ) +
    geom_polygon(
      aes(fill = initial_cluster),
      alpha = 1,
      color = NA
    ) +
    labs(
      x = NULL,
      y = NULL
    ) +
    scale_x_continuous(
      expand = expansion(mult = c(0, 0)),
      limits = range_x
    ) +
    scale_y_reverse(
      expand = expansion(mult = c(0, 0)),
      limits = c(range_y[2], range_y[1])
    ) +
    scale_fill_manual(
      values = initial_cluster_colors_for_final
    ) +
    coord_fixed(ratio = 1) +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "black", color = NA),
      panel.border = element_rect(fill = "transparent", color = "black"),
      panel.grid = element_blank()
    )
}

add_rectangles_and_scale_bar <- function(plot, coord_area1, coord_area2, scale_bar) {
  plot +
    annotate(
      "rect",
      xmin = coord_area1[1],
      xmax = coord_area1[1] + side_length,
      ymin = coord_area1[2],
      ymax = coord_area1[2] + side_length,
      color = "white",
      fill = NA,
      linewidth = 3
    ) +
    annotate(
      "rect",
      xmin = coord_area1[1],
      xmax = coord_area1[1] + side_length,
      ymin = coord_area1[2],
      ymax = coord_area1[2] + side_length,
      color = "black",
      fill = NA,
      linewidth = 2
    ) +
    annotate(
      "rect",
      xmin = coord_area2[1],
      xmax = coord_area2[1] + side_length,
      ymin = coord_area2[2],
      ymax = coord_area2[2] + side_length,
      color = "white",
      fill = NA,
      linewidth = 3
    ) +
    annotate(
      "rect",
      xmin = coord_area2[1],
      xmax = coord_area2[1] + side_length,
      ymin = coord_area2[2],
      ymax = coord_area2[2] + side_length,
      color = "black",
      fill = NA,
      linewidth = 2
    ) +
    annotate(
      "segment",
      x = scale_bar$x[["start"]],
      xend = scale_bar$x[["end"]],
      y = scale_bar$y[["start"]],
      yend = scale_bar$y[["end"]],
      color = "black",
      linewidth = 5
    ) +
    annotate(
      "segment",
      x = scale_bar$x[["start"]],
      xend = scale_bar$x[["end"]],
      y = scale_bar$y[["start"]],
      yend = scale_bar$y[["end"]],
      color = "white",
      linewidth = 2.5
    ) +
    theme(
      plot.margin = unit(c(0, 0, 0, 0), "mm")
    )
}

crop_initial_xenium_map <- function(plot, coord_area) {
  plot +
    coord_fixed(
      ratio = 1,
      xlim = c(coord_area[1], coord_area[1] + side_length),
      ylim = c(coord_area[2] + side_length, coord_area[2])
    ) +
    theme(
      plot.margin = unit(c(0, 0, 0, 0), "mm"),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "black", color = NA),
      panel.border = element_rect(fill = "transparent", color = "black")
    )
}


## Generate and save Xenium maps --------------------------------------------

xenium_mapping_output_dir <- file.path(
  initial_clustering_figure_dir,
  "08_01_04_initial_clustering_xenium_mapping"
)

fs::dir_create(xenium_mapping_output_dir)

for (sample_index in seq_along(tx_number_merge_order)) {
  tx_number <- tx_number_merge_order[[sample_index]]
  tx_id <- paste0("TX5K_", tx_number)
  
  xenium_sample_id <- xenium_sample_id_by_tx_number[[tx_number]]
  stopifnot(!is.null(xenium_sample_id))
  
  xenium_raw_folder <- xenium_raw_folder_by_tx_id[[tx_id]]
  stopifnot(!is.null(xenium_raw_folder))
  
  xenium_raw_dir <- file.path(
    xenium_raw_data_dir,
    xenium_raw_folder
  )
  
  cell_boundary_file <- file.path(
    xenium_raw_dir,
    "cell_boundaries.parquet"
  )
  
  message("Raw directory: ", xenium_raw_dir)
  message("Boundary file: ", cell_boundary_file)
  
  if (!file.exists(cell_boundary_file)) {
    stop(
      "cell_boundaries.parquet was not found:\n",
      cell_boundary_file,
      "\nAvailable files in this directory:\n",
      paste(list.files(xenium_raw_dir), collapse = "\n")
    )
  }
  
  cell_boundary_df <- arrow::read_parquet(cell_boundary_file)
  
  sample_metadata_df <- initial_metadata_df %>%
    dplyr::filter(tx_id == !!tx_id)
  
  stopifnot(nrow(sample_metadata_df) > 0)
  
  boundary_df <- cell_boundary_df %>%
    dplyr::left_join(
      sample_metadata_df %>%
        dplyr::select(
          cell_id_no_prefix,
          initial_cluster
        ),
      by = c("cell_id" = "cell_id_no_prefix")
    ) %>%
    dplyr::mutate(
      initial_cluster = dplyr::case_when(
        is.na(initial_cluster) ~ "Excluded_byQC",
        TRUE ~ initial_cluster
      )
    )
  
  range_x <- range(sample_metadata_df$x, na.rm = TRUE)
  range_y <- range(sample_metadata_df$y, na.rm = TRUE)
  
  image_size <- calculate_image_size(range_x, range_y)
  scale_bar <- calculate_scale_bar(range_x, range_y)
  
  plot_initial_xenium_full <- plot_initial_xenium_map(
    boundary_df = boundary_df,
    range_x = range_x,
    range_y = range_y
  )
  
  plot_initial_xenium_with_rectangles <- add_rectangles_and_scale_bar(
    plot = plot_initial_xenium_full,
    coord_area1 = coordinate_area1[[tx_number]],
    coord_area2 = coordinate_area2[[tx_number]],
    scale_bar = scale_bar
  )
  
  plot_initial_xenium_area1 <- crop_initial_xenium_map(
    plot = plot_initial_xenium_full,
    coord_area = coordinate_area1[[tx_number]]
  )
  
  plot_initial_xenium_area2 <- crop_initial_xenium_map(
    plot = plot_initial_xenium_full,
    coord_area = coordinate_area2[[tx_number]]
  )
  
  ggsave(
    plot = plot_initial_xenium_with_rectangles,
    filename = file.path(
      xenium_mapping_output_dir,
      paste0(
        "initial_clustering_map_",
        stringr::str_to_lower(xenium_sample_id),
        "_whole_tissue.png"
      )
    ),
    width = image_size[["width"]],
    height = image_size[["height"]],
    dpi = figure_dpi,
    bg = figure_bg,
    limitsize = FALSE
  )
  
  ggsave(
    plot = plot_initial_xenium_area1,
    filename = file.path(
      xenium_mapping_output_dir,
      paste0(
        "initial_clustering_map_",
        stringr::str_to_lower(xenium_sample_id),
        "_area1.png"
      )
    ),
    width = 4,
    height = 4,
    dpi = figure_dpi,
    bg = figure_bg
  )
  
  ggsave(
    plot = plot_initial_xenium_area2,
    filename = file.path(
      xenium_mapping_output_dir,
      paste0(
        "initial_clustering_map_",
        stringr::str_to_lower(xenium_sample_id),
        "_area2.png"
      )
    ),
    width = 4,
    height = 4,
    dpi = figure_dpi,
    bg = figure_bg
  )
}
