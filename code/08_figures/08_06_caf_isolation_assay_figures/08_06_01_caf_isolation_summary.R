# CAF isolation summary figure
#
# This script visualizes the clinical characteristics, time to CAF outgrowth,
# and downstream assay usage for 18 paired CAF isolations.
#
# Each CAF pair is shown as one pair-level column. NAT, tumor location, and
# assay usage are pair-level variables. Only the outgrowth row is divided into
# two subcolumns representing 21% O2 and 1% O2 isolation conditions.


## Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
})


## Load settings -------------------------------------------------------------

project_dir <- getwd()

if (basename(project_dir) != "PDAC_CAF_hypoxia_analysis") {
  stop(
    "Please open PDAC_CAF_hypoxia_analysis.Rproj before running this script."
  )
}

source(file.path(project_dir, "code", "00_config", "00_02_paths.R"))


## Define output paths -------------------------------------------------------

caf_isolation_summary_figure_dir <- file.path(
  caf_isolation_figure_dir,
  "08_06_01_caf_isolation_summary_figures"
)

caf_isolation_summary_figure_file <- file.path(
  caf_isolation_summary_figure_dir,
  "caf_isolation_summary.pdf"
)

fs::dir_create(caf_isolation_summary_figure_dir)


## Figure settings -----------------------------------------------------------

pair_numbers <- 1:18

# Horizontal dimensions are defined in shared plot-coordinate units.
# Each CAF pair occupies one pair-level column, and gaps are inserted only
# between adjacent pairs.
width_title <- 5.00
width_each_pair <- 0.92
width_gap_inter_pair <- 0.03
width_padding_left <- 0.15
width_padding_right <- 0.15

# Convert plot-coordinate width to the physical PDF width.
# Adjust this value to scale the full figure horizontally without changing
# the relative widths of the title column, pair columns, or inter-pair gaps.
figure_inches_per_width_unit <- 0.58

width_pair_block <-
  length(pair_numbers) * width_each_pair +
  (length(pair_numbers) - 1) * width_gap_inter_pair

width_total <-
  width_padding_left +
  width_title +
  width_pair_block +
  width_padding_right

figure_width <- width_total * figure_inches_per_width_unit

x_title_left <- width_padding_left
x_title_right <- width_padding_left + width_title
x_title_text <- x_title_left

x_pair_block_left <- x_title_right

pair_layout <- tibble::tibble(
  caf_pair_no = pair_numbers
) %>%
  dplyr::mutate(
    xmin = x_pair_block_left +
      (.data$caf_pair_no - 1) * (width_each_pair + width_gap_inter_pair),
    xmax = .data$xmin + width_each_pair,
    x = (.data$xmin + .data$xmax) / 2
  )


# Outgrowth is the only row split into two subcolumns within each pair.
width_outgrowth_subcolumn <- width_each_pair / 2
outgrowth_x_offset <- width_each_pair / 4

# Row heights are defined in shared plot-coordinate units.
# Changing one value automatically updates all row positions and figure height.
height_pair_no <- 0.58
height_nat <- 0.58
height_location <- 0.58
height_isolation_oxygen <- 0.18
height_outgrowth <- 0.58
height_assay_name <- 0.45
height_legend_title <- 0.20
height_legend_row <- 0.40

height_gap <- 0.02
height_gap_before_legend <- 0.70
height_gap_legend_title <- 0.20

text_size_title <- 7.0
text_size_pair_no <- 7.0
text_size_outgrowth <- 4.5
text_size_assay <- 7.0
text_size_legend_title <- 5.0
text_size_legend_text <- 5.0

icon_size <- 5.0
icon_stroke <- 0.9
legend_icon_size <- 5.0
legend_icon_stroke <- 0.8


# Top and bottom padding are included in the automatic total height.
height_padding_top <- 0.25
height_padding_bottom <- 0.25

# Convert plot-coordinate height to the physical PDF height.
# Adjust this value to scale the full figure vertically without changing
# relative row heights.
figure_inches_per_height_unit <- 0.63


## Vertical layout ------------------------------------------------------------

row_pair_no_top <- height_padding_top
row_pair_no_bottom <- row_pair_no_top + height_pair_no

row_nat_top <- row_pair_no_bottom + height_gap
row_nat_bottom <- row_nat_top + height_nat

row_location_top <- row_nat_bottom + height_gap
row_location_bottom <- row_location_top + height_location

# Keep a gap between Location and the oxygen-condition strip.
row_isolation_oxygen_top <- row_location_bottom + height_gap
row_isolation_oxygen_bottom <-
  row_isolation_oxygen_top + height_isolation_oxygen

# No gap between the oxygen-condition strip and outgrowth tiles.
row_outgrowth_top <- row_isolation_oxygen_bottom
row_outgrowth_bottom <- row_outgrowth_top + height_outgrowth

row_rna_seq_top <- row_outgrowth_bottom + height_gap
row_rna_seq_bottom <- row_rna_seq_top + height_assay_name

row_proliferation_top <- row_rna_seq_bottom + height_gap
row_proliferation_bottom <-
  row_proliferation_top + height_assay_name

row_morphology_top <- row_proliferation_bottom + height_gap
row_morphology_bottom <- row_morphology_top + height_assay_name

row_legend_title_top <- row_morphology_bottom + height_gap_before_legend
row_legend_title_bottom <- row_legend_title_top + height_legend_title
y_legend_title <- (row_legend_title_top + row_legend_title_bottom) / 2

row_legend_1_top <- row_legend_title_bottom + height_gap_legend_title
row_legend_1_bottom <- row_legend_1_top + height_legend_row
y_legend_1 <- (row_legend_1_top + row_legend_1_bottom) / 2

row_legend_2_top <- row_legend_1_bottom + height_gap
row_legend_2_bottom <- row_legend_2_top + height_legend_row
y_legend_2 <- (row_legend_2_top + row_legend_2_bottom) / 2

row_legend_3_top <- row_legend_2_bottom + height_gap
row_legend_3_bottom <- row_legend_3_top + height_legend_row
y_legend_3 <- (row_legend_3_top + row_legend_3_bottom) / 2

height_total <-
  height_padding_top +
  row_legend_3_bottom +
  height_padding_bottom

figure_height <- height_total * figure_inches_per_height_unit

# Convert top-down row positions to the ggplot coordinate system.
to_plot_y <- function(y_from_top) {
  height_total - height_padding_top - y_from_top
}

y_pair_no_plot <-
  to_plot_y(mean(c(row_pair_no_top, row_pair_no_bottom)))

y_nat_plot <-
  to_plot_y(mean(c(row_nat_top, row_nat_bottom)))

y_location_plot <-
  to_plot_y(mean(c(row_location_top, row_location_bottom)))

y_isolation_oxygen_plot <-
  to_plot_y(mean(c(
    row_isolation_oxygen_top,
    row_isolation_oxygen_bottom
  )))

y_outgrowth_plot <-
  to_plot_y(mean(c(row_outgrowth_top, row_outgrowth_bottom)))

# Center the title across the combined oxygen-strip + outgrowth-tile block.
y_outgrowth_group_plot <-
  to_plot_y(mean(c(
    row_isolation_oxygen_top,
    row_outgrowth_bottom
  )))
y_rna_seq_plot <-
  to_plot_y(mean(c(row_rna_seq_top, row_rna_seq_bottom)))

y_proliferation_plot <-
  to_plot_y(mean(c(
    row_proliferation_top,
    row_proliferation_bottom
  )))

y_morphology_plot <-
  to_plot_y(mean(c(
    row_morphology_top,
    row_morphology_bottom
  )))
y_assay_group_plot <- mean(c(y_rna_seq_plot, y_morphology_plot))
y_legend_title_plot <- to_plot_y(y_legend_title)
y_legend_1_plot <- to_plot_y(y_legend_1)
y_legend_2_plot <- to_plot_y(y_legend_2)
y_legend_3_plot <- to_plot_y(y_legend_3)


# Tile heights follow the corresponding row heights.
nat_tile_height <- height_nat
location_tile_height <- height_location
outgrowth_tile_height <- height_outgrowth

nat_colors <- c(
  "None" = "#E5E5E5",
  "GS" = "#DFC0A0",
  "GnP" = "#F4A62A"
)

location_colors <- c(
  "Ph" = "#4C5AA8",
  "Pb" = "#AFC8D8",
  "Pt" = "#9BC9B8"
)

days_color_function <- scales::col_numeric(
  palette = c("#F7F7F7", "#D7C7DF", "#A77DB8", "#6B237F"),
  domain = c(1, 15)
)

oxygen_colors <- c(
  "21" = "#E6A6A6",
  "1" = "#9EAEE0"
)

cell_border_color <- "gray30"
cell_border_linewidth <- 0.25

pair_tile_width <- width_each_pair
outgrowth_tile_width <- width_outgrowth_subcolumn


## Helper functions ----------------------------------------------------------

read_required_csv <- function(file, description) {
  if (!file.exists(file)) {
    stop(description, " not found: ", file, call. = FALSE)
  }

  readr::read_csv(
    file,
    show_col_types = FALSE
  )
}

check_required_columns <- function(data, required_columns, table_name) {
  missing_columns <- setdiff(required_columns, colnames(data))

  if (length(missing_columns) > 0) {
    stop(
      "Missing required columns in ",
      table_name,
      ": ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(NULL)
}

format_legend_label <- function(label, n) {
  paste0(label, " (n = ", n, ")")
}


## Read and validate CAF isolation metadata ---------------------------------

caf_isolation_metadata <- read_required_csv(
  caf_isolation_metadata_file,
  "CAF isolation metadata file"
)

required_columns <- c(
  "caf_pair_no",
  "neoadjuvant_therapy",
  "tumor_location",
  "isolation_oxygen_percent",
  "days_to_outgrowth",
  "used_for_rna_seq",
  "used_for_proliferation_assay",
  "used_for_morphology_assay"
)

check_required_columns(
  caf_isolation_metadata,
  required_columns,
  "CAF isolation metadata"
)

caf_isolation_metadata <- caf_isolation_metadata %>%
  dplyr::select(dplyr::all_of(required_columns)) %>%
  dplyr::arrange(
    .data$caf_pair_no,
    dplyr::desc(.data$isolation_oxygen_percent)
  )

stopifnot(!anyNA(caf_isolation_metadata))
stopifnot(setequal(unique(caf_isolation_metadata$caf_pair_no), pair_numbers))
stopifnot(setequal(unique(caf_isolation_metadata$isolation_oxygen_percent), c(1, 21)))
stopifnot(nrow(caf_isolation_metadata) == 36)


## Prepare plotting data -----------------------------------------------------

pair_level_metadata <- caf_isolation_metadata %>%
  dplyr::group_by(.data$caf_pair_no) %>%
  dplyr::summarise(
    neoadjuvant_therapy = dplyr::first(.data$neoadjuvant_therapy),
    tumor_location = dplyr::first(.data$tumor_location),
    used_for_rna_seq = dplyr::first(.data$used_for_rna_seq),
    used_for_proliferation_assay = dplyr::first(
      .data$used_for_proliferation_assay
    ),
    used_for_morphology_assay = dplyr::first(
      .data$used_for_morphology_assay
    ),
    .groups = "drop"
  )

pair_no_rect_data <- pair_layout %>%
  dplyr::transmute(
    caf_pair_no,
    xmin,
    xmax,
    ymin = to_plot_y(row_pair_no_bottom),
    ymax = to_plot_y(row_pair_no_top)
  )

nat_plot_data <- pair_level_metadata %>%
  dplyr::left_join(pair_layout, by = "caf_pair_no") %>%
  dplyr::transmute(
    caf_pair_no = .data$caf_pair_no,
    x = .data$x,
    y = y_nat_plot,
    fill = unname(nat_colors[.data$neoadjuvant_therapy])
  )

location_plot_data <- pair_level_metadata %>%
  dplyr::left_join(pair_layout, by = "caf_pair_no") %>%
  dplyr::transmute(
    caf_pair_no = .data$caf_pair_no,
    x = .data$x,
    y = y_location_plot,
    fill = unname(location_colors[.data$tumor_location])
  )

outgrowth_plot_data <- caf_isolation_metadata %>%
  dplyr::left_join(pair_layout, by = "caf_pair_no") %>%
  dplyr::mutate(
    oxygen_key = as.character(.data$isolation_oxygen_percent),
    x = .data$x + dplyr::if_else(
      .data$isolation_oxygen_percent == 21,
      -outgrowth_x_offset,
      outgrowth_x_offset
    )
  ) %>%
  dplyr::transmute(
    caf_pair_no = .data$caf_pair_no,
    x = .data$x,
    y = y_outgrowth_plot,
    isolation_oxygen_percent = .data$isolation_oxygen_percent,
    days_to_outgrowth = .data$days_to_outgrowth,
    fill = days_color_function(.data$days_to_outgrowth),
    oxygen_fill = unname(oxygen_colors[.data$oxygen_key])
  )

isolation_oxygen_plot_data <- outgrowth_plot_data %>%
  dplyr::transmute(
    caf_pair_no = .data$caf_pair_no,
    x = .data$x,
    y = y_isolation_oxygen_plot,
    isolation_oxygen_percent = .data$isolation_oxygen_percent,
    fill = .data$oxygen_fill
  )

assay_plot_data <- pair_level_metadata %>%
  dplyr::left_join(pair_layout, by = "caf_pair_no") %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(c(
      "used_for_rna_seq",
      "used_for_proliferation_assay",
      "used_for_morphology_assay"
    )),
    names_to = "assay",
    values_to = "used"
  ) %>%
  dplyr::filter(.data$used) %>%
  dplyr::mutate(
    y = dplyr::recode(
      .data$assay,
      "used_for_rna_seq" = y_rna_seq_plot,
      "used_for_proliferation_assay" = y_proliferation_plot,
      "used_for_morphology_assay" = y_morphology_plot
    ),
    shape = dplyr::recode(
      .data$assay,
      "used_for_rna_seq" = 16,
      "used_for_proliferation_assay" = 15,
      "used_for_morphology_assay" = 0
    )
  )

assay_rect_data <- pair_layout %>%
  dplyr::transmute(
    caf_pair_no,
    xmin,
    xmax,
    ymin = to_plot_y(row_morphology_bottom),
    ymax = to_plot_y(row_rna_seq_top)
  )


## Prepare legend data -------------------------------------------------------

nat_counts <- pair_level_metadata %>%
  dplyr::count(.data$neoadjuvant_therapy) %>%
  tibble::deframe()

location_counts <- pair_level_metadata %>%
  dplyr::count(.data$tumor_location) %>%
  tibble::deframe()

assay_counts <- pair_level_metadata %>%
  dplyr::summarise(
    `RNA-seq` = sum(.data$used_for_rna_seq),
    `Proliferation assay` = sum(.data$used_for_proliferation_assay),
    `Morphology analysis` = sum(.data$used_for_morphology_assay)
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "label",
    values_to = "n"
  )

# Legend positions are defined relative to the full plot width so that they
# remain aligned when the pair-column width or inter-pair gap is changed.
legend_nat_x <- x_title_right
legend_location_x <- width_total * 0.36
legend_oxygen_x <- width_total * 0.49
legend_days_xmin <- width_total * 0.63
legend_days_xmax <- width_total * 0.75
legend_assay_x <- width_total * 0.80

legend_tile_height <- height_legend_row

legend_tile_width <-
  legend_tile_height *
  figure_inches_per_height_unit /
  figure_inches_per_width_unit

legend_text_offset <- legend_tile_width / 2 + 0.18

nat_legend_data <- tibble::tibble(
  x = legend_nat_x,
  y = c(y_legend_1_plot, y_legend_2_plot, y_legend_3_plot),
  fill = unname(nat_colors[c("None", "GS", "GnP")]),
  label = c(
    format_legend_label("None", nat_counts[["None"]]),
    format_legend_label("GS", nat_counts[["GS"]]),
    format_legend_label("GnP", nat_counts[["GnP"]])
  )
)

location_legend_data <- tibble::tibble(
  x = legend_location_x,
  y = c(y_legend_1_plot, y_legend_2_plot, y_legend_3_plot),
  fill = unname(location_colors[c("Ph", "Pb", "Pt")]),
  label = c(
    format_legend_label("Ph", location_counts[["Ph"]]),
    format_legend_label("Pb", location_counts[["Pb"]]),
    format_legend_label("Pt", location_counts[["Pt"]])
  )
)

oxygen_legend_data <- tibble::tibble(
  x = legend_oxygen_x,
  y = c(y_legend_1_plot, y_legend_2_plot),
  fill = unname(oxygen_colors[c("21", "1")]),
  label = c("21% O2", "1% O2")
)

days_legend_breaks <-
  seq(
    legend_days_xmin,
    legend_days_xmax,
    length.out = 16
  )

days_legend_data <- tibble::tibble(
  days_to_outgrowth = 1:15,
  xmin = days_legend_breaks[1:15],
  xmax = days_legend_breaks[2:16],
  ymin = y_legend_1_plot - height_legend_row / 2,
  ymax = y_legend_1_plot + height_legend_row / 2,
  fill = days_color_function(1:15)
)

assay_legend_data <- assay_counts %>%
  dplyr::mutate(
    x = legend_assay_x,
    y = c(y_legend_1_plot, y_legend_2_plot, y_legend_3_plot),
    shape = c(16, 15, 0),
    legend_label = paste0(.data$label, " (n = ", .data$n, ")")
  )


## Create CAF isolation summary figure --------------------------------------

caf_isolation_summary_plot <- ggplot() +
  geom_rect(
    data = pair_no_rect_data,
    aes(
      xmin = .data$xmin,
      xmax = .data$xmax,
      ymin = .data$ymin,
      ymax = .data$ymax
    ),
    inherit.aes = FALSE,
    fill = NA,
    color = cell_border_color,
    linewidth = cell_border_linewidth
  ) +
  geom_tile(
    data = nat_plot_data,
    aes(
      x = .data$x,
      y = .data$y,
      fill = .data$fill
    ),
    width = pair_tile_width,
    height = nat_tile_height,
    color = cell_border_color,
    linewidth = cell_border_linewidth
  ) +
  geom_tile(
    data = location_plot_data,
    aes(
      x = .data$x,
      y = .data$y,
      fill = .data$fill
    ),
    width = pair_tile_width,
    height = location_tile_height,
    color = cell_border_color,
    linewidth = cell_border_linewidth
  ) +
  geom_tile(
    data = isolation_oxygen_plot_data,
    aes(
      x = .data$x,
      y = .data$y,
      fill = .data$fill
    ),
    width = outgrowth_tile_width,
    height = height_isolation_oxygen,
    color = cell_border_color,
    linewidth = cell_border_linewidth
  ) +
  geom_tile(
    data = outgrowth_plot_data,
    aes(
      x = .data$x,
      y = .data$y,
      fill = .data$fill
    ),
    width = outgrowth_tile_width,
    height = outgrowth_tile_height,
    color = cell_border_color,
    linewidth = cell_border_linewidth
  ) +
  geom_text(
    data = outgrowth_plot_data,
    aes(
      x = .data$x,
      y = .data$y,
      label = .data$days_to_outgrowth
    ),
    size = text_size_outgrowth,
    color = "black"
  ) +
  geom_rect(
    data = assay_rect_data,
    aes(
      xmin = .data$xmin,
      xmax = .data$xmax,
      ymin = .data$ymin,
      ymax = .data$ymax
    ),
    inherit.aes = FALSE,
    fill = "#FAFAFA",
    color = cell_border_color,
    linewidth = cell_border_linewidth
  ) +
  geom_point(
    data = assay_plot_data,
    aes(
      x = .data$x,
      y = .data$y,
      shape = factor(.data$shape)
    ),
    size = icon_size,
    stroke = icon_stroke,
    color = "black",
    fill = "white"
  ) +
  annotate(
    "text",
    x = pair_layout$x,
    y = y_pair_no_plot,
    label = pair_layout$caf_pair_no,
    size = text_size_pair_no
  ) +
  annotate(
    "text",
    x = x_title_text,
    y = c(
      y_pair_no_plot,
      y_nat_plot,
      y_location_plot,
      y_outgrowth_group_plot,
      y_assay_group_plot
    ),
    label = c(
      "CAF pair No.",
      "NAT",
      "Location",
      "Days to outgrowth",
      "Assay"
    ),
    hjust = 0,
    fontface = "bold",
    size = text_size_title
  ) +
  annotate(
    "text",
    x = x_title_right - 0.15,
    y = c(y_rna_seq_plot, y_proliferation_plot, y_morphology_plot),
    label = c(
      "RNA-seq",
      "Proliferation",
      "Morphology"
    ),
    hjust = 1,
    size = text_size_assay
  ) +
  geom_tile(
    data = dplyr::bind_rows(
      nat_legend_data,
      location_legend_data,
      oxygen_legend_data
    ),
    aes(
      x = .data$x,
      y = .data$y,
      fill = .data$fill
    ),
    width = legend_tile_width,
    height = legend_tile_height,
    color = "black",
    linewidth = 0.25
  ) +
  geom_text(
    data = dplyr::bind_rows(
      nat_legend_data,
      location_legend_data,
      oxygen_legend_data
    ),
    aes(x = .data$x + legend_text_offset, y = .data$y, label = .data$label),
    hjust = 0,
    size = text_size_legend_text
  ) +
  annotate(
    "text",
    x = c(
      legend_nat_x - legend_tile_width / 2,
      legend_location_x - legend_tile_width / 2,
      legend_oxygen_x - legend_tile_width / 2,
      legend_days_xmin,
      legend_assay_x - 0.15
    ),
    y = y_legend_title_plot,
    label = c(
      "NAT",
      "Location",
      "Isolation oxygen",
      "Days to outgrowth",
      "Assay"
    ),
    hjust = 0,
    fontface = "bold",
    size = text_size_legend_title
  ) +
  geom_rect(
    data = days_legend_data,
    aes(
      xmin = .data$xmin,
      xmax = .data$xmax,
      ymin = .data$ymin,
      ymax = .data$ymax,
      fill = .data$fill
    ),
    color = "grey40",
    linewidth = 0.15
  ) +
  annotate(
    "text",
    x = c(legend_days_xmin, legend_days_xmax),
    y = y_legend_2_plot,
    label = c("1", "15"),
    hjust = c(0, 1),
    size = text_size_legend_text
  ) +
  annotate(
    "text",
    x = c(legend_days_xmin, legend_days_xmax),
    y = y_legend_3_plot,
    label = c("(faster)", "(slower)"),
    hjust = c(0, 1),
    size = text_size_legend_text
  ) +
  geom_point(
    data = assay_legend_data,
    aes(x = .data$x, y = .data$y, shape = factor(.data$shape)),
    size = legend_icon_size,
    stroke = legend_icon_stroke,
    color = "black",
    fill = "white"
  ) +
  geom_text(
    data = assay_legend_data,
    aes(
      x = .data$x + legend_text_offset,
      y = .data$y,
      label = .data$legend_label
    ),
    hjust = 0,
    size = text_size_legend_text
  ) +
  scale_fill_identity() +
  scale_shape_manual(
    values = c("0" = 0, "15" = 15, "16" = 16),
    guide = "none"
  ) +
  scale_x_continuous(
    limits = c(0, width_total),
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    limits = c(0, height_total),
    expand = expansion(mult = 0)
  ) +
  coord_cartesian(clip = "off") +
  theme_void(base_size = 10) +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.margin = margin(8, 8, 8, 8)
  )


## Save figure ---------------------------------------------------------------

ggsave(
  filename = caf_isolation_summary_figure_file,
  plot = caf_isolation_summary_plot,
  device = cairo_pdf,
  width = figure_width,
  height = figure_height,
  units = "in",
  bg = "transparent"
)

message("Saved CAF isolation summary figure: ", caf_isolation_summary_figure_file)
