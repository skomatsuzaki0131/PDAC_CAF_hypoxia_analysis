
## Color palettes ------------------------------------------------------------

cell_type_colors <- c(
  "Acinar" = "#ff00ff",
  "Ductal_like_acinar" = "#ddbcff",
  "Ductal" = "#ddbcff",
  "Normal_ductal" = "#FFFFCC",
  "PanIN1" = "#FFFF00",
  "PanIN2" = "#FD8D3C",
  "PDAC" = "#ff0000",
  "Islet" = "#BF812D",
  "CAF" = "#00ff00",
  "Mural" = "#990000",
  "Endothelial" = "#990000",
  "Lymphoid" = "#00bfff",
  "Lymph_T" = "#00bfff",
  "Lymph_B" = "#00bfff",
  "Plasma" = "#00bfff",
  "Myeloid" = "#8300FF",
  "Mast" = "#8300FF",
  "Nerve" = "#0000ff",
  "Proliferating_Immune" = "#999999",
  "Unclassified" = "#999999"
)

initial_cluster_colors <- c(
  "0" = "#ff00ff",
  "13" = "#ff00ff",
  "24" = "#ff00ff",
  "33" = "#ff00ff",
  "12" = "#ddbcff",
  "1" = "#ddbcff",
  "17" = "#ddbcff",
  "30" = "#ddbcff",
  "23" = "#ddbcff",
  "10" = "#ddbcff",
  "27" = "#ddbcff",
  "20" = "#ddbcff",
  "14" = "#ffff00",
  "9" = "#00ff00",
  "6" = "#00ff00",
  "7" = "#00ff00",
  "11" = "#990000",
  "34" = "#ff0000",
  "3" = "#ff0000",
  "2" = "#00bfff",
  "16" = "#00bfff",
  "15" = "#00bfff",
  "18" = "#00bfff",
  "37" = "#00bfff",
  "5" = "#ff8800",
  "4" = "#ff8800",
  "31" = "#ff8800",
  "26" = "#ff8800",
  "29" = "#0000ff",
  "25" = "#ABABAB",
  "22" = "#ABABAB",
  "28" = "#ABABAB",
  "36" = "#ABABAB",
  "35" = "#ABABAB",
  "32" = "#ABABAB",
  "8" = "#ABABAB",
  "19" = "#ABABAB",
  "21" = "#ABABAB"
)

caf_subcluster_colors_number_named <- c(
  "4" = "#72B7B2",
  "8" = "#C00000",
  "0" = "#4C78A8",
  "1" = "#54A24B",
  "2" = "#B279A2",
  "5" = "#E6AF2E",
  "3" = "#9D755D",
  "7" = "#9C9EDE",
  "6" = "#F58518"
)

caf_subcluster_colors_label_named <- setNames(
  caf_subcluster_colors_number_named,
  paste0("CAF-", names(caf_subcluster_colors_number_named))
)


diverging_colors <- c("#A6D3C5", "#F7F7F7", "#C05A9E")


## Hypoxia mapping ------------------------------------------------------------

hypoxia_roi_color <- "#253494"
normoxia_roi_color <- "#f03b20"
all_caf_color <- cell_type_colors[["CAF"]]


## Niche analysis -------------------------------------------------------------

niche_colors <- c(
  "niche_1"  = "#ffd000",
  "niche_2"  = "#1f78b4",
  "niche_3"  = "#a6cee3",
  "niche_4"  = "#33a02c",
  "niche_5"  = "#b15928",
  "niche_6"  = "#b2df8a",
  "niche_7"  = "#cab2d6",
  "niche_8"  = "#6a3d9a",
  "niche_9"  = "#e31a1c",
  "niche_10" = "#ff7f00"
)

normal_acinar_highlight_colors <- c(
  "Other niches" = "gray80",
  "N1. Normal acinar niche" = "#EF3B2C"
)

## CAF isolation --------------------------------------------------------------

caf_culture_colors <- c(
  "H_CAF" = "#91BFFA", 
  "N_CAF" = "#FEA0A0"
)

caf_culture_colors_4_groups <- c(
  "hypocaf_in_normoxia" = "#7FBFFF",
  "hypocaf_in_hypoxia" = "#607DF0",
  "normocaf_in_normoxia" = "#FF7F7F",
  "normocaf_in_hypomoxia" = "#FFA3A3"
)

nat_colors <- c(
  "None" = "#E5E5E5",
  "GS" = "#E0C3A9",
  "GnP" = "#EEA93C"
)

location_colors <- c(
  "Ph" = "#424F9E",
  "Pb" = "#A2BDD4",
  "Pt" = "#9BC8B6"
)
