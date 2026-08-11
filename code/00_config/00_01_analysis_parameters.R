

## Figure settings ------------------------------------------------------------

figure_dpi <- 500
figure_bg <- "transparent"


## Analysis parameters -------------------------------------------------------

xenium_sample_id_by_tx_number <- c(
  "01" = "Xenium_01",
  "02" = "Xenium_02",
  "11" = "Xenium_03",
  "15" = "Xenium_04",
  "16" = "Xenium_05",
  "19" = "Xenium_06"
)

xenium_sample_label_by_tx_id <- c(
  "TX5K_01" = "Xenium_01",
  "TX5K_02" = "Xenium_02",
  "TX5K_11" = "Xenium_03",
  "TX5K_15" = "Xenium_04",
  "TX5K_16" = "Xenium_05",
  "TX5K_19" = "Xenium_06"
)

xenium_sample_number_display_order <- sprintf("%02d", 1:6)
xenium_sample_id_display_order <- paste0("Xenium_", xenium_sample_number_display_order)
tx_number_merge_order <- c("02", "19", "11", "16", "01", "15")
tx_id_merge_order <- paste0("TX5K_", tx_number_merge_order)


## Parameters for 01_xenium_preprocessing_integration_initial_clustering ------

qc_nfeature_rna <- c(100, 900)
qc_ncount_rna <- c(100, 1800)
qc_information <- paste0(
  "Countable_mag1",
  "_qc_nFeature_RNA:", paste(qc_nfeature_rna, collapse = "~"),
  "_qc_nCount_RNA:", paste(qc_ncount_rna, collapse = "~")
)
qc_information_filename <- gsub(
  pattern = ":",
  replacement = "",
  x = qc_information,
  fixed = TRUE
)

Dim1 <- 20
Res1 <- 1.0


## Parameters for 02_caf_subclustering ----------------------------------------

Dim2 <- 30
Res2 <- 0.5

caf_subcluster_number_order <- c(4, 8, 0, 1, 2, 5, 3, 7, 6)
caf_subcluster_display_order <- paste0("CAF-", caf_subcluster_number_order)

proliferating_caf_subcluster <- "CAF-8"

## Parameters for 03_hypoxia_mapping ------------------------------------------

grid_length_um <- 200
min_cells_per_grid <- 10


roi_side_length_um <- 400

hypoxia_roi_coordinates <- list(
  "01" = c(9400, 1200),
  "02" = c(1400, 9400),
  "11" = c(3400, 3600),
  "15" = c(2200, 7400),
  "16" = c(5200, 3600),
  "19" = c(7600, 4000)
)

normoxia_roi_coordinates <- list(
  "01" = c(3800, 2000),
  "02" = c(1600, 7600),
  "11" = c(2600, 3000),
  "15" = c(4600, 9200),
  "16" = c(1800, 2000),
  "19" = c(8400, 1200)
)


## Parameters for 04_niche_analysis ------------------------------------------- 

niche_n_neighbors <- 40
niche_pc_use <- 15
niche_k <- 10
niche_random_seed <- 123


## Parameters for 05_pseudotime_analysis --------------------------------------

pseudotime_excluded_caf_subcluster <- 6
pseudotime_n_hvgs <- 1500
pseudotime_root_state <- 1
pseudotime_downsample_per_caf_subcluster <- 500
pseudotime_random_seed <- 1234


## Parameters for 06_caf_isolation_assays and 07_bulk_rnaseq ------------------

isolation_oxygen_plot_labels <- c(
  "21" = "'21%'~O[2]",
  "1" = "'1%'~O[2]"
)


## Parameters for 07_bulk_rnaseq ----------------------------------------------

bulk_rnaseq_comparison_plot_labels <- c(
  "Original" = "Isolation",
  "Under21" = "21% O2",
  "Under1" = "1% O2"
)

method <- "SpearmanAverage"
SetName <- "Set1_3pair"

