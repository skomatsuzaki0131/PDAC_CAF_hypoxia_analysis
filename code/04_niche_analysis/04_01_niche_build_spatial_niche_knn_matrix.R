# Build spatial niche kNN matrix
#
# This script builds a k-nearest-neighbor-based spatial niche composition
# matrix for downstream spatial niche clustering.
#
# Important design choice:
# The kNN composition matrix is built from combined cell labels:
# CAF subcluster labels for CAFs (CAF-0 to CAF-8) and initial clustering labels
# for non-CAF cells (Clust.0, Clust.1, ...), after excluding unclassified
# initial clusters.


## Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(Matrix)
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


## Define input and output paths --------------------------------------------

niche_table_dir <- file.path(xenium_table_dir, "niche_analysis")
niche_object_dir <- file.path(xenium_object_dir, "niche_analysis")

fs::dir_create(niche_table_dir)
fs::dir_create(niche_object_dir)


## Define parameters ---------------------------------------------------------

# BuildNicheAssay uses the spatial field-of-view name stored in the
# coordinate-preserving Xenium objects generated during preprocessing.
xenium_fov_name <- "fov"


## Helper functions ----------------------------------------------------------

# The coordinate-preserving Xenium objects are generated in
# 01_01_xenium_preprocessing_integration_initial_clustering.R
# Their cell names are local Xenium cell IDs, whereas the public metadata table
# stores cell IDs prefixed with Xenium sample IDs. Therefore, local_cell_id is
# used to match metadata rows to each coordinate object.
read_xenium_coordinate_object <- function(xenium_sample_id) {
  xenium_coordinate_object_file <- file.path(
    xenium_coordinate_object_dir,
    paste0("xenium_coordinate_object_", xenium_sample_id, ".rds")
  )

  if (!file.exists(xenium_coordinate_object_file)) {
    stop(
      "Xenium coordinate object not found for ",
      xenium_sample_id,
      ". Expected file:\n",
      xenium_coordinate_object_file,
      call. = FALSE
    )
  }

  readRDS(xenium_coordinate_object_file)
}

prepare_initial_metadata_for_niche <- function(metadata_table) {
  required_columns <- c(
    "cell_id",
    "xenium_sample_id",
    "tx_id",
    "initial_cluster",
    "initial_cell_type"
  )

  missing_columns <- setdiff(required_columns, colnames(metadata_table))

  if (length(missing_columns) > 0) {
    stop(
      "Missing required columns in initial metadata: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  metadata_table %>%
    dplyr::mutate(
      cell_id = as.character(.data$cell_id),
      xenium_sample_id = as.character(.data$xenium_sample_id),
      tx_id = as.character(.data$tx_id),
      tx_number = stringr::str_remove(.data$tx_id, "^TX5K_?"),
      local_cell_id = stringr::str_remove(
        .data$cell_id,
        paste0("^", .data$xenium_sample_id, "_")
      ),
      full_id = paste0(
        stringr::str_remove(.data$tx_id, "_"),
        "_",
        .data$local_cell_id
      ),
      initial_cluster = as.character(.data$initial_cluster),
      initial_cell_type = as.character(.data$initial_cell_type)
    )
}

bind_sparse_matrices_by_column <- function(matrix_list) {
  if (length(matrix_list) == 0) {
    stop("No matrices were supplied.", call. = FALSE)
  }

  all_column_names <- matrix_list %>%
    purrr::map(colnames) %>%
    unlist(use.names = FALSE) %>%
    unique()

  aligned_matrices <- purrr::map(
    matrix_list,
    function(x) {
      missing_columns <- setdiff(all_column_names, colnames(x))

      if (length(missing_columns) > 0) {
        zero_matrix <- Matrix::Matrix(
          0,
          nrow = nrow(x),
          ncol = length(missing_columns),
          sparse = TRUE,
          dimnames = list(rownames(x), missing_columns)
        )
        x <- Matrix::cbind2(x, zero_matrix)
      }

      x[, all_column_names, drop = FALSE]
    }
  )

  Reduce(Matrix::rbind2, aligned_matrices)
}


## Read metadata and prepare cell labels -------------------------------------

if (!file.exists(initial_metadata_file)) {
  stop("Initial clustering metadata file not found: ", initial_metadata_file, call. = FALSE)
}

initial_metadata_table <- readr::read_csv(
  initial_metadata_file,
  show_col_types = FALSE
)

if (!file.exists(caf_subcluster_metadata_file)) {
  stop("CAF subcluster metadata file not found: ", caf_subcluster_metadata_file, call. = FALSE)
}

caf_subcluster_metadata_table <- readr::read_csv(
  caf_subcluster_metadata_file,
  show_col_types = FALSE
)

required_caf_metadata_columns <- c(
  "cell_id",
  "xenium_sample_id",
  "tx_id",
  "caf_subcluster"
)

missing_caf_metadata_columns <- setdiff(
  required_caf_metadata_columns,
  colnames(caf_subcluster_metadata_table)
)

if (length(missing_caf_metadata_columns) > 0) {
  stop(
    "Missing required columns in CAF subcluster metadata: ",
    paste(missing_caf_metadata_columns, collapse = ", "),
    call. = FALSE
  )
}

caf_subcluster_metadata_table <- caf_subcluster_metadata_table %>%
  dplyr::mutate(
    cell_id = as.character(.data$cell_id),
    xenium_sample_id = as.character(.data$xenium_sample_id),
    tx_id = as.character(.data$tx_id),
    local_cell_id = stringr::str_remove(
      .data$cell_id,
      paste0("^", .data$xenium_sample_id, "_")
    ),
    full_id = paste0(
      stringr::str_remove(.data$tx_id, "_"),
      "_",
      .data$local_cell_id
    ),
    caf_subcluster = as.character(.data$caf_subcluster)
  ) %>%
  dplyr::select(
    full_id,
    caf_subcluster
  ) %>%
  dplyr::distinct()

if (any(duplicated(caf_subcluster_metadata_table$full_id))) {
  stop(
    "Duplicated full_id values were detected in CAF subcluster metadata.",
    call. = FALSE
  )
}

niche_label_table <- initial_metadata_table %>%
  prepare_initial_metadata_for_niche() %>%
  dplyr::filter(.data$initial_cell_type != "Unclassified") %>%
  dplyr::left_join(
    caf_subcluster_metadata_table,
    by = "full_id"
  ) %>%
  dplyr::mutate(
    label_for_niche_assay = dplyr::if_else(
      !is.na(.data$caf_subcluster),
      .data$caf_subcluster,
      paste0("Clust.", .data$initial_cluster)
    )
  ) %>%
  dplyr::select(
    full_id,
    tx_id,
    tx_number,
    xenium_sample_id,
    cell_id,
    local_cell_id,
    initial_cluster,
    initial_cell_type,
    caf_subcluster,
    label_for_niche_assay
  ) %>%
  dplyr::distinct()

if (any(is.na(niche_label_table$tx_number))) {
  stop("Some tx_id values do not contain a TX5K sample number.", call. = FALSE)
}

if (any(is.na(niche_label_table$local_cell_id))) {
  stop("Some local cell IDs could not be derived from cell_id.", call. = FALSE)
}

if (any(duplicated(niche_label_table$full_id))) {
  duplicated_full_ids <- niche_label_table$full_id[duplicated(niche_label_table$full_id)]
  stop(
    "Duplicated full_id values were detected, for example: ",
    paste(head(duplicated_full_ids), collapse = ", "),
    call. = FALSE
  )
}

readr::write_csv(
  niche_label_table,
  file = file.path(
    niche_table_dir,
    "spatial_niche_input_cell_labels.csv"
  )
)


## Build kNN count matrix ----------------------------------------------------

knn_count_matrix_list <- list()

missing_tx_numbers <- setdiff(
  as.character(tx_number_merge_order),
  unique(niche_label_table$tx_number)
)

if (length(missing_tx_numbers) > 0) {
  stop(
    "Some TX numbers in tx_number_merge_order were not found in niche_label_table: ",
    paste(missing_tx_numbers, collapse = ", "),
    call. = FALSE
  )
}

tx_numbers_to_process <- as.character(tx_number_merge_order)

for (tx_number_i in tx_numbers_to_process) {
  sample_label_table <- niche_label_table %>%
    dplyr::filter(.data$tx_number == tx_number_i)
  
  xenium_sample_id <- unique(sample_label_table$xenium_sample_id)
  
  if (length(xenium_sample_id) != 1) {
    stop(
      "Expected one Xenium sample ID for TX5K_",
      tx_number_i,
      ", but found: ",
      paste(xenium_sample_id, collapse = ", "),
      call. = FALSE
    )
  }
  
  message("Building spatial niche kNN matrix for ", xenium_sample_id)
  
  sample_obj <- read_xenium_coordinate_object(xenium_sample_id)
  
  missing_cells <- setdiff(
    sample_label_table$local_cell_id,
    colnames(sample_obj)
  )
  
  if (length(missing_cells) > 0) {
    stop(
      "Some cells in the metadata were not found in the Xenium object for TX5K_",
      tx_number_i,
      ". Number of missing cells: ",
      length(missing_cells),
      call. = FALSE
    )
  }
  
  sample_obj <- subset(
    sample_obj,
    cells = sample_label_table$local_cell_id
  )
  
  sample_label_table <- sample_label_table %>%
    dplyr::mutate(
      local_cell_id = factor(.data$local_cell_id, levels = colnames(sample_obj))
    ) %>%
    dplyr::arrange(.data$local_cell_id)
  
  if (!identical(colnames(sample_obj), as.character(sample_label_table$local_cell_id))) {
    stop("Cell order mismatch for TX5K_", tx_number_i, call. = FALSE)
  }
  
  sample_obj$label_for_niche_assay <- sample_label_table$label_for_niche_assay
  
  set.seed(niche_random_seed)
  
  sample_obj <- Seurat::BuildNicheAssay(
    object = sample_obj,
    fov = xenium_fov_name,
    group.by = "label_for_niche_assay",
    neighbors = niche_n_neighbors
  )
  
  knn_count_matrix_sample <- t(sample_obj[["niche"]]@counts)
  rownames(knn_count_matrix_sample) <- sample_label_table$full_id
  
  knn_count_matrix_list[[tx_number_i]] <- knn_count_matrix_sample
  
}

knn_count_matrix_all <- bind_sparse_matrices_by_column(knn_count_matrix_list)


## Save outputs --------------------------------------------------------------

saveRDS(
  knn_count_matrix_all,
  file = file.path(
    niche_object_dir,
    "spatial_niche_knn_count_matrix.rds"
  )
)

message("Saved spatial niche kNN matrix with ", nrow(knn_count_matrix_all), " cells and ", ncol(knn_count_matrix_all), " labels.")
