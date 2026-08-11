# Prepare Winter hypoxia ssGSEA scores
#
# This script calculates cell-level Winter hypoxia ssGSEA scores separately
# for each Xenium sample and standardizes the scores within each initial cell
# type in each sample.


## Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(GSVA)
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
source(file.path(project_dir, "code", "00_config", "00_03_genes_and_marker_definitions.R"))
source(file.path(project_dir, "code", "00_config", "00_04_annotation_definitions.R"))


## Define input and output paths --------------------------------------------

hypoxia_mapping_table_dir <- file.path(
  xenium_table_dir,
  "hypoxia_mapping"
)

fs::dir_create(hypoxia_mapping_table_dir)


## Helper functions ----------------------------------------------------------

read_required_csv <- function(file, description) {
  if (!file.exists(file)) {
    stop(description, " not found: ", file, call. = FALSE)
  }
  
  readr::read_csv(file, show_col_types = FALSE)
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

prepare_initial_metadata_for_hypoxia_mapping <- function(metadata_table) {
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
      initial_cell_type = as.character(.data$initial_cell_type),
      initial_cluster_label = paste0("Ini-", .data$initial_cluster)
    )
}

prepare_caf_subcluster_metadata_for_hypoxia_mapping <- function(metadata_table) {
  metadata_table %>%
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
}

add_object_cell_ids <- function(sample_label_table, sample_obj) {
  local_matches <- sample_label_table$local_cell_id %in% colnames(sample_obj)
  full_id_matches <- sample_label_table$full_id %in% colnames(sample_obj)
  
  if (sum(local_matches) >= sum(full_id_matches) && any(local_matches)) {
    sample_label_table %>%
      dplyr::mutate(object_cell_id = .data$local_cell_id)
  } else if (any(full_id_matches)) {
    sample_label_table %>%
      dplyr::mutate(object_cell_id = .data$full_id)
  } else {
    stop(
      "No cells in the metadata table matched the Seurat object cell names.",
      call. = FALSE
    )
  }
}

scale_numeric_vector <- function(x) {
  as.numeric(scale(x))
}


## Read and prepare cell metadata -------------------------------------------

initial_metadata_table <- read_required_csv(
  initial_metadata_file,
  "Initial clustering metadata table"
)

caf_subcluster_metadata_table <- read_required_csv(
  caf_subcluster_metadata_file,
  "CAF subcluster metadata table"
)

check_required_columns(
  initial_metadata_table,
  c(
    "cell_id",
    "xenium_sample_id",
    "tx_id",
    "x",
    "y",
    "initial_cluster",
    "initial_cell_type"
  ),
  "initial clustering metadata table"
)

check_required_columns(
  caf_subcluster_metadata_table,
  c(
    "cell_id",
    "xenium_sample_id",
    "tx_id",
    "caf_subcluster"
  ),
  "CAF subcluster metadata table"
)

initial_metadata_table <- prepare_initial_metadata_for_hypoxia_mapping(
  initial_metadata_table
)

caf_subcluster_metadata_table <- prepare_caf_subcluster_metadata_for_hypoxia_mapping(
  caf_subcluster_metadata_table
)

if (any(duplicated(initial_metadata_table$full_id))) {
  stop("Duplicated full_id values were detected in initial metadata.", call. = FALSE)
}

if (any(duplicated(caf_subcluster_metadata_table$full_id))) {
  stop("Duplicated full_id values were detected in CAF subcluster metadata.", call. = FALSE)
}

cell_label_table <- initial_metadata_table %>%
  dplyr::filter(.data$initial_cell_type != "Unclassified") %>%
  dplyr::left_join(
    caf_subcluster_metadata_table,
    by = "full_id"
  ) %>%
  dplyr::mutate(
    SubClust = dplyr::if_else(
      is.na(.data$caf_subcluster),
      .data$initial_cell_type,
      .data$caf_subcluster
    )
  )


## Define Winter hypoxia gene set -------------------------------------------

winter_hypoxia_gene_set <- list(
  WinterOrig = winter_hypoxia_genes
)


## Calculate Winter hypoxia ssGSEA scores by sample --------------------------

for (tx_number_i in tx_number_merge_order) {
  xenium_sample_id_i <- xenium_sample_id_by_tx_number[[tx_number_i]]
  
  if (is.null(xenium_sample_id_i)) {
    stop(
      "No Xenium sample ID was defined for TX number: ",
      tx_number_i,
      call. = FALSE
    )
  }
  
  message("Preparing Winter hypoxia ssGSEA scores for ", xenium_sample_id_i)
  
  sample_obj <- read_xenium_coordinate_object(xenium_sample_id_i)
  
  sample_label_table <- cell_label_table %>%
    dplyr::filter(.data$tx_number == tx_number_i) %>%
    add_object_cell_ids(sample_obj) %>%
    dplyr::filter(.data$object_cell_id %in% colnames(sample_obj))
  
  if (nrow(sample_label_table) == 0) {
    stop("No analyzable cells were found for TX5K_", tx_number_i, call. = FALSE)
  }
  
  sample_obj <- subset(
    sample_obj,
    cells = sample_label_table$object_cell_id
  )
  
  sample_label_table <- sample_label_table %>%
    dplyr::slice(match(colnames(sample_obj), .data$object_cell_id))
  
  if (!identical(colnames(sample_obj), sample_label_table$object_cell_id)) {
    stop(
      "Cell order mismatch after subsetting sample object for TX5K_",
      tx_number_i,
      call. = FALSE
    )
  }
  
  sample_obj$IniClust <- sample_label_table$initial_cluster_label
  sample_obj$Annot <- sample_label_table$initial_cell_type
  sample_obj$SubClust <- sample_label_table$SubClust
  
  sample_obj <- Seurat::NormalizeData(
    sample_obj,
    assay = "Xenium",
    verbose = FALSE
  )
  
  expression_matrix <- Seurat::GetAssayData(
    sample_obj,
    assay = "Xenium",
    layer = "data"
  ) %>%
    as.matrix()
  
  winter_hypoxia_gene_set_filtered <- list(
    WinterOrig = intersect(winter_hypoxia_gene_set$WinterOrig, rownames(expression_matrix))
  )
  
  if (length(winter_hypoxia_gene_set_filtered$WinterOrig) == 0) {
    stop(
      "No Winter hypoxia genes were found in the expression matrix for TX5K_",
      tx_number_i,
      call. = FALSE
    )
  }
  
  ssgsea_param <- GSVA::ssgseaParam(
    exprData = expression_matrix,
    geneSets = winter_hypoxia_gene_set_filtered,
    minSize = 1,
    maxSize = Inf,
    normalize = TRUE
  )
  
  ssgsea_matrix <- GSVA::gsva(
    ssgsea_param,
    verbose = TRUE
  )
  
  ssgsea_score_table <- as.data.frame(t(ssgsea_matrix)) %>%
    tibble::rownames_to_column(var = "object_cell_id")
  
  colnames(ssgsea_score_table)[-1] <- paste0(
    colnames(ssgsea_score_table)[-1],
    "_ssGSEA"
  )
  
  hypoxia_ssgsea_table <- sample_label_table %>%
    dplyr::select(
      sample = tx_number,
      xenium_sample_id,
      tx_id,
      full_id,
      cell_id = local_cell_id,
      object_cell_id,
      IniClust = initial_cluster_label,
      Annot = initial_cell_type,
      SubClust,
      X = x,
      Y = y
    ) %>%
    dplyr::left_join(
      ssgsea_score_table,
      by = "object_cell_id"
    )
  
  if (any(is.na(hypoxia_ssgsea_table$WinterOrig_ssGSEA))) {
    stop(
      "Some cells are missing Winter ssGSEA scores for TX5K_",
      tx_number_i,
      call. = FALSE
    )
  }
  
  hypoxia_score_table <-
    hypoxia_score_table %>%
    dplyr::group_by(sample, Annot) %>%
    dplyr::mutate(
      Winter_z_celltype = as.numeric(scale(WinterOrig_ssGSEA))
    ) %>%
    dplyr::ungroup()
  
  readr::write_csv(
    hypoxia_ssgsea_normalized_table,
    file = file.path(
      hypoxia_mapping_table_dir,
      paste0(
        "winter_hypoxia_ssgsea_scores_",
        stringr::str_to_lower(xenium_sample_id_i),
        ".csv"
      )
    )
  )
}

message(
  "Saved Winter hypoxia ssGSEA score tables to: ",
  hypoxia_mapping_table_dir
)

