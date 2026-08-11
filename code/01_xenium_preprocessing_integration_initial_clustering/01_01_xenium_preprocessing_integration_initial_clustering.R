# Xenium preprocessing, integration, initial clustering, and cell type annotation
#
# This script contains the organized analysis code used for Xenium preprocessing,
# QC filtering, integration, initial clustering, and initial cell type annotation
# of the integrated Xenium spatial transcriptomic dataset.


## Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(harmony)
  library(future)
})


## Future settings -----------------------------------------------------------

future::plan("sequential")
options(future.globals.maxSize = 10 * 1024^3)


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


## Xenium preprocessing, integration, and initial clustering ------------------

# Xenium output files were loaded using LoadXenium. Cells with zero Xenium
# counts were removed, and the analyzed tissue region was defined by
# sample-specific in silico trimming. Trimmed count matrices were converted to
# Seurat objects, QC-filtered per sample, normalized, merged, integrated using
# Harmony, and initially clustered.

xenium_sample_mapping_file <- file.path(
  xenium_input_dir,
  "metadata",
  "xenium_sample_mapping_metadata.csv"
)

xenium_sample_mapping_metadata <- read.csv(
  xenium_sample_mapping_file,
  stringsAsFactors = FALSE
)

mapping_included <- subset(
  xenium_sample_mapping_metadata,
  included_in_analysis
)

mapping_included <- mapping_included[
  match(tx_id_merge_order, mapping_included$tx_id),
]

stopifnot(identical(mapping_included$tx_id, tx_id_merge_order))

fs::dir_create(xenium_object_dir)
fs::dir_create(xenium_table_dir)


## Create QC-filtered Seurat objects for each Xenium sample ------------------

create_qc_filtered_seurat <- function(
    mapping_row,
    nfeature_range,
    ncount_range
) {
  
  tx_id <- mapping_row$tx_id
  xenium_sample_id <- mapping_row$xenium_sample_id
  raw_data_folder <- mapping_row$raw_data_folder
  
  x_min <- mapping_row$x_min
  x_max <- mapping_row$x_max
  y_min <- mapping_row$y_min
  y_max <- mapping_row$y_max
  
  message("Loading ", xenium_sample_id, " from ", raw_data_folder)
  
  xenium_object <- LoadXenium(
    data.dir = file.path(xenium_raw_data_dir, raw_data_folder),
    fov = "fov"
  )
  
  ## Add spatial coordinates to metadata --------------------------------------
  
  xy_data <- as.data.frame(
    cbind(
      cell = xenium_object@images$fov$centroids@cells,
      xenium_object@images$fov$centroids@coords
    )
  )
  
  rownames(xy_data) <- xy_data$cell
  
  xy_data$X <- as.numeric(xy_data$x)
  xy_data$Y <- as.numeric(xy_data$y)
  
  xenium_object$X <- xy_data[colnames(xenium_object), "X"]
  xenium_object$Y <- xy_data[colnames(xenium_object), "Y"]
  
  xenium_object <- subset(
    xenium_object,
    subset = nCount_Xenium > 0
  )
  
  xenium_object_trimmed <- subset(
    xenium_object,
    subset = X > x_min & X < x_max & Y > y_min & Y < y_max
  )
  
  counts <- GetAssayData(
    object = xenium_object_trimmed,
    assay = "Xenium",
    layer = "counts"
  )
  
  seurat_object <- CreateSeuratObject(
    counts = counts,
    assay = "RNA",
    project = xenium_sample_id,
    min.cells = 3,
    min.features = 0
  )
  
  seurat_object$X <- xenium_object_trimmed@meta.data[colnames(seurat_object), "X"]
  seurat_object$Y <- xenium_object_trimmed@meta.data[colnames(seurat_object), "Y"]
  seurat_object$tx_id <- tx_id
  seurat_object$xenium_sample_id <- xenium_sample_id
  seurat_object$raw_data_folder <- raw_data_folder
  
  n_cells_after_trimming <- ncol(xenium_object_trimmed)
  n_cells_before_qc <- ncol(seurat_object)
  
  seurat_object <- subset(
    seurat_object,
    subset =
      nFeature_RNA > nfeature_range[1] &
      nFeature_RNA < nfeature_range[2] &
      nCount_RNA > ncount_range[1] &
      nCount_RNA < ncount_range[2]
  )
  
  xenium_object_trimmed_qc <- subset(
    xenium_object_trimmed,
    cells = colnames(seurat_object)
  )
  
  fs::dir_create(xenium_coordinate_object_dir)
  
  saveRDS(
    xenium_object_trimmed_qc,
    file = file.path(
      xenium_coordinate_object_dir,
      paste0(
        "xenium_coordinate_object_",
        stringr::str_to_lower(xenium_sample_id),
        ".rds"
      )
    )
  )
  
  seurat_object <- NormalizeData(
    seurat_object,
    verbose = FALSE
  )
  
  message(
    xenium_sample_id,
    ": ",
    format(n_cells_after_trimming, big.mark = ","),
    " cells after trimming; ",
    format(n_cells_before_qc, big.mark = ","),
    " cells before QC; ",
    format(ncol(seurat_object), big.mark = ","),
    " cells after QC"
  )
  
  return(seurat_object)
}

xenium_sample_seurat_object_list <- lapply(
  seq_len(nrow(mapping_included)),
  function(i) {
    create_qc_filtered_seurat(
      mapping_row = mapping_included[i, ],
      nfeature_range = qc_nfeature_rna,
      ncount_range = qc_ncount_rna
    )
  }
)

stopifnot(identical(
  unname(vapply(
    xenium_sample_seurat_object_list, 
    function(x) unique(x$tx_id), 
    character(1))),
  tx_id_merge_order
))

names(xenium_sample_seurat_object_list) <- vapply(
  xenium_sample_seurat_object_list, 
  function(x) unique(x$xenium_sample_id), 
  character(1))

xenium_sample_seurat_object_list <- lapply(
  xenium_sample_seurat_object_list,
  function(obj) {
    sample_id <- unique(obj$xenium_sample_id)
    stopifnot(length(sample_id) == 1)
    
    RenameCells(
      obj,
      new.names = paste(
        sample_id,
        colnames(obj),
        sep = "_"
      )
    )
  }
)

stopifnot(length(xenium_sample_seurat_object_list) == nrow(mapping_included))
stopifnot(identical(names(xenium_sample_seurat_object_list), mapping_included$xenium_sample_id))

raw_initial_object <- Reduce(
  function(x, y) merge(x, y = y, add.cell.ids = NULL),
  xenium_sample_seurat_object_list
)

message("Merged cells: ", format(ncol(raw_initial_object), big.mark = ","))

table(raw_initial_object$xenium_sample_id)
table(raw_initial_object$tx_id)

qc_summary <- mapping_included %>%
  dplyr::select(
    tx_id,
    xenium_sample_id,
    raw_data_folder,
    x_min,
    x_max,
    y_min,
    y_max
  ) %>%
  dplyr::mutate(
    n_cells_after_qc = vapply(xenium_sample_seurat_object_list, ncol, numeric(1))
  )

write.csv(
  qc_summary,
  file = file.path(xenium_table_dir, "xenium_qc_summary.csv"),
  row.names = FALSE
)


## Variable feature selection, scaling, and PCA ------------------------------

raw_initial_object <- raw_initial_object %>%
  FindVariableFeatures(
    selection.method = "vst",
    nfeatures = 2000
  )

raw_initial_object <- ScaleData(
  raw_initial_object,
  features = VariableFeatures(raw_initial_object)
)

set.seed(123)

raw_initial_object <- RunPCA(
  raw_initial_object,
  features = VariableFeatures(raw_initial_object),
  npcs = 50,
  verbose = FALSE
)


## Harmony integration -------------------------------------------------------

set.seed(123)

raw_initial_object <- RunHarmony(
  raw_initial_object,
  group.by.vars = "xenium_sample_id",
  reduction.use = "pca",
  dims.use = 1:50,
  assay.use = "RNA"
)

saveRDS(
  raw_initial_object,
  file = file.path(
    xenium_object_dir,
    "xenium_initial_harmony_raw.rds"
  )
)


## Initial clustering --------------------------------------------------------

set.seed(123)

raw_initial_object <- RunUMAP(
  raw_initial_object,
  reduction = "harmony",
  dims = 1:Dim1
) %>%
  FindNeighbors(
    reduction = "harmony",
    dims = 1:Dim1
  )

set.seed(123)

raw_initial_object <- FindClusters(
  raw_initial_object,
  resolution = Res1
)

## Save raw reconstructed object ---------------------------------------------

saveRDS(
  raw_initial_object,
  file = file.path(
    xenium_object_dir,
    "xenium_initial_clustering_raw.rds"
  )
)

## Load manuscript object for downstream analyses ----------------------------

# Note on reproducibility:
# Raw reprocessing reproduced the manuscript analysis through normalized data,
# variable feature selection, scaling, PCA, and Harmony embedding generation.
# However, FindNeighbors generated an SNN graph that differed from the
# manuscript object. Because clustering was reproducible when the manuscript
# SNN graph was used, the manuscript initial clustering object is used as the
# input for downstream analyses.
#
# The manuscript initial clustering object should be downloaded from GEO and
# placed in inputs/xenium/objects/ before running this section.

manuscript_initial_object <- readRDS(
  file.path(
    xenium_input_dir,
    "objects",
    "xenium_initial_clustering_manuscript.rds"
  )
)


## Export genes retained in the downstream Xenium object ---------------------

xenium_gene_table <- tibble::tibble(
  gene_symbol = rownames(
    manuscript_initial_object[["RNA"]]
  )
)

if (nrow(xenium_gene_table) == 0) {
  stop(
    "No genes were found in the RNA assay of the Xenium object.",
    call. = FALSE
  )
}

if (anyNA(xenium_gene_table$gene_symbol)) {
  stop(
    "Missing gene symbols were found in the Xenium gene list.",
    call. = FALSE
  )
}

if (anyDuplicated(xenium_gene_table$gene_symbol)) {
  stop(
    "Duplicated gene symbols were found in the Xenium gene list.",
    call. = FALSE
  )
}

readr::write_csv(
  xenium_gene_table,
  file.path(
    xenium_table_dir,
    "xenium_gene_symbols.csv"
  )
)


## Rename cell IDs to public sample IDs --------------------------------------

tx_cellnames <- colnames(manuscript_initial_object)

tx_prefix_to_xenium_sample_id <- stats::setNames(
  object = unname(xenium_sample_id_by_tx_number),
  nm = paste0("TX5K", names(xenium_sample_id_by_tx_number))
)

xenium_cellnames <- tx_cellnames

for (tx_prefix in names(tx_prefix_to_xenium_sample_id)) {
  xenium_cellnames <- sub(
    pattern = paste0("^", tx_prefix, "_"),
    replacement = paste0(tx_prefix_to_xenium_sample_id[[tx_prefix]], "_"),
    x = xenium_cellnames
  )
}

manuscript_initial_object <- RenameCells(
  manuscript_initial_object,
  new.names = xenium_cellnames
)

stopifnot(all(grepl("^Xenium_", colnames(manuscript_initial_object))))


## Initial cell type annotation ----------------------------------------------

manuscript_initial_object$initial_cell_type <- unname(
  initial_cluster_annotation[as.character(manuscript_initial_object$seurat_clusters)]
)

stopifnot(!any(is.na(manuscript_initial_object$initial_cell_type)))

manuscript_initial_object$initial_cell_type <- factor(
  manuscript_initial_object$initial_cell_type,
  levels = names(initial_cluster_annotation_list)
)

table(manuscript_initial_object$initial_cell_type, useNA = "ifany")


## Standardize metadata for downstream analyses ------------------------------

manuscript_initial_object$cell_id <- colnames(manuscript_initial_object)

manuscript_initial_object$xenium_sample_number <- stringr::str_extract(
  manuscript_initial_object$cell_id,
  "(?<=Xenium_)\\d+"
)

manuscript_initial_object$xenium_sample_id <- paste0(
  "Xenium_",
  manuscript_initial_object$xenium_sample_number
)

manuscript_initial_object$tx_id <- manuscript_initial_object$orig.ident

stopifnot(setequal(
  unique(manuscript_initial_object$tx_id),
  tx_id_merge_order
))

manuscript_initial_object$initial_cluster <- as.character(
  manuscript_initial_object$seurat_clusters
)

stopifnot(!any(is.na(manuscript_initial_object$xenium_sample_id)))
stopifnot(!any(is.na(manuscript_initial_object$initial_cell_type)))
stopifnot(all(manuscript_initial_object$cell_id == colnames(manuscript_initial_object)))


## Create initial clustering metadata table ----------------------------------

initial_clustering_metadata <- manuscript_initial_object@meta.data %>% 
  tibble::rownames_to_column(var = "cell_id_from_rownames") %>% 
  dplyr::mutate(
    umap_1 = unname(Embeddings(manuscript_initial_object, "umap")[, 1]),
    umap_2 = unname(Embeddings(manuscript_initial_object, "umap")[, 2])
  ) %>% 
  dplyr::select(
    cell_id = cell_id_from_rownames,
    xenium_sample_id,
    xenium_sample_number,
    tx_id,
    x = X,
    y = Y,
    initial_cluster,
    initial_cell_type,
    umap_1,
    umap_2
  )


## Save annotated initial clustering object ----------------------------------

initial_clustering_table_dir <- file.path(
  xenium_table_dir,
  "initial_clustering"
)

fs::dir_create(initial_clustering_table_dir)

saveRDS(
  manuscript_initial_object,
  file = file.path(
    xenium_object_dir,
    "xenium_initial_clustering_annotated.rds"
  )
)

write.csv(
  initial_clustering_metadata,
  file = file.path(
    initial_clustering_table_dir,
    "xenium_initial_clustering_metadata.csv"
  ),
  row.names = FALSE
)
