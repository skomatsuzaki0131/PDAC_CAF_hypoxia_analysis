# CAF extraction, preprocessing, and subclustering
#
# This script extracts CAFs from the annotated initial Xenium object and
# performs CAF subclustering using Harmony embeddings generated during
# the initial all-cell integration step. Harmony is not rerun after CAF
# extraction. Normalized CAF expression data are prepared after clustering
# for downstream analyses.

## Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(future)
  library(fs)
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

caf_subclustering_table_dir <- file.path(xenium_table_dir, "caf_subclustering")

fs::dir_create(xenium_object_dir)

fs::dir_create(caf_subclustering_table_dir)


## Extract CAFs from the annotated initial object ----------------------------

CombinedObj <- readRDS(initial_object_file)

stopifnot("initial_cell_type" %in% colnames(CombinedObj@meta.data))
stopifnot(all(grepl("^Xenium_", colnames(CombinedObj))))

CAFObj <- subset(
  CombinedObj,
  subset = initial_cell_type == "CAF"
)

if (ncol(CAFObj) == 0) {
  stop("No CAF cells were extracted. Check initial_cell_type annotation.")
}

message("Extracted CAFs: ", format(ncol(CAFObj), big.mark = ","))
print(table(CAFObj$initial_cell_type, useNA = "ifany"))

saveRDS(
  CAFObj,
  file = file.path(
    xenium_object_dir,
    "xenium_caf.rds"
  )
)


## CAF subclustering ---------------------------------------------------------

# CAF subclustering was performed using the Harmony embeddings generated
# during the initial all-cell integration step. Harmony was not rerun after
# CAF extraction.

stopifnot("harmony" %in% Reductions(CAFObj))

set.seed(123)

CAFObj <- RunUMAP(
  CAFObj,
  reduction = "harmony",
  dims = 1:Dim2
) %>%
  FindNeighbors(
    reduction = "harmony",
    dims = 1:Dim2
  )

CAFObj <- FindClusters(
  CAFObj,
  resolution = Res2, 
  random.seed = 123
)

stopifnot("RNA_snn_res.0.5" %in% colnames(CAFObj@meta.data))

CAFObj$caf_subcluster <- paste0(
  "CAF-",
  as.character(CAFObj$RNA_snn_res.0.5)
)

CAFObj$caf_subcluster <- factor(
  CAFObj$caf_subcluster,
  levels = paste0("CAF-", caf_subcluster_order)
)

stopifnot(!any(is.na(CAFObj$caf_subcluster)))

print(table(CAFObj$caf_subcluster, useNA = "ifany"))

CAFObj <- CAFObj %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(
    selection.method = "vst",
    nfeatures = 2000
  )

CAFObj <- ScaleData(
  CAFObj,
  features = VariableFeatures(CAFObj)
)

set.seed(123)

CAFObj <- RunPCA(
  CAFObj,
  features = VariableFeatures(CAFObj),
  npcs = 50,
  verbose = FALSE
)

saveRDS(
  CAFObj,
  file = file.path(
    xenium_object_dir,
    "xenium_caf_subclustering_raw.rds"
  )
)


## Standardize CAF metadata --------------------------------------------------

CAFObj$cell_id <- colnames(CAFObj)

CAFObj$xenium_sample_number <- stringr::str_extract(
  CAFObj$cell_id,
  "(?<=Xenium_)\\d+"
)

CAFObj$xenium_sample_id <- paste0(
  "Xenium_",
  CAFObj$xenium_sample_number
)

CAFObj$tx_id <- CAFObj$orig.ident

stopifnot(!any(is.na(CAFObj$caf_subcluster)))
stopifnot(!any(is.na(CAFObj$xenium_sample_id)))
stopifnot(all(CAFObj$cell_id == colnames(CAFObj)))


## Save annotated CAF subclustering object -----------------------------------

saveRDS(
  CAFObj,
  file = file.path(
    xenium_object_dir,
    "xenium_caf_subclustering_annotated.rds"
  )
)

## Create CAF subclustering metadata table -----------------------------------

initial_metadata_df <- read.csv(
  initial_metadata_file,
  stringsAsFactors = FALSE
)

caf_embedding_df <- tibble::tibble(
  cell_id = colnames(CAFObj),
  caf_subcluster = as.character(CAFObj$caf_subcluster),
  caf_umap_1 = unname(Embeddings(CAFObj, "umap")[, 1]),
  caf_umap_2 = unname(Embeddings(CAFObj, "umap")[, 2])
)

caf_metadata_df <- initial_metadata_df %>%
  dplyr::filter(cell_id %in% caf_embedding_df$cell_id) %>%
  dplyr::left_join(
    caf_embedding_df,
    by = "cell_id"
  ) %>%
  dplyr::select(
    cell_id,
    xenium_sample_id,
    xenium_sample_number,
    tx_id,
    x,
    y,
    initial_cluster,
    initial_cell_type,
    caf_subcluster,
    caf_umap_1,
    caf_umap_2
  )

stopifnot(nrow(caf_metadata_df) == ncol(CAFObj))
stopifnot(!any(is.na(caf_metadata_df$caf_subcluster)))

write.csv(
  caf_metadata_df,
  file = file.path(
    caf_subclustering_table_dir,
    "xenium_caf_subclustering_metadata.csv"
  ),
  row.names = FALSE
)
