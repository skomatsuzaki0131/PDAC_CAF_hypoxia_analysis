## Project root --------------------------------------------------------------

project_dir <- getwd()

if (basename(project_dir) != "PDAC_CAF_hypoxia_analysis") {
  stop(
    "Please open PDAC_CAF_hypoxia_analysis.Rproj before running this script."
  )
}


## Input directories ---------------------------------------------------------

input_dir <- file.path(project_dir, "inputs")

gene_set_input_dir <- file.path(input_dir, "gene_sets")
xenium_input_dir <- file.path(input_dir, "xenium")
caf_isolation_input_dir <- file.path(input_dir, "caf_isolation")
bulk_rnaseq_input_dir <- file.path(input_dir, "bulk_rnaseq")

xenium_raw_folder_by_tx_id <- c(
  "TX5K_01" = "TX5K_01",
  "TX5K_02" = "TX5K_02",
  "TX5K_11" = "TX5K_11",
  "TX5K_15" = "TX5K_15_16",
  "TX5K_16" = "TX5K_15_16",
  "TX5K_19" = "TX5K_18_19_20_22"
)


## External large data -------------------------------------------------------

# Raw Xenium output directories are not included in this GitHub repository.
# Local paths are defined in 00_00_local_paths.R, which is excluded from Git.

local_paths_file <- file.path(
  project_dir,
  "code",
  "00_config",
  "00_00_local_paths.R"
)

if (file.exists(local_paths_file)) {
  source(local_paths_file)
}

if (!exists("xenium_raw_data_dir")) {
  stop(
    paste0(
      "xenium_raw_data_dir is not defined.\n",
      "Create code/00_config/00_00_local_paths.R and define the local ",
      "directory containing the downloaded Xenium raw data."
    ),
    call. = FALSE
  )
}

if (!dir.exists(xenium_raw_data_dir)) {
  stop(
    paste0(
      "The Xenium raw data directory was not found: ",
      xenium_raw_data_dir
    ),
    call. = FALSE
  )
}

## Output directories --------------------------------------------------------

output_dir <- file.path(project_dir, "outputs")

## Xenium output directories -------------------------------------------------

xenium_output_dir <- file.path(output_dir, "xenium_outputs")
xenium_object_dir <- file.path(xenium_output_dir, "objects")
xenium_table_dir <- file.path(xenium_output_dir, "tables")
xenium_figure_dir <- file.path(xenium_output_dir, "figures")

xenium_coordinate_object_dir <- file.path(
  xenium_object_dir,
  "xenium_objects_with_coordinates"
)

initial_clustering_figure_dir <- file.path(
  xenium_figure_dir,
  "08_01_initial_clustering_figures"
)

caf_subclustering_figure_dir <- file.path(
  xenium_figure_dir,
  "08_02_caf_subclustering_figures"
)

hypoxia_mapping_figure_dir <- file.path(
  xenium_figure_dir,
  "08_03_hypoxia_mapping_figures"
)

niche_analysis_figure_dir <- file.path(
  xenium_figure_dir,
  "08_04_niche_analysis_figures"
)

pseudotime_analysis_figure_dir <- file.path(
  xenium_figure_dir,
  "08_05_pseudotime_analysis_figures"
)


## CAF isolation assay output directories ------------------------------------

caf_isolation_output_dir <- file.path(output_dir, "caf_isolation_outputs")
caf_isolation_table_dir <- file.path(caf_isolation_output_dir, "tables")
caf_isolation_figure_dir <- file.path(caf_isolation_output_dir, "figures")


## Bulk RNA-seq output directories -------------------------------------------

bulk_rnaseq_output_dir <- file.path(output_dir, "bulk_rnaseq_outputs")
bulk_rnaseq_object_dir <- file.path(bulk_rnaseq_output_dir, "objects")
bulk_rnaseq_table_dir <- file.path(bulk_rnaseq_output_dir, "tables")
bulk_rnaseq_figure_dir <- file.path(bulk_rnaseq_output_dir, "figures")


## Frequently used files -----------------------------------------------------

initial_metadata_file <- file.path(
  xenium_table_dir,
  "initial_clustering",
  "xenium_initial_clustering_metadata.csv"
)

initial_object_file <- file.path(
  xenium_object_dir,
  "xenium_initial_clustering_annotated.rds"
)

caf_object_file <- file.path(
  xenium_object_dir,
  "xenium_caf_subclustering_annotated.rds"
)

caf_subcluster_metadata_file <- file.path(
  xenium_table_dir,
  "caf_subclustering",
  "xenium_caf_subclustering_metadata.csv"
)

caf_marker_table_file <- file.path(
  xenium_table_dir,
  "caf_subclustering",
  "xenium_caf_subcluster_markers_all.csv"
)

caf8_signature_genes_table_file <- file.path(
  gene_set_input_dir,
  "xenium_caf8_signature_genes.csv"
)

caf_isolation_metadata_file <- file.path(
  caf_isolation_input_dir,
  "metadata",
  "caf_isolation_metadata.csv"
)

bulk_rnaseq_tpm_file <- file.path(
  bulk_rnaseq_input_dir,
  "expression_data",
  "tpm.csv"
)