# Xenium gene symbol-to-Entrez ID mapping
#
# This script creates a gene identifier mapping table for genes retained in
# the Xenium RNA assay used for downstream analyses. Gene symbols are mapped
# to Entrez IDs using org.Hs.eg.db, with manually curated mappings applied only
# to genes that remain unmapped. The established TEC mapping to Entrez ID 7006
# is retained, while the secondary TEC mapping to Entrez ID 100124696 is
# removed.


## Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})


## Load settings -------------------------------------------------------------

project_dir <- getwd()

if (basename(project_dir) != "PDAC_CAF_hypoxia_analysis") {
  stop(
    "Please open PDAC_CAF_hypoxia_analysis.Rproj before running this script."
  )
}

source(file.path(project_dir, "code", "00_config", "00_02_paths.R"))

fs::dir_create(xenium_table_dir)


## Read Xenium gene symbols --------------------------------------------------

xenium_gene_table <- read.csv(
  file.path(
    xenium_table_dir,
    "xenium_gene_symbols.csv"
  ),
  stringsAsFactors = FALSE
)

required_gene_table_columns <- "gene_symbol"
missing_gene_table_columns <- setdiff(
  required_gene_table_columns,
  colnames(xenium_gene_table)
)

if (length(missing_gene_table_columns) > 0) {
  stop(
    paste0(
      "The Xenium gene table is missing the following required column: ",
      paste(missing_gene_table_columns, collapse = ", "),
      "."
    ),
    call. = FALSE
  )
}

if (nrow(xenium_gene_table) == 0) {
  stop(
    "The Xenium gene table contains no genes.",
    call. = FALSE
  )
}

if (anyNA(xenium_gene_table$gene_symbol)) {
  stop(
    "Missing gene symbols were found in the Xenium gene table.",
    call. = FALSE
  )
}

if (any(xenium_gene_table$gene_symbol == "")) {
  stop(
    "Empty gene symbols were found in the Xenium gene table.",
    call. = FALSE
  )
}

if (anyDuplicated(xenium_gene_table$gene_symbol)) {
  stop(
    "Duplicated gene symbols were found in the Xenium gene table.",
    call. = FALSE
  )
}

xenium_gene_symbols <- xenium_gene_table$gene_symbol


## Define manually curated Entrez ID mappings -------------------------------

# These manually curated mappings are used only when org.Hs.eg.db does not
# return an Entrez ID for a gene retained in the Xenium analysis object.
manual_entrez_id_by_gene_symbol <- c(
  "H3F3B" = "3021",
  "ARNTL" = "406",
  "H2AFX" = "3014",
  "TMEM173" = "340061",
  "WARS" = "7453",
  "SPATA5L1" = "79029",
  "SPATA5" = "166378",
  "EPRS" = "2058",
  "SLC9A3R2" = "9351",
  "CBWD2" = "150472",
  "HSPB11" = "51668",
  "FAM126B" = "285172",
  "SLC9A3R1" = "9368",
  "GBA" = "2629",
  "PHB" = "5245",
  "EEF1AKNMT" = "51603",
  "CBSL" = "102724560",
  "CCDC113" = "29070",
  "TDGF1" = "6997",
  "GPR1" = "2825",
  "CYHR1" = "50626",
  "CD3EAP" = "10849",
  "DUSP13" = "128854680",
  "DDX58" = "23586",
  "BVES" = "11149",
  "ZBED9" = "114821",
  "ACPP" = "55"
)


## Map gene symbols to Entrez IDs using org.Hs.eg.db -------------------------

orgdb_gene_identifier_mapping <- clusterProfiler::bitr(
  xenium_gene_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db,
  drop = FALSE
) %>%
  tibble::as_tibble() %>%
  dplyr::mutate(
    ENTREZID = as.character(.data$ENTREZID)
  ) %>%
  dplyr::distinct(
    .data$SYMBOL,
    .data$ENTREZID
  )


## Validate and resolve duplicated database mappings ------------------------

# Confirm that TEC is the only symbol with more than one database-derived
# Entrez ID before removing its secondary mapping.
duplicated_orgdb_gene_symbols <- orgdb_gene_identifier_mapping %>%
  dplyr::filter(!is.na(.data$ENTREZID)) %>%
  dplyr::count(
    .data$SYMBOL,
    name = "number_of_entrez_ids"
  ) %>%
  dplyr::filter(.data$number_of_entrez_ids > 1)

expected_duplicated_orgdb_gene_symbols <- tibble::tibble(
  SYMBOL = "TEC",
  number_of_entrez_ids = 2L
)

if (!identical(
  duplicated_orgdb_gene_symbols,
  expected_duplicated_orgdb_gene_symbols
)) {
  stop(
    paste0(
      "Unexpected one-to-many SYMBOL-to-Entrez ID mappings were found. ",
      "Expected only TEC with two Entrez IDs."
    ),
    call. = FALSE
  )
}

tec_entrez_ids <- orgdb_gene_identifier_mapping %>%
  dplyr::filter(
    .data$SYMBOL == "TEC",
    !is.na(.data$ENTREZID)
  ) %>%
  dplyr::pull(.data$ENTREZID) %>%
  sort()

expected_tec_entrez_ids <- sort(c("7006", "100124696"))

if (!identical(tec_entrez_ids, expected_tec_entrez_ids)) {
  stop(
    paste0(
      "Unexpected Entrez IDs were returned for TEC. ",
      "Expected Entrez IDs 7006 and 100124696."
    ),
    call. = FALSE
  )
}

# Retain the established TEC mapping to Entrez ID 7006 and remove the
# secondary database entry excluded in the original Xenium GSEA workflow.
orgdb_gene_identifier_mapping <- orgdb_gene_identifier_mapping %>%
  dplyr::filter(
    !(
      .data$SYMBOL == "TEC" &
        .data$ENTREZID == "100124696"
    )
  )

duplicated_orgdb_gene_symbols_after_tec_removal <-
  orgdb_gene_identifier_mapping %>%
  dplyr::filter(!is.na(.data$ENTREZID)) %>%
  dplyr::count(
    .data$SYMBOL,
    name = "number_of_entrez_ids"
  ) %>%
  dplyr::filter(.data$number_of_entrez_ids > 1)

if (nrow(duplicated_orgdb_gene_symbols_after_tec_removal) > 0) {
  stop(
    "Multiple Entrez IDs remain after removing the secondary TEC mapping.",
    call. = FALSE
  )
}


## Add manually curated mappings --------------------------------------------

# Preserve every Xenium gene, and use manually curated IDs only for genes that
# remain unmapped after database conversion.
xenium_gene_identifier_mapping <- xenium_gene_table %>%
  dplyr::left_join(
    orgdb_gene_identifier_mapping %>%
      dplyr::rename(
        gene_symbol = SYMBOL,
        entrez_id = ENTREZID
      ),
    by = "gene_symbol"
  ) %>%
  dplyr::mutate(
    mapping_source = dplyr::if_else(
      !is.na(.data$entrez_id),
      "org.Hs.eg.db",
      NA_character_
    ),
    manually_curated_entrez_id = unname(
      manual_entrez_id_by_gene_symbol[.data$gene_symbol]
    ),
    entrez_id = dplyr::coalesce(
      .data$entrez_id,
      .data$manually_curated_entrez_id
    ),
    mapping_source = dplyr::case_when(
      !is.na(.data$mapping_source) ~ .data$mapping_source,
      !is.na(.data$manually_curated_entrez_id) ~ "manual",
      TRUE ~ "unmapped"
    )
  ) %>%
  dplyr::select(
    gene_symbol,
    entrez_id,
    mapping_source
  )


## Validate final mapping table ---------------------------------------------

if (nrow(xenium_gene_identifier_mapping) != nrow(xenium_gene_table)) {
  stop(
    "The number of rows changed while creating the Xenium gene mapping table.",
    call. = FALSE
  )
}

if (anyDuplicated(xenium_gene_identifier_mapping$gene_symbol)) {
  stop(
    "Duplicated gene symbols remain in the final Xenium gene mapping table.",
    call. = FALSE
  )
}

if (!identical(
  xenium_gene_identifier_mapping$gene_symbol,
  xenium_gene_table$gene_symbol
)) {
  stop(
    "Gene order changed while creating the Xenium gene mapping table.",
    call. = FALSE
  )
}


## Save Xenium gene identifier mapping --------------------------------------

write.csv(
  xenium_gene_identifier_mapping,
  file = file.path(
    xenium_table_dir,
    "xenium_gene_symbol_entrez_mapping.csv"
  ),
  row.names = FALSE,
  na = ""
)
