## Gene sets and marker definitions -----------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})


## Xenium panel genes --------------------------------------------------------

all_genes_table <- readr::read_csv(
  file = file.path(
    gene_set_input_dir, 
    "xenium_5k_panel_add_entrez.csv"
  ),
  show_col_types = FALSE
) |>
  dplyr::mutate(
    ENTREZID = as.character(ENTREZID)
  )

all_genes_symbol <- all_genes_table$SYMBOL
all_genes_entrez <- all_genes_table$ENTREZID


## Marker and gene set table -------------------------------------------------

marker_table <- read.csv(
  file = file.path(
    gene_set_input_dir, 
    "xenium_gene_sets.csv"),
  header = TRUE
)

myCAF_xenium_panel <- marker_table$gene[marker_table$set_name == "myCAF"]

iCAF_xenium_panel <- marker_table$gene[marker_table$set_name == "iCAF"]

buffa_hypoxia_genes <- marker_table |>
  dplyr::filter(
    set_name == "Buffa" &
      gene %in% all_genes_symbol
  ) |> 
  dplyr::pull(gene)

winter_hypoxia_genes <- marker_table |>
  dplyr::filter(
    set_name == "Winter" &
      gene %in% all_genes_symbol
    ) |>
  dplyr::pull(gene) 


## Human Hallmark gene sets --------------------------------------------------

get_hallmark_term2gene <- function() {
  msigdbr::msigdbr(
    species = "Homo sapiens",
    category = "H"
  ) |>
    dplyr::select(
      term = gs_name,
      gene = entrez_gene
    ) |>
    dplyr::filter(!is.na(.data$gene)) |>
    dplyr::mutate(gene = as.character(.data$gene)) |>
    dplyr::distinct(.data$term, .data$gene)
}
