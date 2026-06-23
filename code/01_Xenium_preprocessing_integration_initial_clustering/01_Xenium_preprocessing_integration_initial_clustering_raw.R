# Xenium preprocessing, integration, initial clustering, and cell type annotation
#
# This script contains the final unrefactored code used for preprocessing,
# QC filtering, integration, initial clustering, and initial cell type annotation
# of the integrated Xenium spatial transcriptomic dataset.
#
# Note: This script will be refactored before publication.


## Load packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(Seurat)
  library(harmony)
  library(presto)
  library(patchwork)
  library(cowplot)
  library(ggpubr)
  library(ggrepel)
  library(ggraph)
  library(RColorBrewer)
  library(gplots)
  library(scales)
  library(beepr)
})


## Input paths ---------------------------------------------------------------

# TODO before publication:
# Replace local absolute paths with project-relative paths or instructions
# for where each input file should be placed.

dir_dataset <- "/Volumes/PortableSSD/[4]myR/[1]dataset"
dir_xenium_raw <- "/Volumes/Extreme SSD/Data"
dir_xenium_data <- "/Volumes/Extreme SSD/Analysis/Data"
dir_analysis <- "/Volumes/Extreme SSD/Analysis/Data/IntegAnalysis"

DirRNAseq <- file.path(dir_dataset, "RNAseq")
DirGO <- file.path(dir_dataset, "RNAseq", "GO")
DirHallmark <- file.path(dir_dataset, "RNAseq", "Hallmark")

DirOutputGO <- "/Volumes/PortableSSD/[3]Graduate school/[3]RNA-seq/analysis/GeneOntonogy"
DirOutputMSigDB <- "/Volumes/PortableSSD/[3]Graduate school/[3]RNA-seq/analysis/MSigDB"


## Gene sets and marker definitions -----------------------------------------

Vec.AllGenes <- read.csv(
  file = file.path(dir_dataset, "Genes_Xenium5k.csv"),
  header = FALSE
)[, 1]

DF.AllGenes_addEntrez <- read.csv(
  file = file.path(dir_dataset, "Genes_Xenium5k_addEntrez.csv"),
  header = TRUE
)

DF.Markers <- read.csv(
  file = file.path(dir_dataset, "MyGeneSet.csv"),
  header = TRUE
) %>%
  dplyr::filter(!Gene %in% c("CD3D", "COL1A2"))

Vec.MarkerGenes_CellType <- subset(
  DF.Markers,
  subset = GeneType == "CellType" &
    Panel5k == "Included" &
    Panel5kQlt != "Duplicated" &
    Panel5kQlt != "MayBePoor"
)$Gene

Vec.CAFmarkers <- subset(
  DF.Markers,
  subset = GeneType == "CAFmarker"
)$Gene

Vec.BuffaOrig <- intersect(
  subset(DF.Markers, subset = Classification1 == "Buffa")$Gene,
  Vec.AllGenes
)

Vec.WinterOrig <- intersect(
  subset(DF.Markers, subset = Classification1 == "WinterCore")$Gene,
  Vec.AllGenes
)

Vec.MoffittActivated <- intersect(
  subset(DF.Markers, subset = Classification1 == "Activated")$Gene,
  Vec.AllGenes
)

Vec.MoffittNormal <- intersect(
  subset(DF.Markers, subset = Classification1 == "Normal")$Gene,
  Vec.AllGenes
)

Vec.Classical <- subset(
  DF.Markers,
  Classification1 == "Classical"
)$Gene

Vec.Basallike <- subset(
  DF.Markers,
  Classification1 == "Basal_like"
)$Gene

Vec.Mix300.NormoCAF <- subset(
  DF.Markers,
  subset = GeneType == "Mix3" & Classification1 == "NormoCAF"
)$Gene

Vec.Mix300.HypoCAF <- subset(
  DF.Markers,
  subset = GeneType == "Mix3" & Classification1 == "HypoCAF"
)$Gene


## Analysis parameters -------------------------------------------------------

TXnumInteg <- c("02", "19", "11", "16", "01", "15")
NumOfSamples <- 6
IntegArgo <- "Harmony"

Mag <- 1
nFeatRNA <- c(100, 900)
nCountRNA <- c(100, 1800)

Dim1 <- 20
Res1 <- 1.0

Dim2 <- 30
Res2 <- 0.5

Num.GridLen <- 200
Num.MinCells <- 10

prolifCAFclust <- "8"

n.neighbors <- 40
k.niche <- 10
pc_use <- 15

Unclassified <- "Exclude"
Labels <- "seurat_clusters"

method <- "SpearmanAverage"

SetName <- "Set1_3pair"
DataSetName <- "rmCAF6"
n_HVGs <- 1000


## Derived labels and output paths ------------------------------------------

QCInfo <- paste0(
  "Countable_mag", Mag,
  "_nFeatRNA:", paste(nFeatRNA, collapse = "~"),
  "_nCountRNA:", paste(nCountRNA, collapse = "~")
)

QCInfo.FileName <- str_replace_all(
  QCInfo,
  pattern = ":",
  replace = ""
)

DirInteg <- paste0(
  dir_analysis, "/",
  "[", NumOfSamples, "case(", paste(TXnumInteg[seq_len(NumOfSamples)], collapse = ","), ")]",
  "_nFeatRNA:", paste(nFeatRNA, collapse = "~"),
  "_nCountRNA:", paste(nCountRNA, collapse = "~"),
  "/"
)


## Initial cluster annotation settings ---------------------------------------

List.Annot <- list(
  Acinar = c(0, 13, 24),
  Ductal_like_acinar = c(33, 1, 12, 17, 30),
  Normal_ductal = 23,
  PanIN = c(10, 27),
  PDAC = 20,
  Islet = 14,
  CAF = c(9, 6, 7),
  Mural = 11,
  Endothelial = c(34, 3),
  Lymph_T = c(2, 16),
  Lymph_B = 15,
  Plasma = c(18, 37),
  Myeloid = c(5, 4, 31),
  Mast = 26,
  Nerve = 29,
  Unclassified = c(8, 19, 21, 22, 25, 28, 32, 35, 36)
)

Vec.Annot <- unlist(
  lapply(names(List.Annot), function(ct) {
    setNames(
      rep(ct, length(List.Annot[[ct]])),
      List.Annot[[ct]]
    )
  })
)

Vec.XeniumID <- c(
  "01" = "Xenium_01",
  "02" = "Xenium_02",
  "11" = "Xenium_03",
  "15" = "Xenium_04",
  "16" = "Xenium_05",
  "19" = "Xenium_06"
)


## Color palettes ------------------------------------------------------------

Vec.colorcode.Final <- c(
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

Vec.OrderedCnumToColor <- c(
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

Vec.SubclustOrder <- c(4, 8, 0, 1, 2, 5, 3, 7, 6)

cols2 <- c(
  "0" = "#4C78A8",
  "1" = "#54A24B",
  "2" = "#B279A2",
  "3" = "#9D755D",
  "4" = "#72B7B2",
  "5" = "#E6AF2E",
  "6" = "#F58518",
  "7" = "#9C9EDE",
  "8" = "#C00000"
)

Vec.Cols.Diverging <- c("#A6D3C5", "#F7F7F7", "#C05A9E")


## Helper functions ----------------------------------------------------------

clean_gs_label_ForFig <- function(x) {
  x %>%
    str_remove("^(HALLMARK_|KEGG_)") %>%
    str_replace_all("_", " ") %>%
    str_squish() %>%
    str_to_lower() %>%
    str_to_sentence() %>%
    str_replace_all("\\bTnfa\\b", "TNF-α") %>%
    str_replace_all("\\bIl6 jak stat3\\b", "IL-6/JAK/STAT3") %>%
    str_replace_all("\\bTgf beta\\b", "TGF-β") %>%
    str_replace_all("\\bMycaf\\b", "myCAF") %>%
    str_replace_all("\\bKras\\b", "KRAS") %>%
    str_replace_all("\\bP53\\b", "p53") %>%
    str_replace_all("\\bIcaf\\b", "iCAF") %>%
    str_replace_all("\\bDna\\b", "DNA") %>%
    str_replace_all("\\bUv\\b", "UV") %>%
    str_replace_all("\\bMyc\\b", "MYC") %>%
    str_replace_all("\\bWnt beta catenin\\b", "Wnt/β-catenin") %>%
    str_replace_all("\\bdn\\b", "down") %>%
    str_replace_all("\\bnfkb\\b", "NF-κB") %>%
    str_replace_all("\\bIl2 stat5\\b", "IL-2/STAT5") %>%
    str_replace_all("\\bG2m\\b", "G2/M") %>%
    str_replace_all("\\bE2f\\b", "E2F") %>%
    str_replace_all("\\bMtorc1\\b", "mTORC1") %>%
    str_replace_all("\\bPi3k akt mtor\\b", "PI3K-AKT-mTOR") %>%
    str_replace_all("targets v", "targets V") %>%
    str_replace_all("orig\\b", " hypoxia") %>%
    str_replace_all("\\biCAF signature\\b", "iCAF") %>%
    str_replace_all("\\biCAF\\b", "iCAF signature") %>%
    str_replace_all("\\bmyCAF signature\\b", "myCAF") %>%
    str_replace_all("\\bmyCAF\\b", "myCAF signature") %>%
    str_replace_all("\\bInterferon alpha\\b", "IFN-α") %>%
    str_replace_all("\\bInterferon gamma\\b", "IFN-γ") %>%
    str_replace_all("Epithelial mesenchymal transition", "EMT")
}

clean_gs_label_ForTable <- function(x) {
  x %>%
    stringi::stri_trans_general("Any-ASCII") %>%
    str_replace_all("[–—−]", "-") %>%
    str_remove("^(HALLMARK_|KEGG_)") %>%
    str_replace_all("_", " ") %>%
    str_squish() %>%
    str_to_lower() %>%
    str_to_sentence() %>%
    str_replace_all("\\bTnfa\\b", "TNF-alpha") %>%
    str_replace_all("\\bIl6 jak stat3\\b", "IL-6/JAK/STAT3") %>%
    str_replace_all("\\bTgf beta\\b", "TGF-beta") %>%
    str_replace_all("\\bMycaf\\b", "myCAF") %>%
    str_replace_all("\\bKras\\b", "KRAS") %>%
    str_replace_all("\\bP53\\b", "p53") %>%
    str_replace_all("\\bIcaf\\b", "iCAF") %>%
    str_replace_all("\\bDna\\b", "DNA") %>%
    str_replace_all("\\bUv\\b", "UV") %>%
    str_replace_all("\\bMyc\\b", "MYC") %>%
    str_replace_all("\\bWnt beta catenin\\b", "Wnt/beta-catenin") %>%
    str_replace_all("\\bdn\\b", "down") %>%
    str_replace_all("\\bnfkb\\b", "NF-kappaB") %>%
    str_replace_all("\\bIl2 stat5\\b", "IL-2/STAT5") %>%
    str_replace_all("\\bG2m\\b", "G2/M") %>%
    str_replace_all("\\bE2f\\b", "E2F") %>%
    str_replace_all("\\bMtorc1\\b", "mTORC1") %>%
    str_replace_all("\\bPi3k akt mtor\\b", "PI3K-AKT-mTOR") %>%
    str_replace_all("targets v", "targets V") %>%
    str_replace_all("orig\\b", " hypoxia") %>%
    str_replace_all("\\biCAF signature\\b", "iCAF") %>%
    str_replace_all("\\biCAF\\b", "iCAF signature") %>%
    str_replace_all("\\bmyCAF signature\\b", "myCAF") %>%
    str_replace_all("\\bmyCAF\\b", "myCAF signature") %>%
    str_replace_all("Epithelial mesenchymal transition", "EMT")
}


## Xenium preprocessing and initial clustering -------------------------------

# Xenium output files were loaded using LoadXenium. Cells with zero Xenium
# counts were removed, and the analyzed tissue region was defined by
# sample-specific in silico trimming. Trimmed count matrices were converted to
# Seurat objects, QC-filtered per sample, normalized, merged, integrated using
# Harmony, and initially clustered.

sample_roi <- tibble::tribble(
  ~TX,  ~Xmin,  ~Xmax,  ~Ymin,  ~Ymax,
  "01",   1029,  10269,    880,   5236,
  "02",    997,   5228,    837,  10791,
  "11",    640,   4345,    746,   4826,
  "15",    735,   7050,   7055,  10991,
  "16",    800,   9550,    900,   5700,
  "19",   5941,   9587,    313,   9303
)

sample_roi <- sample_roi %>%
  dplyr::filter(TX %in% TXnumInteg[seq_len(NumOfSamples)]) %>%
  dplyr::mutate(TX = factor(TX, levels = TXnumInteg[seq_len(NumOfSamples)])) %>%
  dplyr::arrange(TX) %>%
  dplyr::mutate(TX = as.character(TX))


## Create QC-filtered Seurat objects for each Xenium sample ------------------

create_qc_filtered_seurat <- function(TX, Xmin, Xmax, Ymin, Ymax,
                                      Mag, nFeatRNA, nCountRNA) {
  xen_obj <- LoadXenium(
    file.path(dir_xenium_raw, paste0("TX5K_", TX)),
    fov = "fov"
  )
  
  xen_obj <- subset(
    xen_obj,
    subset = nCount_Xenium > 0
  )
  
  xen_obj_trimmed <- subset(
    xen_obj,
    subset = X > Xmin & X < Xmax & Y > Ymin & Y < Ymax
  )
  
  counts <- GetAssayData(
    object = xen_obj_trimmed,
    assay = "Xenium",
    layer = "counts"
  )
  
  seu_obj <- CreateSeuratObject(
    counts = counts,
    assay = "RNA",
    project = paste0("TX5K_", TX),
    min.cells = 3,
    min.features = 0
  )
  
  seu_obj$X <- xen_obj_trimmed@meta.data[colnames(seu_obj), "X"]
  seu_obj$Y <- xen_obj_trimmed@meta.data[colnames(seu_obj), "Y"]
  
  seu_obj <- subset(
    seu_obj,
    subset =
      nFeature_RNA > nFeatRNA[1] &
      nFeature_RNA < nFeatRNA[2] &
      nCount_RNA > nCountRNA[1] &
      nCount_RNA < nCountRNA[2]
  )
  
  seu_obj <- NormalizeData(
    seu_obj,
    verbose = FALSE
  )
  
  output_dir <- file.path(
    dir_xenium_data,
    paste0("TX5K_", TX),
    "Objects"
  )
  
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  saveRDS(
    seu_obj,
    file = file.path(
      output_dir,
      paste0(
        "[SeuObj][TX5K_", TX, "]_Countable_Mag", Mag,
        "_nFeat", nFeatRNA[1], "-", nFeatRNA[2],
        "_nCount", nCountRNA[1], "-", nCountRNA[2],
        ".rds"
      )
    )
  )
  
  message(
    "TX5K_", TX,
    ": ", format(ncol(xen_obj_trimmed), big.mark = ","),
    " cells after trimming; ",
    format(ncol(seu_obj), big.mark = ","),
    " cells after QC"
  )
}

for (i in seq_len(nrow(sample_roi))) {
  create_qc_filtered_seurat(
    TX = sample_roi$TX[i],
    Xmin = sample_roi$Xmin[i],
    Xmax = sample_roi$Xmax[i],
    Ymin = sample_roi$Ymin[i],
    Ymax = sample_roi$Ymax[i],
    Mag = Mag,
    nFeatRNA = nFeatRNA,
    nCountRNA = nCountRNA
  )
}


## Read and merge QC-filtered Seurat objects --------------------------------

read_xenium_seurat <- function(TX, Mag, nFeatRNA, nCountRNA) {
  obj <- readRDS(
    file.path(
      dir_xenium_data,
      paste0("TX5K_", TX),
      "Objects",
      paste0(
        "[SeuObj][TX5K_", TX, "]_Countable_Mag", Mag,
        "_nFeat", nFeatRNA[1], "-", nFeatRNA[2],
        "_nCount", nCountRNA[1], "-", nCountRNA[2],
        ".rds"
      )
    )
  )
  
  colnames(obj) <- paste(
    str_remove(obj$orig.ident, pattern = "_"),
    colnames(obj),
    sep = "_"
  )
  
  return(obj)
}

ObjList <- lapply(
  TXnumInteg[seq_len(NumOfSamples)],
  function(TX) {
    read_xenium_seurat(
      TX = TX,
      Mag = Mag,
      nFeatRNA = nFeatRNA,
      nCountRNA = nCountRNA
    )
  }
)

CombinedObj <- Reduce(
  function(x, y) merge(x, y = y, add.cell.ids = NULL),
  ObjList
)

message("Merged cells: ", format(ncol(CombinedObj), big.mark = ","))


## Variable feature selection, scaling, and PCA ------------------------------

CombinedObj <- CombinedObj %>%
  FindVariableFeatures(
    selection.method = "vst",
    nfeatures = 2000
  )

CombinedObj <- ScaleData(
  CombinedObj,
  features = VariableFeatures(CombinedObj)
)

set.seed(123)

CombinedObj <- RunPCA(
  CombinedObj,
  features = VariableFeatures(CombinedObj),
  npcs = 50,
  verbose = FALSE
)


## Harmony integration -------------------------------------------------------

set.seed(123)

CombinedObj <- RunHarmony(
  CombinedObj,
  group.by.vars = "orig.ident",
  reduction.use = "pca",
  dims.use = 1:50,
  assay.use = "RNA"
)


## Initial clustering --------------------------------------------------------

set.seed(123)

CombinedObj <- CombinedObj %>%
  RunUMAP(
    reduction = "harmony",
    dims = 1:Dim1
  ) %>%
  FindNeighbors(
    reduction = "harmony",
    dims = 1:Dim1
  ) %>%
  FindClusters(
    resolution = Res1
  )


## Initial cell type annotation ----------------------------------------------

CombinedObj$InitialCellType <- Vec.Annot[as.character(CombinedObj$seurat_clusters)]

CombinedObj$InitialCellType <- factor(
  CombinedObj$InitialCellType,
  levels = names(List.Annot)
)

table(CombinedObj$InitialCellType, useNA = "ifany")


## Save annotated initial clustering object ----------------------------------

output_dir <- file.path(DirInteg, "Objects")

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  CombinedObj,
  file = file.path(
    output_dir,
    paste0(
      "[SeuObj][",
      NumOfSamples, "case(",
      paste(TXnumInteg[seq_len(NumOfSamples)], collapse = ","),
      ")]_",
      QCInfo.FileName,
      "_1stPCA50vf2000_ClustDim",
      Dim1,
      "Res",
      Res1,
      "_annotated.rds"
    )
  )
)