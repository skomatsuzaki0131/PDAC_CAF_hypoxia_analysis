library(dplyr)
library(tidyr)
library(tidyverse)
library(stringr)
library(magrittr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(cowplot)
library(patchwork)

library(edgeR)
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)  # ヒト用annotation
library(EnsDb.Hsapiens.v86)

clean_gs_label_ForFig <- function(x){
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
    str_replace_all("\\biCAF signature\\b", "iCAF")  %>%
    str_replace_all("\\biCAF\\b", "iCAF signature") %>%
    str_replace_all("\\bmyCAF signature\\b", "myCAF") %>%
    str_replace_all("\\bmyCAF\\b", "myCAF signature") %>%
    str_replace_all("\\bInterferon alpha\\b", "IFN-α") %>%
    str_replace_all("\\bInterferon gamma\\b", "IFN-γ") %>%
    str_replace_all("Epithelial mesenchymal transition", "EMT")
}
clean_gs_label_ForTable <- function(x){
  x %>%
    stringi::stri_trans_general("Any-ASCII") %>%  # ←追加①
    str_replace_all("[–—−]", "-") %>%             # ←追加②
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
    str_replace_all("\\biCAF signature\\b", "iCAF")  %>%
    str_replace_all("\\biCAF\\b", "iCAF signature") %>%
    str_replace_all("\\bmyCAF signature\\b", "myCAF") %>%
    str_replace_all("\\bmyCAF\\b", "myCAF signature") %>%
    str_replace_all("Epithelial mesenchymal transition", "EMT")
}

DF.Markers = read.csv(file="/Volumes/PortableSSD/[4]myR/[1]dataset/MyGeneSet.csv", header=T) %>% 
  dplyr::filter(!Gene %in% c("CD3D", "COL1A2"))

Directory    = '/Volumes/PortableSSD/[4]myR/[1]dataset/'
DirRNAseq    = "/Volumes/PortableSSD/[4]myR/[1]dataset/RNAseq/"
DirRNAseqNew = '/Volumes/PortableSSD/[4]myR/[1]dataset/RNAseq_new/'
DirKEGG      = '/Volumes/PortableSSD/[4]myR/[1]dataset/RNAseq/KEGG/'
DirOutput    = '/Volumes/PortableSSD/[3]Graduate school/[3]RNA-seq/analysis/KEGG_pathway/'
DirGO = "/Volumes/PortableSSD/[4]myR/[1]dataset/RNAseq/GO/"
DirOutputGO = "/Volumes/PortableSSD/[3]Graduate school/[3]RNA-seq/analysis/GeneOntonogy/"
DirHallmark = '/Volumes/PortableSSD/[4]myR/[1]dataset/RNAseq/Hallmark/'
DirOutputMSigDB = "/Volumes/PortableSSD/[3]Graduate school/[3]RNA-seq/analysis/MSigDB/"


#----------------------------------------------------
#
#   1. DEG analysis
#
#----------------------------------------------------

# Read total count data (Not normalized)
data_all = read.csv("[1]dataset/rnaseq_count.csv", row.names=1)

# Extract two group data of comparison settings
counts_HvsN <-  as.matrix(data_all[c(1:3,10:12,13:15,22:24,25:27,34:36)])   # 1.Simple comparisons
counts_HNvsN <- as.matrix(data_all[c(4:6,10:12,16:18,22:24,28:30,34:36)])  # 2.In 21% comparison: HN vs N
counts_HvsNH <- as.matrix(data_all[c(1:3,7:9,13:15,19:21,25:27,31:33)])    # 3.In 1% comparison: H vs NH
counts_HvsHN <- as.matrix(data_all[c(1:3,4:6,13:15,16:18,25:27,28:30)])    # 4.Before/After switching(Hypo-isolated): H vs HN
counts_NHvsN <- as.matrix(data_all[c(7:9,10:12,19:21,22:24,31:33,34:36)])  # 5.Before/After switching(Normo-isolated): NH vs N

# sample annotation
Vec.sample_names_HvsN <- colnames(counts_HvsN)
Vec.sample_names_HNvsN <- colnames(counts_HNvsN)
Vec.sample_names_HvsNH <- colnames(counts_HvsNH)
Vec.sample_names_HvsHN <- colnames(counts_HvsHN)
Vec.sample_names_NHvsN <- colnames(counts_NHvsN)

# sample information data
Func.sample_info <- function(sample_name){
  DF.sample_info <- data.frame(
    sample = sample_name,
    culture = c(rep("X0904",times=6), rep("X1115",times=6), rep("X1122",times=6)),
    oxygen = rep(c(rep("Hypoxia",times=3), rep("Normoxia",times=3)), times=3),
    stringsAsFactors = FALSE) %>%
    dplyr::mutate(
      culture = factor(culture),
      contrast_group = 
        oxygen %>% 
        str_replace(pattern="Hypoxia", replacement="PositiveSide") %>%
        str_replace(pattern="Normoxia", replacement="NegativeSide") %>% 
        factor(levels = c("NegativeSide", "PositiveSide"))
    )
  rownames(DF.sample_info) <- DF.sample_info$sample
  stopifnot(all(colnames(counts) == rownames(DF.sample_info)))
  DF.sample_info
}
DF.sample_info_HvsN <- Func.sample_info(Vec.sample_names_HvsN)
DF.sample_info_HNvsN <- Func.sample_info(Vec.sample_names_HNvsN)
DF.sample_info_HvsNH <- Func.sample_info(Vec.sample_names_HvsNH)
DF.sample_info_HvsHN <- Func.sample_info(Vec.sample_names_HvsHN)
DF.sample_info_NHvsN <- Func.sample_info(Vec.sample_names_NHvsN)

# DGEList
Func.DGEList <- function(counts, sample_info){
  DGEListObj <- DGEList(
    counts = counts,
    samples = sample_info
  )
  DGEListObj
}
DGEListObj_HvsN <- Func.DGEList(counts_HvsN, DF.sample_info_HvsN)
DGEListObj_HNvsN <- Func.DGEList(counts_HNvsN, DF.sample_info_HNvsN)
DGEListObj_HvsNH <- Func.DGEList(counts_HvsNH, DF.sample_info_HvsNH)
DGEListObj_HvsHN <- Func.DGEList(counts_HvsHN, DF.sample_info_HvsHN)
DGEListObj_NHvsN <- Func.DGEList(counts_NHvsN, DF.sample_info_NHvsN)

# design matrix; culture + oxygen
Func.design.matrix <- function(sample_info){
  Mt.design <- model.matrix(
    ~ culture + contrast_group,
    data = sample_info
  )
  Mt.design
}

MT.design_HvsN <- Func.design.matrix(DF.sample_info_HvsN)
MT.design_HNvsN <- Func.design.matrix(DF.sample_info_HNvsN)
MT.design_HvsNH <- Func.design.matrix(DF.sample_info_HvsNH)
MT.design_HvsHN <- Func.design.matrix(DF.sample_info_HvsHN)
MT.design_NHvsN <- Func.design.matrix(DF.sample_info_NHvsN)

colnames(MT.design_HvsN)
colnames(MT.design_HNvsN)
colnames(MT.design_HvsNH)
colnames(MT.design_HvsHN)
colnames(MT.design_NHvsN)

# filtering
Func.FilterByExpr <- function(DGEList, design.matrix){
  keep <- filterByExpr(
    DGEList,
    design = design.matrix
  )
  DGEList_keep <- DGEList[keep, , keep.lib.sizes = FALSE]
  
  DGEList_keep
}

DGEListObj_keep_HvsN <- Func.FilterByExpr(DGEListObj_HvsN, MT.design_HvsN)
DGEListObj_keep_HNvsN <- Func.FilterByExpr(DGEListObj_HNvsN, MT.design_HNvsN)
DGEListObj_keep_HvsNH <- Func.FilterByExpr(DGEListObj_HvsNH, MT.design_HvsNH)
DGEListObj_keep_HvsHN <- Func.FilterByExpr(DGEListObj_HvsHN, MT.design_HvsHN)
DGEListObj_keep_NHvsN <- Func.FilterByExpr(DGEListObj_NHvsN, MT.design_NHvsN)

# DEG table
Func.DEGtable <- function(DGEListObj_keep, design.matrix){
  # normalization
  DGEListObj_keep_norm <- 
    calcNormFactors(
      DGEListObj_keep,
      method = "TMM"
    )
  # dispersion
  DGEListObj_keep_norm_dispersion <- 
    estimateDisp(
      DGEListObj_keep_norm,
      design.matrix,
      robust = TRUE
    )
  # GLM fitting
  fit <- 
    glmQLFit(
      DGEListObj_keep_norm_dispersion,
      design.matrix,
      robust = TRUE
    )
  # DEG test
  qlf <- 
    glmQLFTest(
      fit,
      coef = "contrast_groupPositiveSide"
    )
  # DEG table
  res <- topTags(
    qlf,
    n = Inf
  )$table
  res_arranged <- 
    dplyr::arrange(res, desc(logFC))
  # out put
  return(res_arranged)
}
DEG.table_HvsN <- Func.DEGtable(DGEListObj_keep_HvsN, MT.design_HvsN)
DEG.table_HNvsN <- Func.DEGtable(DGEListObj_keep_HNvsN, MT.design_HNvsN)
DEG.table_HvsNH <- Func.DEGtable(DGEListObj_keep_HvsNH, MT.design_HvsNH)
DEG.table_HvsHN <- Func.DEGtable(DGEListObj_keep_HvsHN, MT.design_HvsHN)
DEG.table_NHvsN <- Func.DEGtable(DGEListObj_keep_NHvsN, MT.design_NHvsN)

# save DEG table
write.csv(DEG.table_HvsN,  file=paste0(DirRNAseqNew,"[DEGsTable]_Comp1_HvsN.csv"), row.names=T)
write.csv(DEG.table_HNvsN, file=paste0(DirRNAseqNew,"[DEGsTable]_Comp2_HNvsN.csv"), row.names=T)
write.csv(DEG.table_HvsNH, file=paste0(DirRNAseqNew,"[DEGsTable]_Comp3_HvsNH.csv"), row.names=T)
write.csv(DEG.table_HvsHN, file=paste0(DirRNAseqNew,"[DEGsTable]_Comp4_HvsHN.csv"), row.names=T)
write.csv(DEG.table_NHvsN, file=paste0(DirRNAseqNew,"[DEGsTable]_Comp5_NHvsN.csv"), row.names=T)

Func.AddEntrez = function(Compname){
  # Read DEG table
  DF.DEG_0 <- read.table(file=paste0(DirRNAseqNew,"[DEGsTable]_",Compname,".csv"), header=T, sep=",", stringsAsFactors=F, row.names=1)
  Vec.Symbol <- row.names(DF.DEG_0)
  # Add Entrez id
  DF.Sym_and_En <- 
    clusterProfiler::bitr(
      Vec.Symbol, 
      fromType="SYMBOL", 
      toType="ENTREZID", 
      OrgDb=org.Hs.eg.db, 
      drop=FALSE) #変換できなかった遺伝子を"NA"として残す.残さないとDEGのdata frameとズレる.  
  DF.DuplicatedSym <- 
    DF.Sym_and_En %>% 
    rownames_to_column(var="rownumber") %>% 
    group_by(SYMBOL) %>% 
    dplyr::filter(n()>1) %>% 
    cbind(Duplication = "Duplicated")
  if(nrow(DF.DuplicatedSym)==0){
    DF.DEG_1 = cbind(DF.DEG_0, 
                     EntrezID = DF.Sym_and_En$ENTREZID) #DEGのdata.frameに結合
    write.csv(DF.DEG_1, 
              file=paste0(DirRNAseqNew,"[DEGsTable_addENTREZ]_",Compname,".csv"), row.names=TRUE)
  }else{
    print(DF.DuplicatedSym)}
  if(nrow(DF.DuplicatedSym)==2 & 
     DF.DuplicatedSym$ENTREZID[1]=="7006" &
     DF.DuplicatedSym$ENTREZID[2]=="100124696"){
    DF.Sym_and_En_rm = DF.Sym_and_En[-as.numeric(DF.DuplicatedSym[2, "rownumber"]), ]
    DF.DEG_1 = cbind(DF.DEG_0, 
                     EntrezID = DF.Sym_and_En_rm$ENTREZID) #DEGのdata.frameに結合
    write.csv(DF.DEG_1, 
              file=paste0(DirRNAseqNew,"[DEGsTable_addENTREZ]_",Compname,".csv"), row.names=TRUE)
    print("Removed one of TEC (ENTEREZ 100124696)")
  }
}
Func.AddEntrez(Compname = "Comp1_HvsN") # TECのENTEREZ 100124696はchildhoodのものなので削除
Func.AddEntrez(Compname = "Comp2_HNvsN") # TECのENTEREZ 100124696はchildhoodのものなので削除
Func.AddEntrez(Compname = "Comp3_HvsNH") # TECのENTEREZ 100124696はchildhoodのものなので削除
Func.AddEntrez(Compname = "Comp4_HvsHN") # TECのENTEREZ 100124696はchildhoodのものなので削除
Func.AddEntrez(Compname = "Comp5_NHvsN") # TECのENTEREZ 100124696はchildhoodのものなので削除


#----------------------------------------------------
#
#   2. Volcano plot
#
#----------------------------------------------------

DEG1_HvsN_0 <- read.csv(file=paste0(DirRNAseqNew,"[DEGsTable]_Comp1_HvsN.csv"), header=T, row.names=NULL) %>% dplyr::rename(Gene=X)
DEG2_HNvsN_0 <- read.csv(file=paste0(DirRNAseqNew,"[DEGsTable]_Comp2_HNvsN.csv"), header=T, row.names=NULL) %>% dplyr::rename(Gene=X)
DEG3_HvsNH_0 <- read.csv(file=paste0(DirRNAseqNew,"[DEGsTable]_Comp3_HvsNH.csv"), header=T, row.names=NULL) %>% dplyr::rename(Gene=X)
DF.Binded_0 <- bind_rows(
  data.frame("Setting"="Original", DEG1_HvsN_0),
  data.frame("Setting"="Under21", DEG2_HNvsN_0),
  data.frame("Setting"="Under1", DEG3_HvsNH_0)) %>% 
  dplyr::select(Setting, Gene, logFC, FDR) %>% 
  dplyr::mutate(minuslog10FDR = (-1)*log10(FDR),
                ptcol = case_when(logFC>1 & FDR<0.05 ~ "BLUE",
                                  logFC<(-1) & FDR<0.05 ~ "RED",
                                  TRUE ~ "GRAY")) 
MinFDR = min(DF.Binded_0$FDR)
MaxAbs.logFC = max(abs(DF.Binded_0$logFC))
Genes_Annotate = c(
  "TOP2A","TPX2","FOXM1","CENPA","CDK2","CCNA2","BIRC5",
  "CDKN2A", "CDKN1A")
DF.Binded_Lab = subset(DF.Binded_0, Gene%in%Genes_Annotate)

Func.VolPlotForFig = function(SettingName, RightTextCol){
  DF.VolcPlot = subset(DF.Binded_0, Setting==SettingName) %>% 
    dplyr::mutate(ptcol = factor(ptcol, levels=c("GRAY","RED","BLUE"))) %>% 
    dplyr::arrange(ptcol)
  DF.VolcPlotLab = subset(DF.Binded_Lab, Setting==SettingName)
  NudgeYleft=0.6
  NudgeYright=2.0
  NudgeXright=2.5
  NudgeXleft=6
  Force=8
  BoxPadding=0.6
  PointPadding=0.3
  ggplot(DF.VolcPlot, aes(x=logFC, y=minuslog10FDR)) +
    geom_point(aes(fill=ptcol, color=ptcol), 
               shape=21, size=2, stroke=0.1) +
    geom_hline(yintercept=(-1)*log10(0.05), color="black", linetype="dashed", linewidth=0.2) +
    geom_vline(xintercept=c(-1, 1), color="black", linetype="dashed", linewidth=0.2) +
    # right side
    ggrepel::geom_text_repel(
      data=subset(DF.VolcPlotLab, logFC>0), aes(label=Gene),
      color=RightTextCol, segment.color=RightTextCol,
      nudge_y=NudgeYright, nudge_x=NudgeXright, force=Force,
      force_pull=0, direction="both", xlim=c(3.5, Inf), hjust=0,
      box.padding=BoxPadding, point.padding=PointPadding,
      fontface="bold", max.overlaps=100, seed=1) +
    # left side
    ggrepel::geom_text_repel(
      data=subset(DF.VolcPlotLab, logFC<0), aes(label=Gene),
      color=RightTextCol, segment.color=RightTextCol,
      nudge_y=NudgeYleft, nudge_x=-NudgeXleft, force=Force,
      force_pull=0, direction="y", xlim=c(-Inf, -1.2),
      box.padding=BoxPadding, point.padding=PointPadding,
      fontface="bold", max.overlaps=100, seed=1) +
    labs(x="log2 fold change", y=paste0(SettingName,"\n\n-log10(FDR)")) +
    scale_x_continuous(breaks=c(-5, -1, 1, 5)) +
    scale_y_continuous(expand=expansion(mult=c(0, 0.03)),
                       limits=c(0, (-1.1)*log10(MinFDR)),
                       breaks=c(0,5,10)) +
    scale_color_manual(values=c("BLUE"="black", "RED"="black", "GRAY"="gray")) +
    scale_fill_manual(values=c("BLUE"="#2c7fb8", "RED"="#f03b20", "GRAY"="gray")) +
    coord_cartesian(xlim=c(-MaxAbs.logFC, MaxAbs.logFC),
                    clip="off") +
    theme(plot.margin = margin(4, 22, 4, 22),
          legend.position = "none",
          axis.text.x = element_text(face="bold", color="black"),
          axis.text.y = element_text(face="bold", color="black"),
          axis.title = element_text(face="bold", color="black"),
          plot.background = element_rect(fill="transparent", color=NA),
          panel.background = element_rect(fill="white", color=NA),
          panel.grid = element_blank(),
          panel.border = element_rect(fill=NA, color="black"))
}
VolcForPaper=
  (Func.VolPlotForFig("Original","black")|Func.VolPlotForFig("Original",NA))/
  (Func.VolPlotForFig("Under21","black")|Func.VolPlotForFig("Under21",NA))/
  (Func.VolPlotForFig("Under1","black")|Func.VolPlotForFig("Under1",NA)) &
  theme(
    plot.background = element_rect(fill="transparent", color=NA),
    panel.background = element_rect(fill="white", color=NA))
ggsave(plot=VolcForPaper, file="/Volumes/PortableSSD/[4]myR/VolcNew.png", 
       width=9, 
       height=5, 
       dpi=500, 
       bg="transparent")
write.csv(DF.Binded_Lab, file="Volc_AnnotatedGenesNew.csv", row.names=F)
saveRDS(Func.VolPlotForFig("Original",NA), file = "Fig2B_Volc_Upper_New.rds")
saveRDS(Func.VolPlotForFig("Under21",NA), file = "Fig2B_Volc_Mid_New.rds")
saveRDS(Func.VolPlotForFig("Under1",NA), file = "Fig2B_Volc_Lower_New.rds")

#----------------------------------------------------
#
#   3. KEGG
#
#----------------------------------------------------

library(clusterProfiler)
DEG1_HvsN_0 <- read.csv(file=paste0(DirRNAseqNew,"[DEGsTable_addENTREZ]_Comp1_HvsN.csv"), header=T, row.names=NULL) %>% dplyr::rename(Gene=X)
DEG2_HNvsN_0 <- read.csv(file=paste0(DirRNAseqNew,"[DEGsTable_addENTREZ]_Comp2_HNvsN.csv"), header=T, row.names=NULL) %>% dplyr::rename(Gene=X)
DEG3_HvsNH_0 <- read.csv(file=paste0(DirRNAseqNew,"[DEGsTable_addENTREZ]_Comp3_HvsNH.csv"), header=T, row.names=NULL) %>% dplyr::rename(Gene=X)

Func.EnrichKEGG = function(DEGall, CompName){
  DEGsig = 
    DEGall %>% 
    dplyr::filter((logFC>=1 | logFC<=(-1)) & FDR<0.05)
  DEGsig.UpInHCAF = dplyr::filter(DEGsig, logFC > 0)
  DEGsig.DownInHCAF = dplyr::filter(DEGsig, logFC < 0)
  gene_use.UpInHCAF = DEGsig.UpInHCAF$Entrez %>% unique() %>% na.omit()
  gene_use.DownInHCAF = DEGsig.DownInHCAF$Entrez %>% unique() %>% na.omit()
  bg = DEGall$Entrez %>% as.character() %>% unique() %>% na.omit()
  set.seed(1234)
  ResKEGG.UpInHCAF = 
    clusterProfiler::enrichKEGG(
      gene = gene_use.UpInHCAF, 
      universe = bg,
      organism="hsa", 
      keyType = "ncbi-geneid",
      pAdjustMethod = "BH")
  ResKEGG.DownInHCAF = 
    clusterProfiler::enrichKEGG(
      gene=gene_use.DownInHCAF, 
      universe = bg,
      organism="hsa", 
      keyType = "ncbi-geneid",
      pAdjustMethod = "BH")
  DF.Res = 
    rbind(
      dplyr::mutate(ResKEGG.UpInHCAF@result, Setting=CompName, Direction="UP"),
      dplyr::mutate(ResKEGG.DownInHCAF@result, Setting=CompName, Direction="DOWN"))
  return(DF.Res)
}
DF.ResKEGG = 
  bind_rows(
    Func.EnrichKEGG(DEG1_HvsN_0, "1Simple_HvsN"),
    Func.EnrichKEGG(DEG2_HNvsN_0, "2in21pct_HNvsN"),
    Func.EnrichKEGG(DEG3_HvsNH_0, "3in1pct_HvsNH")) %>% 
  dplyr::mutate(
    GeneRatio_num = 
      sapply(GeneRatio, function(x){
        a = as.numeric(strsplit(x, "/")[[1]][1])
        b = as.numeric(strsplit(x, "/")[[1]][2])
        a / b
      }),
    GeneRatio_num_signed = 
      GeneRatio_num*ifelse(Direction=="UP", 1, -1),
    SigLab = case_when(p.adjust<0.001 ~ "***",
                       p.adjust<0.01 ~ "**",
                       p.adjust<0.05 ~ "*",
                       TRUE ~ "n.s."))
DF.ResKEGG_signif <- 
  DF.ResKEGG %>% 
  dplyr::filter(SigLab != "n.s.") %>% 
  dplyr::select(Description, Setting, Direction) %>% 
  pivot_wider(id_cols=Description, 
              names_from=Setting, 
              values_from = Direction)
Vec.DisplayPathway.new <- 
  c("Cell cycle",
    "Motor proteins",
    "Homologous recombination",
    "DNA replication",
    "HIF-1 signaling pathway",
    "JAK-STAT signaling pathway",
    "Hippo signaling pathway",
    "Cadherin signaling",
    "Rap1 signaling pathway",
    "Calcium signaling pathway"
    )
DF.ResKEGG.new = 
  DF.ResKEGG %>% 
  dplyr::filter(p.adjust < 0.05 & Description %in% Vec.DisplayPathway.new) %>% 
  dplyr::mutate(
    Setting = factor(Setting, levels=c("1Simple_HvsN","2in21pct_HNvsN","3in1pct_HvsNH")),
    Description = factor(Description, levels=Vec.DisplayPathway.new))
MaxAbsGR = max(DF.ResKEGG.new$GeneRatio_num)
TextSize=10
PlotKEGGtile = 
  ggplot(DF.ResKEGG.new, 
         aes(x=Setting, y=Description, fill=GeneRatio_num_signed)) +
  geom_tile(color="gray20") +
  geom_text(aes(label=SigLab),
            fontface="bold", color="black", 
            vjust=0.75,
            size=TextSize*0.45) +
  labs(x=NULL, y=NULL,
       fill="Signed\ngene ratio") +
  scale_x_discrete(
    expand=expansion(mult=c(0,0)),
    labels=c("1Simple_HvsN"="Original","3in21pct_HNvsN"="Under21","2in1pct_HvsNH"="Under1")) +
  scale_y_discrete(expand=expansion(mult=c(0,0))) +
  scale_fill_distiller(
    palette="RdBu", 
    direction=1,
    limits=c(-MaxAbsGR,MaxAbsGR),
    breaks=c(-0.1, 0, 0.1),
    guide=guide_colorbar(
      frame.colour="black",
      ticks.colour="black")) +
  theme(
    axis.text.x = element_text(face="bold", color="black", angle=45, hjust=1, size=TextSize*0.5),
    axis.text.y = element_text(face="bold", color="black", size=TextSize*1.3),
    legend.title = element_text(face="bold", color="black", size=TextSize+2, hjust = 0.5),
    legend.text = element_text(face="bold", color="black", size=TextSize),
    plot.background = element_rect(fill="transparent", color=NA),
    legend.background = element_rect(fill="transparent", color=NA),
    panel.background = element_rect(fill="white", color=NA),
    panel.border = element_rect(fill="transparent", color="black"),
    panel.grid = element_blank()
  )
ggsave(plot=PlotKEGGtile,
       file="KEGG_tile_new.png",
       width=5.8, height=3, dpi=500, bg="transparent")
saveRDS(PlotKEGGtile, file="Fig2C_new.rds")

## Pathview 
# DEGs with significance
DEGsig.1_HvsN = dplyr::filter(DEG1_HvsN_0, (logFC < (-1) | logFC > 1) & FDR < 0.05)
DEGsig.2_HNvsN = dplyr::filter(DEG2_HNvsN_0, (logFC < (-1) | logFC > 1) & FDR < 0.05)
DEGsig.3_HvsNH = dplyr::filter(DEG3_HvsNH_0, (logFC < (-1) | logFC > 1) & FDR < 0.05)
# log2FC vector named with EntrezID
Vec.log2FC.sig_1_HvsN = setNames(DEGsig.1_HvsN$logFC, DEGsig.1_HvsN$EntrezID)
Vec.log2FC.sig_2_HNvsN = setNames(DEGsig.2_HNvsN$logFC, DEGsig.2_HNvsN$EntrezID)
Vec.log2FC.sig_3_HvsNH = setNames(DEGsig.3_HvsNH$logFC, DEGsig.3_HvsNH$EntrezID)

Vec.log2FC.sig.rmNA_1_HvsN = Vec.log2FC.sig_1_HvsN[! is.na(names(Vec.log2FC.sig_1_HvsN))]
Vec.log2FC.sig.rmNA_2_HNvsN = Vec.log2FC.sig_2_HNvsN[! is.na(names(Vec.log2FC.sig_2_HNvsN))]
Vec.log2FC.sig.rmNA_3_HvsNH = Vec.log2FC.sig_3_HvsNH[! is.na(names(Vec.log2FC.sig_3_HvsNH))]
# DEG1~3 および DEG1~5をまとめて表示する
library(KEGGREST)
pathway <- keggGet("hsa04110")[[1]]
genes <- pathway$GENE
Vec.EntrezID.CCgenes <- genes[seq(1, length(genes), 2)]
DF.CellCycleGenes = read.csv("[1]dataset/EntrezID_of_genes_belonging_CellCyclePathway.csv", header=TRUE)
# CellCyclePathwayに属する遺伝子のlogFCをmatrixに抽出
v1 <- Vec.log2FC.sig.rmNA_1_HvsN[Vec.EntrezID.CCgenes]
v2 <- Vec.log2FC.sig.rmNA_2_HNvsN[Vec.EntrezID.CCgenes]
v3 <- Vec.log2FC.sig.rmNA_3_HvsNH[Vec.EntrezID.CCgenes]
MT.three = cbind(v1, v2, v3)
colnames(MT.three) = c("Original", "Under21", "Under1")
rownames(MT.three) = Vec.EntrezID.CCgenes
# generate pathview map Cell cycle
library(pathview)
pathview(gene.data = MT.three, 
         pathway.id = "hsa04110", species = "hsa",
         limit = list(gene=3, cpd=1),
         low = list(gene = "red", cpd = "black"),
         mid = list(gene = "white", cpd = "black"), 
         high = list(gene = "blue", cpd = "black"),
         same.layer = FALSE,
         kegg.native = TRUE,
         plot.col.key = F,
         key.pos = "bottomright")
file.rename(from="hsa04110.pathview.multi.png",
            to=paste0("Pathview_CellCycle_new.png"))
# Supp table
DF.ResKEGG_SuppTable <- 
  DF.ResKEGG %>% 
  dplyr::filter(p.adjust < 0.05) %>% 
  group_by(Setting) %>% 
  arrange(desc(Direction), desc(GeneRatio_num), .by_group = T) %>% 
  ungroup() %>% 
  dplyr::select(Setting, Direction ,ID, Description, GeneRatio, BgRatio, p.adjust, geneID) %>% 
  set_colnames(c("Oxygen condition",
                 "Enrichment in Hypo-CAF",
                 "ID",
                 "Description",
                 "Gene ratio",
                 "Background ratio",
                 "p.adjust",
                 "gene ID"))
write.csv(DF.ResKEGG_SuppTable, "SuppTable_KEGG_new.csv", row.names = F)

#----------------------------------------------------
#
#   4. GO analysis (GSEA)
#
#----------------------------------------------------

DEG1_HvsN_0 <- read.csv(file=paste0(DirRNAseqNew,"[DEGsTable_addENTREZ]_Comp1_HvsN.csv"), header=T, row.names=NULL) %>% dplyr::rename(Gene=X)
DEG2_HNvsN_0 <- read.csv(file=paste0(DirRNAseqNew,"[DEGsTable_addENTREZ]_Comp2_HNvsN.csv"), header=T, row.names=NULL) %>% dplyr::rename(Gene=X)
DEG3_HvsNH_0 <- read.csv(file=paste0(DirRNAseqNew,"[DEGsTable_addENTREZ]_Comp3_HvsNH.csv"), header=T, row.names=NULL) %>% dplyr::rename(Gene=X)

Vec.namedFC_DEG1_HvsN <- setNames(DEG1_HvsN_0$logFC, DEG1_HvsN_0$EntrezID)
Vec.namedFC_DEG2_HNvsN <- setNames(DEG2_HNvsN_0$logFC, DEG2_HNvsN_0$EntrezID)
Vec.namedFC_DEG3_HvsNH <- setNames(DEG3_HvsNH_0$logFC, DEG3_HvsNH_0$EntrezID)

Vec.namedFC_DEG1_HvsN_rmNA <- Vec.namedFC_DEG1_HvsN[! is.na(names(Vec.namedFC_DEG1_HvsN))]
Vec.namedFC_DEG2_HNvsN_rmNA <- Vec.namedFC_DEG2_HNvsN[! is.na(names(Vec.namedFC_DEG2_HNvsN))]
Vec.namedFC_DEG3_HvsNH_rmNA <- Vec.namedFC_DEG3_HvsNH[! is.na(names(Vec.namedFC_DEG3_HvsNH))]

Func.gseGO = function(Vec.FC){
  set.seed(1234)
  clusterProfiler::gseGO(
    geneList=Vec.FC,
    OrgDb=org.Hs.eg.db,
    ont="ALL",   #"BP","CC","MF","ALL"から選択
    minGSSize= 50,
    pvalueCutoff = 0.05,
    verbose=FALSE,
    #nPermSimple = 10000,
    eps=0)
}

Res.gseGO_1_HvsN = Func.gseGO(Vec.namedFC_DEG1_HvsN_rmNA)
Res.gseGO_2_HNvsN = Func.gseGO(Vec.namedFC_DEG2_HNvsN_rmNA)
Res.gseGO_3_HvsNH = Func.gseGO(Vec.namedFC_DEG3_HvsNH_rmNA)

Res.gseGO.simple_1_HvsN = clusterProfiler::simplify(Res.gseGO_1_HvsN)
Res.gseGO.simple_2_HNvsN = clusterProfiler::simplify(Res.gseGO_2_HNvsN)
Res.gseGO.simple_3_HvsNH = clusterProfiler::simplify(Res.gseGO_3_HvsNH)

Func.DataWrangling = function(Obj.gseaRes, SettingName){
  DF_0 <- Obj.gseaRes@result
  DF_1 <- DF_0 %>% 
    dplyr::mutate(Direction=
                    ifelse(NES > 0, "UP", "DOWN") %>% 
                    factor(levels = c("UP","DOWN"))
                  ) %>% 
    group_by(Direction) %>% 
    dplyr::arrange(p.adjust, .by_group=TRUE)
  DF_2 <- 
    data.frame(setting = SettingName,
               DF_1)
}
DF.gseaRes.1_HvsN_0 = Func.DataWrangling(Res.gseGO.simple_1_HvsN, "isolation")
DF.gseaRes.2_HNvsN_0 = Func.DataWrangling(Res.gseGO.simple_2_HNvsN, "under21")
DF.gseaRes.3_HvsNH_0 = Func.DataWrangling(Res.gseGO.simple_3_HvsNH, "under1")

DF.gseaRes.bind_0 = 
  bind_rows(DF.gseaRes.1_HvsN_0,
            DF.gseaRes.2_HNvsN_0,
            DF.gseaRes.3_HvsNH_0) %>% 
  dplyr::mutate(minuslog10adjp = (-1)*log10(p.adjust))

Vec.TermsConstantUp <- 
  intersect(subset(DF.gseaRes.1_HvsN_0, ONTOLOGY=="BP" & NES>0 & p.adjust<0.05)$Description, 
            subset(DF.gseaRes.2_HNvsN_0, ONTOLOGY=="BP" & NES>0 & p.adjust<0.05)$Description) %>% 
  intersect(subset(DF.gseaRes.3_HvsNH_0, ONTOLOGY=="BP" & NES>0 & p.adjust<0.05)$Description)
Vec.TermsConstantDown <- 
  intersect(subset(DF.gseaRes.1_HvsN_0, ONTOLOGY=="BP" & NES<0 & p.adjust<0.05)$Description, 
            subset(DF.gseaRes.2_HNvsN_0, ONTOLOGY=="BP" & NES<0 & p.adjust<0.05)$Description) %>% 
  intersect(subset(DF.gseaRes.3_HvsNH_0, ONTOLOGY=="BP" & NES<0 & p.adjust<0.05)$Description)
Vec.TermsToDisplay_0 = 
  c(Vec.TermsConstantDown,
    "SPACE",
    Vec.TermsConstantUp)
Vec.TermsToDisplay_Final = c(
  "chromosome segregation",
  "nuclear division",
  "regulation of cytokinesis",
  "response to type II interferon",
  "inflammatory response",
  "cytokine-mediated signaling pathway",
  "lipid catabolic process")
DF.gseaRes.bind_1 = 
  DF.gseaRes.bind_0 %>% 
  dplyr::filter(Description %in% Vec.TermsToDisplay_Final) %>% 
  dplyr::mutate(Description = factor(Description, levels=Vec.TermsToDisplay_Final))
Lim.NES = max(abs(DF.gseaRes.bind_1$NES))
Plot.DotGO = 
  ggplot(DF.gseaRes.bind_1, aes(x=setting, y=Description)) +
  geom_point(shape=21, aes(fill=NES, size=minuslog10adjp)) +
  labs(x=NULL, y=NULL, size="-log10(Adj.P)") +
  scale_y_discrete(position = "right") +
  scale_fill_distiller(
    palette="RdBu",
    direction=1,
    limits=c(-Lim.NES, Lim.NES),
    breaks=c(2, 0, -2),
    guide=guide_colorbar(
      order=2,
      frame.colour="black", 
      ticks.colour="black",
      label.position="left")) +
  scale_size_continuous(
    range=c(1,10),
    breaks=c(2,10,30),
    guide=guide_legend(
      order=1,
      label.position="left")) +
  theme(axis.text.y = element_text(face="bold", color="black", size=14),
        axis.text.x = element_text(angle=45, hjust=1, face="bold", color="black", size=14),
        legend.key = element_blank(),
        legend.title = element_text(face="bold", color="black", hjust=0.5),
        legend.text = element_text(face="bold", color="black"),
        legend.box = "verticacl",
        legend.box.just = "right",
        legend.justification = "right",
        legend.background = element_rect(fill="transparent", color=NA),
        plot.background = element_rect(fill="transparent", color=NA),
        panel.background = element_rect(fill="white", color="black"),
        panel.grid = element_blank(),
        panel.border = element_rect(fill=NA, color="black"))
ggsave(plot=Plot.DotGO,
       file="GO_new.png",
       width=6.5, height=3.5, dpi=500, bg="transparent")
saveRDS(Plot.DotGO + 
          labs(size="-log10\n(Adj.P)") +
          scale_size_continuous(range=c(1, 8.5),
                                breaks=c(2,10,30),
                                guide=guide_legend(
                                  order=1,
                                  label.position="left")), 
        "Fig2D_new.rds")
DF.gseaRes.bind_2 = DF.gseaRes.bind_0 %>% 
  #dplyr::select(c(-minuslog10adjp, -rank, -enrichmentScore)) %>% 
  dplyr::select(c(setting, Direction, ONTOLOGY, ID, Description, setSize, NES, p.adjust, pvalue, qvalue,
                  leading_edge, core_enrichment)) %>%
  dplyr::rename(Enrichment_in_HypoCAF = Direction)
write.csv(DF.gseaRes.bind_2,
          file="SupplementalTable_bulk_GSEA(GO)_new.csv",
          row.names=F)

#----------------------------------------------------
#
#   5. Human hallmarks GSEA (Hypo-CAF vs Normo-CAF)
#
#----------------------------------------------------
#library(msigdbr)
library(DOSE)
library(enrichplot)
library(ggdendro)
H = read.gmt("/Volumes/PortableSSD/[4]myR/[1]dataset/h.all.v2024.1.Hs.entrez.gmt")
#CategoryC2 = msigdbr(species="Homo sapiens",  category="C2")
#Additional_0 = 
#  data.frame(term = "Buffa_Hypoxia", symbol = Vec.BuffaOrig) %>% 
#  rbind(data.frame(term="Winter_Hypoxia", symbol=Vec.WinterOrig)) %>% 
#  rbind(data.frame(term="KEGG_Cell_cycle", symbol=CategoryC2 %>% filter(gs_name == "KEGG_CELL_CYCLE") %>% pull(gene_symbol) %>% unique())) %>% 
#  rbind(data.frame(term="myCAF_signature", symbol=DF.Markers$Gene[DF.Markers$Classification1=="myCAF"])) %>% 
#  rbind(data.frame(term="iCAF_signature", symbol=DF.Markers$Gene[DF.Markers$Classification1=="iCAF"]))
#Vec.AdditionalEntrez_0 = clusterProfiler::bitr(
#  Additional_0$symbol, fromType="SYMBOL", toType="ENTREZID", 
#  OrgDb=org.Hs.eg.db, drop=FALSE)
#Vec.AdditionalEntrez_1 = 
#  Vec.AdditionalEntrez_0$ENTREZID %>% 
#  setNames(Vec.AdditionalEntrez_0$SYMBOL)
#Additional_1 = Additional_0 %>% 
#  dplyr::mutate(gene = Vec.AdditionalEntrez_1[symbol])
#Additional_2 = dplyr::select(Additional_1, c(term, gene))
#Elyada_0 = read.csv(file='/Volumes/PortableSSD/[4]myR/[1]dataset/ElyadaSupTable_S22_Orthologs.csv')

# ranked gene set
DF.rankedDEG_1_HvsN <- read.csv(file=paste0(DirRNAseqNew,"[DEGsTable_addENTREZ]_Comp1_HvsN.csv"), header=T, row.names=NULL) %>% dplyr::rename(Gene=X)
DF.rankedDEG_2_HNvsN <- read.csv(file=paste0(DirRNAseqNew,"[DEGsTable_addENTREZ]_Comp2_HNvsN.csv"), header=T, row.names=NULL) %>% dplyr::rename(Gene=X)
DF.rankedDEG_3_HvsNH <- read.csv(file=paste0(DirRNAseqNew,"[DEGsTable_addENTREZ]_Comp3_HvsNH.csv"), header=T, row.names=NULL) %>% dplyr::rename(Gene=X)
DF.rankedDEG_4_HvsHN <- read.csv(file=paste0(DirRNAseqNew,"[DEGsTable_addENTREZ]_Comp4_HvsHN.csv"), header=T, row.names=NULL) %>% dplyr::rename(Gene=X)
DF.rankedDEG_5_NHvsN <- read.csv(file=paste0(DirRNAseqNew,"[DEGsTable_addENTREZ]_Comp5_NHvsN.csv"), header=T, row.names=NULL) %>% dplyr::rename(Gene=X)

Vec.namedFC_1_HvsN <- setNames(DF.rankedDEG_1_HvsN$logFC, DF.rankedDEG_1_HvsN$EntrezID)
Vec.namedFC_2_HNvsN <- setNames(DF.rankedDEG_2_HNvsN$logFC, DF.rankedDEG_2_HNvsN$EntrezID)
Vec.namedFC_3_HvsNH <- setNames(DF.rankedDEG_3_HvsNH$logFC, DF.rankedDEG_3_HvsNH$EntrezID)
Vec.namedFC_4_HvsHN <- setNames(DF.rankedDEG_4_HvsHN$logFC, DF.rankedDEG_4_HvsHN$EntrezID)
Vec.namedFC_5_NHvsN <- setNames(DF.rankedDEG_5_NHvsN$logFC, DF.rankedDEG_5_NHvsN$EntrezID)

Vec.namedFC_1_HvsN_rmNA <- Vec.namedFC_1_HvsN[! is.na(names(Vec.namedFC_1_HvsN))]
Vec.namedFC_2_HNvsN_rmNA <- Vec.namedFC_2_HNvsN[! is.na(names(Vec.namedFC_2_HNvsN))]
Vec.namedFC_3_HvsNH_rmNA <- Vec.namedFC_3_HvsNH[! is.na(names(Vec.namedFC_3_HvsNH))]
Vec.namedFC_4_HvsHN_rmNA <- Vec.namedFC_4_HvsHN[! is.na(names(Vec.namedFC_4_HvsHN))]
Vec.namedFC_5_NHvsN_rmNA <- Vec.namedFC_5_NHvsN[! is.na(names(Vec.namedFC_5_NHvsN))]

# Run GSEA
Func.GSEA = function(GeneList){
  set.seed(1234)
  gseaRes <- 
    GSEA(geneList=GeneList,
       minGSSize=10,
       TERM2GENE=rbind(H),
       pvalueCutoff=1,
       verbose=FALSE, 
       eps=0)
  return(gseaRes)
}
gseaRes_1_HvsN <- Func.GSEA(Vec.namedFC_1_HvsN_rmNA)
gseaRes_2_HNvsN <- Func.GSEA(Vec.namedFC_2_HNvsN_rmNA)
gseaRes_3_HvsNH <- Func.GSEA(Vec.namedFC_3_HvsNH_rmNA)
gseaRes_4_HvsHN <- Func.GSEA(Vec.namedFC_4_HvsHN_rmNA)
gseaRes_5_NHvsN <- Func.GSEA(Vec.namedFC_5_NHvsN_rmNA)
List.gseaRes <- list(
  gseaRes_1_HvsN,
  gseaRes_2_HNvsN,
  gseaRes_3_HvsNH,
  gseaRes_4_HvsHN,
  gseaRes_5_NHvsN)

# transform to data.frame
DF.gseaRes_1_HvsN <- gseaRes_1_HvsN@result %>% dplyr::mutate(Setting="Original")
DF.gseaRes_2_HNvsN <- gseaRes_2_HNvsN@result %>% dplyr::mutate(Setting="Under21")
DF.gseaRes_3_HvsNH <- gseaRes_3_HvsNH@result %>% dplyr::mutate(Setting="Under1")
DF.gseaRes_4_HvsHN <- gseaRes_4_HvsHN@result %>% dplyr::mutate(Setting="SwitchHCAF")
DF.gseaRes_5_NHvsN <- gseaRes_5_NHvsN@result %>% dplyr::mutate(Setting="SwitchNCAF")

# conbined HypoCAF vs NormoCAF
DF.CombinedGseaRes_0 <- 
  dplyr::bind_rows(DF.gseaRes_1_HvsN, DF.gseaRes_2_HNvsN, DF.gseaRes_3_HvsNH) %>% 
  dplyr::select(c(Setting, Description, NES, p.adjust)) %>% 
  dplyr::mutate(Setting = factor(Setting, levels=c("Original", "Under21", "Under1")),
                Signif = case_when(p.adjust<0.05 ~ "*",
                                   TRUE ~ ""))
Range.NES = range(DF.CombinedGseaRes_0$NES)
DF.CombinedGseaRes_1 <-
  DF.CombinedGseaRes_0 %>% 
  pivot_wider(id_cols=Setting,
              names_from = Description, 
              values_from = NES) %>% 
  column_to_rownames(var="Setting") %>% as.matrix()
# hierarchical clustering : euclidean, ward.D2
d_hallmark = dist(t(DF.CombinedGseaRes_1), method="euclidean")
Hclust_hallmark <- hclust(d_hallmark, method = "ward.D2")
MethodLabel = "Hierarchical clustering was performed using Ward’s method (ward.D2) with Euclidean distances."
Vec.HallmarkOrder = colnames(DF.CombinedGseaRes_1)[Hclust_hallmark$order]
DF.CombinedGseaRes_0$Description = factor(DF.CombinedGseaRes_0$Description, levels=Vec.HallmarkOrder)
# dendrogram
dend_hallmark = as.dendrogram(Hclust_hallmark)
dend_data = dendro_data(dend_hallmark)
p_dend = 
  ggplot(dend_data$segments) +
  geom_segment(aes(x=y, y=x, xend=yend, yend=xend)) +
  scale_y_continuous(limits = c(0.5, length(Hclust_hallmark$labels) + 0.5),
                     breaks = seq_along(Hclust_hallmark$labels),
                     labels = Hclust_hallmark$labels,
                     expand = expansion(mult=c(0,0))) +
  scale_x_reverse(expand = expansion(mult=c(0.01, 0.01))) +
  theme_void()
TxtSize=10
p_heatmap =
  ggplot(DF.CombinedGseaRes_0, aes(x=Setting, y=Description, fill=NES)) +
  geom_tile(color="gray50") +
  geom_text(aes(label=case_when(p.adjust<0.001 ~ "***",
                                p.adjust<0.01 ~ "**",
                                p.adjust<0.05 ~ "*",
                                TRUE ~ "")),
            size=TxtSize*0.5, vjust=0.8) +
  labs(x=NULL, y=NULL) +
  scale_x_discrete(expand=expansion(mult=c(0,0)),
                   labels=c("Original"="Oxygen level at isolation",
                            "Under21"="Under 21% O2",
                            "Under1"="Under 1% O2")) +
  scale_y_discrete(expand=expansion(mult=c(0,0)),
                   position="right",
                   labels=
                     setNames(clean_gs_label_ForFig(DF.CombinedGseaRes_0$Description),
                              DF.CombinedGseaRes_0$Description)) +
  scale_fill_distiller(
    palette = "RdBu",
    direction = 1,
    limits=c(-max(abs(Range.NES)), max(abs(Range.NES))),
    breaks=c(-2,0,2),
    guide=guide_colorbar(direction="horizontal",
                         title.position="top",
                         #title.vjust=0.75,
                         title.hjust=0.5,
                         frame.colour="black",
                         ticks.colour="black",
                         reverse=FALSE)) +
  #coord_fixed(ratio=0.5) +
  theme(plot.margin = margin(t=0, r=0, b=0, l=0),
        legend.position = "bottom",
        axis.text.x = element_text(face="bold", color="black", angle=45, hjust=1, size=TxtSize*0.5),
        axis.text.y = element_text(face="bold", color="black", size=TxtSize*0.95),
        legend.title = element_text(face="bold", color="black", size=TxtSize*1.4),
        legend.text = element_text(face="bold", color="black", size=TxtSize*1.1),
        legend.background = element_rect(fill="transparent", color=NA),
        plot.background = element_rect(fill="transparent", color=NA),
        panel.background = element_rect(fill="transparent", color=NA),
        panel.border = element_rect(fill="transparent", color="black"))
Plot.joined = 
  p_dend + p_heatmap +
  plot_layout(widths=c(0.4, 0.2)) &
  theme(
    plot.background = element_rect(fill="transparent", color=NA),
    panel.background = element_rect(fill="transparent", color=NA))
ggsave(plot = Plot.joined,
       width=5.5, height=9,
       file="InVitroGSEA_sup4_new.png", dpi=500, bg="transparent")
## Main figure
Vec.HallmarkDisplay = c(#"HALLMARK_HYPOXIA",
  #"HALLMARK_APOPTOSIS",
  "HALLMARK_MYOGENESIS",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_COMPLEMENT",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_P53_PATHWAY",
  #"HALLMARK_UV_RESPONSE_UP",
  #"HALLMARK_MYC_TARGETS_V2",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_MITOTIC_SPINDLE",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_E2F_TARGETS")
DF.gseaRes_2 = DF.CombinedGseaRes_0 %>% 
  dplyr::filter(Description %in% Vec.HallmarkDisplay) %>% 
  dplyr::mutate(minuslog10adjp = (-1)*log10(p.adjust),
                Description = factor(Description, levels=rev(Vec.HallmarkDisplay)))
Range.NES.dot = range(DF.gseaRes_2$NES)
Range.log10P.dot = range(DF.gseaRes_2$minuslog10adjp)
TsizeMain = 10
p_main =
  ggplot(DF.gseaRes_2, aes(x=Setting, y=Description)) +
  geom_point(aes(size=minuslog10adjp, fill=NES),
             shape=21) +
  labs(x=NULL, y=NULL, size="-log10\n(Adj.P)") +
  scale_size_continuous(
    #limits=c(-log10(0.05), max(Range.log10P.dot)),
    range=c(1.5, 8),
    breaks=c(2, 10, 40),
    guide=guide_legend(
      order=1,
      direction="horizontal",
      title.position="top",
      label.position="bottom")) +
  scale_fill_distiller(
    palette="RdBu",
    direction=1,
    breaks=c(-2,0,2),
    limits=c(-max(abs(Range.NES.dot)),
             max(abs(Range.NES.dot))),
    guide=guide_colorbar(
      order=2,
      direction="horizontal",
      title.position="top",
      reverse=FALSE,
      frame.colour="black",
      ticks.colour="black")) +
  scale_x_discrete(labels=c("Original"="Oxygen level at isolation",
                            "Under21"="Under 21% O2",
                            "Under1"="Under 1% O2"),
                   position="top") +
  scale_y_discrete(
    position="right",
    labels=setNames(clean_gs_label_ForFig(DF.gseaRes_2$Description),
                    DF.gseaRes_2$Description)) +
  theme(axis.text.x = element_text(face="bold", color="black", angle=45, hjust=0, size=TsizeMain*1.2),
        axis.text.y = element_text(face="bold", color="black", size=TsizeMain*1.3),
        legend.key = element_blank(),
        legend.text = element_text(face="bold", color="black", size=TsizeMain*1.2),
        legend.title = element_text(face="bold", color="black", size=TsizeMain*1.2, hjust=0.5),
        legend.background = element_rect(fill="transparent", color=NA),
        plot.background = element_rect(fill="transparent", color=NA),
        panel.background = element_rect(fill="white", color=NA),
        panel.border = element_rect(fill=NA, color="black"))
saveRDS(p_main +
          scale_y_discrete(position = "right",
                           labels=setNames(clean_gs_label_ForFig(DF.gseaRes_2$Description),
                                           DF.gseaRes_2$Description)) +
          theme(legend.box = "vertical",
                legend.box.just = "right",
                legend.justification = "right") +
          scale_size_continuous(
            range=c(1.5, 8),
            breaks=c(2, 10, 40),
            guide=guide_legend(
              order=1,
              direction="vertical",
              title.position="top",
              label.position="left")) +
          scale_fill_distiller(
            palette="RdBu",
            direction=1,
            breaks=c(-2,0,2),
            limits=c(-max(abs(Range.NES.dot)),
                     max(abs(Range.NES.dot))),
            guide=guide_colorbar(
              order=2,
              direction="vertical",
              title.position="top",
              label.position="left",
              reverse=FALSE,
              frame.colour="black",
              ticks.colour="black")), 
        "Fig2Eleft_rightYaxis_new.rds")
ggsave(plot = p_main +
         coord_fixed(ratio=0.8),
       width=6.2, height=5.1,
       file="InVitroGSEA_main2_new.png", dpi=500, bg="transparent")


## Sup table
Func.DataProcess <- function(DataFrame, SettingName){
  DF.Processed = DataFrame %>% 
    dplyr::mutate(
      Setting=SettingName,
      Direction=ifelse(NES>0, "UP", "DOWN")) %>% 
    dplyr::filter(! (Description %in% 
                       c("iCAF_signature", "Buffa_Hypoxia", "Winter_Hypoxia", "KEGG_Cell_cycle"))) %>% 
    dplyr::select(
      c(Setting, Direction, Description, setSize, 
        NES, p.adjust, qvalue, leading_edge, core_enrichment)) %>% 
    dplyr::mutate(Direction = factor(Direction, levels=c("UP","DOWN"))) %>% 
    group_by(Direction) %>% 
    dplyr::arrange(p.adjust, .by_group=TRUE) %>% 
    set_rownames(NULL)
  return(DF.Processed)
}
DF.gseaResSup_HvsN_0 <- 
  bind_rows(Func.DataProcess(DF.gseaRes_1_HvsN, "Original"), 
            Func.DataProcess(DF.gseaRes_2_HNvsN, "Under21"), 
            Func.DataProcess(DF.gseaRes_3_HvsNH, "Under1")) %>% 
  dplyr::rename(In_Hypo_CAF = Direction) %>% 
  dplyr::mutate(Description = clean_gs_label_ForTable(Description))
DF.gseaResSup_Switch_0 <-
  bind_rows(Func.DataProcess(DF.gseaRes_4_HvsHN, "Switch_HCAF"), 
            Func.DataProcess(DF.gseaRes_5_NHvsN, "Switch_NCAF")) %>% 
  dplyr::rename(Upreg_in = Direction) %>% 
  dplyr::mutate(Upreg_in = 
                  Upreg_in %>% 
                  str_replace(pattern="UP", replacement="Under1pct") %>% 
                  str_replace(pattern="DOWN", replacement="Under21pct")
                ) %>%
  dplyr::mutate(Description = clean_gs_label_ForTable(Description))
write.csv(DF.gseaResSup_HvsN_0, file="SupplementalTable_bulk_GSEA(Hallmark)_HvsN_ASCII_new.csv", row.names=FALSE)
write.csv(DF.gseaResSup_Switch_0, file="SupplementalTable_bulk_GSEA(Hallmark)_Switch_ASCII_new.csv", row.names=FALSE)


#----------------------------------------------------
#
#   6. Human hallmarks GSEA (Before vs After switching)
#
#----------------------------------------------------
gseaRes_4_HvsHN <- Func.GSEA(Vec.namedFC_4_HvsHN_rmNA)
gseaRes_5_NHvsN <- Func.GSEA(Vec.namedFC_5_NHvsN_rmNA)

DF.gseaRes_4_HvsHN <- gseaRes_4_HvsHN@result %>% dplyr::mutate(Culture="HypoCAF")
DF.gseaRes_5_NHvsN <- gseaRes_5_NHvsN@result %>% dplyr::mutate(Culture="NormoCAF")

DF.gseaRes_SwitchHypoCAF <-
  rbind(slice_max(DF.gseaRes_4_HvsHN, n=5, order_by=NES),
        slice_min(DF.gseaRes_4_HvsHN, n=5, order_by=NES)) %>% 
  dplyr::arrange(-NES) %>% 
  dplyr::mutate(Description = factor(Description, levels=unique(Description)),
                minusLog10adjP = (-1)*log10(p.adjust))
DF.gseaRes_SwitchNormoCAF <- 
  rbind(slice_max(DF.gseaRes_5_NHvsN, n=5, order_by=NES),
        slice_min(DF.gseaRes_5_NHvsN, n=5, order_by=NES)) %>% 
  dplyr::arrange(-NES) %>% 
  dplyr::mutate(Description = factor(Description, levels=unique(Description)),
                minusLog10adjP = (-1)*log10(p.adjust))
Range.NES.Switch = range(c(DF.gseaRes_SwitchHypoCAF$NES,
                           DF.gseaRes_SwitchNormoCAF$NES))
Range.minuslog10adjP.Switch = 
  range(c(DF.gseaRes_SwitchHypoCAF$minusLog10adjP,
          DF.gseaRes_SwitchNormoCAF$minusLog10adjP))
Common.Scale.Size =
  scale_size_continuous(
    name = "-log10(Adj.P)",
    breaks = c(2, 10, 30),
    limits = c(0, max(Range.minuslog10adjP.Switch)),
    range = c(0, 5.8),
    guide = guide_legend(
      override.aes = list(fill="white", color="black"),
      label.position = "bottom"))
Common.vline = geom_vline(xintercept=0)
Common.Theme = theme(
  axis.title.y = element_blank(),#element_text(face="bold", color="black"),
  axis.text.y = element_text(face="bold", color="black", size=11),
  axis.line.y = element_line(color="black"),
  legend.key = element_blank(),
  legend.text = element_text(face="bold", color="black"),
  legend.title = element_text(face="bold", color="black"),
  legend.background = element_rect(fill="transparent", color=NA),
  plot.background = element_rect(fill="transparent", color=NA),
  panel.background = element_rect(fill="transparent", color=NA),
  panel.grid = element_blank(),
  panel.border = element_rect(fill="transparent", color=NA))
Erase.lgd = theme(legend.position="none")
Plot_NormoCAF =
  ggplot(DF.gseaRes_SwitchNormoCAF, aes(x=NES, y=Description)) +
  geom_segment(aes(x=0, xend=NES, 
                   y=Description, yend=Description)) +
  geom_point(aes(size=minusLog10adjP), shape=21, fill="#f03b20") +
  scale_y_discrete(
    labels=setNames(clean_gs_label_ForFig(DF.gseaRes_SwitchNormoCAF$Description) %>% 
                      str_replace(pattern="Epithelial mesenchymal transition",
                                  replacement="EMT"),
                    DF.gseaRes_SwitchNormoCAF$Description),
    expand = expansion(mult = c(0.1, 0.1)),
    position = "right") +
  Common.vline + 
  scale_x_continuous(
    position="top",
    limits=c(-max(abs(Range.NES.Switch)),max(abs(Range.NES.Switch))),
    expand=expansion(mult=c(0.05, 0.05))) +
  Common.Scale.Size +
  theme(legend.position = "bottom",
        plot.margin = margin(t=0, r=0, b=0.1, l=0, "cm"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()) +
  Common.Theme
Plot_HypoCAF =
  ggplot(DF.gseaRes_SwitchHypoCAF, aes(x=NES, y=Description)) +
  geom_segment(aes(x=0, xend=NES, 
                   y=Description, yend=Description)) +
  geom_point(aes(size=minusLog10adjP), shape=21, fill="#2c7fb8") +
  scale_y_discrete(
    labels=setNames(clean_gs_label_ForFig(DF.gseaRes_SwitchHypoCAF$Description) %>% 
                      str_replace(pattern="Epithelial mesenchymal transition",
                                  replacement="EMT"),
                    DF.gseaRes_SwitchHypoCAF$Description),
    expand = expansion(mult = c(0.1, 0.1)),
    position = "right") +
  Common.vline + 
  scale_x_continuous(
    position="bottom",
    limits=c(-max(abs(Range.NES.Switch)),max(abs(Range.NES.Switch))),
    expand=expansion(mult=c(0.05, 0.05))) + 
  Common.Scale.Size +
  theme(legend.position = "bottom",
        plot.margin = margin(t=0.1, r=0, b=0, l=0, "cm"),
        axis.title.x = element_text(face="bold", color="black"),
        axis.text.x = element_text(face="bold", color="black", size=9),
        axis.line.x = element_line(color="black")) +
  Common.Theme
aligned <-
  cowplot::align_plots(
  Plot_NormoCAF + Erase.lgd,
  Plot_HypoCAF + Erase.lgd,
  align="v", axis="r")
plot_alignedpanels = 
  plot_grid(aligned[[1]], aligned[[2]],
            ncol=1, rel_heights=c(1, 1)) 
plot_final =  
  plot_grid(plot_alignedpanels, get_legend(Plot_NormoCAF), 
            ncol=1, rel_heights=c(5,1)) 
ggsave(plot = plot_final,
       width=4, height=7,
       file="InVitroGSEA_BeforeAfterSwitching2_new.png", dpi=500, bg="transparent")
Common.Scale.Size_Final =
  scale_size_continuous(
    name = "-log10(Adj.P)",
    breaks = c(2, 10, 30),
    limits = c(0, max(Range.minuslog10adjP.Switch)),
    range = c(0, 10),
    guide = guide_legend(
      direction = "horizontal",
      override.aes = list(fill="white", color="black"),
      title.position = "top",
      label.position = "bottom"))
saveRDS(Plot_NormoCAF + Common.Scale.Size_Final, file = "Fig2F_top_new.rds")
saveRDS(Plot_HypoCAF + Common.Scale.Size_Final, file = "Fig2F_bottom_new.rds")
saveRDS(get_legend(Plot_NormoCAF + Common.Scale.Size_Final), file = "Fig2F_legend_new.rds")

#----------------------------------------------------
#
#   7. Human hallmarks GSEA, enrichment plot
#
#----------------------------------------------------
library(enrichplot)
#List.gseaRes = readRDS(paste0(DirOutputMSigDB,"[ListRDS]_gseaHallmark_Results_",SetName,".rds"))
Cols.SenesEnrich = c(HALLMARK_P53_PATHWAY="#d7191c",
                     HALLMARK_INTERFERON_GAMMA_RESPONSE="#fdae61",
                     HALLMARK_INTERFERON_ALPHA_RESPONSE="#006e00",
                     HALLMARK_INFLAMMATORY_RESPONSE="#2b83ba")
HMs = c(HALLMARK_E2F_TARGETS="#D55E00",
        HALLMARK_G2M_CHECKPOINT="#4E79A7",
        HALLMARK_MYC_TARGETS_V2="#E69F00",
        HALLMARK_MYC_TARGETS_V1="#7B3294",
        HALLMARK_MITOTIC_SPINDLE="#B3DE69")

Func.Penrich_Switch = function(CompNum,Comparison){
  List.HMs = list(
    "Prolif.HvsN" = 
      c(HALLMARK_E2F_TARGETS="#D55E00", HALLMARK_G2M_CHECKPOINT="#4E79A7",
        HALLMARK_MITOTIC_SPINDLE="#B3DE69"),
    "Senes.HvsN" = 
      c(HALLMARK_INFLAMMATORY_RESPONSE="#2b83ba",
        HALLMARK_P53_PATHWAY="#d7191c", HALLMARK_INTERFERON_GAMMA_RESPONSE="#fdae61",
        HALLMARK_INTERFERON_ALPHA_RESPONSE="#006e00"),
    "Prolif.BAswitch" =
      c(HALLMARK_E2F_TARGETS="#D55E00", HALLMARK_G2M_CHECKPOINT="#4E79A7",
        HALLMARK_MYC_TARGETS_V2="#E69F00", HALLMARK_MYC_TARGETS_V1="#7B3294"))
  HMs = List.HMs[[Comparison]]
  Filename = case_when(
    CompNum==1 & Comparison=="Prolif.HvsN" ~ "Original_Prolif_new",
    CompNum==1 & Comparison=="Senes.HvsN" ~ "Original_Senes_new",
    CompNum==2 & Comparison=="Prolif.HvsN" ~ "Under1_Prolif_new",
    CompNum==2 & Comparison=="Senes.HvsN" ~ "Under1_Senes_new",
    CompNum==3 & Comparison=="Prolif.HvsN" ~ "Under21_Prolif_new",
    CompNum==3 & Comparison=="Senes.HvsN" ~ "Under21_Senes_new",
    CompNum==4 ~ "SwitchHypoCAF_new",
    CompNum==5 ~ "SwitchNormoCAF_new")
  Xlab = case_when(CompNum%in%c(1,2,3) ~ "Upregulated → Downregulated",
                   CompNum%in%c(4,5) ~ "21% O2 → 1% O2")
  LgdLab = clean_gs_label_ForFig(names(HMs)) %>%
    str_replace(pattern="Interferon alpha", replacement="IFNα") %>% 
    str_replace(pattern="Interferon gamma", replacement="IFNγ") %>% 
    setNames(names(HMs))
  EnrichPlot_0 = 
    enrichplot::gseaplot2(
      List.gseaRes[[CompNum]], 
      geneSetID=names(HMs), 
      subplots=1:2)
  p.res = 
    EnrichPlot_0[[1]] +
    labs(y="Enrichment score (ES)", color=NULL) +
    geom_hline(yintercept=0, color="gray70", alpha=0.5) +
    scale_color_manual(
      values=HMs,
      labels=LgdLab) +
    theme(axis.title.y=element_text(face="bold", hjust=0.5, size=18),
          axis.text.y=element_text(face="bold", size=16),
          legend.text=element_text(face="bold", size=15),
          legend.background = element_rect(fill="transparent", color=NA),
          legend.key = element_blank(),
          plot.background = element_rect(fill="transparent", color=NA),
          panel.background = element_rect(fill="transparent", color=NA),
          panel.grid.major.x=element_blank(),
          panel.grid.minor.x=element_blank())
  p.res$layers[[1]]$aes_params$size = 1.3
  p.res$data$Description = 
    factor(p.res$data$Description, levels=rev(names(HMs)))
  p.gene = EnrichPlot_0[[2]] + 
    labs(x=paste0("Gene rank\n(",Xlab,")"), color=NULL) +
    geom_hline(yintercept=seq_along(HMs), color="gray80") +
    scale_color_manual(values=HMs,
                       labels=LgdLab) +
    theme(legend.position="none",
          axis.title.x = element_text(face="bold", size=18),
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          #legend.text=element_text(face="bold", size=15),
          #legend.background = element_rect(fill="transparent", color=NA),
          #legend.key = element_blank(),
          plot.background = element_rect(fill="transparent", color=NA),
          panel.background = element_rect(fill="transparent", color=NA))
  EnrichPlot_Panels = 
    (p.res / p.gene) +
    plot_layout(heights = c(1.8, 0.3))
  EnrichPlot_PanelsAndLgd =
    (EnrichPlot_Panels | guide_area()) +
    plot_layout(
      guides="collect", 
      widths=c(1.0, 1.0)) &
    theme(
      plot.background = element_rect(fill="transparent", color = NA))
  ggsave(plot=EnrichPlot_PanelsAndLgd,
         file=paste0("EnrichmentPlot-",Filename,".png"),
         width=8.5, height=5, bg="transparent")
  # --- ここから追加：cowplotで扱いやすいRDS群 ---
  Erase.LgdLab = theme(legend.position = "none",
                       axis.title.x = element_blank(),
                       axis.title.y = element_blank(),
                       axis.text.x = element_blank(),
                       axis.text.y = element_blank(),
                       axis.ticks.x = element_blank(),
                       plot.margin = margin(0, 0, 0, 0))
  
  p.res_nlgd  = p.res  + Erase.LgdLab
  p.gene_nlgd = p.gene + Erase.LgdLab
  
  lgd = cowplot::get_legend( #get_plot_component("guide-box-right",...
    p.res + 
      guides(color = guide_legend(ncol = 1)) +
      theme(
        legend.position = "right",
        legend.justification = "left",
        legend.box.just = "left",
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0),
        plot.background = element_rect(fill="transparent", color=NA),
        legend.background = element_rect(fill="transparent", color=NA),
        legend.key = element_blank()
      ))
  panel_cow = cowplot::plot_grid(
    p.res_nlgd,
    p.gene_nlgd,
    ncol = 1,
    rel_heights = c(1.8, 0.3),
    align = "v",
    axis = "lr"
  )
  # --- 保存 ---
  saveRDS(p.res,       file = paste0("Fig2E_", Comparison, "_", CompNum, "_res_new.rds"))
  saveRDS(p.gene,      file = paste0("Fig2E_", Comparison, "_", CompNum, "_gene_new.rds"))
  saveRDS(p.res_nlgd,  file = paste0("Fig2E_", Comparison, "_", CompNum, "_res_nlgd_new.rds"))
  saveRDS(p.gene_nlgd, file = paste0("Fig2E_", Comparison, "_", CompNum, "_gene_nlgd_new.rds"))
  saveRDS(lgd,         file = paste0("Fig2E_", Comparison, "_", CompNum, "_legend_new.rds"))
  saveRDS(panel_cow,   file = paste0("Fig2E_", Comparison, "_", CompNum, "_panel_cow_new.rds"))
}
Func.Penrich_Switch(CompNum=1, Comparison="Prolif.HvsN")
Func.Penrich_Switch(CompNum=1, Comparison="Senes.HvsN")
Func.Penrich_Switch(CompNum=2, Comparison="Prolif.HvsN")
Func.Penrich_Switch(CompNum=2, Comparison="Senes.HvsN")
Func.Penrich_Switch(CompNum=3, Comparison="Prolif.HvsN")
Func.Penrich_Switch(CompNum=3, Comparison="Senes.HvsN")
Func.Penrich_Switch(CompNum=4, Comparison="Prolif.BAswitch")
Func.Penrich_Switch(CompNum=5, Comparison="Prolif.BAswitch")

#----------------------------------------------------
#
#   8. Alined Bulk RNA seq figures
#
#----------------------------------------------------
library(grid)
library(gridExtra)

NoLgd <- theme(legend.position="none")
RtLgd <- theme(legend.position="right",
               legend.text = element_text(color="black", face="bold", size=12),
               legend.title = element_text(color="black", face="bold", size=13))
ZeroMgn <- theme(plot.margin=margin(0, 0, 0, 0))
AxText.X <- theme(axis.text.x = element_blank(),
                  axis.title.x = element_blank())
AxText.Y <- theme(axis.text.y = element_text(color = "black", face = "bold", size = 16))
EnrichTheme <- theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  axis.text.x  = element_blank())
ThemeEraseYaxis <- theme(
  axis.text.y  = element_blank(),
  axis.title.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.line.y  = element_blank())
ThemeEraseXYaxistext <- theme(
  plot.title = element_blank(),
  axis.text.x  = element_blank(),
  axis.text.y  = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y = element_blank())
ThemeG <- theme(
  axis.text.y  = element_blank(),
  axis.title.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.line.y  = element_blank())
center_legend_in_block <- function(lgd) {
  arrangeGrob(
    sp,
    lgd,
    sp,
    ncol = 1,
    heights = unit.c(
      unit(1, "null"),
      grobHeight(lgd),
      unit(1, "null")))
}
extract_axis <- function(p, LgdSide) {  #LgdSide: "l"or"r
  gb <- ggplotGrob(p)
  idx <- which(gb$layout$name == paste0("axis-",LgdSide))
  gb$grobs[[idx]]}
center_grob_in_block <- function(gb) {
  arrangeGrob(
    sp,
    gb,
    sp,
    ncol = 1,
    heights = unit.c(
      unit(1, "null"),
      grobHeight(gb),
      unit(1, "null") ) )
}

#+#+#+#+#+# read panels and legends #+#+#+#+#+#
p2B.top <- readRDS("Fig2B_Volc_Upper_new.rds") + ThemeEraseXYaxistext + ZeroMgn
p2B.mid <- readRDS("Fig2B_Volc_Mid_new.rds") + ThemeEraseXYaxistext + ZeroMgn
p2B.bottom <-  readRDS("Fig2B_Volc_Lower_new.rds") + ThemeEraseXYaxistext + ZeroMgn
p2C <- readRDS("Fig2C_new.rds") + ZeroMgn + AxText.Y + AxText.X
p2D <- readRDS("Fig2D_new.rds") + ZeroMgn + AxText.Y + AxText.X
p2E <- readRDS("Fig2Eleft_rightYaxis_new.rds") + ZeroMgn + AxText.Y + AxText.X + 
  scale_x_discrete(position = "bottom")
p2F.sen1 <- readRDS("Fig2E_Senes.HvsN_1_panel_cow_new.rds")
p2F.sen2 <- readRDS("Fig2E_Senes.HvsN_2_panel_cow_new.rds")
p2F.sen3 <- readRDS("Fig2E_Senes.HvsN_3_panel_cow_new.rds")
p2F.pro1 <- readRDS("Fig2E_Prolif.HvsN_1_panel_cow_new.rds")
p2F.pro2 <- readRDS("Fig2E_Prolif.HvsN_2_panel_cow_new.rds")
p2F.pro3 <- readRDS("Fig2E_Prolif.HvsN_3_panel_cow_new.rds")
lgd.sen <- readRDS("Fig2E_Senes.HvsN_1_legend_new.rds")
lgd.pro <- readRDS("Fig2E_Prolif.HvsN_1_legend_new.rds")
p2G.top <- readRDS("Fig2F_top_new.rds") + NoLgd + ZeroMgn + AxText.Y
p2G.bottom <- readRDS("Fig2F_bottom_new.rds") + NoLgd + ZeroMgn + AxText.Y


#+#+#+#+#+# Grob #+#+#+#+#+#
sp <- grid::rectGrob(gp = gpar(col = NA))
p2B.top_gb = ggplotGrob(p2B.top + ZeroMgn + NoLgd)
p2B.mid_gb =  ggplotGrob(p2B.mid + ZeroMgn + NoLgd)
p2B.bottom_gb = ggplotGrob(p2B.bottom + ZeroMgn + NoLgd)
p2C_gb <- ggplotGrob(p2C + NoLgd + AxText.X)
p2D_gb <- ggplotGrob(p2D + NoLgd + AxText.X)
p2E_gb <- ggplotGrob(p2E + NoLgd + AxText.X)
p2C_panel_gb <- ggplotGrob(p2C + ThemeEraseYaxis + NoLgd)
p2C_axis_gb <- extract_axis(p2C + ZeroMgn + NoLgd, "l")
p2D_panel_gb <- ggplotGrob(p2D + ThemeEraseYaxis + NoLgd)
p2D_axis_gb <- extract_axis(p2D + ZeroMgn + NoLgd, "r")
p2E_panel_gb <- ggplotGrob(p2E + ThemeEraseYaxis + NoLgd)
p2E_axis_gb <- extract_axis(p2E + ZeroMgn + NoLgd, "r")
CommonHeight.p2C <- unit.pmax(p2C_panel_gb$heights, p2C_axis_gb$heights)
CommonHeight.p2D <- unit.pmax(p2D_panel_gb$heights, p2D_axis_gb$heights)
CommonHeight.p2E <- unit.pmax(p2E_panel_gb$heights, p2E_axis_gb$heights)
p2C_panel_gb$heights <- CommonHeight.p2C
p2C_axis_gb$heights <- CommonHeight.p2C
p2D_panel_gb$heights <- CommonHeight.p2D
p2D_axis_gb$heights <- CommonHeight.p2D
p2E_panel_gb$heights <- CommonHeight.p2E
p2E_axis_gb$heights <- CommonHeight.p2E
#CommonWidth.CDE <- unit.pmax(p2C_gb$widths, p2D_gb$widths, p2E_gb$widths)
#p2C_gb$widths <- CommonWidth.CDE
#p2D_gb$widths <- CommonWidth.CDE
#p2E_gb$widths <- CommonWidth.CDE
p2F.upper <- cowplot::plot_grid(
  p2F.sen1, sp, p2F.sen3, sp, p2F.sen2,  # original, 21%, 1%
  ncol = 5,
  align = "h",
  axis = "tb",
  rel_widths = c(3, 1, 3, 1, 3) )
p2F.lower <- cowplot::plot_grid(
  p2F.pro1, sp, p2F.pro3, sp, p2F.pro2,  # original, 21%, 1%
  ncol = 5,
  align = "h",
  axis = "tb",
  rel_widths = c(3, 1, 3, 1, 3) )
p2F <- cowplot::plot_grid(
  p2F.upper,
  sp,
  p2F.lower,
  ncol = 1,
  align = "v",
  axis = "lr",
  rel_heights = c(5, 2, 5))
p2F_gb <- ggplotGrob(p2F + ZeroMgn)
p2G.top_panel_gb <- ggplotGrob(p2G.top + AxText.X + ThemeG)
p2G.top_axis_gb <- extract_axis(p2G.top + ZeroMgn, "r")
p2G.bottom_panel_gb <- ggplotGrob(p2G.bottom + AxText.X + ThemeG)
p2G.bottom_axis_gb <- extract_axis(p2G.bottom + ZeroMgn, "r")

lgd.2C <- cowplot::get_legend(p2C + RtLgd + ZeroMgn)
lgd.2D <- cowplot::get_legend(p2D + RtLgd + ZeroMgn)
lgd.2E <- cowplot::get_legend(p2E + ZeroMgn + RtLgd)
lgd.2G <- cowplot::get_legend(
  p2G.top + RtLgd + ZeroMgn +
    theme(legend.background = element_rect(fill = "transparent", color = NA),
          legend.key = element_blank()))
lgd.2C.center_gb <- center_legend_in_block(lgd.2C)
lgd.2D.center_gb <- center_legend_in_block(lgd.2D)
lgd.2E.center_gb <- center_legend_in_block(lgd.2E)

#+#+#+#+#+#   align parameters   #+#+#+#+#+#
TotalHeight = 11 # inch
gap <- 24/25.4   # ≒0.79 inch
Height_CDE <- TotalHeight - 2*gap # inch
Height_B <- 6.3*gap
Height_C <- 3.036745 # 10 lines, inch
Height_D <- 2.125722 # 7 lines, inch
Height_E <- 3.947769 # 13 lines, inch
PanelWidth.CDE <- 1.0 
AxisWidth.CDE <- 4.2
Height_G <-  Height_C + Height_D + gap # inch
gap_inG <- 1/25.4
Height_G_each <- (Height_G - gap_inG) / 2
Height_F <- (Height_B + 0.5*gap + Height_C) - (Height_G + 0.7*gap) # inch

#+#+#+#+#+#   align grobs left   #+#+#+#+#+#
p2B_gb <- arrangeGrob(
  sp, p2B.top_gb, 
  sp, sp,
  sp, p2B.mid_gb, 
  sp, sp,
  sp, p2B.bottom_gb,
  ncol = 2, nrow = 5,
  heights = unit.c(
    unit(Height_B*3/11, "in"), unit(Height_B*1/11, "in"),
    unit(Height_B*3/11, "in"), unit(Height_B*1/11, "in"),
    unit(Height_B*3/11, "in")),
  widths = unit.c(unit((PanelWidth.CDE + AxisWidth.CDE)*0.21, "in"),
                  unit((PanelWidth.CDE + AxisWidth.CDE)*0.79, "in")))
p2C_gb <- arrangeGrob(
  p2C_axis_gb, p2C_panel_gb,
  ncol = 2,
  widths = unit.c(
    unit(AxisWidth.CDE, "in"), unit(PanelWidth.CDE, "in")))
left_plotcol <- arrangeGrob(
  sp,
  p2B_gb,
  sp,
  p2C_gb,
  ncol = 1, 
  heights = unit.c(
    unit(gap + Height_D + gap + Height_E + gap, "in"),
    unit(Height_B, "in"), 
    unit(0.5*gap, "in"), 
    unit(Height_C, "in")))
left_legendcol <- arrangeGrob(
  sp,
  sp,
  sp,
  lgd.2C.center_gb, 
  ncol = 1,
  heights = unit.c(
    unit(gap + Height_D + gap + Height_E + gap, "in"),
    unit(Height_B, "in"), 
    unit(0.5*gap, "in"), 
    unit(Height_C, "in")))

#+#+#+#+#+#   align grobs right   #+#+#+#+#+#
right_legendcol <- arrangeGrob(
  arrangeGrob(sp, lgd.2D.center_gb, sp, ncol = 3, widths = unit.c(unit(0.10, "in"), unit(1.2, "in"), unit(0.30, "in"))),
  arrangeGrob(sp, lgd.2E.center_gb, sp, ncol = 3, widths = unit.c(unit(0.10, "in"), unit(1.2, "in"), unit(0.30, "in"))),
  sp,
  sp,
  ncol = 1, nrow = 4,
  heights = unit.c(
    unit(gap + Height_D + gap, "in"),
    unit(Height_E, "in"), 
    unit(gap, "in"), 
    unit(Height_F + 0.7*gap + Height_G, "in")))
p2DE_gb <- arrangeGrob(
  p2D_panel_gb, p2D_axis_gb, 
  sp,           sp,
  p2E_panel_gb, p2E_axis_gb, 
  ncol = 2, nrow = 3,
  heights = unit.c( 
    unit(Height_D, "in"), 
    unit(gap, "in"),
    unit(Height_E, "in")),
  widths = unit.c(unit(PanelWidth.CDE, "in"), unit(AxisWidth.CDE, "in")))
p2FG_gb <- arrangeGrob(
  p2F_gb,              sp, 
  sp,                  sp,
  p2G.top_panel_gb,    sp,
  sp,                  sp,
  p2G.bottom_panel_gb, sp,
  ncol = 2, nrow = 5,
  heights = unit.c(
    unit(Height_F, "in"), 
    unit(0.7*gap, "in"),
    unit(Height_G_each, "in"), 
    unit(gap_inG, "in"),
    unit(Height_G_each, "in")),
  widths = unit.c(unit(4.1, "in"), unit(1.1, "in")))
right_plotcol <- arrangeGrob(
  sp,
  p2DE_gb,
  sp,
  p2FG_gb,
  ncol = 1, nrow = 4,
  heights = unit.c(
    unit(gap, "in"),
    unit(Height_D + gap + Height_E, "in"),
    unit(gap, "in"),
    unit(Height_F + 0.7*gap + Height_G, "in")))
lgd.2G.center <- arrangeGrob(
  sp,
  lgd.2G,
  sp,
  ncol = 1,
  heights = unit.c(
    unit(1, "null"),
    grobHeight(lgd.2G),
    unit(1, "null") ) )
p2G.top_axis.center <- center_grob_in_block(p2G.top_axis_gb)
p2G.bottom_axis.center <- center_grob_in_block(p2G.bottom_axis_gb)
lgd.2F <- arrangeGrob(
  lgd.sen, 
  sp,
  lgd.pro,
  ncol = 1,
  heights = unit.c(
    unit(5, "null"),
    unit(2, "in"),
    unit(5, "null") ),
  vp = viewport(just = c("left", "center") ) )
right_sidecol <- arrangeGrob(
  arrangeGrob(sp, lgd.2G.center, ncol = 2, widths = unit.c(unit(0.3, "in"), unit(0.9, "in"))),
  lgd.2F,
  sp,
  arrangeGrob(sp, p2G.top_axis.center, ncol = 2, widths = unit.c(unit(0.2, "in"), unit(1.0, "in"))),
  sp,
  arrangeGrob(sp, p2G.bottom_axis.center, ncol = 2, widths = unit.c(unit(0.2, "in"), unit(1.0, "in"))),
  ncol = 1,
  heights = unit.c(
    unit(gap + Height_D + gap + Height_E + gap, "in"), 
    unit(Height_F, "in"),
    unit(0.7*gap, "in"),
    unit(Height_G_each, "in"), 
    unit(gap_inG, "in"),
    unit(Height_G_each, "in") ) )



#+#+#+#+#+#   finalize   #+#+#+#+#+#
png("aligned_plot_gap24_v3_new.png", 
    width = 22, 
    height = (Height_B + 0.5*gap + Height_C) + (gap + Height_E + gap + Height_D + gap), 
    units = "in", 
    res = 500, bg="transparent")
grid.arrange(left_plotcol, left_legendcol, right_legendcol, right_plotcol, right_sidecol, 
             ncol = 5,
             widths = unit.c(unit((AxisWidth.CDE + PanelWidth.CDE), "in"), 
                             unit(1.3, "in"),
                             unit(1.6, "in"),
                             unit(5.2, "in"),
                             unit(1.2, "in")))
dev.off()

#----------------------------------------------------
#
#   9. myCAF, iCAF, apCAF signatures (GSEA)
#
#----------------------------------------------------
library(clusterProfiler)
TermToGene.Ely_0 <- read.csv(file="[1]dataset/ElyadaSupTable_S22_Orthologs.csv", row.names=1)
TermToGene.Ely_1 <- 
  data.frame(term = TermToGene.Ely_0$Subtype,
             gene = TermToGene.Ely_0$ENTREZID)
class(TermToGene.Ely_1$gene)
TermToGene.Ely_2 = na.omit(TermToGene.Ely_1)

# 0. Run GSEA
Func.GSEA.Ely = function(GeneList){
  set.seed(1234)
  GSEA(geneList=GeneList,
       minGSSize=10,
       TERM2GENE=TermToGene.Ely_2,
       pvalueCutoff=1,
       verbose=FALSE, 
       eps=0)
}
gseaRes_1_HvsN = Func.GSEA.Ely(Vec.namedFC_1_HvsN_rmNA)
gseaRes_2_HNvsN = Func.GSEA.Ely(Vec.namedFC_2_HNvsN_rmNA)
gseaRes_3_HvsNH = Func.GSEA.Ely(Vec.namedFC_3_HvsNH_rmNA)
gseaRes_4_HvsHN = Func.GSEA.Ely(Vec.namedFC_4_HvsHN_rmNA)
gseaRes_5_NHvsN = Func.GSEA.Ely(Vec.namedFC_5_NHvsN_rmNA)

DF.Res.Elyada_0 = 
  bind_rows(
    gseaRes_1_HvsN@result %>% cbind(Condition = "Original"),
    gseaRes_2_HNvsN@result %>% cbind(Condition = "Under1"),
    gseaRes_3_HvsNH@result %>% cbind(Condition = "Under21"))　%>% 
  dplyr::mutate(
    Description = factor(Description, levels=c("myCAF","iCAF","apCAF")),
    Condition = factor(Condition, levels=c("Original","Under21","Under1")),
    Signif = case_when(p.adjust<0.001 ~ "***",
                       p.adjust<0.01 ~ "**",
                       p.adjust<0.05 ~ "*",
                       TRUE ~ "n.s."),
    Enrichment_in_HypoCAF = ifelse(NES>0, "UP", "DOWN"))
DF.Res.Elyada_1 = 
  DF.Res.Elyada_0 %>% 
  dplyr::select(c(Condition, Enrichment_in_HypoCAF,Description,
                  setSize, NES, p.adjust, qvalue, leading_edge, core_enrichment)) %>% 
  dplyr::mutate(Description = factor(Description, levels=c("myCAF", "iCAF","apCAF")),
                Condition = factor(Condition, levels=c("Original", "Under21", "Under1"))) %>% 
  dplyr::arrange(Condition) %>% 
  group_by(Condition) %>% 
  dplyr::arrange(Description, .by_group = TRUE)　%>% 
  ungroup()
write.csv(DF.Res.Elyada_1,
          file="SuppTable_GSEA_Elyada_new.csv",
          row.names=F)
Range.NES = range(DF.Res.Elyada_0$NES)
Plot.Elyada = 
  ggplot(DF.Res.Elyada_0, aes(y=Description, x=NES)) +
  geom_segment(aes(yend=Description), xend=0, linewidth=0.35) +
  geom_point(aes(fill=NES), shape=21, size=9.5, stroke=0.5) +
  geom_text(aes(label=Signif,
                x= NES + ifelse(NES>0, 0.35, -0.35),
                hjust = ifelse(NES>0, 0, 1),
                vjust = ifelse(p.adjust<0.05, 0.8, 0.5)),
            fontface="bold", size=5.5) +
  geom_vline(xintercept=0, linewidth = 0.8) +
  labs(y=NULL) +
  scale_x_continuous(limits=c(-1.3*(max(abs(Range.NES))),1.3*max(abs(Range.NES))),
                     #expand=expansion(mult=c(0.1, 0.1)),
                     breaks=c(-2, -1, 0, 1, 2),
                     sec.axis=dup_axis(labels=NULL, name=NULL)) +
  scale_y_discrete(sec.axis=dup_axis(labels=NULL, name=NULL)) +
  scale_fill_distiller(palette="RdBu", 
                       direction=1,
                       limits=c(-(max(abs(Range.NES))),max(abs(Range.NES))),
                       breaks=c(-2, -1, 0, 1, 2),
                       guide=guide_colorbar(direction="vertical",
                                            title.position="top",
                                            title.hjust=0,
                                            frame.colour="black",
                                            ticks.colour="black")) +
  facet_grid(Condition ~ ., switch = "y") +
  theme(panel.spacing = unit(5, "mm"),
        axis.text.x = element_text(face="bold", color="black", size=12),
        axis.text.y = element_text(face="bold", color="black", size=15),
        axis.title = element_text(face="bold", color="black", size=18),
        axis.ticks.x.top = element_blank(),
        axis.ticks.y.right = element_blank(),
        legend.title = element_text(face="bold", color="black", size=15),
        legend.text = element_text(face="bold", color="black", size=12),
        plot.background = element_rect(fill="transparent", color=NA),
        strip.placement = "outside",
        strip.background = element_rect(fill="transparent", color=NA),
        legend.background = element_rect(fill="transparent", color=NA),
        panel.background = element_rect(fill="transparent", color=NA),
        panel.grid = element_blank(),
        panel.border = element_rect(fill="transparent", color=NA),
        axis.line = element_line(color="black"))
ggsave(plot = Plot.Elyada,
       file = "Elyada_new.png",
       width=6, height=6, dpi=500, bg="transparent")

Plot.Elyada.Tile <- 
  ggplot(DF.Res.Elyada_0,
         aes(x = Condition, y = Description)) +
  geom_tile(aes(fill = NES), color = "gray50") +
  geom_text(aes(label = Signif,
                vjust = ifelse(p.adjust < 0.05, 0.9, 0.4),
                size = ifelse(p.adjust < 0.05, 8.0, 7.3)), 
            fontface = "bold") +
  labs(x = NULL, y = NULL) +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  scale_y_discrete(expand = expansion(mult = c(0, 0))) +
  scale_fill_distiller(
    palette="RdBu",
    direction=1,
    limits=c(-(max(abs(Range.NES))),max(abs(Range.NES))),
    breaks=c(1, 0, -1),
    guide = guide_colorbar(
      title.position="top",
      title.hjust=0.5,
      frame.colour="black", 
      ticks.colour="black",
      label.position="right")) +
  scale_size_identity(guide = "none") +
  theme(
    axis.text.x = element_text(face = "bold", color = "black", size = 10),
    axis.text.y = element_text(face = "bold", color = "black", size = 15),
    legend.title = element_text(face="bold", color="black", size=12),
    legend.text = element_text(face="bold", color="black", size=10),
    plot.background = element_rect(fill="transparent", color=NA),
    legend.background = element_rect(fill="transparent", color=NA),
    panel.background = element_rect(fill="transparent", color=NA),
    panel.grid = element_blank(),
    panel.border = element_rect(fill="transparent", color = "black"))
ggsave(plot = Plot.Elyada.Tile,
       file = "Elyada_tile_new.png",
       width = 4, height = 3, dpi=500, bg="transparent")

#----------------------------------------------------








DEG.table_HvsN[c("CDK1","CDC45","BIRC5","CCNA2","CENPA","TOP2A","FOXM1","CCNB1","TPX2"),]
DEG.table_HvsN[c("CDKN1A","CDKN2A"),]
DEG.table_HvsNH[c("CDK1","CDC45","BIRC5","CCNA2","CENPA","TOP2A","FOXM1","CCNB1","TPX2"),]
DEG.table_HvsNH[c("CDKN1A","CDKN2A"),]
DEG.table_HNvsN[c("CDK1","CDC45","BIRC5","CCNA2","CENPA","TOP2A","FOXM1","CCNB1","TPX2"),]
DEG.table_HNvsN[c("CDKN1A","CDKN2A"),]
DEG.table_HvsHN[c("CDK1","CDC45","BIRC5","CCNA2","CENPA","TOP2A","FOXM1","CCNB1","TPX2"),]
DEG.table_HvsHN[c("CDKN1A","CDKN2A"),]
DEG.table_NHvsN[c("CDK1","CDC45","BIRC5","CCNA2","CENPA","TOP2A","FOXM1","CCNB1","TPX2"),]
DEG.table_NHvsN[c("CDKN1A","CDKN2A"),]




