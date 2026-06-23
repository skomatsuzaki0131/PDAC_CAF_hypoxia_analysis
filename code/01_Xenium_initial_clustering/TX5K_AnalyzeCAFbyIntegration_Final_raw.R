# Xenium initial clustering and cell type annotation
# This script contains the final unrefactored code used for preprocessing,
# QC filtering, integration, clustering, UMAP visualization, and cell type annotation
# of the integrated Xenium spatial transcriptomic dataset.
#
# Note: This script will be refactored before publication.


library(dplyr)
library(ggplot2)
library(ggraph)
library(RColorBrewer)
library(gplots)
library(ggrepel)
library(scales)
library(tibble)
library(stringr)
library(magrittr)
library(ggpubr)
library(patchwork)
library(tidyverse)
library(cowplot)
library(patchwork)
library(Seurat)
library(presto)
library(beepr) #beep(expr=NULL, sound=3)
library(fs)

#Define the celltype marker genes
Vec.AllGenes = read.csv(file='/Volumes/PortableSSD/[4]myR/[1]dataset/Genes_Xenium5k.csv', header=F)[ , 1]
DF.AllGenes_addEntrez = read.csv(file='/Volumes/PortableSSD/[4]myR/[1]dataset/Genes_Xenium5k_addEntrez.csv', header=T)
DF.Markers = read.csv(file="/Volumes/PortableSSD/[4]myR/[1]dataset/MyGeneSet.csv", header=T) %>% 
  dplyr::filter(!Gene %in% c("CD3D", "COL1A2"))
Vec.MarkerGenes_CellType = subset(DF.Markers, subset=GeneType=="CellType" &
                                    Panel5k=="Included" &
                                    Panel5kQlt!="Duplicated"&
                                    Panel5kQlt!="MayBePoor")$Gene
Vec.CAFmarkers = subset(DF.Markers, subset=GeneType=="CAFmarker")$Gene
Vec.BuffaOrig = intersect(subset(DF.Markers, subset=Classification1=="Buffa")$Gene,
                          Vec.AllGenes)
Vec.WinterOrig = intersect(subset(DF.Markers, subset=Classification1=="WinterCore")$Gene,
                           Vec.AllGenes)
Vec.MoffittActivated = intersect(subset(DF.Markers, subset=Classification1=="Activated")$Gene,
                                 Vec.AllGenes)
Vec.MoffittNormal = intersect(subset(DF.Markers, subset=Classification1=="Normal")$Gene,
                              Vec.AllGenes)
Vec.Classical = subset(DF.Markers, Classification1=="Classical")$Gene
Vec.Basallike = subset(DF.Markers, Classification1=="Basal_like")$Gene
Vec.Mix300.NormoCAF = subset(DF.Markers, subset=GeneType=="Mix3"&Classification1=="NormoCAF")$Gene
Vec.Mix300.HypoCAF = subset(DF.Markers, subset=GeneType=="Mix3"&Classification1=="HypoCAF")$Gene
Vec.colorcode.Final = 
  c("Acinar"="#ff00ff", "Ductal_like_acinar"="#ddbcff", "Ductal"="#ddbcff",
    "Normal_ductal"="#FFFFCC","PanIN1"="#FFFF00","PanIN2"="#FD8D3C","PDAC"="#ff0000",
    "Islet"="#BF812D",
    "CAF"="#00ff00","Mural"="#990000","Endothelial"="#990000",
    "Lymphoid"="#00bfff", "Lymph_T"="#00bfff","Lymph_B"="#00bfff","Plasma"="#00bfff",
    "Myeloid"="#8300FF", "Mast"="#8300FF", "Nerve"="#0000ff", 
    "Proliferating_Immune"="#999999", "Unclassified"="#999999")
TXnumInteg = c("02","19","11","16","01","15"#,
               #"04"
               #"14"
)
NumOfSamples = 6
IntegArgo = "Harmony"  
Mag=1
nFeatRNA=c(100, 900)
nCountRNA=c(100, 1800)
QCInfo = paste0("Countable_mag",Mag,"_nFeatRNA:",paste(nFeatRNA,collapse="~"),"_nCountRNA:",paste(nCountRNA,collapse="~"))
QCInfo.FileName = str_replace_all(QCInfo, pattern=":", replace="")
DirInteg = 
  paste0('/Volumes/Extreme SSD/Analysis/Data/IntegAnalysis/',
         "[",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
         "_nFeatRNA:",paste(nFeatRNA,collapse="~"),"_nCountRNA:",paste(nCountRNA,collapse="~"),"/")
Dim1=20
Res1=1.0
Dim2=30
Res2=0.5
Num.GridLen = 200
Num.MinCells = 10
prolifCAFclust="8"
n.neighbors=40
k.niche = 10
pc_use = 15
Unclassified="Exclude"
Labels="seurat_clusters"
method="SpearmanAverage"
DirRNAseq = "/Volumes/PortableSSD/[4]myR/[1]dataset/RNAseq/"
DirGO = "/Volumes/PortableSSD/[4]myR/[1]dataset/RNAseq/GO/"
DirOutputGO = "/Volumes/PortableSSD/[3]Graduate school/[3]RNA-seq/analysis/GeneOntonogy/"
DirHallmark = '/Volumes/PortableSSD/[4]myR/[1]dataset/RNAseq/Hallmark/'
DirOutputMSigDB = "/Volumes/PortableSSD/[3]Graduate school/[3]RNA-seq/analysis/MSigDB/"
SetName="Set1_3pair"
DataSetName = "rmCAF6"
n_HVGs = 1000
List.Annot = 
  list(Acinar=c(0,13,24),    Ductal_like_acinar=c(33,1,12,17,30),
       Normal_ductal=23,     PanIN=c(10,27),
       PDAC=c(20),           Islet=c(14),
       CAF=c(9,6,7),         Mural=c(11),
       Endothelial=c(34,3),  Lymph_T=c(2,16),
       Lymph_B=c(15),        Plasma=c(18,37),
       Myeloid=c(5,4,31),    Mast=c(26),
       Nerve=c(29),          Unclassified=c(8, 19, 21, 22, 25, 28, 32, 35, 36) )
Vec.Annot <- unlist(
  lapply(names(List.Annot), function(ct){
    setNames(rep(ct, length(List.Annot[[ct]])),
             List.Annot[[ct]])
  }))
Vec.OrderedCnumToColor =
  c("0"="#ff00ff", "13"="#ff00ff", "24"="#ff00ff", "33"="#ff00ff",
    "12"="#ddbcff", "1"="#ddbcff", "17"="#ddbcff", "30"="#ddbcff", 
    "23"="#ddbcff", "10"="#ddbcff", "27"="#ddbcff", "20"="#ddbcff",
    "14"="#ffff00",
    "9"="#00ff00", "6"="#00ff00", "7"="#00ff00",
    "11"="#990000",
    "34"="#ff0000", "3"="#ff0000",
    "2"="#00bfff", "16"="#00bfff", "15"="#00bfff", "18"="#00bfff", "37"="#00bfff",
    "5"="#ff8800", "4"="#ff8800", "31"="#ff8800", "26"="#ff8800",
    "29"="#0000ff",
    "25"="#ABABAB", "22"="#ABABAB", "28"="#ABABAB",
    "36"="#ABABAB", "35"="#ABABAB", "32"="#ABABAB", "8"="#ABABAB", "19"="#ABABAB", "21"="#ABABAB")
Vec.XeniumID = c("01"="Xenium_01", "02"="Xenium_02",
                 "11"="Xenium_03", "15"="Xenium_04",
                 "16"="Xenium_05", "19"="Xenium_06")
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
Vec.SubclustOrder <- c(4,8,0,1,2,5,3,7,6)
cols2 = 
  c("0"="#4C78A8",  "1"="#54A24B",   "2"="#B279A2",  
    "3"="#9D755D",  "4"="#72B7B2",   "5"="#E6AF2E", 
    "6"="#F58518",  "7"="#9C9EDE",   "8"="#C00000")
Vec.Cols.Diverging <- c("#A6D3C5", "#F7F7F7", "#C05A9E")

#################################################
##      Vitro
##  0. Outgrowth                   ####
library(ggsignif)
library(exactRankTests)
#Data.OG_0 = read.csv('/Volumes/PortableSSD/[4]myR/[1]dataset/CAFlist.csv', header=TRUE)
#Data.OG_1 = 
#  Data.OG_0 %>% t() %>% as.data.frame() %>% 
#  dplyr::mutate(
#    O2 = str_replace(O2, pattern="N", replacement="21pct"),
#    O2 = str_replace(O2, pattern="H", replacement="1pct"),
#    Outgrow = as.integer(Outgrow)) %>% 
#  set_rownames(NULL)
#Data.OG_2 = 
#  Data.OG_1 %>% 
#  pivot_wider(names_from = O2, 
#              values_from = Outgrow,
#              names_glue = "Days_under_{O2}") %>% 
#  dplyr::mutate(delta_days = Days_under_1pct - Days_under_21pct) %>% 
#  dplyr::select(c(Patient, Days_under_21pct, Days_under_1pct, delta_days, Site, NAT)) %>% 
#  rename(CAF_pair_No = Patient)
#write.csv(Data.OG_2, file="SuppTable_OutGrowth_A.csv", row.names=F)
Data.OG_2 <- read.csv(file="SuppTable_OutGrowth_A.csv")
Data.OG_Comp.O2 = 
  Data.OG_2 %>% 
  dplyr::mutate(Days_under_21pct = Days_under_21pct,
                Days_under_1pct = Days_under_1pct) %>% 
  summarise(
    median_21 = median(Days_under_21pct),
    Q1_21 = quantile(Days_under_21pct, 0.25),
    Q3_21 = quantile(Days_under_21pct, 0.75),
    median_1 = median(Days_under_1pct),
    Q1_1 = quantile(Days_under_1pct, 0.25),
    Q3_1 = quantile(Days_under_1pct, 0.75)) %>% 
  dplyr::mutate(
    label_21 = paste0(median_21," (",Q1_21,"-",Q3_21,")"),
    label_1 = paste0(median_1," (",Q1_1,"-",Q3_1,")"),
    p_value = wilcox.exact(x=Data.OG_2$Days_under_21pct,
                           y=Data.OG_2$Days_under_1pct,
                           paired=T)$p.value) #p=0.001953)
write.csv(Data.OG_Comp.O2, file="SuppTable_OutGrowth_B_CompO2.csv", row.names=F)

Data.OG_1 =
  Data.OG_2 %>%
  pivot_longer(
    cols = c(Days_under_21pct, Days_under_1pct),
    names_to = "O2",
    values_to = "Outgrow") %>%
  mutate(
    O2 = recode(O2,
                Days_under_21pct = "21pct",
                Days_under_1pct  = "1pct")) %>%
  rename(Patient = CAF_pair_No) %>%
  select(Patient, O2, Outgrow, Site, NAT)
Data.OG_Comp.Site = 
  Data.OG_1 %>% 
  group_by(Site) %>% 
  summarise(
    n = n(),
    median = median(Outgrow),
    Q1 = quantile(Outgrow, 0.25),
    Q3 = quantile(Outgrow, 0.75)) %>% 
  dplyr::mutate(
    label_1 = paste0(Site," (n=",n,")"),
    label_2 = paste0(median," (",Q1,"-",Q3,")"),
    p_value = kruskal.test(Outgrow ~ Site, data=Data.OG_1)$p.value ,
    Site = factor(Site, levels=c("Pt","Pb","Ph"))) %>% 
  dplyr::arrange(Site)
write.csv(Data.OG_Comp.Site, file="SuppTable_OutGrowth_B_CompSite.csv", row.names=F)
Data.OG_Comp.NAT = 
  Data.OG_1 %>% 
  group_by(NAT) %>% 
  summarise(
    n = n(),
    median = median(Outgrow),
    Q1 = quantile(Outgrow, 0.25),
    Q3 = quantile(Outgrow, 0.75)) %>% 
  dplyr::mutate(
    label_1 = paste0(NAT," (n=",n,")"),
    label_2 = paste0(median," (",Q1,"-",Q3,")"),
    p_value = kruskal.test(Outgrow ~ NAT, data=Data.OG_1)$p.value ,
    NAT = factor(NAT, levels=c("None","GS","GnP"))) %>% 
  dplyr::arrange(NAT)
write.csv(Data.OG_Comp.NAT, file="SuppTable_OutGrowth_B_Comp.NAT.csv", row.names=F)
Data.OG_3 = 
  Data.OG_1 %>% 
  dplyr::mutate(Oxygen = str_replace(O2, pattern="21pct", replacement="N_CAF"),
                Oxygen = str_replace(Oxygen, pattern="1pct", replacement="H_CAF"), 
                Oxygen = factor(Oxygen, levels = c("N_CAF", "H_CAF")),
                Site = factor(Site, levels = c("Ph","Pb","Pt")),
                NAT = factor(NAT, levels = c("None", "GS", "GnP")))
ScaleYcontinu = 
  scale_y_continuous(
    limits=c(0, max(Data.OG_1$Outgrow)+4),
    expand=expansion(mult=c(0.0)))
TextSize=10
CommonLabs = labs(x=NULL, y="Days to outgrowth", fill=NULL)
CommonTheme = 
  theme(
    legend.position="none",
    axis.text.x = element_text(face = "bold", color = "black", size=TextSize*3.0),
    axis.text.y = element_text(face = "bold", color = "black", size=TextSize*3.0),
    axis.title = element_text(face = "bold", color = "black", size=TextSize*3.0),
    legend.text = element_text(face = "bold", color = "black"),
    legend.title = element_text(face = "bold", color = "black"),
    plot.background = element_rect(fill="transparent", color=NA),
    panel.background = element_rect(fill="white", color=NA),
    panel.grid = element_blank(),
    panel.border = element_rect(fill=NA, color="black"))
Plot.CAFoutgrowth.O2 = 
  ggplot(Data.OG_3, aes(x=Oxygen, y=Outgrow, fill=Oxygen)) +
  geom_boxplot(width = 0.45, fill = "white", staplewidth = 0.50, linewidth = 0.8, outliers = F) +
  geom_line(aes(group = Patient), color = "gray60", linewidth=0.4) +
  geom_point(shape=21, 
             size=TextSize*0.4, 
             position=position_jitter(width=0.10, height=0)) +
  geom_signif(comparisons = list(c("H_CAF","N_CAF")),
              map_signif_level = FALSE,
              y_position = 16, 
              annotation=paste0("**"),
              textsize=TextSize*1.5, fontface="bold", vjust=0.3, size=1.5) +
  labs(subtitle=NULL, caption=NULL) + CommonLabs + #paste0("p=",p.value_Oxygen_1," (Wilcoxon signed-rank test)"),
  scale_x_discrete(labels=c(H_CAF="Hypo-\nCAFs", N_CAF="Normo-\nCAFs")) +
  coord_fixed(ratio = 1 / 7.0) +
  scale_fill_manual(
    values=c(H_CAF="#91BFFA", N_CAF="#FEA0A0"))　+
  ScaleYcontinu + CommonTheme + 
  theme(axis.text.x = element_text(color="transparent"))
ggsave(plot=Plot.CAFoutgrowth.O2,
       file="[3]Figures/[Figures]CAFoutgrowth_Plot4.png",
       width=7, height=7, dpi=400, bg="transparent")
Plot.CAFoutgrowth.Site=
  ggplot(Data.OG_3, aes(y=Outgrow, x=Site)) +
  geom_boxplot(width = 0.6, fill = "white", staplewidth = 0.6, linewidth = 0.8, outlier.shape = NA) +
  #geom_dotplot(aes(fill=Site),
  #             binaxis="y", stackdir="center",
  #             binwidth=0.6,
  #             alpha=0.6) +
  geom_point(aes(fill=Site), shape=21, 
             size=TextSize*0.4, alpha=1,
             position=position_jitter(width=0.30, height=0, seed=1236)) +
  geom_text(label="Kruskal–Wallis, p = 0.722", 
            x=3.5, y=18, hjust=1, size=TextSize*0.8, fontface="bold") +
  labs(subtitle="Tumor Location") + CommonLabs + 
  scale_fill_manual(values=c(Ph = "#424F9E",
                             Pb = "#A2BDD4",
                             Pt = "#9BC8B6")) + 
  coord_fixed(ratio = 3 / 15) +
  ScaleYcontinu + CommonTheme
Plot.CAFoutgrowth.NAT =
  ggplot(Data.OG_3, aes(y=Outgrow, x=NAT)) +
  geom_text(label="Kruskal–Wallis, p = 0.269", 
            x=Inf, y=18, hjust=1.05, size=TextSize*0.8, fontface="bold") +
  geom_boxplot(width = 0.6, fill = "white", staplewidth = 0.6, linewidth = 0.8, outlier.shape = NA) +
  #geom_dotplot(aes(fill=NAT), binaxis="y", stackdir="center",
  #             binwidth=0.6,
  #             alpha=0.6) +
  geom_point(aes(fill=NAT), shape=21, 
             size=TextSize*0.4, alpha=1,
             position=position_jitter(width=0.30, height=0, seed=1235)) +
  labs(subtitle="NAT") + CommonLabs +
  scale_fill_manual(values=c(None = "#E5E5E5",
                             GS = "#E0C3A9",
                             GnP = "#EEA93C")) +
  coord_fixed(ratio = 3 / 15) +
  ScaleYcontinu + CommonTheme
Plot.Supp = 
  (Plot.CAFoutgrowth.Site | Plot.CAFoutgrowth.NAT) &
  theme(plot.background = element_rect(fill="transparent", color=NA),
        panel.background = element_rect(fill="white", color=NA))
ggsave(plot=Plot.Supp,
       file="CAFoutgrowth_Supp2.png",
       width=16, height=7, dpi=400)
##  1. Proliferation assay         ####
data_all=read.csv(file="/Volumes/PortableSSD/[4]myR/[1]dataset/proliferation_all.csv", header=T, sep=",", stringsAsFactors=F) %>% 
  rowwise %>% 
  mutate(log2nl = log(nl, 2)) %>% ungroup() %>%
  mutate(Day = rep(c("Day1","Day1","Day1","Day2","Day2","Day2",
                     "Day3","Day3","Day3","Day4","Day4","Day4",
                     "Day5","Day5","Day5","Day6","Day6","Day6",
                     "Day7","Day7","Day7"), times=12),
         oxygen = factor(oxygen, levels=c("hn","h","n","nh")))
Func.ProPlot_2 = function(SAMPLE){
  data_all_2 = 
    dplyr::filter(data_all, patient==SAMPLE)
  PlotData = 
  ggplot(data_all_2, aes(x=time, y=log2nl, color=oxygen)) +
    stat_summary(geom="line", fun="mean",
                 alpha=0.9, 
                 linewidth=1.3) +
    stat_summary(geom="errorbar", 
                 fun.data=function(x) {
                   m <- mean(x, na.rm = TRUE)
                   se <- sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
                   data.frame(y = m, ymin = m - se, ymax = m + se)
                 },
                 alpha=0.55,
                 linewidth=0.65,
                 width=0.35) + 
    stat_summary(geom="point",
                 fun="mean",
                 size=3.0,
                 alpha=0.8) +
    labs(x="Time after seeding (h)", y=expression(bold(""*log[2]*"(NL)") ), color=NULL,
         subtitle=SAMPLE) +
    guides(color = guide_legend(
      override.aes = list(
        linewidth = 1.2,
        shape = 16,
        size = 3))) +
    scale_x_continuous(limits=c(0, 218),
                       breaks=c(24, 72, 120, 168),
                       minor_breaks=NULL, 
                       expand=expansion(mult = c(0,0))) +
    scale_y_continuous(limits=c(0, 4.5),
                       breaks=c(0, 1, 2, 3, 4),
                       oob=scales::squish,
                       expand=expansion(mult=c(0,0))) +
    scale_color_manual(values=c(h = "#607DF0",
                                hn = "#7FBFFF",
                                nh = "#FFA3A3", 
                                n = "#FF7F7F"),
                       labels = c(h="H-CAF", n="N-CAF",
                                  nh="N-CAF\n  under 1% O2",
                                  hn="H-CAF\n under  21% O2")) +
    coord_fixed(ratio=40) +
    theme(text = element_text(face="bold", color="black", size=15),
          axis.text.x = element_text(face="bold", color="black", size=15),
          axis.text.y = element_text(face="bold", color="black", size=15),
          legend.key = element_blank(),
          legend.background = element_rect(fill="transparent", color=NA),
          plot.background = element_rect(fill="transparent", color=NA),
          panel.background = element_rect(fill="transparent", color=NA),
          panel.grid = element_blank(), 
          axis.line = element_line(color="black")) 
}
NL = theme(legend.position="none")
Plot = 
  ggarrange(Func.ProPlot_2(SAMPLE="23-1122") + NL,
            Func.ProPlot_2(SAMPLE="24-0110") + NL,
            Func.ProPlot_2(SAMPLE="24-0207") + NL,
            ggpubr::get_legend(
              Func.ProPlot_2(SAMPLE="23-1122") +
                theme(legend.position = c(0, 0.5), 
                      legend.justification = c(0, 0.5))),
            ncol=4, widths=c(3,3,3,2), nrow=1)
ggsave(plot=Plot,
       file="ProlifAssay.png", #paste0(Dir.OutProlife, "[Figure][Proliferation]_4group_3.png"),
       width=15, height=5, dpi=400, bg="transparent")

##  2. KEGG                        ####
library(clusterProfiler)
Func.EnrichKEGG = function(CompSetting){
  Comp = case_when(
    CompSetting==1 ~ "1Simple_HvsN",
    CompSetting==2 ~ "2in1pct_HvsNH",
    CompSetting==3 ~ "3in21pct_HNvsN")
  DEGall = 
    read.csv(file=paste0(DirRNAseq,"[DEGsTable_addENTREZ]_",SetName,"_",Comp,".csv"))
  DEGsig = 
    DEGall %>% 
    dplyr::filter((logFC>=1 | logFC<=(-1)) & FDR<0.05)
  DEGsig.UpInHCAF = dplyr::filter(DEGsig, logFC>0)
  DEGsig.UpInNCAF = dplyr::filter(DEGsig, logFC<0)
  gene_use.UpInHCAF = DEGsig.UpInHCAF$Entrez %>% unique() %>% na.omit()
  gene_use.UpInNCAF = DEGsig.UpInNCAF$Entrez %>% unique() %>% na.omit()
  bg = DEGall$Entrez %>% as.character() %>% unique() %>% na.omit()
  ResKEGG.UpInHCAF = 
    enrichKEGG(gene = gene_use.UpInHCAF, 
               universe = bg,
               organism="hsa", 
               keyType = "ncbi-geneid",
               pAdjustMethod = "BH")
  ResKEGG.UpInNCAF = 
    enrichKEGG(gene=gene_use.UpInNCAF, 
               universe = bg,
               organism="hsa", 
               keyType = "ncbi-geneid",
               pAdjustMethod = "BH")
  DF.Res = 
    rbind(
      dplyr::mutate(ResKEGG.UpInHCAF@result, Setting=Comp, Direction="UP"),
      dplyr::mutate(ResKEGG.UpInNCAF@result, Setting=Comp, Direction="DOWN"))
  return(DF.Res)
}
DF.ResKEGG = 
  bind_rows(
    Func.EnrichKEGG(1),
    Func.EnrichKEGG(2),
    Func.EnrichKEGG(3))
Vec.KEGG.Core10 <- c(
  "Cell cycle",
  "DNA replication",
  "Homologous recombination",
  "Fanconi anemia pathway",
  "HIF-1 signaling pathway",
  "Glycolysis / Gluconeogenesis",
  "ECM-receptor interaction",
  "PI3K-Akt signaling pathway",
  "Integrin signaling",
  "Cytokine-cytokine receptor interaction")
DF.ResKEGG.sig = 
  dplyr::filter(DF.ResKEGG, p.adjust<0.05 & Description%in%Vec.KEGG.Core10) %>% 
  dplyr::mutate(
    GeneRatio_num = 
      sapply(GeneRatio, function(x){
        a = as.numeric(strsplit(x, "/")[[1]][1])
        b = as.numeric(strsplit(x, "/")[[1]][2])
        a / b
      }),
    GeneRatio_num_signed = 
      GeneRatio_num*ifelse(Direction=="UP", 1, -1),
    Setting = factor(Setting, levels=c("1Simple_HvsN","3in21pct_HNvsN","2in1pct_HvsNH")),
    Description = factor(Description, levels=Vec.KEGG.Core10),
    SigLab = case_when(p.adjust<0.001 ~ "***",
                       p.adjust<0.01 ~ "**",
                       p.adjust<0.05 ~ "*",
                       TRUE ~ "n.s.")
  )
MaxAbsGR = max(DF.ResKEGG.sig$GeneRatio_num)
TextSize=10
PlotKEGGtile = 
  ggplot(DF.ResKEGG.sig, 
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
    breaks=c(-0.2, 0, 0.2),
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
       file="KEGG_tile.png",
       width=5.8, height=3, dpi=500, bg="transparent")
saveRDS(PlotKEGGtile, file="Fig2C.rds")
## Pathview 
  # read significant DEGs data
  DEGsig.1_HvsN = read.csv(paste0(DirRNAseq,"[DEGsTable_addENTREZ_Signif]_",SetName,"_1simple_HvsN.csv"), sep=",", header = TRUE, row.names = 1)
  DEGsig.2_HvsNH = read.csv(paste0(DirRNAseq,"[DEGsTable_addENTREZ_Signif]_",SetName,"_2in1pct_HvsNH.csv"), sep=",", header = TRUE, row.names = 1)
  DEGsig.3_HNvsN = read.csv(paste0(DirRNAseq,"[DEGsTable_addENTREZ_Signif]_",SetName,"_3in21pct_HNvsN.csv"), sep=",", header = TRUE, row.names = 1)
  # generate log2FC vector named with EntrezID
  Vec.log2FC.1_HvsN = setNames(DEGsig.1_HvsN$logFC, DEGsig.1_HvsN$EntrezID)
  Vec.log2FC.2_HvsNH = setNames(DEGsig.2_HvsNH$logFC, DEGsig.2_HvsNH$EntrezID)
  Vec.log2FC.3_HNvsN = setNames(DEGsig.3_HNvsN$logFC, DEGsig.3_HNvsN$EntrezID)
  # DEG1~3 および DEG1~5をまとめて表示する
  library(KEGGREST)
  pathway <- keggGet("hsa04110")[[1]]
  genes <- pathway$GENE
  Vec.EntrezID.CCgenes <- genes[seq(1, length(genes), 2)]
  # CellCyclePathwayに属する遺伝子のlogFCをmatrixに抽出
  v1=Vec.log2FC.1_HvsN[Vec.EntrezID.CCgenes]
  v2=Vec.log2FC.2_HvsNH[Vec.EntrezID.CCgenes]
  v3=Vec.log2FC.3_HNvsN[Vec.EntrezID.CCgenes]
  m_three = cbind(v1,v3,v2)
  colnames(m_three) = c("Original", "Under21", "Under1")
  rownames(m_three) = DF.CellCycleGenes$Entrez_ID
  # generate pathview map Cell cycle
  library(pathview)
  pathview(gene.data = m_three, pathway.id = "hsa04110", species = "hsa",
           limit = list(gene=3, cpd=1),
           low = list(gene = "red", cpd = "black"),
           mid = list(gene = "white", cpd = "black"), 
           high = list(gene = "blue", cpd = "black"),
           same.layer = FALSE,
           kegg.native = TRUE,
           plot.col.key = F,
           key.pos = "bottomright")
  file.rename(from="hsa04110.pathview.multi.png",
              to=paste0("[Figure][KEGG_PathwayMapWithSignifExp_3comp]_CellCycle_",SetName,".png"))
  filesstrings::file.move(files=paste0("[Figure][KEGG_PathwayMapWithSignifExp_3comp]_CellCycle_",SetName,".png"),
                          destinations="/Volumes/PortableSSD/[3]Graduate school/[3]RNA-seq",
                          overwrite=TRUE)
  pathview(gene.data = m_five, pathway.id = "hsa04110", species = "hsa",
           limit = list(gene=2, cpd=1),
           low = list(gene = "blue", cpd = "black"),
           mid = list(gene = "gray", cpd = "black"), 
           high = list(gene = "red", cpd = "black"),
           same.layer = FALSE,
           kegg.native = TRUE)


##  3. GO                          ####
  DirOutputGO = "/Volumes/PortableSSD/[3]Graduate school/[3]RNA-seq/analysis/GeneOntonogy/"
  List.gseaRes = readRDS(paste0(DirOutputGO,"[RDSdata][gseGO][ListOfEnrichResult]_Set1_3pair_Simplified.rds"))
  Func.DataWrangl = function(DF){
    DF %>% 
      dplyr::mutate(Direction=ifelse(NES>0,"UP","DOWN"),
                    Direction=factor(Direction, levels=c("UP","DOWN"))) %>% 
      group_by(Direction) %>% 
      dplyr::arrange(p.adjust, .by_group=TRUE)
  }
  DF.gseaRes.1_HvsN_1 = cbind(setting="isolation", List.gseaRes[[1]]@result) %>% Func.DataWrangl()
  DF.gseaRes.2_HvsNH_1 = cbind(setting="under1", List.gseaRes[[2]]@result) %>% Func.DataWrangl()
  DF.gseaRes.3_HNvsN_1 = cbind(setting="under21", List.gseaRes[[3]]@result) %>% Func.DataWrangl()
  DF.gseaRes.bind_0 = 
    bind_rows(DF.gseaRes.1_HvsN_1,
              DF.gseaRes.3_HNvsN_1,
              DF.gseaRes.2_HvsNH_1) %>% 
    dplyr::mutate(minuslog10adjp = (-1)*log10(p.adjust))
  Vec.TermsToDisplay_0 = 
    c(intersect(subset(DF.gseaRes.1_HvsN_1, ONTOLOGY=="BP"&NES<0)$Description, 
                subset(DF.gseaRes.2_HvsNH_1, ONTOLOGY=="BP"&NES<0)$Description) %>% 
      intersect(subset(DF.gseaRes.3_HNvsN_1, ONTOLOGY=="BP"&NES<0)$Description),
      intersect(subset(DF.gseaRes.1_HvsN_1, ONTOLOGY=="BP"&NES>0)$Description, 
                subset(DF.gseaRes.2_HvsNH_1, ONTOLOGY=="BP"&NES>0)$Description) %>% 
      intersect(subset(DF.gseaRes.3_HNvsN_1, ONTOLOGY=="BP"&NES>0)$Description))
  Vec.TermsToDisplay_Final = c(
    "chromosome segregation",
    "nuclear division",
    "regulation of chromosome separation",
    "inflammatory response",
    "collagen metabolic process",
    "lipid catabolic process",
    "inorganic ion transmembrane transport")
  DF.gseaRes.bind_1 = DF.gseaRes.bind_0 %>% 
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
         file=paste0("[Figure][gseGO(GSEA)][DotPlot_Final].png"),
         width=6.5, height=3.5, dpi=500, bg="transparent")
  saveRDS(Plot.DotGO + 
            labs(size="-log10\n(Adj.P)") +
            scale_size_continuous(range=c(1, 8.5),
                                  breaks=c(2,10,30),
                                  guide=guide_legend(
                                    order=1,
                                    label.position="left")), 
          "Fig2D.rds")
  DF.gseaRes.bind_2 = DF.gseaRes.bind_0 %>% 
    #dplyr::select(c(-minuslog10adjp, -rank, -enrichmentScore)) %>% 
    dplyr::select(c(setting, Direction, ONTOLOGY, ID, Description, setSize, NES, p.adjust, pvalue, qvalue,
                    leading_edge, core_enrichment)) %>%
    dplyr::rename(Enrichment_in_HypoCAF = Direction)
  write.csv(DF.gseaRes.bind_2,
            file="SupplementalTable_bulk_GSEA(GO).csv",
            row.names=F)

##  4. GSEA 3 settings             ####
  library(ggdendro)
  List.gseaRes = readRDS(paste0(DirOutputMSigDB,"[ListRDS]_gseaHallmark_Results_",SetName,".rds"))
  gseaRes_1_HvsN = List.gseaRes[[1]]@result %>% dplyr::mutate(Setting="Original")
  gseaRes_2_HvsNH = List.gseaRes[[2]]@result %>% dplyr::mutate(Setting="Under1")
  gseaRes_3_HNvsN = List.gseaRes[[3]]@result %>% dplyr::mutate(Setting="Under21")
  #gseaRes_4_HvsHN = List.gseaRes[[4]]
  #gseaRes_5_NHvsN = List.gseaRes[[5]]
  DF.gseaRes_0 = 
    dplyr::bind_rows(gseaRes_1_HvsN, gseaRes_3_HNvsN, gseaRes_2_HvsNH) %>% 
    dplyr::select(c(Setting, Description, NES, p.adjust)) %>% 
    dplyr::mutate(Setting = factor(Setting, levels=c("Original", "Under21", "Under1")),
                  Signif = case_when(p.adjust<0.05 ~ "*",
                                     TRUE ~ "")) %>% 
    dplyr::filter(!(Description %in% c("iCAF_signature",
                                       "KEGG_Cell_cycle",
                                       "Buffa_Hypoxia",
                                       "Winter_Hypoxia")))
  Range.NES = range(DF.gseaRes_0$NES)
  DF.gseaRes_1 = DF.gseaRes_0 %>% 
    pivot_wider(id_cols=Setting,
                names_from = Description, 
                values_from = NES) %>% 
    column_to_rownames(var="Setting") %>% as.matrix()
  # hierarchical clustering : euclidean, ward.D2
  d_hallmark = dist(t(DF.gseaRes_1), method="euclidean")
  Hclust_hallmark <- hclust(d_hallmark, method = "ward.D2")
  MethodLabel = "Hierarchical clustering was performed using Ward’s method (ward.D2) with Euclidean distances."
  Vec.HallmarkOrder = colnames(DF.gseaRes_1)[Hclust_hallmark$order]
  DF.gseaRes_0$Description = factor(DF.gseaRes_0$Description, levels=Vec.HallmarkOrder)
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
    ggplot(DF.gseaRes_0, aes(x=Setting, y=Description, fill=NES)) +
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
                       setNames(clean_gs_label_ForFig(DF.gseaRes_0$Description),
                                DF.gseaRes_0$Description)) +
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
          panel.background = element_rect(fill="transpanrent", color=NA),
          panel.border = element_rect(fill="transpanrent", color="black"))
  Plot.joined = 
    p_dend + p_heatmap +
    plot_layout(widths=c(0.4, 0.2)) &
    theme(
      plot.background = element_rect(fill="transparent", color=NA),
      panel.background = element_rect(fill="transparent", color=NA))
  ggsave(plot = Plot.joined,
         width=5.5, height=9,
         file="InVitroGSEA_sup4.png", dpi=500, bg="transparent")
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
    DF.gseaRes_2 = DF.gseaRes_0 %>% 
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
                  direction="vertical",
                  title.position="top",
                  label.position="left",
                  reverse=FALSE,
                  frame.colour="black",
                  ticks.colour="black")), 
            "Fig2Eleft_rightYaxis.rds")
    ggsave(plot = p_main +
             coord_fixed(ratio=0.8),
           width=6.2, height=5.1,
           file="InVitroGSEA_main2.png", dpi=500, bg="transparent")
##  5. GSEA before/after switching ####
List.gseaRes = readRDS(paste0(DirOutputMSigDB,"[ListRDS]_gseaHallmark_Results_",SetName,".rds"))
gseaRes_4_HvsHN = List.gseaRes[[4]]@result %>% 
  dplyr::mutate(Culture="HypoCAF")　%>% 
  dplyr::filter(! Description %in% c("iCAF_signature","Buffa_Hypoxia",
                                     "Winter_Hypoxia","KEGG_Cell_cycle"))
gseaRes_5_NHvsN = List.gseaRes[[5]]@result %>% 
  dplyr::mutate(Culture="NormoCAF")　%>% 
  dplyr::filter(! Description %in% c("iCAF_signature","Buffa_Hypoxia",
                                     "Winter_Hypoxia","KEGG_Cell_cycle"))
gseaRes_BAscitch_HypoCAF = 
  rbind(slice_max(gseaRes_4_HvsHN, n=5, order_by=NES),
        slice_min(gseaRes_4_HvsHN, n=5, order_by=NES)) %>% 
  dplyr::arrange(-NES) %>% 
  dplyr::mutate(Description = factor(Description, levels=unique(Description)),
                minusLog10adjP = (-1)*log10(p.adjust))
gseaRes_BAscitch_NormoCAF = 
  rbind(slice_max(gseaRes_5_NHvsN, n=5, order_by=NES),
        slice_min(gseaRes_5_NHvsN, n=5, order_by=NES)) %>% 
  dplyr::arrange(-NES) %>% 
  dplyr::mutate(Description = factor(Description, levels=unique(Description)),
                minusLog10adjP = (-1)*log10(p.adjust))
Range.NES.BAswitch = range(c(gseaRes_BAscitch_HypoCAF$NES,
                             gseaRes_BAscitch_NormoCAF$NES))
Range.minuslog10adjP.BAswitch = 
  range(c(gseaRes_BAscitch_HypoCAF$minusLog10adjP,
          gseaRes_BAscitch_NormoCAF$minusLog10adjP))
Common.Scale.Size =
  scale_size_continuous(
    name = "-log10(Adj.P)",
    breaks = c(2, 10, 30),
    limits = c(0, max(Range.minuslog10adjP.BAswitch)),
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
#Erase.lgd = theme(legend.position="none")
Plot_NormoCAF =
  ggplot(gseaRes_BAscitch_NormoCAF, aes(x=NES, y=Description)) +
  geom_segment(aes(x=0, xend=NES, 
                   y=Description, yend=Description)) +
  geom_point(aes(size=minusLog10adjP), shape=21, fill="#f03b20") +
  scale_y_discrete(
    labels=setNames(clean_gs_label_ForFig(gseaRes_BAscitch_NormoCAF$Description) %>% 
                      str_replace(pattern="Epithelial mesenchymal transition",
                                  replacement="EMT"),
                    gseaRes_BAscitch_NormoCAF$Description),
    expand = expansion(mult = c(0.1, 0.1)),
    position = "right") +
  Common.vline + 
  scale_x_continuous(
    position="top",
    limits=c(-max(abs(Range.NES.BAswitch)),max(abs(Range.NES.BAswitch))),
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
  ggplot(gseaRes_BAscitch_HypoCAF, aes(x=NES, y=Description)) +
  geom_segment(aes(x=0, xend=NES, 
                   y=Description, yend=Description)) +
  geom_point(aes(size=minusLog10adjP), shape=21, fill="#2c7fb8") +
  scale_y_discrete(
    labels=setNames(clean_gs_label_ForFig(gseaRes_BAscitch_HypoCAF$Description) %>% 
                      str_replace(pattern="Epithelial mesenchymal transition",
                                  replacement="EMT"),
                    gseaRes_BAscitch_HypoCAF$Description),
    expand = expansion(mult = c(0.1, 0.1)),
    position = "right") +
  Common.vline + 
  scale_x_continuous(
    position="bottom",
    limits=c(-max(abs(Range.NES.BAswitch)),max(abs(Range.NES.BAswitch))),
    expand=expansion(mult=c(0.05, 0.05))) + 
  Common.Scale.Size +
  theme(legend.position = "bottom",
        plot.margin = margin(t=0.1, r=0, b=0, l=0, "cm"),
        axis.title.x = element_text(face="bold", color="black"),
        axis.text.x = element_text(face="bold", color="black", size=9),
        axis.line.x = element_line(color="black")) +
  Common.Theme
aligned = cowplot::align_plots(
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
       file="InVitroGSEA_BeforeAfterSwitching2.png", dpi=500, bg="transparent")
Common.Scale.Size_Final =
  scale_size_continuous(
    name = "-log10(Adj.P)",
    breaks = c(2, 10, 30),
    limits = c(0, max(Range.minuslog10adjP.BAswitch)),
    range = c(0, 10),
    guide = guide_legend(
      direction = "horizontal",
      override.aes = list(fill="white", color="black"),
      title.position = "top",
      label.position = "bottom"))
saveRDS(Plot_NormoCAF + Common.Scale.Size_Final, file = "Fig2F_top.rds")
saveRDS(Plot_HypoCAF + Common.Scale.Size_Final, file = "Fig2F_bottom.rds")
saveRDS(get_legend(Plot_NormoCAF + Common.Scale.Size_Final), file = "Fig2F_legend.rds")
## . Enrichment plot
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
    CompNum==1 & Comparison=="Prolif.HvsN" ~ "Original_Prolif",
    CompNum==1 & Comparison=="Senes.HvsN" ~ "Original_Senes",
    CompNum==2 & Comparison=="Prolif.HvsN" ~ "Under1_Prolif",
    CompNum==2 & Comparison=="Senes.HvsN" ~ "Under1_Senes",
    CompNum==3 & Comparison=="Prolif.HvsN" ~ "Under21_Prolif",
    CompNum==3 & Comparison=="Senes.HvsN" ~ "Under21_Senes",
    CompNum==4 ~ "SwitchHypoCAF",
    CompNum==5 ~ "SwitchNormoCAF")
  Xlab = case_when(CompNum%in%c(1,2,3) ~ "Upregulated → Downregulated",
                   CompNum%in%c(4,5) ~ "21% O2 → 1% O2")
  LgdLab = clean_gs_label_ForFig(names(HMs)) %>%
    str_replace(pattern="Interferon alpha", replacement="IFNα") %>% 
    str_replace(pattern="Interferon gamma", replacement="IFNγ") %>% 
    setNames(names(HMs))
  EnrichPlot_0 = 
    enrichplot::gseaplot2(List.gseaRes[[CompNum]], 
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

lgd = cowplot::get_legend(
  p.res + theme(
    legend.position = "right",
    legend.justification = "left",
    legend.box.just = "left",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    plot.background = element_rect(fill="transparent", color=NA),
    legend.background = element_rect(fill="transparent", color=NA),
    legend.key = element_blank()
  )
)

panel_cow = cowplot::plot_grid(
  p.res_nlgd,
  p.gene_nlgd,
  ncol = 1,
  rel_heights = c(1.8, 0.3),
  align = "v",
  axis = "lr"
)

# 必要なら完成版もcowplotで保存
full_cow = cowplot::plot_grid(
  panel_cow,
  lgd,
  ncol = 2,
  rel_widths = c(1, 0.22),
  align = "h",
  axis = "tb"
)

# --- 保存 ---
saveRDS(p.res,       file = paste0("Fig2E_", Comparison, "_", CompNum, "_res.rds"))
saveRDS(p.gene,      file = paste0("Fig2E_", Comparison, "_", CompNum, "_gene.rds"))
saveRDS(p.res_nlgd,  file = paste0("Fig2E_", Comparison, "_", CompNum, "_res_nlgd.rds"))
saveRDS(p.gene_nlgd, file = paste0("Fig2E_", Comparison, "_", CompNum, "_gene_nlgd.rds"))
saveRDS(lgd,         file = paste0("Fig2E_", Comparison, "_", CompNum, "_legend.rds"))
saveRDS(panel_cow,   file = paste0("Fig2E_", Comparison, "_", CompNum, "_panel_cow.rds"))
saveRDS(full_cow,    file = paste0("Fig2E_", Comparison, "_", CompNum, "_full_cow.rds"))

}
Func.Penrich_Switch(CompNum=1, Comparison="Prolif.HvsN")
Func.Penrich_Switch(CompNum=1, Comparison="Senes.HvsN")
Func.Penrich_Switch(CompNum=2, Comparison="Prolif.HvsN")
Func.Penrich_Switch(CompNum=2, Comparison="Senes.HvsN")
Func.Penrich_Switch(CompNum=3, Comparison="Prolif.HvsN")
Func.Penrich_Switch(CompNum=3, Comparison="Senes.HvsN")
Func.Penrich_Switch(CompNum=4, Comparison="Prolif.BAswitch")
Func.Penrich_Switch(CompNum=5, Comparison="Prolif.BAswitch")

## 5. GSEA(Hallmark) sup table
List.gseaRes = readRDS(paste0(DirOutputMSigDB,"[ListRDS]_gseaHallmark_Results_",SetName,".rds"))
Func.DataProcess = function(DataFrame, SettingName){
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
DF.gseaRes_1_HvsN = List.gseaRes[[1]]@result %>% Func.DataProcess("Original")
DF.gseaRes_2_HvsNH = List.gseaRes[[2]]@result %>% Func.DataProcess("Under1")
DF.gseaRes_3_HNvsN = List.gseaRes[[3]]@result %>% Func.DataProcess("Under21")
DF.gseaRes_4_HvsHN = List.gseaRes[[4]]@result %>% Func.DataProcess("Switch_HCAF")
DF.gseaRes_5_NHvsN = List.gseaRes[[5]]@result %>% Func.DataProcess("Switch_NCAF")
DF.gseaResSup_HvsN_0 = bind_rows(DF.gseaRes_1_HvsN, 
                                 DF.gseaRes_3_HNvsN, 
                                 DF.gseaRes_2_HvsNH) %>% 
  dplyr::rename(In_Hypo_CAF = Direction) %>% 
  dplyr::mutate(Description = clean_gs_label_ForTable(Description))
DF.gseaResSup_Switch_0 = bind_rows(DF.gseaRes_4_HvsHN, 
                                   DF.gseaRes_5_NHvsN) %>% 
  dplyr::rename(Upreg_in = Direction) %>% 
  dplyr::mutate(Upreg_in = str_replace(Upreg_in, pattern="UP", replacement="Under1pct") ) %>%
  dplyr::mutate(Upreg_in = str_replace(Upreg_in, pattern="DOWN", replacement="Under21pct") ) %>% 
  dplyr::mutate(Description = clean_gs_label_ForTable(Description))
write.csv(DF.gseaResSup_HvsN_0, file="SupplementalTable_bulk_GSEA(Hallmark)_HvsN_ASCII.csv", row.names=FALSE)
write.csv(DF.gseaResSup_Switch_0, file="SupplementalTable_bulk_GSEA(Hallmark)_Switch_ASCII.csv", row.names=FALSE)
str(DF.gseaResSup_HvsN_0)
str(DF.gseaResSup_Switch_0)


##  6. Main figure 2               ####
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
p2B.top <- readRDS("Fig2B_Volc_Upper.rds") + ThemeEraseXYaxistext + ZeroMgn
p2B.mid <- readRDS("Fig2B_Volc_Mid.rds") + ThemeEraseXYaxistext + ZeroMgn
p2B.bottom <-  readRDS("Fig2B_Volc_Lower.rds") + ThemeEraseXYaxistext + ZeroMgn
p2C <- readRDS("Fig2C.rds") + ZeroMgn + AxText.Y + AxText.X
p2D <- readRDS("Fig2D.rds") + ZeroMgn + AxText.Y + AxText.X
p2E <- readRDS("Fig2Eleft_rightYaxis.rds") + ZeroMgn + AxText.Y + AxText.X + 
  scale_x_discrete(position = "bottom")
p2F.sen1 <- readRDS("Fig2E_Senes.HvsN_1_panel_cow.rds")
p2F.sen2 <- readRDS("Fig2E_Senes.HvsN_2_panel_cow.rds")
p2F.sen3 <- readRDS("Fig2E_Senes.HvsN_3_panel_cow.rds")
p2F.pro1 <- readRDS("Fig2E_Prolif.HvsN_1_panel_cow.rds")
p2F.pro2 <- readRDS("Fig2E_Prolif.HvsN_2_panel_cow.rds")
p2F.pro3 <- readRDS("Fig2E_Prolif.HvsN_3_panel_cow.rds")
lgd.sen <- readRDS("Fig2E_Senes.HvsN_1_legend.rds")
lgd.pro <- readRDS("Fig2E_Prolif.HvsN_1_legend.rds")
p2G.top <- readRDS("Fig2F_top.rds") + NoLgd + ZeroMgn + AxText.Y
p2G.bottom <- readRDS("Fig2F_bottom.rds") + NoLgd + ZeroMgn + AxText.Y


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
Height_F <- (Height_B + 0.70*gap + Height_C) - (Height_G + 0.95*gap) # inch

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
    unit(0.70*gap, "in"), 
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
    unit(0.70*gap, "in"), 
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
    unit(Height_F + 0.95*gap + Height_G, "in")))
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
    unit(0.95*gap, "in"),
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
    unit(Height_F + 0.95*gap + Height_G, "in")))
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
    unit(0.95*gap, "in"),
    unit(Height_G_each, "in"), 
    unit(gap_inG, "in"),
    unit(Height_G_each, "in") ) )



#+#+#+#+#+#   finalize   #+#+#+#+#+#
png("aligned_plot_gap24_v4.png", 
    width = 22, 
    height = (Height_B + 0.70*gap + Height_C) + (gap + Height_E + gap + Height_D + gap), 
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






##  7. CAF signatures (CPM)        ####
DF.CPM = read.csv("[1]dataset/rnaseq_allCPM.csv", row.names=1)
DF.Elyada = read.csv(file="[1]dataset/ElyadaSupTable_S22_Orthologs.csv", row.names=1)
DF.logCPM = log2(DF.CPM+1)
DF.Z_logCPM = t(scale(t(DF.logCPM)))
# define gene sets
Vec.myCAFsignat.Xen = DF.Markers$Gene[DF.Markers$Classification1=="myCAF"]
Vec.iCAFsignat.Xen = DF.Markers$Gene[DF.Markers$Classification1=="iCAF"]
Vec.apCAFsignat.Xen = c("CD74","HLA-DRA")
Vec.myCAFsignat.Ely = subset(DF.Elyada, Subtype=="myCAF")$Ortholog_Human
Vec.iCAFsignat.Ely = subset(DF.Elyada, Subtype=="iCAF")$Ortholog_Human
Vec.apCAFsignat.Ely = subset(DF.Elyada, Subtype=="apCAF")$Ortholog_Human
Vec.myCAFsignat.Ely = Vec.myCAFsignat.Ely[!Vec.myCAFsignat.Ely=="H19"]
Vec.iCAFsignat.Ely = Vec.iCAFsignat.Ely[!Vec.iCAFsignat.Ely=="LY6S"]
Vec.apCAFsignat.Ely[Vec.apCAFsignat.Ely=="NHERF1"] = "SLC9A3R1"
# calculate (each 36 sample)
Vec.Score.myCAF.Xen_0 = apply(DF.Z_logCPM[Vec.myCAFsignat.Xen, ], 2, mean, na.rm=TRUE)
Vec.Score.iCAF.Xen_0 = apply(DF.Z_logCPM[Vec.iCAFsignat.Xen, ], 2, mean, na.rm=TRUE)
Vec.Score.apCAF.Xen_0 = apply(DF.Z_logCPM[Vec.apCAFsignat.Xen, ], 2, mean, na.rm=TRUE)
Vec.Score.myCAF.Ely_0 = apply(DF.Z_logCPM[Vec.myCAFsignat.Ely, ], 2, mean, na.rm=TRUE)
Vec.Score.iCAF.Ely_0 = apply((DF.Z_logCPM[Vec.iCAFsignat.Ely, ]), 2, mean, na.rm=TRUE)
Vec.Score.apCAF.Ely_0 = apply(DF.Z_logCPM[Vec.apCAFsignat.Ely, ], 2, mean, na.rm=TRUE)
# bind
DF.Scores_0 = 
  data.frame(Vec.Score.myCAF.Xen_0,
             Vec.Score.iCAF.Xen_0,
             Vec.Score.apCAF.Xen_0,
             Vec.Score.myCAF.Ely_0,
             Vec.Score.iCAF.Ely_0,
             Vec.Score.apCAF.Ely_0) %>% 
  dplyr::mutate(Sample = str_replace(colnames(DF.Z_logCPM), 
                                     pattern="_1|_2|_3",
                                     replace=""))%>% 
  tidyr::separate(Sample, into=c("Pair","Condition"), sep="_", remove=FALSE)
DF.Scores_1_ByCondition = DF.Scores_0 %>%   # mean of replicated
  group_by(Condition) %>% 
  summarise(
    myCAF.Xen = mean(Vec.Score.myCAF.Xen_0, na.rm=TRUE),
    iCAF.Xen = mean(Vec.Score.iCAF.Xen_0, na.rm=TRUE),
    apCAF.Xen = mean(Vec.Score.apCAF.Xen_0, na.rm=TRUE),
    myCAF.Ely = mean(Vec.Score.myCAF.Ely_0, na.rm=TRUE),
    iCAF.Ely = mean(Vec.Score.iCAF.Ely_0, na.rm=TRUE),
    apCAF.Ely = mean(Vec.Score.apCAF.Ely_0, na.rm=TRUE)) %>% 
  dplyr::mutate(Condition = factor(Condition, levels=c("N","NH","gap","HN","H")))
CommonTheme = 
  theme(axis.text.x = element_text(face="bold", color="black"),
        axis.text.y = element_text(face="bold", color="black"),
        axis.title.y = element_text(face="bold", color="black"),
        legend.text = element_text(face="bold", color="black"),
        plot.background = element_rect(fill="transparent", color=NA),
        legend.background = element_rect(fill="transparent", color=NA),
        panel.background = element_rect(fill="white", color=NA),
        panel.border = element_rect(fill=NA, color="black"),
        panel.grid = element_blank())
Func.Dot.Agg = function(SignatName){
  DF = DF.Scores_1_ByCondition[ , c("Condition", SignatName)] %>% 
    set_colnames(c("Condition", "Signature"))
  DF.ForLine_N = 
    subset(DF.Scores_1_ByCondition, Condition %in% c("N","NH"))[ , c("Condition", SignatName)] %>% 
    set_colnames(c("Condition", "Signature"))
  DF.ForLine_H = 
    subset(DF.Scores_1_ByCondition, Condition %in% c("HN","H"))[ , c("Condition", SignatName)] %>% 
    set_colnames(c("Condition", "Signature"))
  P = ggplot(DF, aes(x=Condition, y=Signature)) + 
    geom_line(data=DF.ForLine_N, aes(group=1), color="black") + 
    geom_line(data=DF.ForLine_H, aes(group=1), color="black") + 
    geom_point(aes(fill=Condition), shape=21, size=PtSize) + 
    labs(y=paste0(SignatName,"\n\nSignature score\n(mean Z-scored log2CPM)"),
         x=NULL, fill=NULL) +
    scale_fill_manual(values=c("N"="#f03b20","NH"="#f03b20","HN"="#2c7fb8","H"="#2c7fb8")) + 
    CommonLabs + 
    scale_x_discrete(
      drop = FALSE, 
      breaks = c("N","NH","HN","H"),
      labels = c("21%\n\nN","1%\n\nNH"," 21%\n\nHN","1%\n\nH")) +
    CommonTheme
  return(P)
}
P.Agg_1 = Func.Dot.Agg("myCAF.Xen")
P.Agg_2 = Func.Dot.Agg("iCAF.Xen")
P.Agg_3 = Func.Dot.Agg("apCAF.Xen")
P.Agg_4 = Func.Dot.Agg("myCAF.Ely")
P.Agg_5 = Func.Dot.Agg("iCAF.Ely")
P.Agg_6 = Func.Dot.Agg("apCAF.Ely")
( P.Agg_1 | P.Agg_2 | P.Agg_3 )/( P.Agg_4 | P.Agg_5 | P.Agg_6 ) +
  plot_layout(guides = "collect")


DF.Scores_2_BySample = DF.Scores_0 %>% 
  group_by(Sample) %>% 
  summarise(
    myCAF.Xen = mean(Vec.Score.myCAF.Xen_0, na.rm=TRUE),
    iCAF.Xen = mean(Vec.Score.iCAF.Xen_0, na.rm=TRUE),
    apCAF.Xen = mean(Vec.Score.apCAF.Xen_0, na.rm=TRUE),
    myCAF.Ely = mean(Vec.Score.myCAF.Ely_0, na.rm=TRUE),
    iCAF.Ely = mean(Vec.Score.iCAF.Ely_0, na.rm=TRUE),
    apCAF.Ely = mean(Vec.Score.apCAF.Ely_0, na.rm=TRUE)) %>% 
  tidyr::separate(Sample, into=c("Pair","Condition"), sep="_", remove=FALSE) %>% 
  dplyr::mutate(Condition = factor(Condition, levels=c("N","NH","gap","HN","H")))
Func.Dot.BySample = function(SignatName, PtSize){
  DF = DF.Scores_2_BySample[ , c("Sample","Pair","Condition", SignatName)] %>% 
    set_colnames(c("Sample","Pair","Condition", "Signature"))
  DF.ForLine_N = subset(DF.Scores_2_BySample, Condition %in% c("N","NH"))[ , c("Sample","Pair","Condition", SignatName)] %>% 
    set_colnames(c("Sample","Pair","Condition", "Signature"))
  DF.ForLine_H = subset(DF.Scores_2_BySample, Condition %in% c("HN","H"))[ , c("Sample","Pair","Condition", SignatName)] %>% 
    set_colnames(c("Sample","Pair","Condition", "Signature"))
  P = 
    ggplot(DF, aes(x=Condition, y=Signature, fill=Pair)) + 
    geom_line(data=DF.ForLine_N, aes(group=Pair), color="black") + 
    geom_line(data=DF.ForLine_H, aes(group=Pair), color="black") + 
    geom_point(shape=21, size=PtSize) + 
    labs(y=paste0(SignatName,"\n\nSignature score\n(mean Z-scored log2CPM)"),
         x=NULL, fill=NULL) +
    scale_fill_manual(
      values = c("X0904"="#BFE8BF", "X1115"="#FFF4B8", "X1122"="#D8CFF0"),
      labels = c("X0904"="Pair #2", "X1115"="Pair #8", "X1122"="Pair #9")) + 
    scale_x_discrete(
      drop = FALSE, 
      breaks = c("N","NH","HN","H"),
      labels = c("21%\n\nN","1%\n\nNH"," 21%\n\nHN","1%\n\nH")) + 
    CommonTheme
  return(P)
}
P.ByS.1 = Func.Dot.BySample("myCAF.Xen", PtSize=4)
P.ByS.2 = Func.Dot.BySample("iCAF.Xen", PtSize=4)
P.ByS.3 = Func.Dot.BySample("apCAF.Xen", PtSize=4)
P.ByS.4 = Func.Dot.BySample("myCAF.Ely", PtSize=4)
P.ByS.5 = Func.Dot.BySample("iCAF.Ely", PtSize=4)
P.ByS.6 = Func.Dot.BySample("apCAF.Ely", PtSize=4)
( P.ByS.1 | P.ByS.2 | P.ByS.3 )/( P.ByS.4 | P5 | P6 ) +
  plot_layout(guides = "collect")

DF.logCPM[c("PDGFRA","PDGFRB","DES","FN1","ISLR"), c(1:3, 10:12, 13:15, 22:24, 25:27, 34:36)]


##  8. Expression heatmap          ####
library(gplots)
library(circlize)
library(ComplexHeatmap)
library(colorRamp2)
library(grid)
library(gridExtra)

DF.TPM.4comp = read.csv("/Volumes/PortableSSD/[4]myR/[1]dataset/rnaseq_allTPM.csv", header=TRUE, row.names=1)
DF.TPM.2comp = DF.TPM.4comp[, c(1:3, 10:12, 13:15, 22:24, 25:27, 34:36)]
DF.logTPM.2comp = log2(DF.TPM.2comp + 1)
# Select high-SD genes
Func.SelectHighSD = function(DataFrame, Ngenes){
  NamedVec.SD = apply(DataFrame, 1, sd, na.rm = TRUE)
  Vec.OrderIndexBySD = order(NamedVec.SD, decreasing = TRUE)
  DF.OrderedBySD = DataFrame[Vec.OrderIndexBySD, , drop=FALSE]
  DF_OrderedBySD_topNgenes = slice_head(DF.OrderedBySD, n=Ngenes)
  MT_Z_SDtopNgenes = scale(t(DF_OrderedBySD_topNgenes)) %>% as.matrix() %>% t()
  if(sd(MT_Z_SDtopNgenes[1, ])==1){
    return(MT_Z_SDtopNgenes)
  }
}
MT.Z.2comp_n2000 = Func.SelectHighSD(DF.logTPM.2comp, Ngenes=2000)
MT.Z.2comp_n2000_BySample <- MT.Z.2comp_n2000

# Parameter
FontSize <- 16
# Size (inch)  
Height.Plot <- 7.018373   # HeatmapBody + BottomColorBar*2
Height.BottomColBar <- 0.236  # inch
Height.HeatmapBody <- Height.Plot - 2*Height.BottomColBar
Height.Legend = unit(1.968504, "in")

Width.Plot <- 5.2 # PanelCのAxisText+Panel（= ClustTree + Heatmapbody）
Width.HeatmapBody <- 4.108
Width.ColumnGap <- 0.0315
Width.Tree <- Width.Plot - Width.HeatmapBody
Width.RightAnnot <- 1.44 # = Linkage + Text
Width.Link <- 0.30

ColorScale = colorRamp2(
  c(-2, 0, 2),
  c("#98CCBE", "#F7F7F7", "#C05A9E"))
Colors.Oxy = c(
  N = "#E8B7C4",
  H = "#B9C6E6")
Colors.Sample = c(
  Pair2 = "#C8E3C8",
  Pair8 = "#F3E8B4",
  Pair9 = "#D6CBE8")


Func.Index.2comp = function(Gene){ 
  Idx = which(rownames(MT.Z.2comp_n2000)==Gene)
  if(length(Idx) == 0) return(NA_integer_)
  Idx
}
Vec.ProSenesGenes_0 = c("TOP2A","BIRC5","FOXM1","CDK1","CCNA2","TPX2","CCNB1","CENPA","CDC45",
                        "CDC25","ZNF77",#,"EFNB1","EBP"
                        "CDKN2A","CDKN1A")
Vec.ProSenesGenes.2comp = sapply(Vec.ProSenesGenes_0, Func.Index.2comp)
Vec.ProSenesGenes.2comp
Vec.ProSenesGenes.2comp_rmNA = Vec.ProSenesGenes.2comp[!is.na(Vec.ProSenesGenes.2comp)]
RowAnnot.2comp = 
  rowAnnotation(
    foo = anno_mark(
      at = unname(Vec.ProSenesGenes.2comp_rmNA),
      labels = names(Vec.ProSenesGenes.2comp_rmNA),
      labels_gp = gpar(col="black", fontface="bold", fontsize=FontSize),
      link_width = unit(Width.Link, "in")),
    width = unit(Width.RightAnnot, "in"))
Deco.2comp_BySample = function(){
  FontSize.Deco = FontSize
  for(i in 1:3){
    decorate_annotation("CAF_pair", 
                        slice = i, {
                          grid.rect(x = unit(0.5, "npc"),
                                    y = unit(0.5, "npc"),
                                    width  = unit(1, "npc"),
                                    height = unit(1, "npc"),
                                    gp = gpar(fill = NA, col = "black", lwd = 0.7))                    
                        })
    decorate_annotation("Oxygen", slice = i, {
      grid.rect(
        x = unit(0.25, "npc"),
        y = unit(0.50, "npc"),
        width  = unit(0.5, "npc"),
        height = unit(1, "npc"),
        gp = gpar(fill = NA, col = "black", lwd = 0.7))
      grid.rect(
        x = unit(0.75, "npc"),
        y = unit(0.50, "npc"),
        width  = unit(0.5, "npc"),
        height = unit(1, "npc"),
        gp = gpar(fill = NA, col = "black", lwd = 0.7))
    })
  }
}
# Plot
Func.ComplexHeatmap = function(
    FileName, GeneClustK){
  png(file=paste0("heat_",FileName,".png"),
      width = 6.5, height = 10, units = "in", res = 500, bg="transparent")
  set.seed(1234)
  # Color Bars
  ColorBarBottom.2comp_BySample = 
    HeatmapAnnotation(
      CAF_pair = anno_simple(
        c(rep("Pair2",times=6),rep("Pair8",times=6),rep("Pair9",times=6)),
        col = Colors.Sample,
        border = FALSE,
        height = unit(Height.BottomColBar, "in")),
      Oxygen = anno_simple(
        rep(c("N","N","N","H","H","H"),times=3),
        col = Colors.Oxy,
        border = FALSE,
        height = unit(Height.BottomColBar, "in")),
      show_legend = FALSE, annotation_name_side = "left", 
      annotation_name_gp = gpar(fontface="bold", fontsize=FontSize),
      annotation_label = c("",""))
  # Heatmap body
  hm1 = 
    ComplexHeatmap::Heatmap(
      MT.Z.2comp_n2000_BySample,
      name = "main",
      show_row_names = FALSE,
      show_column_names = FALSE,
      column_title = NULL, column_title_side = "bottom",
      row_title = NULL,
      row_title_gp = gpar(fontface = "bold"),
      row_title_rot = 90,
      row_km = GeneClustK, border = FALSE,
      na_col = "black", 
      col = ColorScale,
      show_heatmap_legend = FALSE,
      bottom_annotation = ColorBarBottom.2comp_BySample,
      right_annotation = RowAnnot.2comp,
      cluster_rows = TRUE,
      clustering_distance_rows = "euclidean",
      clustering_method_rows = "ward.D2", # Default="complete"
      show_row_dend = TRUE,
      row_dend_width = unit(Width.Tree, "in"),
      row_dend_gp = gpar(col = "grey30", lwd = 0.6),
      cluster_columns = FALSE,
      column_split=c(rep("Pair2", 6), rep("Pair8", 6), rep("Pair9", 6)),
      column_gap=unit(Width.ColumnGap, "in"),
      width = unit(Width.HeatmapBody, "in"),
      height = unit(Height.HeatmapBody, "in"))
  # Legend; color scale bar
  Lgd1 = Legend(title="Z score (±2 clip)", 
                title_position="topcenter",
                at=c(-2,  0 ,  2),
                labels=c("-2",  "0" ,  "2"),
                col_fun=ColorScale, border="black", 
                title_gp = gpar(fontface="bold", fontsize=FontSize*0.8),
                labels_gp = gpar(fontface="bold", fontsize=FontSize*0.7),
                direction = "horizontal",
                legend_width  = unit(1.20, "in"),
                legend_height = unit(0.23, "in"))
  Legends = packLegend(Lgd1, direction="horizontal")
  # draw
  draw(hm1, 
       annotation_legend_list = Legends,
       annotation_legend_side = "bottom")
  Deco.2comp_BySample()
  dev.off()
}
Func.ComplexHeatmap(
  FileName="2compBySample_SpearmanWard_final_k=6",
  GeneClustK = 6)


##  9. GSEA myCAF/iCAF/CAF-8 score ####

library(clusterProfiler)

# DEGs from Hypo-/Normo-CAFs
DirRNAseqNew = '/Volumes/PortableSSD/[4]myR/[1]dataset/RNAseq_new/'

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


# Term to gene: major 3 subtype signature
TermToGene.Ely_0 <- read.csv(file="/Volumes/PortableSSD/[4]myR/[1]dataset/ElyadaSupTable_S22_Orthologs.csv", row.names=1)
TermToGene.Ely_1 <- 
  data.frame(
    term = TermToGene.Ely_0$Subtype,
    gene = as.character(TermToGene.Ely_0$ENTREZID))
TermToGene.Ely_2 <- na.omit(TermToGene.Ely_1)

# Term to gene: CAF-8
Vec.SymbolToEntrez <- 
  setNames(as.character(DF.AllGenes_addEntrez$ENTREZID),
           DF.AllGenes_addEntrez$SYMBOL)
DF.DEGsCAF8_0 <-  read.csv(
  file = paste0(DirInteg, "DEG_table/[SubclustsOfCAFs_DEGs_All][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2,".csv"), 
  row.names=1) %>% 
  dplyr::filter(cluster == "8") %>% 
  dplyr::mutate(entrez = Vec.SymbolToEntrez[gene])

DF.DEGsCAF8_1 <- DF.DEGsCAF8_0 %>% 
  dplyr::filter(
    avg_log2FC > 1 &
      p_val_adj < 0.05 &
      pct.1 > 0.2)
DF.DEGsCAF8_1$gene
DF.DEGsCAF8_2 <- 
  DF.DEGsCAF8_1 %>% 
  dplyr::select(gene, avg_log2FC, pct.1, p_val_adj) %>% 
  dplyr::mutate(pct.1 = 100*pct.1) %>% 
  set_colnames(c("Gene symbol", "Avarage log2 fold change", "Percent expression in CAF-8 (%)", "Adjusted P value"))
write.csv(DF.DEGsCAF8_2,
          file = "SuppTable_CAF8_marker_genes.csv",
          row.names = FALSE)
TermToGene.CAF8 <- 
  data.frame(
    term = "CAF8",
    gene = dplyr::pull(DF.DEGsCAF8_1, entrez),
    row.names = NULL)

# Run GSEA
Func.GSEA.Ely = function(GeneList){
  set.seed(1234)
  GSEA(geneList = GeneList,
       minGSSize = 10,
       TERM2GENE = rbind(
         TermToGene.Ely_2,
         TermToGene.CAF8),
       pvalueCutoff = 1,
       verbose = FALSE, 
       eps = 0)
}
gseaRes_1_HvsN = Func.GSEA.Ely(Vec.namedFC_1_HvsN_rmNA)
gseaRes_2_HNvsN = Func.GSEA.Ely(Vec.namedFC_2_HNvsN_rmNA)
gseaRes_3_HvsNH = Func.GSEA.Ely(Vec.namedFC_3_HvsNH_rmNA)
gseaRes_4_HvsHN = Func.GSEA.Ely(Vec.namedFC_4_HvsHN_rmNA)
gseaRes_5_NHvsN = Func.GSEA.Ely(Vec.namedFC_5_NHvsN_rmNA)

# Result data
DF.Res.CAFsubtype_0 = 
  bind_rows(
    gseaRes_1_HvsN@result %>% cbind(Condition = "Original"),
    gseaRes_2_HNvsN@result %>% cbind(Condition = "Under1"),
    gseaRes_3_HvsNH@result %>% cbind(Condition = "Under21"))　%>% 
  dplyr::mutate(
    Description = factor(Description, levels=c("CAF8","myCAF","iCAF","apCAF")),
    Condition = factor(Condition, levels=c("Original","Under21","Under1")),
    minuslog10adjP = (-1)*log10(p.adjust),
    Signif = case_when(p.adjust<0.001 ~ "***",
                       p.adjust<0.01 ~ "**",
                       p.adjust<0.05 ~ "*",
                       TRUE ~ "n.s."),
    Enrichment_in_HypoCAF = ifelse(NES>0, "UP", "DOWN"))

DF.Res.CAFsubtype_1 = 
  DF.Res.CAFsubtype_0 %>% 
  dplyr::select(c(Condition, Description, Enrichment_in_HypoCAF,
                  setSize, NES, p.adjust, qvalue, leading_edge, core_enrichment)) %>% 
  dplyr::mutate(Description = factor(Description, levels=c("CAF8", "myCAF", "iCAF","apCAF")),
                Condition = factor(Condition, levels=c("Original", "Under21", "Under1"))) %>% 
  dplyr::arrange(Condition) %>% 
  group_by(Condition) %>% 
  dplyr::arrange(Description, .by_group = TRUE)　%>% 
  ungroup()
write.csv(DF.Res.CAFsubtype_1,
          file="SuppTable_GSEA_CAFsubtype.csv",
          row.names=F)

# Plot
Range.NES = range(DF.Res.CAFsubtype_0$NES)
Plot.CAFsubtype = 
  ggplot(DF.Res.CAFsubtype_0, aes(y=Description, x=NES)) +
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
ggsave(plot = Plot.CAFsubtype,
       file = "CAFsubtype_new.png",
       width=6, height=6, dpi=500, bg="transparent")

Plot.CAFsubtype.Tile <- 
  DF.Res.CAFsubtype_0 %>% 
  mutate(Description = factor(Description, levels = rev(levels(Description)))) %>% 
    ggplot(aes(x = Condition, y = Description)) +
    geom_tile(aes(fill = NES), color = "gray50") +
    geom_text(aes(label = Signif,
                  vjust = ifelse(p.adjust < 0.05, 0.9, 0.4),
                  size = ifelse(p.adjust < 0.05, 8.0, 5.5)), 
              fontface = "bold") +
    labs(x = NULL, y = NULL) +
    scale_x_discrete(expand = expansion(mult = c(0, 0)),
                     labels = function(x) {
                       label_map <- c(
                         Original = "bold('Isolation')",
                         Under21 = "bold('21% ' * O[2])",
                         Under1 = "bold('1% ' * O[2])")
                       parse(text = unname(label_map[x]))
                     }
                     ) +
    scale_y_discrete(expand = expansion(mult = c(0, 0)),
                     labels = c(CAF8 = "CAF-8")) +
    scale_fill_distiller(
      palette="RdBu",
      direction=1,
      limits=c(-(max(abs(Range.NES))),max(abs(Range.NES))),
      breaks=c(2, 0, -2),
      guide = guide_colorbar(
        title.position="top",
        title.hjust=0.5,
        frame.colour="black", 
        ticks.colour="black",
        label.position="right")) +
    scale_size_identity(guide = "none") +
    theme(
      aspect.ratio = 1.2,
      axis.text.x = element_text(face = "bold", color = "black", size = 15,
                                 angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(face = "bold", color = "black", size = 15),
      legend.title = element_text(face="bold", color="black", size=12),
      legend.text = element_text(face="bold", color="black", size=10),
      plot.background = element_rect(fill="transparent", color=NA),
      legend.background = element_rect(fill="transparent", color=NA),
      panel.background = element_rect(fill="transparent", color=NA),
      panel.grid = element_blank(),
      panel.border = element_rect(fill="transparent", color = "black"))
Plot.CAFsubtype.Tile
ggsave(plot = Plot.CAFsubtype.Tile,
       file = "/Volumes/PortableSSD/[4]myR/CAFsubtype_tile_new2.png",
       width = 5, height = 5, dpi=500, bg="transparent")

##

library(enrichplot)
Func.getGSEAcurve <- function(gseaRes, Condition){
  DF_0 <- enrichplot:::gsInfo(gseaRes, geneSetID = "CAF8")
  DF_1 <- DF_0 %>% 
    dplyr::mutate(
      condition = Condition,
      RankPct = x / max(x) * 100)
  }
DF.CAF8curve <- 
  bind_rows(
    Func.getGSEAcurve(gseaRes = gseaRes_1_HvsN, Condition = "Isolation"),
    Func.getGSEAcurve(gseaRes = gseaRes_2_HNvsN, Condition = "Under21"),
    Func.getGSEAcurve(gseaRes = gseaRes_3_HvsNH, Condition = "Under1")) %>% 
  dplyr::mutate(condition = factor(condition, levels = c("Isolation", "Under21", "Under1")))
Plot.CAF8curve <- 
  DF.CAF8curve %>% 
  ggplot(aes(x = RankPct, y = runningScore)) +
  geom_line(aes(color = condition), 
            linewidth = 1.0) +
  labs(title = "CAF-8 signature",
       y = "Enrichment score", 
       x = "Rank in ordered gene list (%)") +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0)),
    breaks = c(0, 25, 50, 75, 100)) +
  scale_color_manual(
    values = c("#4D4D4D", "#008B8B", "#B2A100"),
    labels = function(x) {
      label_map <- c(
        Isolation = "bold('Isolation')",
        Under21 = "bold('21% ' * O[2])",
        Under1 = "bold('1% ' * O[2])")
      parse(text = unname(label_map[x]))
    }) +
  geom_hline(yintercept = 0, color = "gray50", linewidth = 0.2) + 
  guides(color = guide_legend(
    override.aes = list(linewidth = 1.5),
    keywidth = unit(0.35, "in"))) +
  theme(
    aspect.ratio = 1.0,
    legend.position = c(0.78, 0.78),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    panel.border = element_rect(fill = "transparent", color = NA),
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(face = "bold", color = "black", size = 12),
    axis.text.y = element_text(face = "bold", color = "black", size = 12),
    legend.text = element_text(face = "bold", color = "black", size = 15),
    axis.title = element_text(face = "bold", color = "black", size = 15),
    legend.title = element_blank()
  )
ggsave(ggarrange(Plot.CAFsubtype.Tile,
                 Plot.CAF8curve,
                 ncol = 2, align = "h"),
       filename = "/Volumes/PortableSSD/[4]myR/CAFsubtype_TileAndRunning.png",
       width = 10, height = 5, dpi = 500)

## ↓↓落選↓↓
Plot.CAF8curve.facet <- 
  DF.CAF8curve %>% 
  ggplot(aes(x = RankPct, y = runningScore)) +
  geom_line(#aes(color = condition), 
            linewidth = 1.0) +
  labs(title = "CAF-8",
       y = "Enrichment score", 
       x = "Rank in\nordered gene list (%)") +
  facet_grid(condition ~ ., 
             scales = "free_y",
             axes = "all") +
  geom_hline(yintercept = 0, color = "gray50", linewidth = 0.2) + 
  theme(
    aspect.ratio = 1.0,
    legend.position = "none",
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    panel.border = element_rect(fill = "transparent", color = NA),
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(face = "bold", color = "black"),
    axis.text.y = element_text(face = "bold", color = "black"),
    legend.text = element_text(face = "bold", color = "black"),
    axis.title = element_text(face = "bold", color = "black"),
    legend.title = element_blank()
  )
  Plot.CAF8curve.facet

  
##  10. visualization of CAF-8 signature ####

# CAF8 signature
Vec.CAF8signature <- 
 read.csv(file = "SuppTable_CAF8_marker_genes.csv") %>% 
 dplyr::pull(Gene.symbol)
 
# DEG table of each comparison settings
DirRNAseqNew = '/Volumes/PortableSSD/[4]myR/[1]dataset/RNAseq_new/'
  
DF.rankedDEG_1_HvsN <- 
  read.csv(file=paste0(DirRNAseqNew,"[DEGsTable_addENTREZ]_Comp1_HvsN.csv"), 
           header=T, row.names=NULL) %>% 
  dplyr::rename(Gene = X) %>% 
  dplyr::mutate(Comparison = "Isolation")
DF.rankedDEG_2_HNvsN <- 
  read.csv(file=paste0(DirRNAseqNew,"[DEGsTable_addENTREZ]_Comp2_HNvsN.csv"), 
           header=T, row.names=NULL) %>% 
  dplyr::rename(Gene=X) %>% 
  dplyr::mutate(Comparison = "Under21")
DF.rankedDEG_3_HvsNH <- 
  read.csv(file=paste0(DirRNAseqNew,"[DEGsTable_addENTREZ]_Comp3_HvsNH.csv"), 
           header=T, row.names=NULL) %>% 
  dplyr::rename(Gene=X) %>% 
  dplyr::mutate(Comparison = "Under1")

DF.signature_1_HvsN <- dplyr::filter(DF.rankedDEG_1_HvsN, Gene %in% Vec.CAF8signature)
DF.signature_2_HNvsN <- dplyr::filter(DF.rankedDEG_2_HNvsN, Gene %in% Vec.CAF8signature)
DF.signature_3_HvsNH <- dplyr::filter(DF.rankedDEG_3_HvsNH, Gene %in% Vec.CAF8signature)

DF.log2FC_0 <- 
  bind_rows(DF.signature_1_HvsN,
           DF.signature_2_HNvsN,
           DF.signature_3_HvsNH) %>% 
  dplyr::mutate(
    Gene = factor(Gene, levels = rev(DF.signature_1_HvsN$Gene)),
    Comparison = factor(Comparison, levels = c("Isolation", "Under21", "Under1")),
    Significance = ifelse(FDR<0.05, "Signif", "ns")
  )

CAF_subtype_Lollipop <- 
ggplot(DF.log2FC_0, aes(x = logFC, y = Gene)) +
  geom_hline(yintercept = c(1:25), color = "gray90", linewidth = 0.1) +
  geom_point(aes(color = Significance)) +
  geom_segment(aes(x = 0, xend = logFC, y = Gene, yend = Gene, color = Significance),
               linewidth = 0.6) +
  labs(y = NULL, 
       x = expression(bold(log[2]~"fold change (Hypo-CAF / Normo-CAF)"))) +
  geom_vline(xintercept = 0) +
  scale_x_continuous(limits = c(-2.5, 2.5)) +
  scale_color_manual(
    values = c(Signif = "Black", ns = "gray80"),
    labels = c(Signif = "FDR < 0.05", ns = "n.s.")) +
  facet_grid(
    . ~ Comparison,
    labeller = labeller(
      Comparison = as_labeller(c(
        Isolation = "bold(Isolation)",
        Under21   = "bold(21*'%'~O[2])",
        Under1    = "bold(1*'%'~O[2])"), 
        default = label_parsed))) +
  theme(legend.position = "bottom",
        panel.spacing = unit(5, "mm"),
        axis.text.x = element_text(face="bold", color="black", size=11),
        axis.text.y = element_text(face="bold", color="black", size=11),
        axis.title = element_text(face="bold", color="black", size=15),
        legend.title = element_blank(),
        legend.text = element_text(face="bold", color="black", size=12),
        plot.background = element_rect(fill="transparent", color=NA),
        strip.background = element_rect(fill="transparent", color=NA),
        legend.background = element_rect(fill="transparent", color=NA),
        panel.background = element_rect(fill="transparent", color=NA),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.border = element_rect(fill="transparent", color=NA),
        axis.line = element_line(color="black"))
ggsave(CAF_subtype_Lollipop,
       filename = "/Volumes/PortableSSD/[4]myR/CAFsubtype_Lollipop.png",
       width = 8, height = 5, dpi = 500)


CAF_subtype_HM <- 
DF.log2FC_0 %>% 
  dplyr::mutate(SignifLab = ifelse(FDR<0.05, "", "n.s.")) %>% 
  ggplot(aes(x = Comparison, y = Gene)) +
  geom_tile(aes(fill = logFC), color = "gray50") +
  geom_text(aes(label = SignifLab), color = "gray30", fontface = "bold", size = 3) +
  labs(x = NULL, y = NULL, 
       fill = expression(bold(log[2]~"fold change"))) +
  scale_x_discrete(expand = expansion(mult = c(0, 0)),
                   labels = c(Under21 = "21% O2",
                             Under1 = "1% O2")) +
  scale_y_discrete(expand = expansion(mult = c(0, 0))) +
  scale_fill_gradient2(
    low = "#A6D3C5", 
    mid = "#F7F7F7",
    high = "#C05A9E",
    midpoint = 0,
    limits = c(-2.5, 2.5),
    guide = guide_colorbar(
      #direction="horizontal",
      title.position="top",
      title.hjust=0.5,
      frame.colour="black",
      ticks.colour="black")) +
  theme(
    aspect.ratio = 2.3,
    legend.position = "bottom",
    axis.text.x = element_text(face="bold", color="black", size=10, 
                               angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(face="bold", color="black", size=10),
    legend.title = element_text(face="bold", color="black", size=10),
    legend.text = element_text(face="bold", color="black", size=8),
    plot.background = element_rect(fill="transparent", color=NA),
    strip.background = element_rect(fill="transparent", color=NA),
    legend.background = element_rect(fill="transparent", color=NA),
    panel.background = element_rect(fill="transparent", color=NA),
    panel.border = element_rect(fill="transparent", color="black"))
#CAF_subtype_HM
ggsave(CAF_subtype_HM,
       filename = "/Volumes/PortableSSD/[4]myR/CAFsubtype_Heatmap.png",
       width = 5, height = 6, dpi = 500)
  
#################################################
##
##     Supplemental table
##  0. All genes                  ####
SeuObj.CAF_0 = readRDS(paste0(DirInteg,"Objects/[SeuObj_ExtractedCAF][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                              QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2,".rds") )
DF.Allgenes = data.frame(Number = 1:nrow(SeuObj.CAF_0),
                         Genes = rownames(SeuObj.CAF_0))
write.csv(DF.Allgenes, file="AllGenes.csv", row.names=F)
##  1. Initial clustering summary ####
DataList = readRDS(
  file=paste0(DirInteg,"[MetaDataList][InitialClust][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
              "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".rds"))
DF.Data.Ini_0 = DataList[["MetaData"]]
DF.Data.Ini_1 = 
  table(DF.Data.Ini_0$seurat_clusters) %>% as.data.frame() %>% 
  set_colnames(c("initial_clusters", "cell_number")) %>% 
  dplyr::mutate(pct_of_total = sprintf("'%.1f%%", (cell_number/nrow(DF.Data.Ini_0))*100),
                cell_type = Vec.Annot[as.character(initial_clusters)]) %>% 
  dplyr::mutate(initial_clusters = paste0("C",initial_clusters),
                cell_number = formatC(cell_number, big.mark=","))
DF.Data.Ini_2 =
  table(DF.Data.Ini_0$sample, DF.Data.Ini_0$seurat_clusters) %>% as.data.frame() %>% 
  set_colnames(c("sample", "initial_clusters","cell_number")) %>% 
  dplyr::mutate(cell_type=Vec.Annot[as.character(initial_clusters)]) %>% 
  group_by(sample, cell_type) %>% 
  summarize(cell_number = sum(cell_number), .groups="drop") %>% 
  pivot_wider(names_from=sample,values_from=cell_number) %>% 
  column_to_rownames(var="cell_type")#%>% 
DF.Data.Ini_2.1 = DF.Data.Ini_2[names(List.Annot), ]
for(i in 1:ncol(DF.Data.Ini_2.1)){
  cell_sum = sum(DF.Data.Ini_2.1[,i])
  pct_cell = (DF.Data.Ini_2.1[ ,i]/cell_sum)*100
  DF.Data.Ini_2.1[,i] = 
    paste0(
      formatC(DF.Data.Ini_2.1[,i], big.mark=","),
      " (",sprintf("%.1f",pct_cell),"%)")
}
DF.Data.Ini_2.2 = 
  rbind("Total_cell_in_sample" = paste0("n = ", format(colSums(DF.Data.Ini_2), big.mark=",", trim=TRUE)),
        DF.Data.Ini_2.1) %>% 
  set_colnames(Vec.XeniumID[colnames(DF.Data.Ini_2.1)])
write.csv(DF.Data.Ini_2.2, "/Volumes/PortableSSD/[4]myR/SuppTable_IniClust_CellTypeCompositionAcrossSample.csv", row.names = T)
DF.Data.Ini_3 = 
  data.frame(n_cell = apply(DF.Data.Ini_2, 1, sum)) %>% 
  dplyr::mutate(pct_of_total = n_cell/sum(DF.Data.Ini_2)*100,
                n_cell_chr = formatC(n_cell, big.mark = ","),
                clusts = List.Annot[rownames(.)],
                clusts = sapply(clusts, function(x) paste0("C", sort(x), collapse=", ")))
DF.Data.Ini_4 = 
  DF.Data.Ini_3[names(List.Annot), ] %>% 
  rownames_to_column(var="cell_type") 
write.csv(DF.Data.Ini_4, "/Volumes/PortableSSD/[4]myR/SuppTable_IniClust_summary.csv", row.names = F)
##  2. Sub-clustering summary     ####
SeuObj.CAF_0 = readRDS(paste0(DirInteg,"Objects/[SeuObj_ExtractedCAF][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                              QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2,".rds") )
DF.SummarySub_0 = data.frame(
  "full_id"=colnames(SeuObj.CAF_0),
  "sample"=SeuObj.CAF_0$orig.ident,
  "subcluster"=SeuObj.CAF_0$seurat_clusters)
DF.SummarySub_SubclustsProp = 
  table(DF.SummarySub_0$subcluster) %>% as.data.frame() %>% 
  set_colnames(c("subcluster","cell_number")) %>% 
  dplyr::mutate(pct_of_CAF = sprintf("%.2f%%", (cell_number/nrow(DF.SummarySub_0))*100),
                pct_of_totalQCcells = sprintf("%.2f%%", (cell_number/1096788)*100),
                subcluster=paste0("CAF-",subcluster))
write.csv(DF.SummarySub_SubclustsProp,
          file="SubclustSummary.csv", row.names=F)
DF.SummarySub_BySample_count = 
  table(DF.SummarySub_0$subcluster, DF.SummarySub_0$sample) %>% as.data.frame() %>% 
  pivot_wider(names_from="Var2", values_from="Freq") %>% 
  column_to_rownames(var="Var1") %>% 
  set_rownames(paste0("CAF-",rownames(.))) %>% 
  set_colnames(as.character(Vec.XeniumID))
write.csv(DF.SummarySub_BySample_count,
          file="SubclustSummary_BySample_RawCount.csv", row.names=T)
DF.SummarySub_BySample_prop = 
  DF.SummarySub_BySample_count / rowSums(DF.SummarySub_BySample_count) * 100
write.csv(DF.SummarySub_BySample_prop,
          file="SubclustSummary_BySample_Prop.csv", row.names=T)

##  3. DEG table of sub-clusters  ####
DF.DEG_0 = read.csv(file = '/Volumes/Extreme SSD/Analysis/Data/IntegAnalysis/[6case(02,19,11,16,01,15)]_nFeatRNA:100~900_nCountRNA:100~1800/DEG_table/[SubclustsOfCAFs_DEGs_All][6case(02,19,11,16,01,15)]_Countable_mag1_nFeatRNA100~900_nCountRNA100~1800_1stPCA50vf2000_ClustDim20Res1_2ndPCA50vf2000ClustDim30Res0.5.csv')
DF.DEG_1 = DF.DEG_0 %>% 
  dplyr::mutate(subcluster = paste0("CAF-",cluster)) %>% 
  group_by(cluster) %>% 
  dplyr::arrange(desc(avg_log2FC), .by_group=TRUE) %>% ungroup() %>% 
  dplyr::select(subcluster, gene, avg_log2FC, pct.1, pct.2, p_val ,p_val_adj)
write.csv(DF.DEG_1,
          file="DEG.csv", row.names = F)
##  4. Mean and Z-score of ssGSVA ####
DF.gsvaScore_0 = 
  read.csv(file=paste0(DirInteg,"/[DataTable_ExtractedCAF_ssGSVA][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                       QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".csv"))
DF.gsvaScore_1 = DF.gsvaScore_0 %>% 
  dplyr::select(c(-X, -Y,
                  -apCAF_ssGSVA,
                  -Proliferation_ssGSVA,
                  -REACTOME_Cellular_senescence_ssGSVA)) %>% 
  dplyr::mutate(full_id = paste(DF.gsvaScore_0$sample, DF.gsvaScore_0$cell_id, sep="_")) %>% 
  dplyr::relocate(full_id, .before=sample) %>% 
  set_rownames(.$full_id)
SeuObj.CAF_0 = readRDS(paste0(DirInteg,"Objects/[SeuObj_ExtractedCAF][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                              QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2,".rds") )
DF.Meta_0 = data.frame("full_id"=colnames(SeuObj.CAF_0),
                       "seurat_clusters"=paste0("CAF-", SeuObj.CAF_0$seurat_clusters))
if( identical(DF.Meta_0$full_id, DF.gsvaScore_1$full_id) ){
  DF.Meta_1 = inner_join(DF.Meta_0, DF.gsvaScore_1, by="full_id")
}
DF.ScoreMean_0 = DF.Meta_1 %>% 
  dplyr::select(-c("full_id","sample","cell_id") ) %>% 
  group_by(seurat_clusters) %>% 
  summarise(across(everything(), mean)) %>% 
  column_to_rownames(var="seurat_clusters") %>%
  #set_colnames(str_remove(colnames(.), pattern="_ssGSVA")) %>% 
  #set_colnames(str_remove(colnames(.), pattern="HALLMARK_")) %>% 
  #set_colnames(tools::toTitleCase(tolower(colnames(.)))) %>% 
  dplyr::rename(Buffa_hypoxia_ssGSVA = BuffaOrig_ssGSVA,
                Winter_hypoxia_ssGSVA = WinterOrig_ssGSVA,
                myCAF_signature_ssGSVA = myCAF_ssGSVA,
                iCAF_signature_ssGSVA = iCAF_ssGSVA)
DF.ScoreMean_1 = DF.ScoreMean_0 %>% 
  rownames_to_column(var="Subcluster") %>% 
  pivot_longer(cols=-Subcluster,
               names_to="Gene_set", values_to="Mean_ssGSVA_SCORE")
DF.Zscore_0 = base::scale(DF.ScoreMean_0) %>% as.data.frame() %>% 
  rownames_to_column(var="Subcluster") %>% 
  pivot_longer(cols=-Subcluster,
               names_to="Gene_set", values_to="Z_score_within_gene_set")
DF.ssGSVAscore = 
  inner_join(DF.ScoreMean_1, DF.Zscore_0, by=c("Subcluster","Gene_set"))　%>% 
  dplyr::mutate(Gene_set = 
                  str_remove(Gene_set, pattern="_ssGSVA") %>% 
                  clean_gs_label_ForTable())
write.csv(DF.ssGSVAscore,
          file="SupplementalTable_ssGSVA_ASCII.csv",
          row.names=FALSE)
##  5. CellCycleScoring()         ####
SeuObj.CAF_1 = SeuObj.CAF_0 %>% 
  CellCycleScoring(s.features=cc.genes$s.genes,
                   g2m.features=cc.genes$g2m.genes,
                   set.ident=T#,
                   #seed=1234
                   )
DF.CCscores_0 = data.frame("Subcluster"=SeuObj.CAF_1$seurat_clusters,
                           "S.Score"=SeuObj.CAF_1$S.Score,
                           "G2M.Score"=SeuObj.CAF_1$G2M.Score) %>% 
  group_by(Subcluster) %>% 
  summarize(Mean_S_phase_score=mean(S.Score),
            Mean_G2.M_phase_score=mean(G2M.Score)) %>% 
  dplyr::mutate(Subcluster=paste0("CAF-",Subcluster))
write.csv(DF.CCscores_0,
          file="SupplementalTable_CellCycleScoring.csv",
          row.names=FALSE)
rm(SeuObj.CAF_1)
##  6. GSEA                       ####
DF.GSEA.res_0 = read.csv(
  file=paste0(DirInteg, "/[DataTable_ExtractedCAF_PseudoBulkGSEA][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
              QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2,".csv"),
  row.names=NULL)
DF.GSEA.res_1 <- 
  DF.GSEA.res_0 %>% 
  dplyr::select(c(-Description, -leading_edge, -rank, 
                  -core_enrichment, -Significance, -pvalue, 
                  -enrichmentScore, -qvalue)) %>%  
  dplyr::mutate(Subclust = paste0("CAF-",Subclust),
                ID = clean_gs_label_ForTable(ID)) %>% 
  dplyr::rename(Subcluster = Subclust,
                Direction = Sign,
                Gene_set = ID,
                FDR = p.adjust)
readr::write_excel_csv(DF.GSEA.res_1, "SupplementalTable_GSEA_readr.csv")
write.csv(DF.GSEA.res_1,
          file="SupplementalTable_GSEA.csv",
          row.names=FALSE)
##  7. Gene coverage              ####
library(msigdbr)
# 1. Get gene sets
CategoryH = msigdbr(species = "Homo sapiens", category = "H")
CategoryH_list = 
  CategoryH %>% 
  split(.$gs_name) %>%
  lapply(function(x) unique(x$gene_symbol))
CategoryC2 = msigdbr(species="Homo sapiens",  category="C2")
KEGG_Cell_cycle = 
  CategoryC2 %>% 
  filter(gs_name == "KEGG_CELL_CYCLE") %>% 
  pull(gene_symbol) %>% 
  unique()
DF.CAF8marker <- read.csv(
  file="SuppTable_CAF8_marker_genes.csv", row.names=NULL)
DF.Coverage_0 = 
  data.frame(
    GeneSet=character(),
    OrigGenes=character(),
    NumOfOrigGenes=numeric(),
    Intersect=character(),
    NumOfIntersect=numeric(),
    Coverage=numeric())
Func.Coverage = function(SetName, GeneVec){
  Vec.Intersect = intersect(GeneVec, Vec.AllGenes)
  DF.Add = 
    data.frame(
      GeneSet = SetName,
      OrigGenes = paste(sort(GeneVec), collapse=", "),
      NumOfOrigGenes = length(GeneVec),
      Intersect = paste(sort(Vec.Intersect), collapse=", "),
      NumOfIntersect = length(Vec.Intersect),
      Coverage = length(Vec.Intersect) / length(GeneVec) * 100)
}
DF.Coverage_1 = 
  bind_rows(
    DF.Coverage_0,
    Func.Coverage("Sphase", cc.genes$s.genes),
    Func.Coverage("G2Mphase", cc.genes$g2m.genes),
    Func.Coverage("KEGG_cell_cycle", KEGG_Cell_cycle),
    Func.Coverage("Winter", DF.Markers[DF.Markers$Classification1=="WinterCore", ]$Gene),
    Func.Coverage("Buffa", DF.Markers[DF.Markers$Classification1=="Buffa", ]$Gene),
    Func.Coverage("CAF8", DF.CAF8marker$Gene.symbol),
    Func.Coverage("myCAF", DF.Markers[DF.Markers$Classification1=="myCAF", ]$Gene),
    Func.Coverage("iCAF", DF.Markers[DF.Markers$Classification1=="iCAF", ]$Gene)
  )　%>% 
  dplyr::select(c(GeneSet, NumOfOrigGenes, NumOfIntersect, Coverage, Intersect))
write.csv(DF.Coverage_1, "SuppTable_coverage1.csv", row.names=F)
DF.Coverage_2 = DF.Coverage_0
Vec.Hallmarks = CategoryH_list %>% names()
for(i in 1:length(CategoryH_list)){
  DF.Coverage_2 = 
    rbind(DF.Coverage_2,
          Func.Coverage(Vec.Hallmarks[i], CategoryH_list[[Vec.Hallmarks[i]]]))
}
DF.Coverage_3 = DF.Coverage_2 %>% 
  dplyr::select(c(GeneSet, NumOfOrigGenes, NumOfIntersect, Coverage, Intersect)) %>% 
  dplyr::mutate(GeneSet = clean_gs_label_ForTable(GeneSet))
write.csv(DF.Coverage_3, "SuppTable_coverage2_ASCII.csv", row.names=F)
##
#################################################
##
##     Xenium data Initial
##  1. Dim Plot Initial                             ####
DataList = readRDS(
  file=paste0(DirInteg,"[MetaDataList][InitialClust][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
              "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".rds"))
DF.MetaData_0 = DataList[["MetaData"]]
set.seed(123)
DF.UMAP_1 = dplyr::slice_sample(DF.MetaData_0, prop=1)
DF.UMAP_1_rmUnclassified = 
  DF.UMAP_1 %>% 
  dplyr::filter(!(seurat_clusters %in% List.Annot[["Unclassified"]]))
DF.UMAPcenter = 
  DF.UMAP_1 %>% 
  group_by(seurat_clusters) %>% 
  dplyr::summarize(umap_1 = median(umap_1), umap_2 = median(umap_2)) %>% 
  dplyr::mutate(Label = paste0("C",seurat_clusters,".",
                               Vec.Annot[as.character(seurat_clusters)] ) )
CommonColorManualDimIni = scale_color_manual(values=Vec.OrderedCnumToColor) 
CommonThemeDimIni = theme(plot.margin=unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
                          legend.position="none",
                          #axis.text.x = element_text(face="bold", color="black", size=15),
                          #axis.text.y = element_text(face="bold", color="black", size=15),
                          axis.text.x = element_blank(),
                          axis.text.y = element_blank(),
                          axis.ticks = element_blank(),
                          axis.title = element_text(face="bold", color="black", size=17),
                          plot.background = element_rect(fill="transparent", color=NA),
                          panel.background = element_rect(fill="transparent", color=NA),
                          panel.grid = element_blank(),
                          panel.border = element_rect(fill="transparent", color="black"))
DimPlot1 = 
  ggplot(DF.UMAP_1_rmUnclassified, aes(x=umap_1, y=umap_2)) +
  geom_point(aes(color=seurat_clusters), size=0.3, alpha=0.2) + 
  geom_point(
    data=dplyr::filter(DF.UMAPcenter, !(seurat_clusters %in% List.Annot[["Unclassified"]])), 
    color="black", size=0.6) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  CommonColorManualDimIni + coord_fixed(ratio=1) + CommonThemeDimIni
DimPlot2 = 
  ggplot(DF.UMAP_1, aes(x=umap_1, y=umap_2)) +
  geom_point(aes(color=seurat_clusters), size=0.3, alpha=0.2) + 
  geom_point(
    data=dplyr::filter(DF.UMAPcenter, seurat_clusters %in% List.Annot[["Unclassified"]]), 
    color="black", size=0.6) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  CommonColorManualDimIni + coord_fixed(ratio=1) + CommonThemeDimIni
DimPlot3 = 
  ggplot(DF.UMAP_1, aes(x=umap_1, y=umap_2)) +
  geom_point(aes(color=seurat_clusters), size=0.3, alpha=0.2) + 
  geom_text_repel(
    data=DF.UMAPcenter, 
    aes(label=Label),
    color="black", box.padding=0.3, point.padding=0.2, segment.alpha=0.5) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  CommonColorManualDimIni + coord_fixed(ratio=1) + CommonThemeDimIni
DimPlotCombined <- 
  (DimPlot1 | DimPlot2 | DimPlot3) &
  theme(plot.background = element_rect(fill="transparent", color=NA),
        panel.background = element_rect(fill="transparent", color=NA))
ggsave(
  plot = DimPlotCombined,
  filename = paste0(DirInteg,"[Figure][InitialClust_DimPlot_Final]2.png"), 
  dpi=500, width=18, height=6, bg = "transparent" )

##  2. Feature plot representative cell type marker ####
GeneSet = c("COL1A1", "AMY2A", "SOX9", "KRT19",
            "SCG2", "PECAM1",
            "PTPRC", "CD3E", "CD4", "MZB1", "CD68", "KIT", "MPZ","CD19","CD3E",
            "ITGAX", "CD34")
MT.Exp_0 = 
  readRDS(paste0(DirInteg,"Objects/[SeuObj][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                 QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,
                 "_JoinLayer.rds") ) %>% 
  GetAssayData(layer="data") 
DF.MetaData_0_rmUnclassified = 
  DF.MetaData_0 %>% 
  dplyr::select(c(cell_id, umap_1, umap_2, seurat_clusters)) %>% 
  dplyr::filter(!(seurat_clusters %in% List.Annot[["Unclassified"]]))
for(i in 1:length(GeneSet)){
  TargetGene <- GeneSet[i]
  DF.Exp_1 <- 
    MT.Exp_0[TargetGene, ] %>% as.data.frame() %>% 
    set_colnames("TargetGene") %>% 
    rownames_to_column(var="cell_id")
  DF.MetaData_1 = 
    DF.MetaData_0_rmUnclassified %>% 
    inner_join(DF.Exp_1, by="cell_id")
  #DF.MetaData_1[ DF.MetaData_1$Label=="Unclassified", "TargetGene"] = NA
  set.seed(123)
  DF.MetaData_2 = DF.MetaData_1[sample(nrow(DF.MetaData_1)), ]
  #DF.MetaData_2 = dplyr::arrange(DF.MetaData_1, TargetGene)
  FeatPlot <-  
    ggplot(DF.MetaData_2, aes(x=umap_1, y=umap_2)) +
    geom_point(aes(color=TargetGene), size=0.1, alpha=0.2) +
    labs(title=TargetGene, color=NULL,
         x = "UMAP 1", y = "UMAP 2") +
    scale_color_viridis_c(
      option="magma",
      trans="sqrt",
      breaks = c(0, 2, 4, 6),
      guide=guide_colorbar(
        frame.colour="black",
        ticks.colour="black",
        barwidth=unit(4.5, "mm"))) +
    coord_fixed(ratio=1) +
    theme(legend.position = "right",
          plot.margin = margin(t=1, r=1, b=1, l=1, "mm"),
          legend.margin = margin(t=0, r=0, b=0, l=0, "cm"),
          plot.title = element_text(face="bold", color="black", size = 15),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_text(face="bold", color="black", size = 10),
          legend.text = element_text(face="bold", color="black"),
          plot.background = element_rect(fill="transparent", color=NA),
          legend.background = element_rect(fill="transparent", color=NA),
          panel.background = element_rect(fill="transparent", color=NA),
          panel.grid = element_blank(),
          panel.border = element_rect(fill="transparent", color="black"))
  ggsave(FeatPlot,
         file=paste0("/Volumes/PortableSSD/[4]myR/",TargetGene,".png"),
         width=3, height=3, dpi=500, bg="transparent")
}

##  3. Cluster composition within each sample       ####
DF.ClustCellCount_0 = 
  table(DF.MetaData_0$seurat_clusters) %>% as.data.frame() %>% 
  set_colnames(c("seurat_clusters","n_cell")) %>% 
  dplyr::mutate(prop = (n_cell/nrow(DF.MetaData_0))*100,
                Cell_type = Vec.Annot[as.character(seurat_clusters)])
DF.ClustCellCount.Unclasified = 
  subset(DF.ClustCellCount_0, Cell_type == "Unclassified")
#
DF.ClustCellCountBySample_0 = 
  table(DF.MetaData_0$sample, DF.MetaData_0$seurat_clusters) %>% as.data.frame() %>% 
  set_colnames(c("sample","seurat_clusters","n_cell")) %>% 
  dplyr::mutate(celltype = Vec.Annot[as.character(seurat_clusters)])
Vec.CelltypeToSupplab = c("Acinar"="Normal_acinar",
                          "Ductal_like_acinar"="Ductal",
                          "Normal_ductal"="Ductal",
                          "PanIN"="Ductal",
                          "PDAC"="Ductal",
                          "Islet"="Islet",
                          "CAF"="CAF",
                          "Mural"="EndothelialMural",
                          "Endothelial"="EndothelialMural",
                          "Lymph_B"="Lymphoid",
                          "Lymph_T"="Lymphoid",
                          "Plasma"="Lymphoid",
                          "Myeloid"="Myeloid",
                          "Mast"="Myeloid",
                          "Nerve"="Nerve",
                          "Unclassified"="Unclassified")
DF.ClustCellCountBySample_1 = 
  DF.ClustCellCountBySample_0 %>% 
  dplyr::mutate(Supplab = Vec.CelltypeToSupplab[celltype])
DF.ClustCellCountBySample_2 = 
  DF.ClustCellCountBySample_1 %>% 
  group_by(sample, Supplab) %>% 
  dplyr::summarise(n_cell_Supplab = sum(n_cell),
                   .groups = "drop")
DF.ClustCellCountBySample_2$Supplab = 
  factor(DF.ClustCellCountBySample_2$Supplab, 
         levels=unique(unname(Vec.CelltypeToSupplab)))
DF.CelltypCellCountBySample_1 = 
  DF.CelltypCellCountBySample_0 %>%
  dplyr::mutate(n_total_insample = sum(n_cell),
                prop_insample = (n_cell/n_total_insample)*100)
Plot.Bar = 
  ggplot(DF.ClustCellCountBySample_2,
         aes(x=sample, y=n_cell_Supplab, fill=Supplab)) +
  geom_bar(stat="identity", position="fill", color="gray30", linewidth=0.3) +
  labs(y="Cell type proportion\nwithin sample (%)", fill=NULL,
       x=NULL) + #"Sample ID", ) +
  scale_x_discrete(expand=expansion(mult=c(0.12, 0.12)),
                   labels=c("01"="01", "02"="02", "11"="03",
                            "15"="04", "16"="05", "19"="06")) +
  scale_y_continuous(labels=function(x) x*100,
                     expand=expansion(mult=0)) +
  scale_fill_manual(values=c("Normal_acinar"="#ff00ff",
                             "Ductal"="#ddbcff",
                             "Islet"="#ffff00",
                             "CAF"="#00ff00",
                             "EndothelialMural"="#990000",
                             "Lymphoid"="#00bfff",
                             "Myeloid"="#ff8800",
                             "Nerve"="#0000ff",
                             "Unclassified"="#D0D0D0"),
                    labels=c("Normal_acinar"="Normal acinar",
                             "Ductal"="Ductal",
                             "EndothelialMural"="Endothelial and Mural",
                             "Lymphoid"="Lymphoid",
                             "Myeloid"="Myeloid",
                             "Nerve"="Neuronal")) +
  coord_fixed(ratio=6.5) +
  theme(plot.margin=unit(c(0.5, 0.1, 0.1, 0.1), "cm"),
        legend.position="right",
        legend.key = element_blank(),
        legend.text = element_text(face="bold", size=10),
        axis.text = element_text(face="bold", color="black", size=12),
        axis.title = element_text(face="bold", size=12),
        plot.background = element_rect(fill="transparent", color=NA),
        legend.background = element_rect(fill="transparent", color=NA),
        panel.background = element_rect(fill="white", color=NA),
        panel.border = element_rect(fill=NA, color="black", linewidth=1),
        panel.grid = element_blank())
ggsave(
  plot = Plot.Bar,
  file=paste0(DirInteg,"[Figure][InitialClust_BarPlot2_ForPaper][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
              "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".png"), 
  dpi=500, width=6, height=4, bg = "transparent" )
##  4. Dot Plot Initial                             ####
SeuObjIni_0 <-
  readRDS(paste0(DirInteg,"Objects/[SeuObj][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                 QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_JoinLayer.rds") )
genes_use=c(
  "AMY2A","CPA1",
  "EPCAM","CDH1","SOX9",
  "FGFR3","C6","HNF1B",
  "MUC5B","TFF1","MUC1","KRT19","MUC5AC","CLDN18","CTSE",
  "LAMC2","COL17A1","FAM83A",　#"TNS4","SLC2A1","CDH3","MACC1","DKK1",#"CEACAM5",
  "INS","SCG2",
  "FAP","COL5A1","COL1A1",  #"PDGFRB","ACTA2","TAGLN",
  "MYH11","NOTCH3","RGS5",
  "PECAM1","CD34","FLT4",  #"IL3RA",
  "PTPRC","CD8A","GZMA","CD3E","CD4",
  "CD19",  #"CD79A",
  "CD38","TNFRSF17","MZB1",   #"CD27",
  "CSF1R","CD68","ITGAX",   #"CD14","ITGAM","FCGR3A",
  "LAMP3", "CCR7", "LY75",
  "KIT","CMA1",
  "SCN7A","MPZ","NCAM1")

#Order.Ini=c(
#  0,13,24,33,12,1,17,30,23,10,27,20,
#  14,9,6,7,11,34,3,32,2,8,16,15,18,37,5,4,21,
#  19,31,26,29,25,22,28,36,35)
Order.Ini=c(
  0,13,24,
  1,12,33,
  30,17,23,10,27,20,
  14,
  7,6,9,11,34,3,32,2,8,16,15,37,18,5,4,21,
  19,31,26,29,25,22,28,36,35)

Acinar=c(0,13,24)
Ductal_like_acinar=c(33,1,12,17,30)
Normal_ductal=23
PanIN=c(10,27)
PDAC = 20
Islet=c(14)
CAF=c(9,6,7)
Mural=c(11)
Endothelial=c(34,3)
Lymph_T=c(32,2,8,16)
Lymph_B=c(15)
Plasma=c(18,37)
Myeloid=c(5,4,21,19,31)
Mast=c(26)
Nerve=c(29)
Unclassified=c(36,35,25,22,28)
Vec.OrderedCnumToCT = 
  case_when(
    Order.Ini %in% Acinar ~ paste0("C",Order.Ini,". Acinar cell"),
    Order.Ini %in% Ductal_like_acinar ~ paste0("C",Order.Ini,". Ductal-like acinar cell"),
    Order.Ini %in% Normal_ductal ~ paste0("C",Order.Ini,". Ductal cell (Normal)"),
    Order.Ini %in% PanIN ~ paste0("C",Order.Ini,". Ductal cell (PanIN)"),
    Order.Ini %in% PDAC ~ paste0("C",Order.Ini,". Ductal cell (PDAC)"),
    Order.Ini %in% Islet ~ paste0("C",Order.Ini,". Islet cell"),
    Order.Ini %in% CAF ~ paste0("C",Order.Ini,". CAF"),
    Order.Ini %in% Mural ~ paste0("C",Order.Ini,". Mural cell"),
    Order.Ini %in% Endothelial ~ paste0("C",Order.Ini,". Endothelial cell/Pericyte"),
    Order.Ini %in% Lymph_T ~ paste0("C",Order.Ini,". T Lymphocyte"),
    Order.Ini %in% Lymph_B ~ paste0("C",Order.Ini,". B Lymphocyte"),
    Order.Ini %in% Plasma ~ paste0("C",Order.Ini,". Plasma celll"),
    Order.Ini %in% Myeloid ~ paste0("C",Order.Ini,". Myeloid cell"),
    Order.Ini %in% Mast ~ paste0("C",Order.Ini,". Mast cell"),
    Order.Ini %in% Nerve ~ paste0("C",Order.Ini,". Neuronal cell"),
    Order.Ini %in% Unclassified ~ paste0("C",Order.Ini,". Unclassified"),
    TRUE ~ paste0("C",Order.Ini,". Others") )

Plot.Dot_0 <- 
  SeuObjIni_0 %>% 
  subset(!(subset = seurat_clusters %in% List.Annot[["Unclassified"]])) %>% 
  DotPlot(
    features = genes_use,
    scale = TRUE)
DF.MarkerExp_0 <- Plot.Dot_0$data
DF.MarkerExp_1 <- 
  DF.MarkerExp_0 %>% 
  group_by(features.plot) %>% 
  dplyr::mutate(
    #Ylab = paste0("CAF-", id), 
    NormExp = {
      Range = max(avg.exp) - min(avg.exp)
      if (Range == 0) { rep(0, n = length(unique(DF.DEG_1$cluster))) 
      } else {
        (avg.exp - min(avg.exp)) / Range
      }
    },
    id = factor(id, levels = rev(Order.Ini))
  ) 
TextSize = 10
DotPlot <- 
  DF.MarkerExp_1 %>% 
  ggplot(aes(x = features.plot, y = id)) +
  geom_point(shape = 21,
             aes(size = pct.exp, fill = NormExp)) +
  labs(fill = "Scaled\nexpression\n(0–1 per gene)",
       size = "Percent\nexpressed") +
  scale_y_discrete(
    labels = setNames(Vec.OrderedCnumToCT, Order.Ini)) +
  scale_fill_gradientn(
    breaks = c(0, 0.5, 1),
    colors = 
      c("gray90", "white", "#fff7bc", "#fec44f", "#d95f0e"),
    guide = guide_colorbar(
      frame.colour = "black",
      ticks.colour = "black")) +
  scale_size(
    range = c(1, TextSize),
    guide = guide_legend(
      title.position = "top",
      direction = "vertical")) +
  theme(
    plot.title = element_blank(),
    axis.title = element_blank(),
    legend.position = "right",
    #legend.box.margin = margin(t=70),
    legend.key = element_blank(),
    axis.text.x = element_text(color = "black", face = "bold", size = TextSize*2,
                               angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(color = "black", face = "bold", size = TextSize*2),
    legend.title = element_text(color = "black", face = "bold", size = TextSize*1.5, margin = margin(b = 10)),
    legend.text = element_text(color = "black", face = "bold", size = TextSize*1.3),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = "transparent"),
    panel.background = element_rect(fill = "transparent", color = "transparent"),
    panel.border = element_rect(fill = "transparent", color = "black"),
    panel.grid = element_blank())
ggsave(DotPlot,
       file=paste0(DirInteg,"[Figure][InitialClust_DotPlot_CellTypeAndCAFmarkers_Final][",
                   NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                   "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".png"), 
       dpi=500, width=20, height=11, bg="transparent")

# Old version
Acinar=c(0,13,24)
Ductal_like_acinar=c(33,12,1,17,30)
Normal_ductal=23
PanIN=c(10,27)
PDAC = 20
Islet=c(14)
CAF=c(9,6,7)
Mural=c(11)
Endothelial=c(34,3)
Plasma=c(18,37)
Lymph_B=c(15)
Lymph_T=c(16,2)
Myeloid=c(5,4,31)
Mast=c(26)
Nerve=c(29)
Unclassified=c(21,32,8,19,36,35,25,22,28)
Order.Ini=c(Acinar,Ductal_like_acinar,Normal_ductal,PanIN,PDAC,Islet,CAF,Mural,
            Endothelial,Plasma,Lymph_B,Lymph_T,Myeloid,Mast,Nerve,Unclassified)
Vec.Cell.use = 
  WhichCells(SeuObjIni_0,
             idents = setdiff(unique(Idents(SeuObjIni_0)), Unclassified))
SeuObjIni_1 = subset(SeuObjIni_0, cells=Vec.Cell.use)
saveRDS(SeuObjIni_1, file="20260314.rds")
Vec.OrderedCnumToCT = 
  case_when(
    Order.Ini %in% Acinar ~ paste0("C",Order.Ini,". Acinar cell"),
    Order.Ini %in% Ductal_like_acinar ~ paste0("C",Order.Ini,". Ductal-like acinar cell"),
    Order.Ini %in% Normal_ductal ~ paste0("C",Order.Ini,". Ductal cell (Normal)"),
    Order.Ini %in% PanIN ~ paste0("C",Order.Ini,". Ductal cell (PanIN)"),
    Order.Ini %in% PDAC ~ paste0("C",Order.Ini,". Ductal cell (PDAC)"),
    Order.Ini %in% Islet ~ paste0("C",Order.Ini,". Islet cell"),
    Order.Ini %in% CAF ~ paste0("C",Order.Ini,". CAF"),
    Order.Ini %in% Mural ~ paste0("C",Order.Ini,". Mural cell"),
    Order.Ini %in% Endothelial ~ paste0("C",Order.Ini,". Endothelial cell"),
    Order.Ini %in% Lymph_T ~ paste0("C",Order.Ini,". T Lymphocyte"),
    Order.Ini %in% Lymph_B ~ paste0("C",Order.Ini,". B Lymphocyte"),
    Order.Ini %in% Plasma ~ paste0("C",Order.Ini,". Plasma celll"),
    Order.Ini %in% Myeloid ~ paste0("C",Order.Ini,". Myeloid cell"),
    Order.Ini %in% Mast ~ paste0("C",Order.Ini,". Mast cell"),
    Order.Ini %in% Nerve ~ paste0("C",Order.Ini,". Neuronal cell"),
    Order.Ini %in% Unclassified ~ paste0("C",Order.Ini,". Unclassified"),
    TRUE ~ paste0("C",Order.Ini,". Others") )
Order.Ini_Disp = Order.Ini[1:(length(Order.Ini)-length(Unclassified))]
genes_use=c("AMY2A","CPA1",
            "FGFR3","C6","HNF1B",
            "SOX9","EPCAM","CDH1",
            "TFF1","MUC1","KRT19","CTSE","CLDN18", #"MUC5B","MUC5AC", 
            "LAMC2","COL17A1","FAM83A",　#"TNS4","SLC2A1","CDH3","MACC1","DKK1",#"CEACAM5",
            "INS","SCG2",
            "FAP","COL5A1","COL1A1",  #"PDGFRB","ACTA2","TAGLN",
            "RGS5","NOTCH3","MYH11",
            "CD34","FLT4", #"PECAM1",  #"IL3RA",
            #"PTPRC",
            "MZB1","CD38","TNFRSF17",   #"CD27",
            "CD19",  #"CD79A",
            "CD4","CD3E","CD8A","GZMA",
            "ITGAX","CSF1R","CD68",   #"CD14","ITGAM","FCGR3A",
            "LAMP3", "LY75", #"CCR7", 
            "KIT","CMA1",
            "SCN7A","MPZ","NCAM1")
Idents(SeuObjIni_1) = factor(SeuObjIni_1$seurat_clusters, levels=Order.Ini_Disp)
dot.scale = 8
p1 =
  DotPlot(SeuObjIni_1, features=genes_use) +
  geom_point(aes(size = pct.exp),
             color = "black", shape=21, stroke=0.2)
p2 = 
  p1 +
  scale_size_continuous(range = c(0, dot.scale),
                        breaks = c(1, 40, 80),
                        labels = function(x) paste0(round(x), "%")) +
  #scale_x_discrete(expand=expansion(mult=c(0.013, 0.013))) +
  scale_y_discrete(labels=Vec.OrderedCnumToCT,
                   expand=expansion(mult=c(0.02, 0.02))) +
  labs(x = NULL, y = NULL, size = "Percent\nExprtessed",
       color = "Scaled\nAverage\nExpression") +
  guides(size = guide_legend(title = "Percent\nExpressed"),
         color = guide_colorbar(title = "Scaled\nAverage\nExpression",
                                direction="vertical",
                                frame.colour="black",
                                ticks.colour="black")) +
  theme(axis.text.x = element_text(face = "bold", size = 19, angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(face = "bold", size = 19, color = "black"),
        legend.text = element_text(face = "bold", size = 17, hjust= 0, color = "black"),
        legend.title = element_text(face = "bold", size = 19),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        panel.border = element_rect(fill = "transparent", color = "black"),
        panel.grid = element_blank(),
        axis.line = element_blank()) +
  scale_color_distiller(
    palette="RdBu", 
    limits=c(-max(abs(p1@data$avg.exp.scaled)),max(abs(p1@data$avg.exp.scaled))),
    breaks=c(-2,0,2)) +
  annotate(geom="rect", xmin=0.5, xmax=2.5, ymin=0.5, ymax=3.5, fill="orange", color="orange",alpha=0.2) +
  annotate(geom="rect", xmin=2.5, xmax=5.5, ymin=3.5, ymax=8.5, fill="orange", color="orange",alpha=0.2) +
  annotate(geom="rect", xmin=5.5, xmax=8.5, ymin=3.5, ymax=12.5, fill="orange", color="orange",alpha=0.2) +
  annotate(geom="rect", xmin=8.5, xmax=11.5, ymin=8.5, ymax=12.5, fill="orange", color="orange",alpha=0.2) +
  annotate(geom="rect", xmin=11.5, xmax=13.5, ymin=9.5, ymax=12.5, fill="orange", color="orange",alpha=0.2) +
  annotate(geom="rect", xmin=13.5, xmax=16.5, ymin=11.5, ymax=12.5, fill="orange", color="orange",alpha=0.2) + #PDAC
  annotate(geom="rect", xmin=16.5, xmax=18.5, ymin=12.5, ymax=13.5, fill="orange", color="orange",alpha=0.2) +
  # CAF
  annotate(geom="rect", xmin=18.5, xmax=21.5, ymin=13.5, ymax=16.5, fill="orange", color="orange",alpha=0.2) +
  # Mural
  annotate(geom="rect", xmin=21.5, xmax=24.5, ymin=16.5, ymax=17.5, fill="orange", color="orange",alpha=0.2) +
  # Endo
  annotate(geom="rect", xmin=24.5, xmax=26.5, ymin=17.5, ymax=19.5, fill="orange", color="orange",alpha=0.2) + 
  annotate(geom="rect", xmin=26.5, xmax=29.5, ymin=19.5, ymax=21.5, fill="orange", color="orange",alpha=0.2) + # CD8T
  annotate(geom="rect", xmin=29.5, xmax=30.5, ymin=21.5, ymax=22.5, fill="orange", color="orange",alpha=0.2) + # CD4T
  annotate(geom="rect", xmin=30.5, xmax=32.5, ymin=22.5, ymax=23.5, fill="orange", color="orange",alpha=0.2) + # B
  annotate(geom="rect", xmin=31.5, xmax=34.5, ymin=23.5, ymax=24.5, fill="orange", color="orange",alpha=0.2) + # Plasma
  annotate(geom="rect", xmin=34.5, xmax=37.5, ymin=24.5, ymax=26.5, fill="orange", color="orange",alpha=0.2) + # Myeloid
  annotate(geom="rect", xmin=37.5, xmax=39.5, ymin=26.5, ymax=27.5, fill="orange", color="orange",alpha=0.2) + # Myeloid
  annotate(geom="rect", xmin=39.5, xmax=41.5, ymin=27.5, ymax=28.5, fill="orange", color="orange",alpha=0.2) + # Mast
  annotate(geom="rect", xmin=41.5, xmax=44.5, ymin=28.5, ymax=29.5, fill="orange", color="orange",alpha=0.2) # Neural
ggsave(p2,
       file=paste0(DirInteg,"[Figure][5InitialClust_DotPlot_CellTypeAndCAFmarkers_Final][",
                   NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                   "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".png"), 
       dpi=500, width=21, height=11, bg="transparent")

##
#################################################
##
##     Xenium data Subclustering
##  1. Subcluster composition within each sample ####
SeuObj.CAF_0 = readRDS(
  paste0(DirInteg,"Objects/[SeuObj_ExtractedCAF][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
         QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2,".rds") )
DF.CellTable_0 = 
  table(SeuObj.CAF_0$seurat_clusters, SeuObj.CAF_0$orig.ident) %>% 
  as.data.frame() %>% 
  set_colnames(c("cluster","sample","cell")) %>% 
  dplyr::mutate(cluster = paste0("CAF-",cluster))
DF.CellsBySample_Labs_0 = 
  DF.CellTable_0 %>% 
  group_by(sample) %>% 
  dplyr::mutate(total=sum(cell)) %>% 
  ungroup() %>% rowwise() %>% 
  dplyr::mutate(prop=cell/total,
                PctLab=paste0(formatC(prop*100, format="f", digits=1),"%"),
                FullLab=paste0(formatC(cell, big.mark=","),"\n",
                               "(",PctLab,")")) %>% 
  ungroup() %>% group_by(sample) %>% 
  dplyr::mutate(cumsum=cumsum(prop),
                coord_y=1-cumsum+prop/2)
DF.CellsByClust_Labs_0 = 
  dplyr::arrange(DF.CellTable_0, cluster) %>% 
  group_by(cluster) %>% 
  dplyr::mutate(total=sum(cell)) %>% 
  ungroup() %>% rowwise() %>% 
  dplyr::mutate(prop=cell/total,
                PctLab=paste0(formatC(prop*100, format="f", digits=1),"%"),
                FullLab=paste0(formatC(cell, big.mark=","),"\n",
                               "(",PctLab,")")) %>% 
  ungroup() %>% group_by(cluster) %>% 
  dplyr::mutate(cumsum=cumsum(prop),
                coord_y=1-cumsum+prop/2)
TxtSizeSubCompos = 10
CommonTheme_bar = theme(
  axis.text.x = element_text(face="bold", color="black", size=TxtSizeSubCompos*1.3),
  axis.text.y = element_text(face="bold", color="black", size=TxtSizeSubCompos*1.5),
  axis.title = element_text(face="bold", color="black", size=TxtSizeSubCompos),
  legend.text = element_text(face="bold", color="black", size=TxtSizeSubCompos*1.3),
  legend.title = element_text(face="bold", color="black", size=TxtSizeSubCompos*1.3),
  legend.key = element_blank(),
  plot.background = element_rect(fill="transparent", color=NA),
  legend.background = element_rect(fill="transparent", color=NA),
  panel.background = element_rect(fill="transparent", color=NA),
  panel.grid = element_blank(),
  axis.line = element_line(color = "black"))
P1 = 
  ggplot(DF.CellTable_0, aes(x=sample, y=cell, fill=cluster)) +
  geom_bar(stat="identity", position="fill", color="gray60", linewidth=0.2) +
  geom_text(data=DF.CellsBySample_Labs_0, 
            aes(y=coord_y, 
                #label=case_when(prop>0.08 ~label,
                #                prop>0.03 ~str_replace(label, pattern="\n", replace="  "),
                #                TRUE ~ "")
                label=case_when(prop>0.03 ~ PctLab,
                                TRUE ~ "")),
            lineheight=0.9,
            size=TxtSizeSubCompos*0.6,
            fontface="bold") +
  geom_text(data=dplyr::slice(group_by(DF.CellsBySample_Labs_0, sample), 1), 
            aes(label=formatC(total, big.mark=",")), 
            y=1.03, 
            size=TxtSizeSubCompos*0.5) +
  labs(y="CAF subcluster composition\nwithin each sample (%)", x=NULL, fill=NULL) +
  scale_y_continuous(labels=function(x) x*100,
                     expand=expansion(mult=c(0, 0.05))) +
  scale_x_discrete(labels=setNames(Vec.XeniumID, paste0("TX5K_",names(Vec.XeniumID)))) +
  scale_fill_manual(values=setNames(cols2, paste0("CAF-",names(cols2)))) +
  CommonTheme_bar
P2 = 
  ggplot(DF.CellTable_0, aes(x=cluster, y=cell, fill=sample)) +
  geom_bar(stat="identity", position="fill", color="gray60", linewidth=0.2) +
  geom_text(data=DF.CellsByClust_Labs_0, 
            aes(y=coord_y, 
                label=case_when(prop>0.05 ~ PctLab,
                                TRUE ~ "")),
            lineheight=0.8, 
            size=TxtSizeSubCompos*0.6,
            fontface="bold") +
  geom_text(data=dplyr::slice(group_by(DF.CellsByClust_Labs_0, cluster), 1), 
            aes(label=formatC(total, big.mark=",")), 
            y=1.02, 
            size=TxtSizeSubCompos*0.5) +
  labs(y="Sample composition\nwithin each CAF subcluster (%)", x=NULL, 
       fill="Sample ID") +
  scale_y_continuous(labels=scales::percent_format(),
                     expand=expansion(mult=c(0, 0.05))) +
  scale_fill_manual(labels=setNames(Vec.XeniumID, paste0("TX5K_",names(Vec.XeniumID))),
                    values=c('#6E8FB8','#E5A857','#9FC9C1','#D98585','#8BBF7A','#E9D58A')) +
  CommonTheme_bar
BarPlot = ggarrange(P1,P2,ncol=1, nrow=2, align="hv")
ggsave(BarPlot, 
       file="SubClust_CellComposition.png",
       dpi = 500, width = TxtSizeSubCompos*1.1, height = TxtSizeSubCompos, bg = "transparent") 
# Dot plot
DF.CellsBySample_Labs_1 = 
  DF.CellsBySample_Labs_0 %>% 
  dplyr::select(c(cluster, sample, cell, prop)) %>% 
  dplyr::rename(PropWithinSample = prop)
DF.CellsByClust_Labs_1 = 
  DF.CellsByClust_Labs_0 %>% 
  dplyr::select(c(cluster, sample, cell, prop)) %>% 
  dplyr::rename(PropWithinClust = prop)
DF.CellsTableSub = 
  inner_join(DF.CellsBySample_Labs_1, DF.CellsByClust_Labs_1, 
             by=c("cluster", "sample"))
ggplot(DF.CellsTableSub, aes(x=cluster, y=sample)) +
  geom_point(shape=21, color="black",
             aes(size=PropWithinClust, 
                 fill=PropWithinSample)) +
  labs(size="Cell prop\nwithin clust",
       fill="Cell prop\nwithin sample") +
  scale_fill_gradient(
    low = "#d6e6f2",
    high = "#2c7fb8")
##  2. Dimplot                                   ####
SeuObj.CAF_0 = readRDS(paste0(DirInteg,"Objects/[SeuObj_ExtractedCAF][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                              QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2,".rds") )

CommonTitle.Sub = paste0(NumOfSamples,"case Integration(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")\n",
                         QCInfo,"\n1stClust:PCA50vf2000_ClustDim",Dim1,"Res",formatC(Res1,digits=2,format="f"),
                         "\n2ndClust:PCA50vf2000_ClustDim",Dim2,"Res",formatC(Res2,digits=2,format="f"))
NumOfSubclust = length(unique(SeuObj.CAF_0$seurat_clusters))
CommonCaption.Sub = paste0("Total ",format(ncol(SeuObj.CAF_0), big.mark=",", scientific=F), " CAFs, ",NumOfSubclust, " clusters.")
CommonCAFLab = 
  paste0("CAF-", 0:(NumOfSubclust-1)) %>% 
  setNames(as.character(0:(NumOfSubclust-1)))
DF.UMAP_0 = Embeddings(SeuObj.CAF_0, "umap") %>% as.data.frame() %>%
  rownames_to_column(var = "full_id") %>% 
  dplyr::mutate("cell_id"=str_remove(colnames(SeuObj.CAF_0), pattern="TX5K.*_"),
                "seurat_clusters" = SeuObj.CAF_0$seurat_clusters,
                "sample" = str_remove(SeuObj.CAF_0$orig.ident, pattern="TX5K_"))
set.seed(123)
DF.UMAP_1 = dplyr::slice_sample(DF.UMAP_0, prop=1)
DF.UMAPcenter = group_by(DF.UMAP_1, seurat_clusters) %>% 
  dplyr::summarize(umap_1 = median(umap_1), umap_2 = median(umap_2)) %>% 
  dplyr::mutate(label = CommonCAFLab[seurat_clusters])
DimPlot_ScaleX = scale_x_continuous(limits = range(DF.UMAP_0$umap_1), breaks=c(-10,-5,0,5,10))
DimPlot_ScaleY = scale_y_continuous(limits = range(DF.UMAP_0$umap_2), breaks=c(-10,-5,0,5,10))
Col.Qualitative = c("#ebac23","#b80058","#008cf9","#1f78b4",
                    "#6a3d9a","#006e00","#33a02c",
                    "#b24502","#ff7f00","#d163e6",
                    "#ff9287","#878500","#ffff99","#cab2d6")
cols2 <- c(
  "0"="#4C78A8",  # blue
  "1"="#54A24B",  # green
  "2"="#B279A2",  # purple
  "3"="#9D755D",  # brown
  "4"="#72B7B2",  # teal
  "5"="#E6AF2E",  # softer yellow
  "6"="#F58518",  # orange (少し落ち着かせる)
  "7"="#9C9EDE",  # lavender blue
  "8"="#C00000"   # 強調
)
DimPlot1 = 
  ggplot(DF.UMAP_1, aes(x=umap_1, y=umap_2)) +
  geom_point(data=dplyr::filter(DF.UMAP_1, seurat_clusters != "8"),
             aes(color=seurat_clusters), 
             size=0.4, alpha=0.6) +
  geom_point(data=dplyr::filter(DF.UMAP_1, seurat_clusters == "8"),
             color="#C00000", 
             size=0.6, alpha=0.5) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  DimPlot_ScaleX + DimPlot_ScaleY +
  scale_color_manual(values=cols2)+
  coord_fixed(ratio=1) + 
  theme(plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid = element_blank(),
        legend.position="none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(face="bold", size=18),
        axis.line = element_line(color="black")) # + facet_wrap(. ~ seurat_clusters)
ggsave(
  plot=DimPlot1,
  file=paste0(DirInteg,"[Figure][SubClust_DimPlot2_ForPaper_AllSubclust][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
              "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,
              "_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"), 
  dpi=400, width=6, height=6, bg = "transparent") 
DimPlot.Split.list = list()
for(i in 0:(NumOfSubclust-1)){
  Plot.Split = 
    ggplot(subset(DF.UMAP_1, subset=seurat_clusters==i), 
           aes(x=umap_1, y=umap_2)) +
    geom_point(data=DF.UMAP_1, color="gray90", size=0.4, alpha=0.6) +
    geom_point(color=cols2[i+1], size=0.4, alpha=0.6) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    DimPlot_ScaleX + DimPlot_ScaleY +
    coord_fixed(ratio=1) + 
    theme(plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent", color = NA),
          panel.grid = element_blank(),
          legend.position="none",
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_text(face = "bold"),
          axis.line = element_line(color="black"))
  DimPlot.Split.list[[paste0("CAF_",i)]] = Plot.Split
}
DimPlot2 <- 
  wrap_plots(DimPlot.Split.list,
             ncol = 5,
             nrow = 2) &
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background  = element_rect(fill = "transparent", color = NA))
ggsave(
  plot=DimPlot2,
  file=paste0(DirInteg,"[Figure][SubClust_DimPlot2_ForPaper_SingleSubclust][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
              "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,
              "_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,"2.png"),
  dpi=400, width=13.3, height=5.3, bg="transparent" ) 
##  3. volcano                                   ####
DF.DEGs_0 = read.csv(
  file = paste0(DirInteg, "DEG_table/[SubclustsOfCAFs_DEGs_All][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2,".csv"), 
  row.names=1) %>% dplyr::select(-p_val)
DF.DEGs_1 = DF.DEGs_0 %>% rowwise() %>% 
  dplyr::mutate(minuslog10AdjP = (-1)*log10(p_val_adj),
                minuslog10AdjP_Volc = pmin(minuslog10AdjP, 300),
                colors_by_p = 
                  ifelse(p_val_adj<0.05 & (avg_log2FC>1 | avg_log2FC<(-1) ),
                         paste0("CAF-",cluster), "ns") )
DF.DEGs_1_CAF8 = DF.DEGs_1 %>% subset(cluster==8)
CommonTheme.Volc = theme(legend.position="none",
                         axis.text = element_text(face="bold", color="black"),
                         #axis.text.x = element_text(color="black"),
                         #axis.text.y = element_text(color="black"),
                         axis.title = element_text(face="bold", color="black"),
                         axis.ticks = element_line(color="black"),
                         strip.text = element_text(face = "bold"),
                         strip.background = element_rect(fill=NA, color=NA),
                         plot.background = element_rect(fill="transparent", color=NA),
                         panel.background = element_rect(fill="transparent", color=NA),
                         #panel.grid = element_line(color = "gray90"),
                         panel.grid.major.x = element_blank(),
                         panel.grid.major.y = element_blank(),
                         panel.grid.minor.x = element_blank(),
                         panel.grid.minor.y = element_blank(),
                         panel.border = element_rect(fill="transparent", color = "black"))
BaseLine.h = geom_hline(yintercept=0, color="black")
BaseLine.v = geom_vline(xintercept=0, color="black")
Plot.Volc.Single = 
  ggplot(DF.DEGs_1_CAF8,aes(x=avg_log2FC, y=minuslog10AdjP_Volc)) +
  geom_point(data=subset(DF.DEGs_1_CAF8, subset=colors_by_p=="ns"), 
             color="gray70", size=0.72) +
  geom_point(data=subset(DF.DEGs_1_CAF8, subset=colors_by_p!="ns"), 
             aes(fill=colors_by_p), 
             shape=21, size=1.44, color="black", stroke=0.05) +
  labs(x="log2(Fold change)", y="-log10(Adj.p-value)") +
  scale_y_continuous(limits=c(0, 300),
                     breaks=c(0,100,200,300),
                     labels=c(0,100,200,300))+
  scale_x_continuous(breaks=c(-4, 0, 4),
                     limits=c(-max(abs(DF.DEGs_1_CAF8$avg_log2FC)),
                              max(abs(DF.DEGs_1_CAF8$avg_log2FC)))) +
  scale_fill_manual(values=c("CAF-8"=cols2[["8"]])) +
  BaseLine.h + BaseLine.v + CommonTheme.Volc +
  theme(aspect.ratio = 4.4/4.82)
ggsave(Plot.Volc.Single,
       file=paste0(DirInteg, "/[Figure][SubclustOfCAF_Volcano_Single_ForPaper][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                   QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2,".png"), 
       width=3, height=4.5, dpi=600)
Plot.Volc.Muptiple = 
  ggplot(DF.DEGs_1,aes(x=avg_log2FC, y=minuslog10AdjP_Volc)) +
  geom_point(data=subset(DF.DEGs_1, subset=colors_by_p=="ns"), 
             color="gray70", size=0.3) +
  geom_point(data=subset(DF.DEGs_1, subset=colors_by_p!="ns"), 
             aes(fill=colors_by_p), 
             shape=21, size=0.9, stroke=0.05, color="black") +
  labs(x="log2(Fold change)", y="-log10(Adj.p-value)") +
  scale_y_continuous(limits=c(0, 300))+
  scale_x_continuous(breaks=c(-6, -3, 0, 3, 6)) +
  scale_fill_manual(values=
                      setNames(cols2, paste0("CAF-",names(cols2)))) +
  facet_grid( cluster ~ .) +
  BaseLine.h + BaseLine.v + CommonTheme.Volc
ggsave(Plot.Volc.Muptiple,
       file=paste0(DirInteg, "/[Figure][SubclustOfCAF_Volcano_Multiple_ForPaper][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                   QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2,".png"), 
       width=4.5, height=4.5, dpi=600)

##  4. Cluster-level GSEA                        ####

library(org.Hs.eg.db)
library(msigdbr)
library(clusterProfiler)
# gene sets
H = clusterProfiler::read.gmt("/Volumes/PortableSSD/[4]myR/[1]dataset/h.all.v2024.1.Hs.entrez.gmt")
CategoryC2 = msigdbr(species="Homo sapiens",  category="C2")
Additional_0 = 
  data.frame(term = "Buffa_Hypoxia", symbol = Vec.BuffaOrig) %>% 
  rbind(data.frame(term="Winter_Hypoxia", symbol=Vec.WinterOrig)) %>% 
  rbind(data.frame(term="KEGG_Cell_cycle", symbol=CategoryC2 %>% filter(gs_name == "KEGG_CELL_CYCLE") %>% pull(gene_symbol) %>% unique())) %>% 
  rbind(data.frame(term="myCAF_signature", symbol=DF.Markers$Gene[DF.Markers$Classification1=="myCAF"])) %>% 
  rbind(data.frame(term="iCAF_signature", symbol=DF.Markers$Gene[DF.Markers$Classification1=="iCAF"]))
Vec.AdditionalEntrez_0 = clusterProfiler::bitr(
  Additional_0$symbol, fromType="SYMBOL", toType="ENTREZID", 
  OrgDb=org.Hs.eg.db, drop=FALSE)
Vec.AdditionalEntrez_1 = 
  Vec.AdditionalEntrez_0$ENTREZID %>% 
  setNames(Vec.AdditionalEntrez_0$SYMBOL)
Additional_1 = Additional_0 %>% 
  dplyr::mutate(gene = Vec.AdditionalEntrez_1[symbol])
Additional_2 = dplyr::select(Additional_1, c(term, gene))
# read DEG table
DEG_0 = read.csv(file = paste0(DirInteg, "DEG_table/[SubclustsOfCAFs_DEGs_All][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                               QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2,".csv"), 
                 row.names=1)
DF.Sym_and_En_0 = clusterProfiler::bitr(DEG_0$gene, fromType="SYMBOL", toType="ENTREZID", 
                                        OrgDb=org.Hs.eg.db, drop=FALSE)
# manually fill NA cell, failed to map
DF.Sym_and_En_0[DF.Sym_and_En_0$SYMBOL=="H3F3B", "ENTREZID"]=3021
DF.Sym_and_En_0[DF.Sym_and_En_0$SYMBOL=="ARNTL", "ENTREZID"]=406
DF.Sym_and_En_0[DF.Sym_and_En_0$SYMBOL=="H2AFX", "ENTREZID"]=3014
DF.Sym_and_En_0[DF.Sym_and_En_0$SYMBOL=="TMEM173", "ENTREZID"]=340061 # STING1
DF.Sym_and_En_0[DF.Sym_and_En_0$SYMBOL=="WARS", "ENTREZID"]=7453
DF.Sym_and_En_0[DF.Sym_and_En_0$SYMBOL=="SPATA5L1", "ENTREZID"]=79029
DF.Sym_and_En_0[DF.Sym_and_En_0$SYMBOL=="SPATA5", "ENTREZID"]=166378
DF.Sym_and_En_0[DF.Sym_and_En_0$SYMBOL=="EPRS", "ENTREZID"]=2058
DF.Sym_and_En_0[DF.Sym_and_En_0$SYMBOL=="SLC9A3R2", "ENTREZID"]=9351
DF.Sym_and_En_0[DF.Sym_and_En_0$SYMBOL=="CBWD2", "ENTREZID"]=150472
DF.Sym_and_En_0[DF.Sym_and_En_0$SYMBOL=="HSPB11", "ENTREZID"]=51668
DF.Sym_and_En_0[DF.Sym_and_En_0$SYMBOL=="FAM126B", "ENTREZID"]=285172
DF.Sym_and_En_0[DF.Sym_and_En_0$SYMBOL=="SLC9A3R1", "ENTREZID"]=9368
DF.Sym_and_En_0[DF.Sym_and_En_0$SYMBOL=="GBA", "ENTREZID"]=2629
DF.Sym_and_En_0[DF.Sym_and_En_0$SYMBOL=="PHB", "ENTREZID"]=5245
DF.Sym_and_En_0[DF.Sym_and_En_0$SYMBOL=="EEF1AKNMT", "ENTREZID"]=51603
DF.Sym_and_En_0[DF.Sym_and_En_0$SYMBOL=="CBSL", "ENTREZID"]=102724560
DF.Sym_and_En_0[DF.Sym_and_En_0$SYMBOL=="CCDC113", "ENTREZID"]=29070
DF.Sym_and_En_0[DF.Sym_and_En_0$SYMBOL=="TDGF1", "ENTREZID"]=6997
DF.Sym_and_En_0[DF.Sym_and_En_0$SYMBOL=="GPR1", "ENTREZID"]=2825
DF.Sym_and_En_0[DF.Sym_and_En_0$SYMBOL=="CYHR1", "ENTREZID"]=50626
DF.Sym_and_En_0[DF.Sym_and_En_0$SYMBOL=="CD3EAP", "ENTREZID"]=10849
DF.Sym_and_En_0[DF.Sym_and_En_0$SYMBOL=="DUSP13", "ENTREZID"]=128854680  #DUSP13A
DF.Sym_and_En_0[DF.Sym_and_En_0$SYMBOL=="DDX58", "ENTREZID"]=23586
DF.Sym_and_En_0[DF.Sym_and_En_0$SYMBOL=="BVES", "ENTREZID"]=11149
DF.Sym_and_En_0[DF.Sym_and_En_0$SYMBOL=="ZBED9", "ENTREZID"]=114821
DF.Sym_and_En_0[DF.Sym_and_En_0$SYMBOL=="ACPP", "ENTREZID"]=55
DF.Sym_and_En_1 = DF.Sym_and_En_0 %>% 
  dplyr::filter(DF.Sym_and_En_0$ENTREZID != 100124696) %>% # TEC is duplicated
  set_colnames(c("gene", "EntrezID"))
DEG_1 = left_join(DEG_0, DF.Sym_and_En_1, by="gene")
DEG_2 = group_by(DEG_1, cluster) %>% 
  dplyr::arrange(desc(avg_log2FC), .by_group=TRUE)

# Run GSEA() at cluster-level
NumObClust=DEG_0$cluster %>% unique() %>% length()
DF.ResGSEA = data.frame()
List.ResGSEA = list()
for(TargetClust in 0:(NumObClust-1)){
  DEG_TargetClust = subset(DEG_2, subset=cluster==TargetClust)
  GeneList = DEG_TargetClust$avg_log2FC
  names(GeneList) = DEG_TargetClust$EntrezID
  set.seed(1234)
  GeneList = GeneList + rnorm(length(GeneList), mean=0 , sd=1e-10)
  Res.Pos=GSEA(
    geneList=GeneList,
    minGSSize=10,
    TERM2GENE=rbind(H,Additional_2),
    pvalueCutoff=1,
    verbose=FALSE, eps=0,
    scoreType = "pos")
  Res.Neg=GSEA(
    geneList=GeneList,
    minGSSize=10,
    TERM2GENE=rbind(H,Additional_2),
    pvalueCutoff=1,
    verbose=FALSE, eps=0,
    scoreType = "neg")
  DF.ResGSEA_TargetClust_0 = rbind(cbind(Res.Pos@result, "Sign"="UP"), 
                                   cbind(Res.Neg@result, "Sign"="DOWN"))
  DF.ResGSEA_TargetClust_1 = 
    if(nrow(DF.ResGSEA_TargetClust_0) > 0){
      DF.ResGSEA_TargetClust_0 %>% 
        cbind("Subclust"=TargetClust) %>% 
        dplyr::mutate(Significance=
                        case_when(p.adjust<0.001 ~ "***",
                                  p.adjust<0.01 ~ "**",
                                  p.adjust<0.05 ~ "*",
                                  TRUE ~ "ns")) %>% 
        dplyr::relocate(Significance, .before=ID) %>% 
        dplyr::relocate(Sign, .before=Significance) %>%  
        dplyr::relocate(Subclust, .before=Sign) 
    } else { NULL }
  rownames(DF.ResGSEA_TargetClust_1) = NULL
  DF.ResGSEA = rbind(DF.ResGSEA, DF.ResGSEA_TargetClust_1)
  List.ResGSEA[[paste0(as.character(TargetClust),"_UP")]] = Res.Pos
  List.ResGSEA[[paste0(as.character(TargetClust),"_DOWN")]] = Res.Neg
}
write.csv(DF.ResGSEA,
          file=paste0(DirInteg, "/[DataTable_ExtractedCAF_PseudoBulkGSEA][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                      QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2,".csv"),
          row.names=F)
saveRDS(List.ResGSEA,
        file=paste0(DirInteg, "/[DataRDS_ExtractedCAF_PseudoBulkGSEA][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                    QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2,".rds"))

# visualize, Bar plot
library(ggtext)
DF.ResGSEA_0 = read.csv(
  file=paste0(DirInteg, "/[DataTable_ExtractedCAF_PseudoBulkGSEA][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
              QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2,".csv"),
  row.names=NULL) %>% 
  dplyr::mutate(Subclust = case_when(Sign=="UP" ~ paste0("CAF",Subclust,"_UP"),
                                     Sign=="DOWN" ~ paste0("CAF",Subclust,"_DOWN")))
NumOfSubclust = (DF.ResGSEA_0$Subclust %>% unique() %>% length())/2

DF.ResGSEA_0.UP.Top10 = DF.ResGSEA_0 %>% 
  dplyr::mutate(Subclust = str_remove(Subclust, pattern="_.*")) %>% 
  subset(subset=Sign=="UP") %>% 
  group_by(Subclust) %>% 
  slice_max(order_by=NES, n=10) %>% 
  dplyr::select(Subclust, Sign, Significance, ID, NES, p.adjust)
DF.ResGSEA_0.DOWN.Top10 = DF.ResGSEA_0 %>% 
  dplyr::mutate(Subclust = str_remove(Subclust, pattern="_.*")) %>% 
  subset(subset=Sign=="DOWN") %>% 
  group_by(Subclust) %>% 
  slice_min(order_by=NES, n=10) %>% 
  dplyr::select(Subclust, Sign, Significance, ID, NES, p.adjust)
DF.ResGSEA_1 = 
  rbind(DF.ResGSEA_0.UP.Top10, DF.ResGSEA_0.DOWN.Top10) %>% 
  group_by(Subclust) %>% 
  dplyr::arrange(NES, .by_group=TRUE) %>% 
  dplyr::mutate(
    ID_plot = fct_inorder(paste0(as.character(ID), "___", Subclust))) %>% 
  ungroup()
List.DF = DF.ResGSEA_1 %>% 
  group_split(Subclust)

m = as.numeric(quantile(abs(DF.ResGSEA_1$NES), 0.975, na.rm=TRUE))
List.Plot = list()
for(i in 1:NumOfSubclust){
  DF.OneForPlot = List.DF[[i]] %>% 
    dplyr::arrange(NES) %>% 
    dplyr::mutate(ID_plot = fct_inorder(paste0(as.character(ID),"___",Subclust)),
                  NES_fill = ifelse(Significance=="ns", NA, NES),
                  ColumnBorder = ifelse(Significance=="ns", "NS", "Signif"),
                  LabelCore = clean_gs_label_ForFig(as.character(ID)),
                  AxisLabel = ifelse(
                    Significance == "ns",
                    paste0("<span style='color:gray60'>", LabelCore, "</span>"),
                    LabelCore))
  Plot.One = 
    ggplot(DF.OneForPlot, 
           aes(x=NES, y=ID_plot, fill=NES_fill, color=ColumnBorder)) +
    geom_col() +
    labs(title = paste0("CAF-",(i-1)),
         y=NULL, fill="NES") +
    scale_x_continuous(limits = c(-max(abs(DF.ResGSEA_1$NES), na.rm=TRUE),
                                  max(abs(DF.ResGSEA_1$NES), na.rm=TRUE)),
                       breaks = c(-10, -5, 0, 5, 10)) + 
    scale_y_discrete(
      labels = setNames(DF.OneForPlot$AxisLabel,
                        DF.OneForPlot$ID_plot)) +
    #scale_fill_distiller(type="div", palette="RdBu", direction=-1,
    scale_fill_gradient2(low = Vec.Cols.Diverging[1], 
                         mid = Vec.Cols.Diverging[2], 
                         high = Vec.Cols.Diverging[3],
                         limits = c(-m, m),
                         oob=scales::squish,
                         na.value = "gray90",
                         guide=guide_colorbar(#direction="horizontal",
                           #title.position="top",
                           #title.hjust=0.5,
                           frame.colour="black",
                           ticks.colour="black")) +
    scale_color_manual(values=c("Signif"="black",
                                "NS"=NA),
                       guide="none") +
    geom_vline(xintercept=0) +
    theme(text = element_text(face="bold"),
          axis.text.y = ggtext::element_markdown(face = "bold"),
          axis.text.x = element_text(color = "black"),
          axis.text = element_text(color="black"),
          legend.key = element_blank(),
          plot.background = element_rect(fill="transparent", color=NA),
          panel.background = element_rect(fill="white", color=NA),
          panel.grid = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"))
  List.Plot[[i]] = Plot.One
}
Plot.All = 
  (wrap_plots(List.Plot, nrow = 3) +
     plot_layout(guides = "collect")) &
  theme(
    legend.position = "bottom",
    legend.background = element_rect(fill="transparent", color=NA),
    plot.background = element_rect(fill="transparent", color=NA),
    panel.background = element_rect(fill="white", color=NA))
ggsave(plot=Plot.All,
       filename = "subclust_cluster_level_gsea3.png", 
       width=17.5, height=11.0, dpi=500, 
       bg="transparent")


##  5. DotPlot                                   ####
Seu.CAF_0 = readRDS(paste0(DirInteg,"Objects/[SeuObj_ExtractedCAF][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                           QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2,".rds") )
CommonTitle.Sub = paste0(NumOfSamples,"case Integration(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")\n",
                         QCInfo,"\n1stClust:PCA50vf2000_ClustDim",Dim1,"Res",formatC(Res1,digits=2,format="f"),
                         "\n2ndClust:PCA50vf2000_ClustDim",Dim2,"Res",formatC(Res2,digits=2,format="f"))
NumOfSublust = length(unique(Seu.CAF_0$seurat_clusters))
CommonCaption.Sub = paste0("Total ",format(ncol(Seu.CAF_0), big.mark=",", scientific=F), " CAFs, ",NumOfSublust, " clusters.")
Seu.CAF_0$seurat_clusters = factor(Seu.CAF_0$seurat_clusters, levels=c(4,8,0,1,2,5,3,7,6))
Idents(Seu.CAF_0) = "seurat_clusters"
# 3. Dot plot : CAF markers
DotPlotTheme = theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                     axis.text.y = element_text(hjust = 0),
                     axis.text = element_text(face = "bold", size = 19),
                     legend.text = element_text(face = "bold", size = 15),
                     legend.title = element_text(face = "bold", size = 15),
                     axis.title = element_text(face = "bold", size = 17),
                     plot.background = element_rect(fill="transparent", color=NA),
                     panel.background = element_rect(fill="white", color=NA))
DotPlotLgdLab = guides(size = guide_legend(title = "Percent\nExpressed"),
                       color = guide_colorbar(title = "Scaled\nAverage\nExpression",
                                              direction="vertical",
                                              frame.colour="black",
                                              ticks.colour="black"))
DotPlotColScale <- scale_color_gradientn(
  limits=c(2.5, -2.5),
  colors=c(colors = c("#045a8d","#2b8cbe","#74a9cf",
                      "gray90",
                      "#fcbba1","#ef3b2c","#99000d")))
DotPlot <-
  DotPlot(Seu.CAF_0, cluster.idents=FALSE, 
          dot.scale = 10, 
          features=c(#DF.Markers$Gene[DF.Markers$Classification1=="myCAF"],
            #DF.Markers$Gene[DF.Markers$Classification1=="iCAF"],
            "COL5A2","COL10A1","FN1","ACTA2","POSTN","SULF1","THBS2","MYL9","FAP","ADAMTS12","MMP11","MYLK",
            "TPM1","TAGLN","CCN2",
            "BIRC5","TOP2A","CENPF","NUSAP1","PTTG1","STMN1","ENO1",
            #"HOPX",
            "NOTCH1","JAG1","DTX1","WNT5A",
            #"CXCL8","IL6",
            "CXCL12","PDGFRA","DPT",#"IL1R1",
            "PLA2G2A","VEGFA","CCL2","CFH","C4A","C4B",
            "CD74","HLA-DRA","VIM","APOD") ) &
  geom_point(aes(size=pct.exp), shape=21, stroke=0.5, color="gray20") &
  annotate(geom="rect", xmin=0.5, xmax=12.5, ymin=0.5, ymax=4.5, fill="orange", color="orange", alpha=0.03) &
  annotate(geom="rect", xmin=12.5, xmax=15.5, ymin=0.5, ymax=1.5, fill="orange", color="orange", alpha=0.03) &
  annotate(geom="rect", xmin=15.5, xmax=22.5, ymin=1.5, ymax=2.5, fill="orange", color="orange", alpha=0.03) &
  annotate(geom="rect", xmin=22.5, xmax=26.5, ymin=2.5, ymax=3.5, fill="orange", color="orange", alpha=0.03) &
  annotate(geom="rect", xmin=26.5, xmax=29.5, ymin=4.5, ymax=7.5, fill="orange", color="orange", alpha=0.03) &
  annotate(geom="rect", xmin=29.5, xmax=32.5, ymin=4.5, ymax=5.5, fill="orange", color="orange", alpha=0.03) &
  annotate(geom="rect", xmin=32.5, xmax=35.5, ymin=5.5, ymax=6.5, fill="orange", color="orange", alpha=0.03) &
  annotate(geom="rect", xmin=35.5, xmax=37.5, ymin=7.5, ymax=8.5, fill="orange", color="orange", alpha=0.03) &
  annotate(geom="rect", xmin=37.5, xmax=39.5, ymin=8.5, ymax=9.5, fill="orange", color="orange", alpha=0.03) &
  scale_y_discrete(labels=c("0"="CAF-0: ECM remodeling, NOTCH/WNT-high",
                            "1"="CAF-1: ECM remodeling, weak markers",
                            "2"="CAF-2: Inflammatory, PLA2G2A-high",
                            "3"="CAF-3: Inflammatory, weak markers",
                            "4"="CAF-4: ECM remodeling, ACTA2-high",
                            "5"="CAF-5: Inflammatory, CFH/C4-high",
                            "6"="CAF-6: Peri-neural fibroblast",
                            "7"="CAF-7: apCAF",
                            "8"="CAF-8: ECM remodeling, Proliferating"),
                   expand=expansion(mult=0)) &
  scale_x_discrete(expand=expansion(mult=0)) &
  labs(y=NULL, x=NULL) & 
  DotPlotLgdLab & RotatedAxis() & DotPlotTheme & DotPlotColScale
DotPlot
ggsave(DotPlot,
       file=paste0(DirInteg,"[Figure][SubClust_DotPlot_FinalMarkers][",
                   NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                   "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"), 
       dpi=400, width=19, height=5, bg="transparent")
##  6. Pseudo-time analysis                      ####
library(monocle)
library(Seurat)
library(Matrix)
library(igraph)
library(clusterProfiler)
library(org.Hs.eg.db)
library(grid)
library(gridExtra)
library(ggridges)
library(msigdbr)

DataSetName = "rmCAF6"
Seu.Obj <- 
  readRDS(paste0(
    DirInteg,"Objects/[SeuObj_ExtractedCAF][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
    QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2,".rds") ) %>% 
  subset(seurat_clusters != 6)
cds <- readRDS(
  file=paste0(DirInteg,"/Pseudotime/[CellDataSet][",DataSetName,"_HVG",n_HVGs,"][",
              NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
              "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".rds"))
Num.States <- length(unique(pData(cds)$State))
DF.Coord_0 <- data.frame(
  Component1 = reducedDimS(cds)[1, ],
  Component2 = reducedDimS(cds)[2, ])
DF.Coord_1 <- DF.Coord_0 %>% 
  cbind(pData(cds)[rownames(DF.Coord_0), 
                   c("orig.ident", "seurat_clusters", "Pseudotime", "State")]) %>% 
  dplyr::mutate(
    State = paste0("State ",State),
    seurat_clusters = paste0("CAF-", seurat_clusters),
    seurat_clusters = factor(seurat_clusters, levels = paste0("CAF-", Vec.SubclustOrder)))
DF.Node <- 
  t(reducedDimK(cds)) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "node")　%>% 
  dplyr::rename(Component1 = V1, Component2 = V2)
DF.Edge <- 
  igraph::as_data_frame(
    minSpanningTree(cds),
    what = "edges")
DF.EdgeCoord <- 
  DF.Edge %>% 
  left_join(
    DF.Node,
    by = c("from" = "node")) %>% 
  dplyr::rename(x = Component1, y = Component2) %>% 
  left_join(
    DF.Node,
    by = c("to" = "node")) %>% 
  dplyr::rename(xend = Component1, yend = Component2)
CommonTheme <- 
  theme(
    axis.title = element_text(face = "bold", color = "black", size = 10),
    legend.title = element_text(face = "bold", color = "black", size = 10),
    axis.text.x = element_text(face = "bold", color = "black", size = 10),
    axis.text.y = element_text(face = "bold", color = "black", size = 10),
    legend.text = element_text(face = "bold", color = "black", size = 10),
    plot.background = element_rect(fill = NA, color = NA),
    panel.background = element_rect(fill = NA, color = NA),
    legend.background = element_rect(fill = NA, color = NA),
    legend.key = element_blank(),
    panel.border = element_rect(fill = NA, color = "black"),
    panel.grid = element_blank())
CommonLabs <- 
  labs(x = "Component 1", y = "Component 2", color = NULL)
p2 <- 
  rbind(
    dplyr::filter(DF.Coord_1, seurat_clusters != "CAF-8"),
    dplyr::filter(DF.Coord_1, seurat_clusters == "CAF-8")) %>% 
  ggplot(aes(x = Component1, y = Component2)) +
  geom_segment(data = DF.EdgeCoord,
               aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_point(aes(color = seurat_clusters),
             size = 0.8) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_color_manual(
    values = setNames(as.character(cols2), paste0("CAF-", 0:8))) +
  CommonLabs + CommonTheme +
  theme(aspect.ratio = 0.845)
p3 <- 
  ggplot(DF.Coord_1, aes(x = Component1, y = Component2)) +
  geom_segment(data = DF.EdgeCoord,
               aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_point(aes(color = Pseudotime)) + 
  labs(x = "Component 1", y = "Component 2", color = "Pseudotime") + 
  scale_color_viridis(
    option = "inferno",
    guide=guide_colorbar(direction="vertical",
                         title.position="top",
                         title.hjust=0.5,
                         frame.colour="black",
                         ticks.colour="black"),
    breaks = c(0, 4,8)) +
  CommonTheme +
  theme(aspect.ratio = 0.845)
DF.Component_0 <- 
  table(
    State = pData(cds)$State,
    CAF = pData(cds)$seurat_clusters) %>% 
  as.data.frame() %>% 
  dplyr::mutate(State = paste0("State ", State),
                CAF = paste0("CAF-", CAF),
                CAF = factor(CAF, levels = paste0("CAF-", Vec.SubclustOrder)))
p4 <- 
  ggplot(DF.Component_0,
         aes(x = State, y = Freq, fill = CAF))+
  geom_bar(stat="identity", color = "gray50", linewidth = 0.3) +
  labs(x = NULL, y = "Number of CAF cells", fill = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.03)),
                     breaks = c(0, 500, 1000)) +
  scale_fill_manual(values = setNames(as.character(cols2), paste0("CAF-", 0:8))) +
  CommonTheme + 
  theme(axis.text.x = element_text(face = "bold", color = "black", 
                                   angle = 90, hjust = 1, vjust = 0.5, size = 10)) +
  theme(aspect.ratio = 0.615)
p4
Num.WidthsRight <- 8
pAligned <- 
  (p2|p4|p3) +
  plot_layout(widths = c(5, 8, 5)) &
  theme(plot.background = element_rect(fill = NA, color = NA))
#pAligned
ggsave(
  pAligned,
  file=paste0(DirInteg,"/Pseudotime/[BasePlots_final][",DataSetName,"_HVG",n_HVGs,"][",
              NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
              "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"),
  dpi=400, width = 13, height = 5, bg="transparent")
DF.Coord_1 %>% 
  rownames_to_column(var = "cell_id") %>% 
  write.csv(
    file=paste0(DirInteg,"/Pseudotime/[MetaData][",DataSetName,"_HVG",n_HVGs,"][",
                NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".csv"),
    row.names = F)
DF.Component_1 <- 
  DF.Component_0 %>% 
  pivot_wider(names_from = State, values_from = Freq) %>% 
  column_to_rownames(var = "CAF")
DF.Component_2 <-
  DF.Component_1 %>% 
  dplyr::mutate(
    n_total = apply(DF.Component_1, MARGIN = 1, FUN = sum))
write.csv(
  DF.Component_2,
  file=paste0(DirInteg,"/Pseudotime/[ClusterComposition][",DataSetName,"_HVG",n_HVGs,"][",
              NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
              "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".csv"),
  row.names = F)




Func.CellDataSet <- function(SeuObj, DataSetName, n_HVGs){
  # 1. Down sampling
  SeuObj@graphs <- list()
  Idents(SeuObj) <- "seurat_clusters"
  set.seed(1234)
  SeuObj_ds <- subset(
    SeuObj,
    downsample = 500
  )
  SeuObj_ds <- JoinLayers(SeuObj_ds)
  ExprMat <- GetAssayData(
    SeuObj_ds,
    assay = "RNA",
    layer = "counts"
  )
  
  # 2. Cell metadata
  AnnoDF.Phenotype <- new(
    "AnnotatedDataFrame",
    data = SeuObj_ds@meta.data
  )
  
  # 3. Gene metadata
  DF.FeatureData <- data.frame(
    gene_short_name = rownames(ExprMat),
    row.names = rownames(ExprMat))
  AnnoDF.Feature <- new(
    "AnnotatedDataFrame",
    data = DF.FeatureData)
  
  # 4. Monocle2 CellDataSet
  cds <- newCellDataSet(
    ExprMat,
    phenoData = AnnoDF.Phenotype,
    featureData = AnnoDF.Feature,
    expressionFamily = negbinomial.size())
  
  # 5. Size factor / dispersion
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  DispTable <- dispersionTable(cds)
  
  # 6. Define ordering genes
  SeuObj_ds <- FindVariableFeatures(SeuObj_ds, nfeatures = n_HVGs)
  Vec.HVGs <- VariableFeatures(SeuObj_ds)
  OrderingGenes <- intersect(Vec.HVGs, rownames(cds))
  
  cds <- setOrderingFilter(cds, OrderingGenes)
  cds <- reduceDimension(
    cds,
    max_components = 2,
    method = "DDRTree"
  )
  
  cds <- orderCells(cds)
  #Warning messages:
  #1: In dfs(graph = graph, root = root, mode = mode, unreachable = unreachable,  :
  #   Argument `neimode' is deprecated; use `mode' instead
  #2: In dfs(graph = graph, root = root, mode = mode, unreachable = unreachable,  :
  #   Argument `neimode' is deprecated; use `mode' instead
  List.Summary <- list(
    "mean_expression" = summary(DispTable$mean_expression),
    "dispersion_empirical" = summary(DispTable$dispersion_empirical),
    "dispersion_fit" = summary(DispTable$dispersion_fit),
    "NA" = colSums(is.na(DispTable)),
    "TheNumberOfOrderingGenes" = length(OrderingGenes)
  )
  saveRDS(cds,
         file=paste0(DirInteg,"/Pseudotime/[CellDataSet][",DataSetName,"_HVG",n_HVGs,"][",
                     NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                     "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".rds"))
  saveRDS(List.Summary,
          file=paste0(DirInteg,"/Pseudotime/[SummaryList][",DataSetName,"_HVG",n_HVGs,"][",
                      NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                      "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".rds"))
  return(List.Summary)
}
Func.CellDataSet(SeuObj = Seu.CAF_0, DataSetName = "WholeCAF", n_HVGs = 2000)
Func.CellDataSet(SeuObj = Seu.CAF_rmCAF6, DataSetName = "rmCAF6", n_HVGs = 5000)
Func.CellDataSet(SeuObj = Seu.CAF_rmCAF6, DataSetName = "rmCAF6", n_HVGs = 4500)
Func.CellDataSet(SeuObj = Seu.CAF_rmCAF6, DataSetName = "rmCAF6", n_HVGs = 4000)
Func.CellDataSet(SeuObj = Seu.CAF_rmCAF6, DataSetName = "rmCAF6", n_HVGs = 3500)
Func.CellDataSet(SeuObj = Seu.CAF_rmCAF6, DataSetName = "rmCAF6", n_HVGs = 3000)
Func.CellDataSet(SeuObj = Seu.CAF_rmCAF6, DataSetName = "rmCAF6", n_HVGs = 2500)
Func.CellDataSet(SeuObj = Seu.CAF_rmCAF6, DataSetName = "rmCAF6", n_HVGs = 2000)
Func.CellDataSet(SeuObj = Seu.CAF_rmCAF6, DataSetName = "rmCAF6", n_HVGs = 1500)
Func.CellDataSet(SeuObj = Seu.CAF_rmCAF6, DataSetName = "rmCAF6", n_HVGs = 1000)
 
Func.BasePlots <- function(DataSetName, n_HVGs){
  cds <- readRDS(
    file=paste0(DirInteg,"/Pseudotime/[CellDataSet][",DataSetName,"_HVG",n_HVGs,"][",
                NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".rds")
  )
  Num.States <- length(unique(pData(cds)$State))
  DF.Coord_0 <- data.frame(
    Component1 = reducedDimS(cds)[1, ],
    Component2 = reducedDimS(cds)[2, ])
  DF.Coord_1 <- DF.Coord_0 %>% 
    cbind(pData(cds)[rownames(DF.Coord_0), 
                     c("orig.ident", "seurat_clusters", "Pseudotime", "State")]) %>% 
    dplyr::mutate(State = paste0("State ",State),
                  seurat_clusters = paste0("CAF-", seurat_clusters),
                  seurat_clusters = factor(seurat_clusters, levels = paste0("CAF-", c(4,8,0,1,2,5,3,7,6))))
  DF.Node <- 
    t(reducedDimK(cds)) %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "node")　%>% 
    dplyr::rename(Component1 = V1, Component2 = V2)
  DF.Edge <- 
    igraph::as_data_frame(
      minSpanningTree(cds),
      what = "edges")
  DF.EdgeCoord <- 
    DF.Edge %>% 
    left_join(
      DF.Node,
      by = c("from" = "node")) %>% 
    dplyr::rename(x = Component1, y = Component2) %>% 
    left_join(
      DF.Node,
      by = c("to" = "node")) %>% 
    dplyr::rename(xend = Component1, yend = Component2)
  CommonTheme <- 
    theme(
      axis.title = element_text(face = "bold", color = "black", size = 10),
      legend.title = element_text(face = "bold", color = "black", size = 10),
      axis.text.x = element_text(face = "bold", color = "black", size = 10),
      axis.text.y = element_text(face = "bold", color = "black", size = 10),
      legend.text = element_text(face = "bold", color = "black", size = 10),
      plot.background = element_rect(fill = NA, color = NA),
      panel.background = element_rect(fill = NA, color = NA),
      legend.background = element_rect(fill = NA, color = NA),
      legend.key = element_blank(),
      panel.border = element_rect(fill = NA, color = "black"),
      panel.grid = element_blank())
  CommonLabs <- 
    labs(x = "Component 1", y = "Component 2", color = NULL)
  p1 <- 
    ggplot(DF.Coord_1, aes(x = Component1, y = Component2)) +
    geom_segment(data = DF.EdgeCoord,
                 aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_point(aes(color = State)) +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    CommonLabs + CommonTheme
  p2 <- 
    rbind(
      dplyr::filter(DF.Coord_1, seurat_clusters != "CAF-8"),
      dplyr::filter(DF.Coord_1, seurat_clusters == "CAF-8")) %>% 
    ggplot(aes(x = Component1, y = Component2)) +
    geom_segment(data = DF.EdgeCoord,
                 aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_point(aes(color = seurat_clusters),
               size = 0.8) +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    scale_color_manual(
      values = setNames(as.character(cols2), paste0("CAF-", 0:8))) +
    CommonLabs + CommonTheme
  p3 <- 
    ggplot(DF.Coord_1, aes(x = Component1, y = Component2)) +
    geom_segment(data = DF.EdgeCoord,
                 aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_point(aes(color = Pseudotime)) + 
    labs(x = "Component 1", y = "Component 2", color = "Pseudotime") + 
    scale_color_viridis(
      option = "inferno",
      guide=guide_colorbar(direction="vertical",
                           title.position="top",
                           title.hjust=0.5,
                           frame.colour="black",
                           ticks.colour="black")) +
    CommonTheme
  DF.Component_0 <- 
    table(
      State = pData(cds)$State,
      CAF = pData(cds)$seurat_clusters) %>% 
    as.data.frame() %>% 
    dplyr::mutate(State = paste0("State ", State),
                  CAF = paste0("CAF-", CAF))
  p4 <- 
    ggplot(DF.Component_0,
           aes(x = State, y = Freq, fill = CAF))+
    geom_bar(stat="identity", color = "gray50", linewidth = 0.3) +
    labs(x = NULL, y = "Subcluster composition\nwithin each state", fill = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    scale_fill_manual(values = setNames(as.character(cols2), paste0("CAF-", 0:8))) +
    CommonTheme + 
    theme(axis.text.x = element_text(face = "bold", color = "black", angle = 45, hjust = 1, size = 10))
  p5 <- 
    ggplot(DF.Component_0,
           aes(x = CAF, y = Freq, fill = State))+
    geom_bar(stat="identity", color = "gray50", linewidth = 0.3) +
    labs(x = NULL, y = "State composition\nwithin each CAF subcluster", fill = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    CommonTheme + 
    theme(axis.text.x = element_text(face = "bold", color = "black", angle = 45, hjust = 1, size = 10))
  Num.WidthsRight <- 
    case_when(Num.States > 8 ~ 10,
              Num.States > 5 ~ 8,
              Num.States > 4 ~ 5,
              TRUE ~ 3)
  pAligned <- 
    ((p1/p2/p3)|(p4/p5)) +
    plot_layout(widths = c(5, Num.WidthsRight)) &
    theme(
      plot.background = element_rect(fill = NA, color = NA))
  pAligned
  ggsave(
    pAligned,
    file=paste0(DirInteg,"/Pseudotime/[BasePlots][",DataSetName,"_HVG",n_HVGs,"][",
                NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"),
    dpi=400, width= 5 + Num.WidthsRight, height=10, bg="transparent")
  DF.Coord_1 %>% 
    rownames_to_column(var = "cell_id") %>% 
    write.csv(
      file=paste0(DirInteg,"/Pseudotime/[MetaData][",DataSetName,"_HVG",n_HVGs,"][",
                  NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                  "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".csv"),
      row.names = F)
  DF.Component_1 <- 
    DF.Component_0 %>% 
    pivot_wider(names_from = State, values_from = Freq) %>% 
    column_to_rownames(var = "CAF")
  DF.Component_2 <-
    DF.Component_1 %>% 
    dplyr::mutate(
      n_total = apply(DF.Component_1, MARGIN = 1, FUN = sum))
  write.csv(
    DF.Component_2,
    file=paste0(DirInteg,"/Pseudotime/[ClusterComposition][",DataSetName,"_HVG",n_HVGs,"][",
                NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".csv"),
    row.names = F)
}
Func.BasePlots(DataSetName = "WholeCAF", n_HVGs = 2000)
Func.BasePlots(DataSetName = "rmCAF6", n_HVGs = 5000)
Func.BasePlots(DataSetName = "rmCAF6", n_HVGs = 4500)
Func.BasePlots(DataSetName = "rmCAF6", n_HVGs = 4000)
Func.BasePlots(DataSetName = "rmCAF6", n_HVGs = 3500)
Func.BasePlots(DataSetName = "rmCAF6", n_HVGs = 3000)
Func.BasePlots(DataSetName = "rmCAF6", n_HVGs = 2500)
Func.BasePlots(DataSetName = "rmCAF6", n_HVGs = 2000)
Func.BasePlots(DataSetName = "rmCAF6", n_HVGs = 1500)
Func.BasePlots(DataSetName = "rmCAF6", n_HVGs = 1000)

Func.StateFreq <- function(DataSetName, n_HVGs){
  Unclassified = "Exclude"
  Labels = "seurat_clusters"
  n.neighbors = 40
  k.niche = 10
  pc_use = 15
  Vec.NicheAnnot = 
    c("niche_1"="N7.Plasma cell niche",        "niche_2"="N1.Normal acinar niche",
      "niche_3"="N2.Atrophic acinar niche",    "niche_4"="N8.Lymphoid niche",
      "niche_5"="N5.Islet niche",              "niche_6"="N6.Vascular niche",
      "niche_7"="N9.Lymphoid–myeloid niche",   "niche_8"="N3.Normal–PanIN niche",
      "niche_9"="N4.Cancer niche",             "niche_10"="N10.Neuronal niche")
  DF.NicheRes = read.csv(
    file=paste0(
      DirInteg,"/NicheAnalysisData/[KmeansClust_Final][",Unclassified,"Unclassified_Labeled",Labels,"_Neighbor",n.neighbors,
      "_npc",pc_use,"_k",k.niche,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
      "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".csv"),
    row.names=1) %>% 
    dplyr::mutate(nicheLab = Vec.NicheAnnot[niche]) %>% 
    dplyr::mutate(cell_id = rownames(.), .before = everything())
  cds <- readRDS(
    file=paste0(DirInteg,"/Pseudotime/[CellDataSet][",DataSetName,"_HVG",n_HVGs,"][",
                NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".rds")
  )
  Vec.NicheLabs <- 
    DF.NicheRes$nicheLab %>% 
    unique() %>% 
    str_sort(numeric = TRUE)
  DF.NicheAndState_0 <- 
    Biobase::pData(cds) %>% 
    rownames_to_column(var = "cell_id") %>% 
    inner_join(DF.NicheRes, by = "cell_id") %>% 
    dplyr::mutate(State = paste0("State",State))
  DF.NicheAndState_1 <- 
    table(DF.NicheAndState_0$State,
          DF.NicheAndState_0$nicheLab) %>% 
    as.data.frame() %>% 
    set_colnames(c("State", "Niche", "N_cell"))
  CommonTheme <- 
    theme(
      legend.position = "bottom", 
      legend.direction = "vertical",
      legend.justification = "top",
      axis.title = element_text(face = "bold", color = "black", size = 10),
      legend.title = element_text(face = "bold", color = "black", size = 10),
      axis.text.x = element_text(face = "bold", color = "black", size = 10),
      axis.text.y = element_text(face = "bold", color = "black", size = 10),
      legend.text = element_text(face = "bold", color = "black", size = 10),
      plot.background = element_rect(fill = NA, color = NA),
      panel.background = element_rect(fill = NA, color = NA),
      legend.background = element_rect(fill = NA, color = NA),
      legend.key = element_blank(),
      panel.border = element_rect(fill = NA, color = "black"),
      panel.grid = element_blank())
  p1 <- 
    ggplot(DF.NicheAndState_1,
         aes(y = Niche, x = N_cell)) +
    geom_bar(stat = "identity", aes(fill = State), color = "gray40", linewidth = 0.2) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.03))) +
    labs(y = NULL) + CommonTheme
  p2 <- 
    ggplot(DF.NicheAndState_1,
         aes(y = State, x = N_cell)) +
    geom_bar(stat = "identity", aes(fill = Niche), color = "gray40", linewidth = 0.2) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.03))) +
    labs(y = NULL) + CommonTheme
  DF.Coord_0 <- data.frame(
    Component1 = reducedDimS(cds)[1, ],
    Component2 = reducedDimS(cds)[2, ])
  DF.Coord_1 <- 
    DF.Coord_0 %>% 
    cbind(Biobase::pData(cds)[
      rownames(DF.Coord_0), c("orig.ident", "seurat_clusters", "Pseudotime", "State")
      ]) %>% 
    rownames_to_column(var = "cell_id") %>% 
    inner_join(DF.NicheRes, by = "cell_id") %>% 
    ungroup() %>% 
    dplyr::mutate(
      State = paste0("State ",State),
      seurat_clusters = paste0("CAF-", seurat_clusters),
      seurat_clusters = factor(seurat_clusters, levels = paste0("CAF-", c(4,8,0,1,2,5,3,7,6))),
      nicheLab = factor(nicheLab, levels = Vec.NicheLabs)) %>% 
    dplyr::arrange(desc(nicheLab))
  p3 <- 
    ggplot(DF.Coord_1, aes(x = Component1, y = Component2)) +
      geom_point(aes(color = nicheLab)) +
      labs(x = "Component 1", y = "Component 2", color = NULL) +
      guides(color = guide_legend(override.aes = list(size = 5))) +
      scale_color_manual(
        values = setNames(c("red", rep("gray", times = 9)),
                          Vec.NicheLabs)) + CommonTheme
  ggsave(
    plot = (p1|p2|p3) & CommonTheme,
    file=paste0(
      DirInteg,"/Pseudotime/[StateNicheComposition][",DataSetName,"_HVG",n_HVGs,"][",
      NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
      "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"),
    width = 12, height = 6, dpi = 400, bg = "transparent")
  DF.NicheComposition_0 <-
    table(DF.Coord_1$nicheLab, DF.Coord_1$State) %>% 
    as.data.frame() %>% 
    pivot_wider(names_from = Var1, values_from = Freq, id_cols = Var2) %>% 
    column_to_rownames(var = "Var2") %>% t() %>% as.data.frame()
  DF.NicheComposition_1 <- 
    DF.NicheComposition_0 %>% 
    dplyr::mutate(n_total = apply(DF.NicheComposition_0, 1, sum))
  write.csv(
    DF.NicheComposition_1,
    file=paste0(DirInteg,"/Pseudotime/[NicheComposition][",DataSetName,"_HVG",n_HVGs,"][",
                NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".csv"),
    row.names = T)
}
Func.StateFreq(DataSetName = "rmCAF6", n_HVGs = 5000)
Func.StateFreq(DataSetName = "rmCAF6", n_HVGs = 3000)
Func.StateFreq(DataSetName = "rmCAF6", n_HVGs = 2500)
Func.StateFreq(DataSetName = "rmCAF6", n_HVGs = 2000)
Func.StateFreq(DataSetName = "rmCAF6", n_HVGs = 1500)
Func.StateFreq(DataSetName = "rmCAF6", n_HVGs = 1000)

Func.BasePlots_SetRoot <- function(DataSetName, n_HVGs, RootState){
  cds <- readRDS(
    file=paste0(DirInteg,"/Pseudotime/[CellDataSet][",DataSetName,"_HVG",n_HVGs,"][",
                NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".rds")
  )
  cds_rootstate <- orderCells(
    cds,
    root_state = RootState
  )
  DF.Coord_0 <- data.frame(
    Component1 = reducedDimS(cds_rootstate)[1, ],
    Component2 = reducedDimS(cds_rootstate)[2, ])
  DF.Coord_1 <- DF.Coord_0 %>% 
    cbind(pData(cds_rootstate)[rownames(DF.Coord_0), 
                     c("orig.ident", "seurat_clusters", "Pseudotime", "State")]) %>% 
    dplyr::mutate(State = paste0("State ",State),
                  seurat_clusters = paste0("CAF-", seurat_clusters),
                  seurat_clusters = factor(seurat_clusters, levels = paste0("CAF-", c(4,8,0,1,2,5,3,7,6))))
  DF.Node <- 
    t(reducedDimK(cds_rootstate)) %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "node")　%>% 
    dplyr::rename(Component1 = V1, Component2 = V2)
  DF.Edge <- 
    igraph::as_data_frame(
      minSpanningTree(cds_rootstate),
      what = "edges")
  DF.EdgeCoord <- 
    DF.Edge %>% 
    left_join(
      DF.Node,
      by = c("from" = "node")) %>% 
    dplyr::rename(x = Component1, y = Component2) %>% 
    left_join(
      DF.Node,
      by = c("to" = "node")) %>% 
    dplyr::rename(xend = Component1, yend = Component2)
  CommonTheme <- 
    theme(
      axis.title = element_text(face = "bold", color = "black", size = 10),
      legend.title = element_text(face = "bold", color = "black", size = 10),
      axis.text.x = element_text(face = "bold", color = "black", size = 10),
      axis.text.y = element_text(face = "bold", color = "black", size = 10),
      legend.text = element_text(face = "bold", color = "black", size = 10),
      plot.background = element_rect(fill = NA, color = NA),
      panel.background = element_rect(fill = NA, color = NA),
      legend.background = element_rect(fill = NA, color = NA),
      legend.key = element_blank(),
      panel.border = element_rect(fill = NA, color = "black"),
      panel.grid = element_blank())
  CommonLabs <- 
    labs(x = "Component 1", y = "Component 2", color = NULL)
  p1 <- 
    ggplot(DF.Coord_1, aes(x = Component1, y = Component2)) +
    geom_point(aes(color = State)) +
    geom_segment(data = DF.EdgeCoord,
                 aes(x = x, y = y, xend = xend, yend = yend)) +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    CommonLabs + CommonTheme
  p2 <- 
    ggplot(DF.Coord_1, aes(x = Component1, y = Component2)) +
    geom_point(aes(color = seurat_clusters)) +
    geom_segment(data = DF.EdgeCoord,
                 aes(x = x, y = y, xend = xend, yend = yend)) +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    scale_color_manual(
      values = setNames(as.character(cols2), paste0("CAF-", 0:8))) +
    CommonLabs + CommonTheme
  p3 <- 
    ggplot(DF.Coord_1, aes(x = Component1, y = Component2)) +
    geom_point(aes(color = Pseudotime)) + 
    geom_segment(data = DF.EdgeCoord,
                 aes(x = x, y = y, xend = xend, yend = yend)) +
    labs(x = "Component 1", y = "Component 2", color = "Pseudotime") + 
    scale_color_viridis(
      option = "inferno",
      guide=guide_colorbar(direction="vertical",
                           title.position="top",
                           title.hjust=0.5,
                           frame.colour="black",
                           ticks.colour="black")) +
    CommonTheme
  ggsave(
    (p1/p2/p3) &
      theme(
        plot.background = element_rect(fill = NA, color = NA)),
    file=paste0(DirInteg,"/Pseudotime/[BasePlots_Time][",DataSetName,"_HVG",n_HVGs,"_RootState",RootState,"][",
                NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"),
    dpi=400, width=5, height=10, bg="transparent")
  DF.Coord_1 %>% 
    rownames_to_column(var = "cell_id") %>% 
    write.csv(
      file=paste0(DirInteg,"/Pseudotime/[MetaData][",DataSetName,"_HVG",n_HVGs,"_RootState",RootState,"][",
                  NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                  "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".csv"),
      row.names = F)
}
Func.BasePlots_SetRoot(DataSetName = "rmCAF6", n_HVGs = 1500, RootState = 1)

Func.DEGsAlongTime <- function(DataSetName){
  cds <- readRDS(
    file=paste0(DirInteg,"/Pseudotime/[CellDataSet][",DataSetName,"][",
                NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".rds"))
  deg.pseudo <- differentialGeneTest(
    cds,
    fullModelFormulaStr = "~sm.ns(Pseudotime)")
  Genes.pseudo <- 
    rownames(dplyr::filter(deg.pseudo, qval < 0.0001))
  write.table(
    Genes.pseudo,
    file = paste0(DirInteg,"/Pseudotime/[DEGsAlongTime][",DataSetName,"][",
                  NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                  "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".csv"),
    sep = ",", row.names = FALSE, col.names = FALSE)
}
Func.DEGsAlongTime(DataSetName = "WholeCAF")
Func.DEGsAlongTime(DataSetName = "rmCAF6")

Func.DEGsAlongTime_HmAndClust <- function(DataSetName, k.clust){
  # 1. gene heatmap and gene clustering
  cds <- readRDS(
    file=paste0(DirInteg,"/Pseudotime/[CellDataSet][",DataSetName,"][",
                NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".rds")
  )
  Genes.pseudo <- read.csv(
    file = paste0(DirInteg,"/Pseudotime/[DEGsAlongTime][",DataSetName,"][",
                  NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                  "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".csv"),
    header = FALSE)[ , 1]
  set.seed(1234)
  PH <- 
    plot_pseudotime_heatmap(
      cds[Genes.pseudo, ],
      cluster_rows = TRUE,
      hclust_method = "ward.D2",
      num_clusters = k.clust,
      cores = 1,
      show_rownames = FALSE,
      return_heatmap = TRUE)
  Vec.GeneCluster <- cutree(PH$tree_row, k = k.clust)
  DF.GeneCluster <- data.frame(
    gene = names(Vec.GeneCluster),
    cluster = Vec.GeneCluster)
  DF.GOres <- data.frame()
  for(i in 1:k.clust){
    UniverseGenes <- rownames(Seu.CAF_0)
    GeneList <- 
      DF.GeneCluster %>% 
      dplyr::filter(cluster == i) %>% 
      rownames()
    ego <- clusterProfiler::enrichGO(
      gene          = GeneList,
      universe      = UniverseGenes,
      OrgDb         = org.Hs.eg.db,
      keyType       = "SYMBOL",
      ont           = "BP",
      pAdjustMethod = "BH")
    ego.simple <- clusterProfiler::simplify(
      ego,
      cutoff = 0.7,
      by = "p.adjust",
      select_fun = min)
    if(ego.simple@result %>% nrow() >1){
      DF.GOres.target <- 
        ego.simple@result %>% 
        cbind(GeneClust = paste0("GeneGroup", i))
    }else{
      DF.GOres.target <- NULL
    }
    DF.GOres <- 
      rbind(DF.GOres,
            DF.GOres.target)
  }
  DF.GOres.top5 <- 
    DF.GOres %>% 
    group_by(GeneClust) %>% 
    dplyr::arrange(p.adjust) %>% 
    slice_head(n = 5) %>% ungroup() %>% 
    dplyr::mutate(
      minuslog10adjp = (-1)*log10(p.adjust)#,
      #Description = factor(Description, levels = rev(.$Description)),
      #GeneClust = factor(GeneClust, levels = paste0("GeneGroup", k.clust:1))
    )
  DF.GOres.top5.plot <- 
    DF.GOres.top5 %>% 
    dplyr::mutate(GeneClust = factor(GeneClust, levels = paste0("GeneGroup", 1:k.clust))) %>% 
    group_by(GeneClust) %>% 
    dplyr::arrange(desc(p.adjust), .by_group = TRUE) %>% 
    dplyr::mutate(term_order = row_number()) %>% 
    ungroup() %>% 
    dplyr::mutate(
      group_num = as.numeric(GeneClust),
      y_pos = row_number() + (group_num - 1)*0.7,
      significance = ifelse(p.adjust<0.05, "", "n.s."))
  DF.GroupLab <- 
    DF.GOres.top5.plot %>% 
    dplyr::group_by(GeneClust) %>% 
    dplyr::summarise(
      y_mid = mean(y_pos),
      .groups = "drop")
  pBar <- 
    ggplot(DF.GOres.top5.plot, 
           aes(x = minuslog10adjp, y = y_pos)) +
    geom_col(orientation = "y") +
    #geom_text(aes(label = significance), hjust = -0.3) +
    labs(x = "–log10 adj.p", y = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.03))) +
    scale_y_continuous(
      breaks = DF.GOres.top5.plot$y_pos,
      labels = DF.GOres.top5.plot$Description,
      expand = expansion(mult = c(0.01, 0.01)),
      sec.axis = dup_axis(
        breaks = DF.GroupLab$y_mid,
        labels = DF.GroupLab$GeneClust,
        name = NULL)) +
    theme(
      plot.background = element_rect(fill = NA, color = NA),
      panel.background = element_rect(fill = NA, color = NA),
      panel.border = element_rect(fill = NA, color = "black"),
      axis.title.x = element_text(face = "bold", color = "black"),
      axis.text.y.left = element_text(face = "bold", color = "black"),
      axis.text.y.right = element_text(face = "bold", color = "black", angle = -90, hjust = 0.5, vjust = 0.5),
      axis.text.x = element_text(face = "bold", color = "black"),
      axis.ticks.y.right = element_blank())
  # ridge plot
  DF.Coord_1 <- read.csv(
    file=paste0(DirInteg,"/Pseudotime/[MetaData][",DataSetName,"][",
                NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".csv"))
  p.ridge <- 
    ggplot(DF.Coord_1,
           aes(x = Pseudotime, y = seurat_clusters)) +
    geom_density_ridges(scale = 1.5, alpha = 0.6) +
    labs(x = NULL, y = NULL) +
    scale_x_continuous(
      limits = range(DF.Coord_1$Pseudotime),
      expand = expansion(mult = c(0, 0))) +
    theme(
      axis.ticks = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_text(face = "bold", color = "black"),
      plot.background = element_rect(fill = NA, color = NA),
      panel.background = element_rect(fill = NA, color = NA),
      legend.background = element_rect(fill = NA, color = NA),
      legend.key = element_blank(),
      panel.border = element_rect(fill = NA, color = NA))
  # save
  ObjestList <- list(
    "Heatmap" = PH,
    "Barplot" = pBar,
    "GOres" = DF.GOres,
    "GOres.top5" = DF.GOres.top5)
  saveRDS(
    ObjestList,
    file=paste0(DirInteg,"/Pseudotime/[DEGsAlongTime_ObjList][",DataSetName,"_k",k.clust,"][",
                NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".rds"))
}
Func.DEGsAlongTime_HmAndClust(DataSetName = "WholeCAF", k.clust = 6)
Func.DEGsAlongTime_HmAndClust(DataSetName = "rmCAF6", k.clust = 7)

Func.DEGsAlongTime_AlignedFig <- function(DataSetName, k.clust){
  ObjestList <- readRDS(
    file=paste0(DirInteg,"/Pseudotime/[DEGsAlongTime_ObjList][",DataSetName,"_k",k.clust,"][",
                NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".rds"))
  PH <- ObjestList[["Heatmap"]]
  pBar <- ObjestList[["Barplot"]]
  DF.GOres <- ObjestList[["GOres"]]
  DF.GOres.top5 <- ObjestList[["GOres.top5"]]
  
  GetGrob <- function(PH, name){
    PH$gtable$grobs[[which(PH$gtable$layout$name == name)]]
  }
  G.Tree   <- GetGrob(PH, "row_tree")
  G.Anno   <- GetGrob(PH, "row_annotation")
  G.Body   <- GetGrob(PH, "matrix")
  G.Legend <- GetGrob(PH, "legend")
  sp <- nullGrob()
  center_grob_in_block <- function(gb) {
    arrangeGrob(
      nullGrob(), gb, nullGrob(),
      ncol = 1,
      heights = unit.c(
        unit(1, "null"),
        grobHeight(gb),
        unit(1, "null")))
  }
  Legend.Center <- center_grob_in_block(G.Legend)
  Legend.Center.Left <- 
    ggdraw() +
    draw_grob(Legend.Center) + 
    theme(plot.margin = margin(0,0,0,0))
  Legend.Title <- 
    ggplot() +
    annotate("text",
             x = 1,  y = 0.5,
             label = "Row Z-score",
             vjust = -0.5,
             angle = 90,
             fontface = "bold",
             size = 4 ) +
    coord_cartesian(
      xlim = c(0.5, 1),
      ylim = c(0, 1),
      expand = FALSE) +
    theme_void() + 
    theme(plot.margin = margin(0,0,0,0))
  Legend.Block <- cowplot::plot_grid(
    Legend.Title,
    Legend.Center.Left,
    nrow = 1,
    rel_widths = c(0.5, 0.5)
  )
  G.RidgeAxis <- cowplot::get_plot_component(p.ridge, "axis-l")
  G.RidgePanel <- cowplot::get_plot_component(
    p.ridge + theme(axis.text.y = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.title.y = element_blank()),
    "panel")
  Blank <- ggplot() + theme_void()
  # heatmap左側
  G.HeatLeft <- cowplot::plot_grid(
    G.Tree,
    G.Anno,
    nrow = 1,
    rel_widths = c(0.85, 0.15))
  # 列幅を共通化
  Widths <- c(0.15, 0.35, 0.10)
  Ridge.Row <- cowplot::plot_grid(
    G.RidgeAxis,
    G.RidgePanel,
    Blank,
    nrow = 1,
    rel_widths = Widths)
  Heatmap.Row <- cowplot::plot_grid(
    G.HeatLeft,
    G.Body,
    Legend.Block,
    nrow = 1,
    rel_widths = Widths)
  RidgeHeat <- cowplot::plot_grid(
    Ridge.Row,
    Heatmap.Row,
    ncol = 1,
    rel_heights = c(0.55, 1))
  Final <- cowplot::plot_grid(
    RidgeHeat,
    ObjestList[["Barplot"]],
    ncol =2,
    rel_widths = c(0.45, 0.55)
  )
  ggsave(
    Final,
    file=paste0(DirInteg,"/Pseudotime/[Plot_DEGsAlongTime][",DataSetName,"][",
                NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"),
    dpi=400, width=15, height=10, bg="transparent")
}
Func.DEGsAlongTime_AlignedFig(DataSetName = "WholeCAF", k.clust = 6)
Func.DEGsAlongTime_AlignedFig(DataSetName = "rmCAF6", k.clust = 7)

Plot.Heatmap <- 
  wrap_elements(full = ObjestList[["Heatmap"]]$gtable)
(p.ridge/Plot.Heatmap)|ObjestList[["Barplot"]]

Func.DEGsAmongState <- function(SeuObj, DataSetName, n_HVGs){
  # Read data
  DF.Coord_1 <- read.csv(
    file=paste0(DirInteg,"/Pseudotime/[MetaData][",DataSetName,"_HVG",n_HVGs,"][",
                NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".csv"))
  SeuObj@graphs <- list()
  Idents(SeuObj) <- "seurat_clusters"
  set.seed(1234)
  SeuObj_ds <- subset(
    SeuObj,
    downsample = 500
  )
  SeuObj_ds <- JoinLayers(SeuObj_ds)
  if(identical(colnames(SeuObj_ds), DF.Coord_1$cell_id)){
    SeuObj_ds$State <- DF.Coord_1$State
    Idents(SeuObj_ds) <- factor(
      SeuObj_ds$State,
      levels = paste0("State ", 1:length(unique(SeuObj_ds$State))))
  }
  # DEG table among states
  set.seed(1234)
  DF.DEGsAmongState_0 <- FindAllMarkers(
    SeuObj_ds, 
    only.pos = FALSE,
    slot = "data",
    min.pct = 0.00, 
    min.cell.feature = 0, 
    min.cells.group = 0, 
    logfc.threshold = 0, 
    return.thresh = Inf,
    random.seed = 1234,
    test.use = "wilcox")
  DF.DEGsAmongState_1 <- 
    DF.DEGsAmongState_0 %>% 
    dplyr::mutate(entrez = tibble::deframe(DF.AllGenes_addEntrez)[gene]) %>% 
    group_by(cluster) %>% 
    dplyr::arrange(desc(avg_log2FC), .by_group = TRUE) %>% 
    ungroup() %>% 
    dplyr::select(c(cluster, gene, avg_log2FC, pct.1, pct.2, p_val, p_val_adj, entrez)) 
  write.csv(
    DF.DEGsAmongState_1,
    file=paste0(DirInteg,"/Pseudotime/[DataTable_DEGsAmongState][",DataSetName,"_HVG",n_HVGs,"][",
                NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".csv"),
    row.names = F)
  DF.DEGsAmongState_2 <- 
    DF.DEGsAmongState_1 %>%
    dplyr::filter(avg_log2FC > 0) %>% 
    dplyr::mutate(
      PointCol = ifelse(avg_log2FC > 1 & 
                          #pct.1 > 0.4 & 
                          p_val_adj < 0.05, "Colored", "Gray"),
      PointCol = factor(PointCol, levels = c("Gray", "Colored"))) %>% 
    group_by(PointCol) %>% 
    dplyr::arrange(avg_log2FC, .by_group = TRUE)
  DF.DEGsAmongState_3 <- DF.DEGsAmongState_2 %>%
    dplyr::filter(PointCol == "Colored" ) 
  Plot.Volc <- 
  DF.DEGsAmongState_2 %>% 
    ggplot(aes(y = avg_log2FC, x = cluster)) +
    geom_point(
      position = "jitter",
      aes(color = PointCol,
          fill = ifelse(PointCol == "Colored", avg_log2FC, NA),
          shape = PointCol),
      size = 3) +
    labs(y = "log2 fold change", x = NULL, fill = "log 2 fold change") +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.03)),
      breaks = c(1, 4, 8)) +
    scale_shape_manual(
      values = c("Colored" = 21, "Gray" = 19),
      guide = "none") +
    scale_color_manual(
      values = c("Colored" = "black", "Gray" = "gray80"),
      guide = "none") +
    scale_fill_distiller(
      palette = "Oranges",
      na.value = "gray80",
      guide = guide_colorbar(
        title.hjust = 0.5,
        ticks.colour = "black",
        frame.colour = "black")) +
    theme(legend.position = "none",
          axis.text.x = element_text(face = "bold", color = "black", size = 10),
          axis.text.y = element_text(face = "bold", color = "black", size = 10),
          axis.title = element_text(face = "bold", color = "black", size = 10),
          legend.title = element_text(face = "bold", color = "black", size = 10),
          legend.text = element_text(face = "bold", color = "black", size = 10),
          plot.background = element_rect(fill = NA, color = NA),
          panel.background = element_rect(fill = NA, color = NA),
          legend.background = element_rect(fill = NA, color = NA),
          panel.border = element_rect(fill = NA, color = "black"),
          #axis.line = element_line(color = "black"),
          strip.placement = "outside",
          strip.background = element_rect(fill = NA, color = NA),
          panel.grid = element_blank())
  ggsave(
    Plot.Volc,
    file=paste0(DirInteg,"/Pseudotime/[Plot_DEGsAmongState_Volcano][",DataSetName,"_HVG",n_HVGs,"][",
                NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"),
    dpi = 400, width = 5, height = 4, bg = "transparent")
  }
Func.DEGsAmongState(SeuObj = Seu.CAF_0, DataSetName = "WholeCAF")
Func.DEGsAmongState(SeuObj = Seu.CAF_rmCAF6, DataSetName = "rmCAF6", n_HVGs = 3000)
Func.DEGsAmongState(SeuObj = Seu.CAF_rmCAF6, DataSetName = "rmCAF6", n_HVGs = 2500)
Func.DEGsAmongState(SeuObj = Seu.CAF_rmCAF6, DataSetName = "rmCAF6", n_HVGs = 2000)
Func.DEGsAmongState(SeuObj = Seu.CAF_rmCAF6, DataSetName = "rmCAF6", n_HVGs = 1500)

Func.GeneSetAnalAmongState <- function(SeuObj, DataSetName, n_HVGs){
  # 1. Read DEGs table
  DF.DEGsAmongState_1 <- 
    read.csv(
    file=paste0(DirInteg,"/Pseudotime/[DataTable_DEGsAmongState][",DataSetName,"_HVG",n_HVGs,"][",
                NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".csv"),
    row.names = NULL)
  NumOfState <- 
    DF.DEGsAmongState_1$cluster %>% 
    unique() %>% length()
  # 2. Hallmarks GSEA
  # 2-1. read Hallmarks gene set
  H = clusterProfiler::read.gmt("/Volumes/PortableSSD/[4]myR/[1]dataset/h.all.v2024.1.Hs.entrez.gmt") %>% 
    dplyr::filter(gene %in% DF.AllGenes_addEntrez$ENTREZID)
  CategoryC2 = msigdbr(species="Homo sapiens",  category="C2")
  Additional_0 = 
    data.frame(term = "Buffa_Hypoxia", symbol = Vec.BuffaOrig) %>% 
    rbind(data.frame(term="Winter_Hypoxia", symbol=Vec.WinterOrig)) %>% 
    rbind(data.frame(term="KEGG_Cell_cycle", symbol=CategoryC2 %>% filter(gs_name == "KEGG_CELL_CYCLE") %>% pull(gene_symbol) %>% unique())) %>% 
    rbind(data.frame(term="myCAF_signature", symbol=DF.Markers$Gene[DF.Markers$Classification1=="myCAF"])) %>% 
    rbind(data.frame(term="iCAF_signature", symbol=DF.Markers$Gene[DF.Markers$Classification1=="iCAF"]))
  Vec.AdditionalEntrez_0 = clusterProfiler::bitr(
    Additional_0$symbol, 
    fromType="SYMBOL", 
    toType="ENTREZID", 
    OrgDb=org.Hs.eg.db, 
    drop=FALSE)
  Vec.AdditionalEntrez_1 = 
    Vec.AdditionalEntrez_0$ENTREZID %>% 
    setNames(Vec.AdditionalEntrez_0$SYMBOL)
  Additional_1 = Additional_0 %>% 
    dplyr::mutate(gene = Vec.AdditionalEntrez_1[symbol])
  Additional_2 = dplyr::select(Additional_1, c(term, gene))
  # 2-2. Run GSEA()
  DF.gseaRes_0 <- data.frame()
  for(i in 1:NumOfState){
    Vec.RankedLog2FC <- 
      DF.DEGsAmongState_1 %>% 
      dplyr::filter(cluster == paste0("State ", i))　%>% 
      dplyr::select(entrez, avg_log2FC) %>% 
      tibble::deframe()
    gseaRes = clusterProfiler::GSEA(
      geneList = Vec.RankedLog2FC,
      minGSSize = 10,
      TERM2GENE = rbind(H,Additional_2),
      pvalueCutoff = 1,
      verbose = FALSE, 
      eps = 0)
    DF.gseaRes_target <- 
      cbind(State = paste0("State",i),
            gseaRes@result) %>% 
      dplyr::mutate(Direction = ifelse(NES > 0, "UP", "DOWN")) %>% 
      relocate(Direction, .after = State)
    DF.gseaRes_0 <- rbind(
      DF.gseaRes_0,
      DF.gseaRes_target
    )
  }
  DF.gseaRes_1 <- 
    DF.gseaRes_0 %>% 
    dplyr::filter(NES > 0) %>% 
    slice_max(order_by = NES, by = State, n = 5) %>% 
    dplyr::mutate(
      ID_lab = paste(ID, State, sep = "_"),
      ID_lab = factor(ID_lab, levels = rev(ID_lab)),
      minuslog10adjP = (-1)*log10(p.adjust))
  PlotBar.GSEA.hallmarks <- 
    ggplot(DF.gseaRes_1, aes(x = NES, y = ID_lab)) +
    geom_bar(stat = "identity",
             color = "gray10", linewidth = 0.2,
             aes(fill = ifelse(p.adjust<0.05, minuslog10adjP, NA))) +
    labs(fill = "–log10(adjusted P value)", y = NULL, x = "NES") +
    scale_x_continuous(
      expand = expansion(mult = c(0, 0.03)),
      breaks = c(0, 1, 2)) +
    scale_y_discrete(label = function(x) 
      str_remove(x, pattern = "_State.") %>% 
        clean_gs_label_ForFig()) +
    scale_fill_distiller(
      palette = "Oranges",
      guide = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        ticks.colour = "black",
        frame.colour = "black")) +
    facet_grid(State ~ ., scales = "free_y", space = "free_y", switch = "y") +
    theme(legend.position = "bottom",
          legend.justification = "center",
          axis.text.x = element_text(face = "bold", color = "black", size = 10),
          axis.text.y = element_text(face = "bold", color = "black", size = 10),
          axis.title = element_text(face = "bold", color = "black", size = 10),
          legend.title = element_text(face = "bold", color = "black", size = 10),
          legend.text = element_text(face = "bold", color = "black", size = 10),
          plot.background = element_rect(fill = NA, color = NA),
          panel.background = element_rect(fill = NA, color = NA),
          legend.background = element_rect(fill = NA, color = NA),
          #panel.border = element_rect(fill = NA, color = "black"),
          axis.line = element_line(color = "black"),
          strip.placement = "outside",
          strip.background = element_rect(fill = NA, color = NA),
          panel.grid = element_blank())
  ggsave(
    PlotBar.GSEA.hallmarks,
    file=paste0(DirInteg,"/Pseudotime/[Plot_DEGsAmongState_GSEA(GSEA)][",DataSetName,"_HVG",n_HVGs,"][",
                NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"),
    dpi=400, width=3.5, height=6.0, bg="transparent")
  DF.gseaRes_0 %>% 
    dplyr::select(c(State, Direction, Description, setSize, NES, p.adjust, qvalue, leading_edge, core_enrichment)) %>% 
    dplyr::mutate(Description = clean_gs_label_ForTable(Description)) %>%
    group_by(State) %>% 
    dplyr::arrange(desc(NES), .by_group = TRUE) %>% 
    write.csv(
      file=paste0(DirInteg,"/Pseudotime/[DataTable_GseaHallmarkAmongState][",DataSetName,"_HVG",n_HVGs,"][",
                  NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                  "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".csv"),
      row.names = F)
  # 3. GO analysis GSEA
  DF.gsegoRes_0 <- data.frame()
  for(i in 1:NumOfState){
    Vec.RankedLog2FC <- 
      DF.DEGsAmongState_1 %>% 
      dplyr::filter(cluster == paste0("State ", i))　%>% 
      dplyr::select(entrez, avg_log2FC) %>% 
      tibble::deframe()
    set.seed(1234)
    gsegoRes = clusterProfiler::gseGO(
      geneList=Vec.RankedLog2FC,
      OrgDb=org.Hs.eg.db,
      ont="ALL",   #"BP","CC","MF","ALL"から選択
      minGSSize= 50,
      pvalueCutoff = 0.05,
      verbose=FALSE,
      nPermSimple = 10000,
      eps=0)
    if(nrow(gsegoRes@result) > 0){
    DF.gsegoRes_target <- 
      cbind(State = paste0("State",i),
            gsegoRes@result) %>% 
      dplyr::arrange(desc(NES)) %>% 
      dplyr::mutate(Direction = ifelse(NES > 0, "UP", "DOWN")) %>% 
      relocate(Direction, .after = State)
    DF.gsegoRes_0 <- rbind(
      DF.gsegoRes_0,
      DF.gsegoRes_target)}
    }
  DF.gsegoRes_1 <- 
    DF.gsegoRes_0 %>% 
    dplyr::filter(NES > 0) %>% 
    slice_max(order_by = NES, by = State, n = 5) %>% 
    dplyr::mutate(
      ID_lab = paste(Description, State, sep = "_"),
      ID_lab = factor(ID_lab, levels = rev(ID_lab)),
      minuslog10adjP = (-1)*log10(p.adjust))
  PlotBar.GSEA.GO <- 
    ggplot(DF.gsegoRes_1, aes(x = NES, y = ID_lab)) +
    geom_bar(stat = "identity",
             color = "gray10", linewidth = 0.2,
             aes(fill = ifelse(p.adjust<0.05, minuslog10adjP, NA))) +
    labs(fill = "–log10\n(adjusted P value)", y = NULL, x = "NES") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.03))) +
    scale_y_discrete(label = function(x) 
      str_remove(x, pattern = "_State.") %>% 
        clean_gs_label_ForFig()) +
    scale_fill_distiller(
      palette = "Oranges",
      guide = guide_colorbar(
        title.hjust = 0.5,
        ticks.colour = "black",
        frame.colour = "black")) +
    facet_grid(State ~ ., scales = "free_y", space = "free_y", switch = "y") +
    theme(axis.text.x = element_text(face = "bold", color = "black", size = 10),
          axis.text.y = element_text(face = "bold", color = "black", size = 10),
          axis.title = element_text(face = "bold", color = "black", size = 10),
          legend.title = element_text(face = "bold", color = "black", size = 10),
          legend.text = element_text(face = "bold", color = "black", size = 10),
          plot.background = element_rect(fill = NA, color = NA),
          panel.background = element_rect(fill = NA, color = NA),
          legend.background = element_rect(fill = NA, color = NA),
          #panel.border = element_rect(fill = NA, color = "black"),
          axis.line = element_line(color = "black"),
          strip.placement = "outside",
          strip.background = element_rect(fill = NA, color = NA),
          panel.grid = element_blank())
  ggsave(
    PlotBar.GSEA.GO,
    file=paste0(DirInteg,"/Pseudotime/[Plot_DEGsAmongState_GSEA(GO)][",DataSetName,"_HVG",n_HVGs,"][",
                NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"),
    dpi=400, width=8, height=8, bg="transparent")
  # 4. GO analysis (ORA)
  DF.GO.ORA.res <- data.frame()
  for(i in 1:NumOfState){
    UniverseGenes <- Vec.AllGenes
    GeneList <- 
      DF.DEGsAmongState_1 %>% 
      dplyr::filter(
        cluster == paste0("State ",i) &
        avg_log2FC > 1 &
        p_val_adj < 0.05) %>% 
      dplyr::pull(gene)
    if(length(GeneList) >0){
    ego <- clusterProfiler::enrichGO(
      gene          = GeneList,
      universe      = UniverseGenes,
      OrgDb         = org.Hs.eg.db,
      keyType       = "SYMBOL",
      ont           = "BP",
      pAdjustMethod = "BH")
    ego.simple <- clusterProfiler::simplify(
      ego,
      cutoff = 0.7,
      by = "p.adjust",
      select_fun = min)
      if(ego.simple@result %>% nrow() >1){
        DF.GO.ORA.res.target <- 
          ego.simple@result %>% 
          dplyr::mutate(State = paste0("State ", i), .before = ID)
        DF.GO.ORA.res <- 
          rbind(DF.GO.ORA.res,
                DF.GO.ORA.res.target)
      }
    }
  }
  DF.GO.ORA.res.top5 <- 
    DF.GO.ORA.res %>% 
    group_by(State) %>% 
    dplyr::arrange(p.adjust) %>% 
    slice_head(n = 5) %>% ungroup() %>% 
    dplyr::mutate(
      minuslog10adjp = (-1)*log10(p.adjust)#,
      #Description = factor(Description, levels = rev(.$Description)),
      #GeneClust = factor(GeneClust, levels = paste0("GeneGroup", k.clust:1))
    )
  DF.GO.ORA.res.top5.plot <- 
    DF.GO.ORA.res.top5 %>% 
    dplyr::mutate(State = factor(State, levels = paste0("State ", 1:NumOfState))) %>% 
    group_by(State) %>% 
    dplyr::arrange(desc(p.adjust), .by_group = TRUE) %>% 
    dplyr::mutate(term_order = row_number()) %>% 
    ungroup() %>% 
    dplyr::mutate(
      group_num = as.numeric(State),
      y_pos = row_number() + (group_num - 1)*0.7,
      significance = ifelse(p.adjust<0.05, "", "n.s."))
  DF.GroupLab <- 
    DF.GO.ORA.res.top5.plot %>% 
    dplyr::group_by(State) %>% 
    dplyr::summarise(
      y_mid = mean(y_pos),
      .groups = "drop")
  #pBar <- 
    ggplot(DF.GO.ORA.res.top5.plot, 
           aes(x = minuslog10adjp, y = y_pos)) +
    geom_col(orientation = "y") +
    #geom_text(aes(label = significance), hjust = -0.3) +
    labs(x = "–log10 (adjusted P value)", y = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.03))) +
    scale_y_continuous(
      breaks = DF.GO.ORA.res.top5.plot$y_pos,
      labels = DF.GO.ORA.res.top5.plot$Description,
      expand = expansion(mult = c(0.01, 0.01)),
      sec.axis = dup_axis(
        breaks = DF.GroupLab$y_mid,
        labels = DF.GroupLab$State,
        name = NULL)) +
    theme(
      plot.background = element_rect(fill = NA, color = NA),
      panel.background = element_rect(fill = NA, color = NA),
      panel.border = element_rect(fill = NA, color = "black"),
      axis.title.x = element_text(face = "bold", color = "black"),
      axis.text.y.left = element_text(face = "bold", color = "black"),
      axis.text.y.right = element_text(face = "bold", color = "black", angle = -90, hjust = 0.5, vjust = 0.5),
      axis.text.x = element_text(face = "bold", color = "black"),
      axis.ticks.y.right = element_blank())
}
GeneSetAnalAmongState(SeuObj = Seu.CAF_rmCAF6, DataSetName = "rmCAF6", n_HVGs = 1500)








##  7. ssGSEA score (insert)                     ####
Dim1=20
Res1=1.0
Dim2=30
Res2=0.5
method="EuclideanWardd2"
MidPoint.Proliferation=0.65
MidPoint.KEGG_CellCycle=0.7
MidPoint.E2F=0.50
MidPoint.G2M=0.65
MidPoint.MITOTIC_SPINDLE=0.75
MidPoint.MYC_TARGETS_V1=0.65
MidPoint.MYC_TARGETS_V2=0.65
MidPoint.BuffaOrig=0.70
MidPoint.WinterOrig=0.70
MidPoint.HYPOXIA=0.6
MidPoint.myCAF=0.70
MidPoint.iCAF=0.55
MidPoint.apCAF=0.999
  SeuObj.CAF_0 = readRDS(paste0(DirInteg,"Objects/[SeuObj_ExtractedCAF][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                                QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2,".rds") )
  DF.gsvaScore_0 = read.csv(file=paste0(DirInteg,"/[DataTable_ExtractedCAF_ssGSVA][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                                        QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".csv"))
  rownames(DF.gsvaScore_0) <- paste(
    DF.gsvaScore_0$sample,
    DF.gsvaScore_0$cell_id, 
    sep="_")
  DF.gsvaScore_0 = DF.gsvaScore_0 %>% 
    dplyr::mutate(full_id = rownames(DF.gsvaScore_0)) %>% 
    dplyr::select(c(-Proliferation_ssGSVA, 
                    -REACTOME_Cellular_senescence_ssGSVA,
                    -apCAF_ssGSVA))
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
  # 1.Heatmap of all pathways
  DF.Meta_0 = 
    data.frame("full_id"=colnames(SeuObj.CAF_0),
               "seurat_clusters"=paste0("CAF-", SeuObj.CAF_0$seurat_clusters),
               SeuObj.CAF_0@reductions[["umap"]]@cell.embeddings)
  if( identical(DF.Meta_0$full_id, DF.gsvaScore_0$full_id) ){
    DF.Meta_1 = inner_join(DF.Meta_0, DF.gsvaScore_0, by="full_id")
  }
  DF.ScoreMean_0 = 
    DF.Meta_1 %>% 
    dplyr::select(-c("full_id","umap_1","umap_2",
                     "sample","cell_id","X","Y") ) %>% 
    group_by(seurat_clusters) %>% 
    summarise(across(everything(), mean))
  MT.ScoreMean_0 = as.matrix(DF.ScoreMean_0[ , -1])
  rownames(MT.ScoreMean_0) = DF.ScoreMean_0$seurat_clusters
  DF.Zscore_0 = base::scale(MT.ScoreMean_0) %>% as.data.frame()
  colnames(DF.Zscore_0) = str_remove(colnames(DF.Zscore_0), pattern="_ssGSVA")
  DF.Zscore_1 = rownames_to_column(DF.Zscore_0, var="seurat_clusters")
  DF.Zscore_2 = pivot_longer(DF.Zscore_1,
                             -seurat_clusters,
                             names_to = "pathway",
                             values_to = "zscore")
  # hierarchical clustering : pearson/euclidean
  method="SpearmanAverage"
  d_pathway = dist(t(DF.Zscore_0), method="euclidean")
  class(d_pathway) = "dist"
  Hclust_pathway <- hclust(d_pathway, method = "ward.D2")
  d_clusters = dist(DF.Zscore_0, method="euclidean")
  class(d_clusters) = "dist"
  Hclust_clusters <- hclust(d_clusters, method = "ward.D2")
  MethodLabel = "Hierarchical clustering was performed using Ward’s method (ward.D2) with Euclidean distances."
  
  Vec.HclustOrder.Pathway = colnames(DF.Zscore_0)[Hclust_pathway$order]
  DF.Zscore_2$pathway = factor(DF.Zscore_2$pathway,
                               levels=Vec.HclustOrder.Pathway)
  levels(DF.Zscore_2$pathway) = clean_gs_label_ForFig(levels(DF.Zscore_2$pathway))
  Vec.HclustOrder.Subclust = rownames(DF.Zscore_0)[Hclust_clusters$order]
  DF.Zscore_2$seurat_clusters = factor(DF.Zscore_2$seurat_clusters,
                                       levels=paste0("CAF-", rev(Vec.SubclustOrder))) 
  # dendrogram
  library(ggdendro)
  dend_pathway = as.dendrogram(Hclust_pathway)
  dend_data = dendro_data(dend_pathway)
  p_dend = 
    ggplot(dend_data$segments) +
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
    scale_x_continuous(limits = c(0.5, length(Hclust_pathway$labels) + 0.5),
                       breaks = seq_along(Hclust_pathway$labels),
                       labels = Hclust_pathway$labels,
                       expand = expansion(mult=c(0,0))) +
    coord_fixed(ratio=0.2) + 
    #scale_y_reverse() +
    theme_void()
  RectLineWidth=0.8
  R1 <- annotate(geom="rect", xmin=1.5, xmax=3.5, ymin=5.5, ymax=6.5, color="red", fill=NA, linewidth=RectLineWidth)
  R2 <- annotate(geom="rect", xmin=3.5, xmax=8.5, ymin=0.5, ymax=2.5, color="red", fill=NA, linewidth=RectLineWidth)
  R3 <- annotate(geom="rect", xmin=10.5, xmax=15.5, ymin=7.5, ymax=9.5, color="red", fill=NA, linewidth=RectLineWidth)
  R4 <- annotate(geom="rect", xmin=15.5, xmax=21.5, ymin=5.5, ymax=6.5, color="red", fill=NA, linewidth=RectLineWidth)
  R5 <- annotate(geom="rect", xmin=17.5, xmax=21.5, ymin=7.5, ymax=8.5, color="red", fill=NA, linewidth=RectLineWidth)
  R6 <- annotate(geom="rect", xmin=21.5, xmax=24.5, ymin=7.5, ymax=8.5, color="red", fill=NA, linewidth=RectLineWidth)
  R7 <- annotate(geom="rect", xmin=21.5, xmax=31.5, ymin=2.5, ymax=3.5, color="red", fill=NA, linewidth=RectLineWidth)
  R8 <- annotate(geom="rect", xmin=31.5, xmax=32.5, ymin=2.5, ymax=5.5, color="red", fill=NA, linewidth=RectLineWidth)
  R9 <- annotate(geom="rect", xmin=34.5, xmax=35.5, ymin=2.5, ymax=5.5, color="red", fill=NA, linewidth=RectLineWidth)
  R10 <- annotate(geom="rect", xmin=36.5, xmax=37.5, ymin=2.5, ymax=5.5, color="red", fill=NA, linewidth=RectLineWidth)
  R11 <- annotate(geom="rect", xmin=41.5, xmax=42.5, ymin=2.5, ymax=3.5, color="red", fill=NA, linewidth=RectLineWidth)
  R12 <- annotate(geom="rect", xmin=43.5, xmax=44.5, ymin=2.5, ymax=5.5, color="red", fill=NA, linewidth=RectLineWidth)
  R13 <- annotate(geom="rect", xmin=50.5, xmax=52.5, ymin=3.5, ymax=4.5, color="red", fill=NA, linewidth=RectLineWidth)
  p_heat = 
    ggplot(DF.Zscore_2, 
           aes(y=seurat_clusters, x=pathway, fill=zscore)) +
    geom_tile(color="gray50") +
    labs(x=NULL, y=NULL, fill="Scaled mean ssGSEA ES") +
    scale_x_discrete(expand=expansion(mult=c(0,0))) +
    scale_y_discrete(expand=expansion(mult=c(0,0)),
                     position="left") +
    scale_fill_gradientn(limits=c(max(abs(range(DF.Zscore_2$zscore))),
                                  -max(abs(range(DF.Zscore_2$zscore)))),
                         breaks=c(-2, 0, 2),
                         colors=
                           #c("#5B6F9E", "#F7F7F7", "#B8C76A"),
                           #c("#6A6FA3", "#F7F7F7", "#4E9A91"),
                           #c("#A6D3C5", "#F7F7F7", "#B36A9C"),
                           c("#A6D3C5", "#F7F7F7", "#C05A9E"),
                           #c("#A6D3C5", "#F7F7F7", "#C23A8B"),
                           #c("#4F5D75", "#F7F7F7", "#2A9D8F"),
                           #c("#5B7FB0", "#F7F7F7", "#7A8F3A"),
                           #c("#6F8DB8", "#F7F7F7", "#A88D55"),
                           #c("#5B7FB0", "#F7F7F7", "#9C7A3A"),
                           #c("#2C3E75", "#F7F7F7", "#3A9C8C"),
                           #c("#1F9E89", "#F7F7F7", "#7B3294"),
                           #c("#6F8DB8", "#F7F7F7", "#7A3E9D"),
                           #c("#5B7FB0", "#F7F7F7", "#7A3E9D"),
                           #c("#3B6EA8", "#F7F7F7", "#7A3E9D"),
                           #c("#4dac26","#b8e186","#f7f7f7","#f1b6da","#d01c8b"),
                         guide = guide_colorbar(
                           #direction="horizontal",
                           title.position="top",
                           title.hjust=0.5,
                           frame.colour="black",
                           ticks.colour="black")) +
    coord_fixed(ratio=1.15) + 
    #R1 + R2 + R3 + R4 + R5 + R6 + R7 + R8 + R9 + R10 + R11 + R12 + R13 +
    theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
          legend.position="bottom",
          axis.text = element_text(face="bold", color="black"),
          axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5, color="black",
                                     size=12),
          axis.text.y = element_text(color="black",
                                     size=12),
          legend.text = element_text(face = "bold", color="black"),
          legend.title = element_text(face = "bold"),
          legend.background = element_rect(fill="transparent", color=NA),
          plot.background = element_rect(fill="transparent", color=NA),
          panel.background = element_rect(fill="transparent", color=NA),
          panel.border = element_rect(fill="transparent", color = "black"),
          panel.grid = element_blank())
  #p_heat
  p_heatdend <- 
    p_dend + p_heat + 
    plot_layout(ncol = 1) +
    plot_annotation(caption = MethodLabel) &
    theme(plot.background = element_rect(fill="transparent", color=NA),
          panel.background = element_rect(fill="transparent", color=NA))
  #p_heatdend
  ggsave(plot = p_heatdend,
         file=paste0(DirInteg,"[Figure][SubClust_ssGSVAscores_Heatmap_ForPaper2][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                     "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,
                     "_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"), 
         dpi=500, width=12, height=8, bg = "transparent" )
  Disp = c("myCAF signature", "EMT", "Hedgehog signaling", "TGF-β signaling", "Hypoxia", 
           "Buffa hypoxia", "Winter hypoxia", "Glycolysis",
           "E2F targets", "G2/M checkpoint", "MYC targets V2", 
           "Notch signaling", "Wnt/β-catenin signaling",
           "iCAF signature", "IL-2/STAT5 signaling", "IL-6/JAK/STAT3 signaling", "TNF-α signaling via NF-κB",
           "Inflammatory response", "IFN-α response", 
           "Heme metabolism", "Peroxisome")
  DF.Zscore_2$pathway %>% unique()
  
  DF.Zscore_3 <- 
    DF.Zscore_2 %>% 
    #dplyr::mutate(pathway = clean_gs_label_ForFig(pathway)) %>% 
    dplyr::filter(pathway %in% Disp) %>% 
    dplyr::mutate(
      pathway = factor(pathway, levels = Disp),
      seurat_clusters = factor(seurat_clusters, 
                               levels = paste0("CAF-", rev(Vec.SubclustOrder))))
  p_heat_Presentation <- 
    DF.Zscore_3 %>% 
    ggplot(aes(y=seurat_clusters, x=pathway, fill=zscore)) +
    geom_tile(color="gray50") +
    labs(x=NULL, y=NULL, fill="Scaled mean ssGSEA ES") +
    scale_x_discrete(expand=expansion(mult=c(0,0))) +
    scale_y_discrete(expand=expansion(mult=c(0,0)),
                     position="left") +
    scale_fill_gradientn(limits=c(max(abs(range(DF.Zscore_2$zscore))),
                                  -max(abs(range(DF.Zscore_2$zscore)))),
                         breaks=c(-2, 0, 2),
                         colors=
                           #c("#A6D3C5", "#F7F7F7", "#7A3E9D"),
                           c("#A6D3C5", "#F7F7F7", "#C05A9E"),
                           #c("#5B7FB0", "#F7F7F7", "#7A3E9D"),
                           #c("#5B6F9E", "#F7F7F7", "#B8C76A"),
                           #c("#4dac26","#b8e186","#f7f7f7","#f1b6da","#d01c8b"),
                         guide=guide_colorbar(#direction="horizontal",
                           title.position="top",
                           title.hjust=0.5,
                           frame.colour="black",
                           ticks.colour="black",
                           barwidth = 15,
                           barheight = 1)) +
    coord_fixed(ratio=0.9) + # 1.15) +
    theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
          legend.position = "bottom",
          axis.text = element_text(face="bold", color="black"),
          axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5, 
                                     color="black", size=21),
          axis.text.y = element_text(color="black", size=21),
          legend.text = element_text(face = "bold", color="black", size = 15),
          legend.title = element_text(face = "bold", color="black", size = 12),
          legend.background = element_rect(fill="transparent", color=NA),
          plot.background = element_rect(fill="transparent", color=NA),
          panel.background = element_rect(fill="transparent", color=NA),
          panel.border = element_rect(fill="transparent", color = "black"),
          axis.line = element_blank(),
          panel.grid = element_blank())
  ggsave(plot = p_heat_Presentation,
         file="ssGSVA_heatmap_presentatiion_2.png", 
         dpi=500, width=10, height=8, bg = "transparent" )

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
  # 3. Violin plot : ssGSVA scores of sub-clusters
  if ( identical(colnames(SeuObj.CAF_0), rownames(DF.gsvaScore_0)) ){ 
    SeuObj.CAF_0@meta.data = cbind(SeuObj.CAF_0@meta.data, DF.gsvaScore_0) 
    DF.gsvaScore_2 = data.frame(SeuObj.CAF_0@reductions[["umap"]]@cell.embeddings[,c("umap_1","umap_2")],
                                "nCount_RNA"=SeuObj.CAF_0$nCount_RNA,
                                "nFeature_RNA"=SeuObj.CAF_0$nFeature_RNA,
                                "subclust"=SeuObj.CAF_0$seurat_clusters,
                                DF.gsvaScore_0)
  }

  Func.Violin = function(ScoreName){
    DF = 
      DF.gsvaScore_2[,c("subclust", ScoreName)] %>% 
      set_colnames(c("subclust", "Target")) %>% 
      dplyr::mutate(subclust = factor(subclust, levels = Vec.SubclustOrder))
    NamedVec.colors = setNames(cols2, as.character(0:8))
    #DF.Average = summarise(group_by(DF, subclust),
    #                       Mean = mean(Target)) %>% 
    #  dplyr::arrange(Mean)
    #DF$subclust = factor(DF$subclust, levels=as.character(DF.Average$subclust))
    ggplot(DF, aes(x=subclust, y=Target, fill=subclust)) + 
      geom_violin(trim=TRUE, #adjust=2
                  scale = "width",
                  width = 0.6) + 
      geom_boxplot(width=0.2, fill="white", outlier.shape = NA) +
      labs(x=NULL, y=NULL,
           subtitle=
             str_remove(ScoreName, pattern="_ssGSVA") %>% 
             str_replace(pattern="HALLMARK_", replacement="HALLMARK\n") %>% 
             str_replace(pattern="KEGG_", replacement="KEGG\n") %>% 
             str_replace(pattern="REACTOME_", replacement="REACTOME\n") %>% 
             str_replace(pattern="Orig", replacement="_original") %>% 
             str_replace(pattern="MYC_TARGETS_V", replacement="MYC_v") %>% 
             str_replace(pattern="INTERFERON", replacement="IFN")) +
      scale_x_discrete(labels=
                         setNames(paste0("CAF-",(0:8)),as.character(0:8))) +
      scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8)) +
      scale_fill_manual(values=NamedVec.colors)+ 
      coord_cartesian(ylim = c(NA,
                               quantile(DF$Target, 0.999, na.rm = TRUE))) +
      theme(legend.position="none",
            axis.text.x = element_text(face="bold", color="black", size=12, angle=45, hjust=1),
            axis.text.y = element_text(face="bold", color="black", size=10),
            #axis.title = element_text(face = "bold"),
            #legend.text = element_text(face = "bold"),
            #legend.title = element_text(face = "bold"),
            plot.background = element_rect(fill="transparent", color=NA),
            panel.background = element_rect(fill="transparent", color=NA),
            panel.grid = element_blank(),
            axis.line = element_line(color = "black"))
  }
  VlnPlot.GSVAbySubclust = 
    ggarrange(
      Func.Violin(ScoreName="myCAF_ssGSVA"),
      Func.Violin(ScoreName="iCAF_ssGSVA"),
      ncol=2, nrow=1, align="hv")
  ggsave(VlnPlot.GSVAbySubclust,
         file=paste0(DirInteg,"[Figure][SubClust_ssGSVAscores_Violin_ForPaper2][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                     "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"), 
         dpi=400, width=8, height=3, bg = "transparent" )
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
  # 4. Feature plot 
  #Vec.ColorScale = c("#000004","#000004","#000004", "#1b0c41", "#51116c", "#8c2964", "#c84646", "#f9c629", "#f9c629", "#f9c629")
  Vec.ColorScale = c("#253494", "#2c7fb8","#41b6c4", "#a1dab4", "#ffffcc", "#fecc5c", "#fd8d3c","#f03b20", "#bd0026")
  Func.FeatPlotGSVA = function(ColumnName, MidPoint){
    DF.Score_0 <-
      DF.gsvaScore_2[ , c("umap_1","umap_2",paste0(ColumnName,"_ssGSVA"))] %>% 
      set_colnames(c("umap_1","umap_2","TargetSignature")) %>% 
      dplyr::arrange(TargetSignature)
    Plot <-
      ggplot(DF.Score_0, 
             aes(x=umap_1, y=umap_2, color=TargetSignature)) + 
      geom_point(size=0.1) + 
      labs(x = "UMAP 1", y = "UMAP 2",
           color="ssGSEA\nscore") + 
      coord_fixed(ratio=1) + 
      theme(axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.title = element_text(face="bold", color="black", size=18),
            legend.text = element_text(face="bold", color="black", size=15),
            legend.title = element_text(face = "bold", color="black", size=15),
            legend.box.margin = margin(t=0, r=0, b=0, l=0),
            plot.background = element_rect(fill="transparent", color=NA),
            legend.background = element_rect(fill="transparent", color=NA),
            panel.background = element_rect(fill="transparent", color=NA),
            panel.grid = element_blank(),
            axis.line = element_line(color = "black"),
            plot.margin = margin(0,1,0,1)) + 
      scale_color_gradientn(
        colors=Vec.ColorScale, 
      #scale_color_viridis(
      #  option="viridis",
        values=scales::rescale(
          quantile(range(DF.Score_0$TargetSignature), 
                   probs=c(quantile(c(0, MidPoint),probs=c(0, 0.25, 0.5, 0.75)),
                           quantile(c(MidPoint, 1),probs=c(0, 0.25, 0.5, 0.75, 1.0))))),
        limits=range(DF.Score_0$TargetSignature),
        breaks=c(0.0, 0.2, 0.4, 0.6),
        guide=guide_colorbar(frame.colour="black", ticks.colour="black"))
    return(Plot)
  }
  FeatPlot.gsvaScore = 
    ggarrange(
      Func.FeatPlotGSVA("myCAF", MidPoint.myCAF),
      Func.FeatPlotGSVA("iCAF", MidPoint.iCAF),
      ncol=2, nrow=1, align="hv")
  ggsave(plot=FeatPlot.gsvaScore,
         file=paste0(DirInteg,"[Figure][SubClust_ssGSVAscores_FeaturePlot_ForPaper][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                     "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,"2.png"), 
         dpi=500, width=10, height=5, bg = "transparent" )
  #  FeatPlot.gsvaScore.Byclust = ggarrange(
  #    FeatPlot.Proliferation + facet_wrap(subclust ~ ., nrow=1),
  #    FeatPlot.BuffaOrig + facet_wrap(subclust ~ ., nrow=1),
  #    FeatPlot.WinterOrig + facet_wrap(subclust ~ ., nrow=1),
  #    FeatPlot.myCAF + facet_wrap(subclust ~ ., nrow=1),
  #    FeatPlot.iCAF + facet_wrap(subclust ~ ., nrow=1),
  #    FeatPlot.apCAF + facet_wrap(subclust ~ ., nrow=1),
  #    ncol=1, align="hv") %>% 
  #    annotate_figure(top=CommonTitle.Sub, 
  #                    bottom=paste0(CommonCaption.Sub, "\nMidpoints of scores are shifted; Proliferation:",MidPoint.Proliferation,
  #                                  ", BuffaOriginal:",MidPoint.BuffaOrig,", WinterOriginal:",MidPoint.WinterOrig,
  #                                  ", myCAF:",MidPoint.myCAF,", iCAF:",MidPoint.iCAF,", apCAF:",MidPoint.apCAF) )
  #  ggsave(plot=FeatPlot.gsvaScore.Byclust,
  #         file=paste0(DirInteg,"[Figure][SubClust_ssGSVAscores_FeaturePlot_BySignature][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
  #                     "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"), 
  #         dpi=300, width=16, height=12, bg = "white" )
  
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
  # 5. Corralation
  #library(ggcorrplot)
  #ggcorrplot(corr=MT.Correlation_0, tl.cex=3, type="upper", show.diag=TRUE)
  library(corrplot)
  Vec.SignatureNames = c("BuffaOrig_ssGSVA","WinterOrig_ssGSVA","HALLMARK_HYPOXIA_ssGSVA",
                         "myCAF_ssGSVA","iCAF_ssGSVA","apCAF_ssGSVA",
                         "Proliferation_ssGSVA","HALLMARK_E2F_TARGETS_ssGSVA","HALLMARK_G2M_CHECKPOINT_ssGSVA",
                         "HALLMARK_MITOTIC_SPINDLE_ssGSVA",
                         "HALLMARK_MYC_TARGETS_V1_ssGSVA","HALLMARK_MYC_TARGETS_V2_ssGSVA",
                         "HALLMARK_P53_PATHWAY_ssGSVA",
                         "HALLMARK_INTERFERON_ALPHA_RESPONSE_ssGSVA",
                         "HALLMARK_INTERFERON_GAMMA_RESPONSE_ssGSVA",
                         "KEGG_Cell_cycle_ssGSVA","REACTOME_Cellular_senescence_ssGSVA")
  DF.gsvaScore_1 = DF.gsvaScore_0[, Vec.SignatureNames]
  MT.Correlation_0 = cor(DF.gsvaScore_1)
  colors=brewer.pal(n=10, name="RdBu")
  colors=rev(c(colors[1:4], "#ffffff", "#ffffff", colors[7:10]))
  corrplot(MT.Correlation_0, tl.cex=0.5, ,method="pie", type="upper", tl.srt=45,
           col=colors)
  Vec.ShowSignatures = c("BuffaOrig_ssGSVA","WinterOrig_ssGSVA","HALLMARK_HYPOXIA_ssGSVA",
                         #"KEGG_CellCycle_ssGSVA",
                         "HALLMARK_G2M_CHECKPOINT_ssGSVA","HALLMARK_E2F_TARGETS_ssGSVA","HALLMARK_MITOTIC_SPINDLE_ssGSVA",
                         "HALLMARK_MYC_TARGETS_V1_ssGSVA","HALLMARK_MYC_TARGETS_V2_ssGSVA"#,
                         #"myCAF_ssGSVA","iCAF_ssGSVA","apCAF_ssGSVA"
  )
  MT.Correlation_1 = MT.Correlation_0[Vec.ShowSignatures,Vec.ShowSignatures]
  png(paste0(DirInteg,"[Figure][SubClust_ssGSVAscores_Corrplot][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
             "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"),
      width=3000, height=2000, res=400)
  corrplot(MT.Correlation_1, method="pie", type="upper", 
           tl.cex=0.8, 
           tl.srt=45, 
           cl.length = 11,   # 例：-1〜1 を 0.2刻み相当の段数
           cl.cex = 0.8,
           cl.ratio = 0.25,
           col=colors)
  dev.off()

Func.gsvaScore.Sub(Dim1=20, Res1=1.0, Dim2=30, Res2=0.5,
                   method=="EuclideanWardd2",  # or method="SpearmanAverage"
                   MidPoint.Proliferation=0.65,
                   MidPoint.KEGG_CellCycle=0.7,
                   MidPoint.E2F=0.50,
                   MidPoint.G2M=0.65,
                   MidPoint.MITOTIC_SPINDLE=0.75,
                   MidPoint.MYC_TARGETS_V1=0.65,
                   MidPoint.MYC_TARGETS_V2=0.65,
                   MidPoint.BuffaOrig=0.70,
                   MidPoint.WinterOrig=0.70,
                   MidPoint.HYPOXIA=0.6,
                   MidPoint.myCAF=0.70,
                   MidPoint.iCAF=0.6,
                   MidPoint.apCAF=0.999)




##
##
##  8. Cell cycle scoring                        ####
SeuObj.CAF_0 = readRDS(paste0(DirInteg,"Objects/[SeuObj_ExtractedCAF][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                              QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2,".rds") )
SeuObj.CAF_1 = SeuObj.CAF_0 %>% 
  CellCycleScoring(
    s.features=intersect(cc.genes$s.genes, Vec.AllGenes),
    g2m.features=intersect(cc.genes$g2m.genes, Vec.AllGenes),
    set.ident=T,
    seed=1234)
DF.CCscores_0 <-  data.frame(
  "subclust"=SeuObj.CAF_1$seurat_clusters,
  "S.Score"=SeuObj.CAF_1$S.Score,
  "G2M.Score"=SeuObj.CAF_1$G2M.Score,
  "phase"=SeuObj.CAF_1$Phase)
DF.MeanScore_0 <- DF.CCscores_0 %>% 
  group_by(subclust) %>% 
  summarise(Mean.S=mean(S.Score),
            Mean.G2M=mean(G2M.Score))
# violin plot
CommonTheme.Vln <- theme(
  legend.position="none",
  axis.text = element_text(face="bold", color="black", size=25),
  axis.title = element_text(face="bold", color="black", size=10),
  legend.text = element_text(face = "bold"),
  legend.title = element_text(face = "bold"),
  plot.background = element_rect(fill="transparent", color=NA),
  panel.background = element_rect(fill="transparent", color=NA),
  panel.grid = element_blank(),
  panel.border = element_rect(fill=NA, color = "black"))
DF.CCscores_1 <- 
  DF.CCscores_0 %>% 
  dplyr::mutate(subclust = factor(subclust, levels = Vec.SubclustOrder))
Vln.S <- 
  DF.CCscores_1 %>% 
  ggplot(aes(x=subclust, y=S.Score, fill=subclust)) + 
  geom_violin(scale="width", linewidth=0.8) + 
  geom_boxplot(width=0.15, fill="white", outliers=F, linewidth=0.8) +
  scale_x_discrete(labels=function(x) paste0("CAF-", x)) +
  scale_fill_manual(values=cols2) +
  labs(x=NULL) + CommonTheme.Vln + theme(axis.text.x = element_text(angle=45, hjust=1))
Vln.G2M <- 
  DF.CCscores_1 %>% 
  ggplot(aes(x=subclust, y=G2M.Score, fill=subclust)) + 
  geom_violin(scale="width", linewidth=0.8) + 
  geom_boxplot(width=0.15, fill="white", outliers=F, linewidth=0.8) +
  scale_x_discrete(labels=function(x) paste0("CAF-", x)) +
  scale_y_continuous(breaks=c(-0.3, 0, 0.3, 0.6)) +
  scale_fill_manual(values=cols2) +
  labs(x=NULL) + CommonTheme.Vln + theme(axis.text.x = element_text(angle=45, hjust=1))
Plot.CCscore.Vln = 
  ggarrange(
    Vln.S, Vln.G2M,
    ncol=1, align="hv")
ggsave(plot=Plot.CCscore.Vln,
       file=paste0(DirInteg,"[Figure][SubClust_CellCycleScores_BySubclust_VlnPlot_ForPaper][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                   "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"), 
       dpi=500, width=8, height=10, bg = "transparent" )
# scatter plot
Plot.Scores.Scatter <- 
  ggplot(DF.CCscores_0, aes(S.Score, G2M.Score)) +
  geom_point(data = subset(DF.CCscores_0, subclust !="8"),
             shape = 16,
             color="grey80", 
             size=2.0, alpha=0.5) +
  geom_point(data=subset(DF.CCscores_0, subclust == "8"), 
             shape=21, 
             fill="#C00000", color="black", 
             size=2.0, stroke=0.05) +
  labs(title=NULL) +
  scale_y_continuous(breaks=c(-0.3, 0, 0.3, 0.6)) +
  theme(aspect.ratio = 1.3,
        legend.position="none",
        axis.text = element_text(face="bold", color="black", size=15),
        axis.title = element_text(face="bold", color="black"),
        plot.background = element_rect(fill="transparent", color=NA),
        panel.background = element_rect(fill="transparent", color=NA),
        panel.border = element_rect(fill=NA, color=NA),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"))
ggsave(plot=Plot.Scores.Scatter,
       file=paste0(DirInteg,"[Figure][SubClust_CellCycleScores_BySubclust_ScatterPlot_ForPaper][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                   "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"), 
       dpi=500, width=4, height=5.5, bg = "transparent" )

##  9. Hypoxia mapping                           ####

# 1. Cell label
List.InitialClustData = readRDS(
  file=paste0(DirInteg,"[MetaDataList][InitialClust][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
              "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".rds"))
DF.ClustData_0 <- List.InitialClustData$MetaData[ , c("cell_id", "sample", "seurat_clusters")]
DF.ClustData_1 <- 
  DF.ClustData_0 %>% 
  dplyr::mutate(
    Annot = Vec.Annot[as.character(seurat_clusters)],
    seurat_clusters = paste0("Ini-",seurat_clusters)) %>% 
  rename(IniClust = seurat_clusters)
SeuObj.CAF_0 = readRDS(
  paste0(DirInteg,"Objects/[SeuObj_ExtractedCAF][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
         QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2,".rds") )
DF.CAFclustData <- data.frame(
  cell_id = colnames(SeuObj.CAF_0),
  SubClust = paste0("CAF-", SeuObj.CAF_0$seurat_clusters)
)
DF.ClustData_2 <- 
  left_join(DF.ClustData_1, DF.CAFclustData, by = "cell_id") %>% 
  dplyr::mutate(SubClust = ifelse(is.na(SubClust), Annot, SubClust))

# 2. Insert cell label into individual sample seurat object
for(i in 1:6){
  TX <- TXnumInteg[i]
  Func.ReadSeuObjForInteg = function(TXnum, Mag, nFeatRNA, nCountRNA, Dim1, Res1){
    TXdl=if(TXnum %in% c("27","28")){ "27_28" }else{ TXnum }
    Directory = paste0("/Volumes/Extreme SSD/Analysis/Data/TX5K_",TXdl,"/")
    Dir0 = "/Volumes/Extreme SSD/Analysis/Data/"
    SeuratObj = readRDS(paste0(Directory,"Objects/[SeuObj][TX5K_",TXnum,"]_Countable_Mag",Mag,
                               "_nFeat",nFeatRNA[1],"-",nFeatRNA[2],
                               "_nCount",nCountRNA[1],"-",nCountRNA[2],".rds") )
    colnames(SeuratObj) = paste(str_remove(SeuratObj$orig.ident, pattern="_"), 
                                colnames(SeuratObj), sep="_")
    return(SeuratObj)}
  SeuObj.AllCelltype_0 <- Func.ReadSeuObjForInteg(TXnum=TX, Mag=1, nFeatRNA=c(100, 900), nCountRNA=c(100, 1800))
  DF.ClustData.indiv.rmUC_0 <- 
    DF.ClustData_2 %>% 
    dplyr::filter(sample == TX) %>% 
    dplyr::filter(Annot != "Unclassified") %>% 
    column_to_rownames(var = "cell_id")
  SeuObj.AllCelltype_1 <- subset(
    SeuObj.AllCelltype_0, 
    cells = rownames(DF.ClustData.indiv.rmUC_0))
  if(identical(colnames(SeuObj.AllCelltype_1), rownames(DF.ClustData.indiv.rmUC_0))){
    SeuObj.AllCelltype_1$IniClust <- as.character(DF.ClustData.indiv.rmUC_0$IniClust)
    SeuObj.AllCelltype_1$Annot <- as.character(DF.ClustData.indiv.rmUC_0$Annot)
    SeuObj.AllCelltype_1$SubClust <- as.character(DF.ClustData.indiv.rmUC_0$SubClust)
    SeuObj.AllCelltype_2 <- 
      NormalizeData(SeuObj.AllCelltype_1, verbose=FALSE)
    saveRDS(SeuObj.AllCelltype_2,
            file=paste0(DirInteg,"HypoxiaMapping/[SeuObj_TX",TX,"][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                        QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".rds"))
  }
}

# 3. ssGSEA
library(msigdbr)
library(GSVA)
# 3-1. Get gene sets
CategoryH <- msigdbr(species = "Homo sapiens", collection = "H")
CategoryH_list = CategoryH %>% split(.$gs_name) %>%
  lapply(function(x) unique(x$gene_symbol))
CategoryC2 = msigdbr(species="Homo sapiens",  category="C2")
Pathways = list(
  BuffaOrig = Vec.BuffaOrig,
  WinterOrig = Vec.WinterOrig,
  HMhypoxia = 
    dplyr::filter(CategoryH, gs_name == "HALLMARK_HYPOXIA") %>% 
    dplyr::pull(gene_symbol),
  BuffaMeta = 
    dplyr::filter(CategoryC2, gs_name == "BUFFA_HYPOXIA_METAGENE") %>% 
    dplyr::pull(gene_symbol),
  WinterMeta = 
    dplyr::filter(CategoryC2, gs_name == "WINTER_HYPOXIA_METAGENE") %>% 
    dplyr::pull(gene_symbol))
Pathways.FilterUnobservedGenes = 
  lapply(Pathways, function(gs) {
    intersect(gs, Vec.AllGenes)
  })

# 3-2.Calculate GSVA
for(i in 1:6){
  TX <- TXnumInteg[i]
  SeuObj.AllCelltype_2 <- readRDS(
    file=paste0(DirInteg,"HypoxiaMapping/[SeuObj_TX",TX,"][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".rds"))
  expr <-  
    GetAssayData(
      SeuObj.AllCelltype_2, assay="RNA", layer="data") %>% 
    as.matrix()
  set.seed(1234)
  param = ssgseaParam(
    exprData = expr, 
    geneSets = Pathways.FilterUnobservedGenes,
    minSize = 1, 
    maxSize = Inf, 
    normalize = TRUE)  # default
  MT.ssGSEAres_0 = gsva(param, verbose=TRUE)
  MT.ssGSEAres_1 <- t(MT.ssGSEAres_0)
  colnames(MT.ssGSEAres_1) = paste0(colnames(MT.ssGSEAres_1), "_ssGSEA")
  DF.ssGSEAres_2 <-
    as.data.frame(MT.ssGSEAres_1) %>% 
    rownames_to_column(var="full_id")
  if(identical(DF.ssGSEAres_2$full_id, 
               colnames(SeuObj.AllCelltype_2))){
    DF.ssGSEAres_3 <- data.frame(
      "sample" = str_remove(DF.ssGSEAres_2$full_id, pattern="_.*"),
      "cell_id" = str_remove(DF.ssGSEAres_2$full_id, pattern="TX5K.*_"),
      "IniClust" = SeuObj.AllCelltype_2$IniClust,
      "Annot" = SeuObj.AllCelltype_2$Annot,
      "SubClust" = SeuObj.AllCelltype_2$SubClust,
      "X" = SeuObj.AllCelltype_2$X,
      "Y" = SeuObj.AllCelltype_2$Y,
      stringsAsFactors = FALSE) %>% 
      cbind(DF.ssGSEAres_2)
  }
  write.csv(
    DF.ssGSEAres_3,
    file=paste0(DirInteg,"HypoxiaMapping/[ssGSEAscoreTable_TX",TX,"][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".csv"),
    row.names=FALSE)
}
DF.ssGSEAres_trans_2 =
  DF.ssGSEAres_4 %>% 
  dplyr::mutate(cluster = factor(cluster, levels = names(Vec.OrderedCnumToColor)))
ggplot(DF.ssGSEAres_trans_2,
       aes(x = cluster, y = BuffaOrig)) +
  geom_violin(aes(fill = annot))
ggplot(DF.ssGSEAres_trans_2,
       aes(x = annot, y = BuffaOrig)) +
  geom_violin(aes(fill = annot))

# 4. Cell lineage normalize / Cluster normalize
for(i in 1:6){
  TX <- TXnumInteg[i]
  DF.ssGSEAres_3 <- 
    read.csv(
      file=paste0(DirInteg,"HypoxiaMapping/[ssGSEAscoreTable_TX",TX,"][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                  QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".csv"))
  DF.ssGSEAres_4.Normalized <- 
    DF.ssGSEAres_3 %>% 
    group_by(Annot) %>% 
    dplyr::mutate(
      X = as.numeric(X),
      Y = as.numeric(Y),
      Winter_z_celltype = as.numeric(scale(WinterOrig_ssGSEA)),
      Buffa_z_celltype  = as.numeric(scale(BuffaOrig_ssGSEA)),
      Hallmark_z_celltype = as.numeric(scale(HMhypoxia_ssGSEA))
    ) %>% 
    ungroup() %>% 
    group_by(IniClust) %>%
    mutate(
      Winter_z_cluster = as.numeric(scale(WinterOrig_ssGSEA)),
      Buffa_z_cluster  = as.numeric(scale(BuffaOrig_ssGSEA)),
      Hallmark_z_cluster = as.numeric(scale(HMhypoxia_ssGSEA))
    ) %>%
    ungroup()
  write.csv(
    DF.ssGSEAres_4.Normalized,
    file=paste0(DirInteg,"HypoxiaMapping/[ssGSEAscoreTable_Normalized_TX",TX,"][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".csv"),
    row.names=FALSE)
}

# 5. Local Hypoxia score derived Non-CAF cells
for(i in 1:6){
  if(i == 1){
    DF.LocalNonCAFscoreOfCAF <- data.frame()
  }
  TX <- TXnumInteg[i]
  DF.ssGSEAres_4.Normalized <- 
    read.csv(
      file=paste0(DirInteg,"HypoxiaMapping/[ssGSEAscoreTable_Normalized_TX",TX,"][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                  QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".csv"))
  DF.ssGSEAres_4.Normalized_Grid <- 
    DF.ssGSEAres_4.Normalized %>% 
    rowwise() %>% 
    dplyr::mutate(
      GridX.edge = floor(X / Num.GridLen),
      GridX.center = GridX.edge + 0.5,
      GridY.edge = floor(Y / Num.GridLen),
      GridY.center = GridY.edge + 0.5,
      Grid_id = paste(GridX.edge, GridY.edge, sep="_"))
  DF.Grid.nonCAF <-  
    DF.ssGSEAres_4.Normalized_Grid %>% 
    dplyr::filter(Annot != "CAF")
  Num.MinCells <- 30
  DF.Grid.nonCAF.sum_0 = 
    DF.Grid.nonCAF %>% 
    group_by(Grid_id, GridX.center, GridY.center) %>% 
    summarise(
      n_cells = n(),
      mean_Winter_celltype = mean(Winter_z_celltype),
      mean_Buffa_celltype = mean(Buffa_z_celltype),
      mean_Hallmark_celltype = mean(Hallmark_z_celltype),
      mean_Winter_cluster = mean(Winter_z_cluster),
      mean_Buffa_cluster = mean(Buffa_z_cluster),
      mean_Hallmark_cluster = mean(Hallmark_z_cluster),
      .groups = "drop") %>% 
    dplyr::mutate(
      mean_Winter_celltype = ifelse(n_cells < Num.MinCells, NA, mean_Winter_celltype),
      mean_Buffa_celltype = ifelse(n_cells < Num.MinCells, NA, mean_Buffa_celltype),
      mean_Hallmark_celltype = ifelse(n_cells < Num.MinCells, NA, mean_Hallmark_celltype),
      mean_Winter_cluster = ifelse(n_cells < Num.MinCells, NA, mean_Winter_cluster),
      mean_Buffa_cluster = ifelse(n_cells < Num.MinCells, NA, mean_Buffa_cluster),
      mean_Hallmark_cluster = ifelse(n_cells < Num.MinCells, NA, mean_Hallmark_cluster) 
    )
  RangeGridX = range(DF.Grid.nonCAF$GridX.edge)
  RangeGridY = range(DF.Grid.nonCAF$GridY.edge)
  Vec.AllGridId = as.vector(outer(seq(RangeGridX[1], RangeGridX[2]), 
                                  seq(RangeGridY[1], RangeGridY[2]), paste, sep="_"))
  Vec.DroppedGridIds = setdiff(Vec.AllGridId, unique(DF.Grid.nonCAF.sum_0$Grid_id))
  DF.naGrid = 
    DF.Grid.nonCAF.sum_0[rep(NA_integer_, length(Vec.DroppedGridIds)), ] %>% 
    dplyr::mutate(Grid_id = Vec.DroppedGridIds,
                  GridX.center = str_remove(Grid_id, pattern="_.*"),
                  GridX.center = as.numeric(GridX.center) + 0.5,
                  GridY.center = str_remove(Grid_id, pattern=".*_"),
                  GridY.center = as.numeric(GridY.center) + 0.5)
  DF.Grid.nonCAF.sum_1 = 
    rbind(DF.Grid.nonCAF.sum_0, DF.naGrid) %>%
    dplyr::mutate(
      GridX.center_um = GridX.center * Num.GridLen,
      GridY.center_um = GridY.center * Num.GridLen)
  RangeX = c(min(DF.ssGSEAres_4.Normalized$X), max(DF.ssGSEAres_4.Normalized$X))
  RangeY = c(min(DF.ssGSEAres_4.Normalized$Y), max(DF.ssGSEAres_4.Normalized$Y))
  
  Pathways <- c("Winter", "Hallmark", "Buffa")
  Normalize <- "celltype" #/ "cluster"
  CircleClust <- prolifCAFclust
  for(j in 1:1){
    Pathway <- Pathways[j]
    DF.Grid.nonCAF.sum_2 <- 
      DF.Grid.nonCAF.sum_1[ , c("GridX.center_um", "GridY.center_um", paste0("mean_", Pathway,"_",Normalize))] %>% 
      set_colnames(c("GridX.center_um", "GridY.center_um", "TargetPathway"))
    Xen.HypoxiaHeatmap = 
      ggplot() +
      geom_tile(
        data=DF.Grid.nonCAF.sum_2,
        aes(x=GridX.center_um, y=GridY.center_um, fill=TargetPathway),
        width = Num.GridLen,
        height = Num.GridLen,
        color = "gray60") + 
      geom_point(
        data=subset(DF.ssGSEAres_4.Normalized, SubClust == paste0("CAF-", CircleClust)),
        aes(x = X, y = Y), color="black",
        size=3.0, shape=21, stroke=1.5, fill=NA) +
      labs(subtitle=paste0("Local mean of Non-CAF derived signature score\n",
                           "TX5K_",TX,", CAF-",CircleClust, "\nPathway:",Pathway,"\nNormalize:Among each ",Normalize),
           fill = NULL, x = NULL, y = NULL) +
      scale_x_continuous(expand=expansion(mult=c(0,0))) +
      scale_y_reverse(expand=expansion(mult=c(0,0))) + 
      coord_cartesian(xlim = RangeX, ylim = RangeY,
                      ratio=1) +
      scale_fill_gradientn(colors=c("#bd0026", "#f03b20", "#fd8d3c", "#fecc5c", "gray95",
                                    "#a1dab4", "#41b6c4", "#2c7fb8", "#253494"),
                           na.value="gray80",
                           limits=c(-max(abs(DF.Grid.nonCAF.sum_2$TargetPathway), na.rm=TRUE), 
                                    max(abs(DF.Grid.nonCAF.sum_2$TargetPathway), na.rm=TRUE)),
                           guide=guide_colorbar(
                             title.hjust=0.5,
                             title.vjust=0.5,
                             frame.colour="black",
                             ticks.colour="black")) +
      theme(legend.background = element_rect(fill = "transparent", color = NA),
            axis.text.x = element_text(face = "bold", color = "black"),
            axis.text.y = element_text(face = "bold", color = "black"),
            legend.text = element_text(face = "bold", color = "black"),
            plot.subtitle = element_text(size = 20, face = "bold"),
            plot.background = element_rect(fill = "transparent", color = NA),
            panel.background = element_rect(fill = "transparent", color = NA),
            panel.grid = element_blank(),
            panel.border = element_rect(fill = "transparent", color = "black") )
    ggsave(
      plot = Xen.HypoxiaHeatmap,
      filename = paste0(DirInteg,"HypoxiaMapping/MapFigure/[HypoxiaMap_TX",TX,"_CAF",CircleClust,"_",Pathway,"_",Normalize,"-norm][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                        QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".png"),
      width = diff(RangeX)/500, height = diff(RangeY)/500, dpi = 100, bg = "transparent")
    Xen.HypoxiaHeatmap.NoLgd <- 
      Xen.HypoxiaHeatmap +
      theme(plot.margin = margin(0, 0, 0, 0),
            plot.subtitle = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "none",
            title = element_blank())
    ggsave(
      plot = Xen.HypoxiaHeatmap.NoLgd,
      filename = paste0(DirInteg,"HypoxiaMapping/MapFigure/[HypoxiaMap(Panel)_TX",TX,"_CAF",CircleClust,"_",Pathway,"_",Normalize,"-norm][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                        QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".png"),
      width = diff(RangeX)/500, height = diff(RangeY)/500, dpi = 500, bg = "transparent")
    Xen.HypoxiaHeatmap.ForLgd <- 
      Xen.HypoxiaHeatmap + 
      labs(fill = NULL) + 
      theme(plot.margin = margin(0, 0, 0, 0), 
            legend.box.margin = margin(0, 0, 0, 0))
    if(Pathway == "Winter"){
      ggsave(
        plot = get_legend(Xen.HypoxiaHeatmap.ForLgd),
        filename = paste0(DirInteg,"HypoxiaMapping/MapFigure/[HypoxiaMap(Legend)_TX",TX,"_CAF",CircleClust,"_",Pathway,"_",Normalize,"-norm][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                          QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".png"),
        width = 1, height = 3, dpi = 500, bg = "transparent")
      saveRDS(Xen.HypoxiaHeatmap.NoLgd,
              file=paste0(DirInteg,"/XeniumView_Sub/",QCInfo.FileName,
                          "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,
                          "/ssGSVA/[Figure][ExtractedCAF_XenView_GSVA_HypoxicAreaAndProlifCAF_TX5K_",TX,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                          "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".rds"))
    }
  }
  DF.CAF_GridScore_0 <- 
    DF.ssGSEAres_4.Normalized_Grid %>% 
    dplyr::filter(Annot == "CAF")
  Vec.GridId.GridScore <- 
    setNames(DF.Grid.nonCAF.sum_1$mean_Winter_celltype, 
             DF.Grid.nonCAF.sum_1$Grid_id)
  DF.CAF_GridScore_1 <- 
    DF.CAF_GridScore_0 %>% 
    dplyr::mutate(GridScore = Vec.GridId.GridScore[Grid_id])
  DF.CAF_GridScore_2 <- 
    DF.CAF_GridScore_1 %>% 
    dplyr::select(full_id, SubClust, GridScore) %>% 
    dplyr::mutate(SubClust = factor(SubClust, levels = paste0("CAF-", Vec.SubclustOrder)))
  DF.LocalNonCAFscoreOfCAF <- 
    rbind(DF.LocalNonCAFscoreOfCAF,
          DF.CAF_GridScore_2)
  if(i == 6){
    write.csv(DF.LocalNonCAFscoreOfCAF,
              file = paste0(DirInteg,"HypoxiaMapping/[ssGSEAscoreTable_Normalized_AllCases][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                            QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".csv"))
  }
}
DF.LocalNonCAF <- 
  DF.LocalNonCAFscoreOfCAF %>% 
  ungroup() %>% 
  dplyr::mutate(sample = str_remove(full_id, pattern = "_.*")) %>% 
  dplyr::filter(!is.na(GridScore)) %>% 
  ggplot(aes(x = SubClust, y = GridScore)) +
  labs(subtitle = paste0("TX",TX), x = NULL, y = "Local non-CAF hypoxia score") +
  geom_violin(aes(fill = SubClust), scale = "width") +
  geom_boxplot(outliers = FALSE, width = 0.3) +
  scale_fill_manual(values = setNames(as.character(cols2), paste0("CAF-", names(cols2)))) +
  theme(plot.background = element_rect(fill = "transparent", color = "transparent"),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "black"),
        panel.grid = element_blank()) +
  facet_wrap(. ~ sample, ncol = 3, scales = "free_y")


#################################################
##
##     Xenium mapping (New)
##  1. Initial clust                         ####

Coordinate.Area1 = list(
  "01" = c(2500, 2900), "02" = c(1300, 9500), "11" = c(2100, 3350), 
  "16" = c(4600, 5200), "15" = c(1250, 10550), "19" = c(7770, 6670))
Coordinate.Area2 <- list(
  "01" = c(6200, 3100), "02" = c(3800, 2500), "11" = c(2650, 3000),
  "16" = c(4800, 3400), "15" = c(5400, 7450), "19" = c(8000, 750))

##  2. Hypoxic / Normoxic area               ####

# define coorinate ( original coordinate !!!)
SideLength = 400
Coordinate.Hypo <- list(
  "01" = c(9400, 1200), 
  "02" = c(1400, 9400), 
  "11" = c(3400, 3600),
  "15" = c(2200, 7400),
  "16" = c(5200, 3600), 
  "19" = c(7600, 4000))
Coordinate.Normo <- list(
  "01" = c(3800, 2000),
  "02" = c(1600, 7600),
  "11" = c(2600, 3000),
  "15" = c(4600, 9200), 
  "16" = c(1800, 2000), 
  "19" = c(8400, 1200))

Func.HypoxiaHeatmap.Rect <-  function(TX){
  Xen.HypoxiaHeatmap = readRDS(
    file=paste0(DirInteg,"/XeniumView_Sub/",QCInfo.FileName,
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,
                "/ssGSVA/[Figure][ExtractedCAF_XenView_GSVA_HypoxicAreaAndProlifCAF_TX5K_",TX,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".rds"))
  Xen.HypoxiaHeatmap@layers$geom_point$aes_params$size = 5
  Xen.HypoxiaHeatmap@layers$geom_point$aes_params$stroke = 3
  LengthX = diff(Xen.HypoxiaHeatmap@coordinates$limits$x)
  LengthY = diff(Xen.HypoxiaHeatmap@coordinates$limits$y)
  if(LengthY > LengthX){ # long
    SaveWidth=9
    SaveHeight=9/LengthX*LengthY
  } else {  # Wide
    SaveWidth=9/LengthY*LengthX
    SaveHeight=9
  }
  ggsave(plot = 
           Xen.HypoxiaHeatmap +
           annotate("rect", 
                    xmin=Coordinate.Hypo[[TX]][1], 
                    xmax=Coordinate.Hypo[[TX]][1] + SideLength,
                    ymin=Coordinate.Hypo[[TX]][2], 
                    ymax=Coordinate.Hypo[[TX]][2] + SideLength, 
                    color="black", fill=NA, linewidth=4) +
           annotate("rect", 
                    xmin=Coordinate.Hypo[[TX]][1], 
                    xmax=Coordinate.Hypo[[TX]][1] + SideLength,
                    ymin=Coordinate.Hypo[[TX]][2], 
                    ymax=Coordinate.Hypo[[TX]][2] + SideLength, 
                    color="#253494", fill=NA, linewidth=3) +
           annotate("rect", 
                    xmin=Coordinate.Normo[[TX]][1], 
                    xmax=Coordinate.Normo[[TX]][1] + SideLength,
                    ymin=Coordinate.Normo[[TX]][2], 
                    ymax=Coordinate.Normo[[TX]][2] + SideLength, 
                    color="black", fill=NA, linewidth=4) +
           annotate("rect", 
                    xmin=Coordinate.Normo[[TX]][1], 
                    xmax=Coordinate.Normo[[TX]][1] + SideLength,
                    ymin=Coordinate.Normo[[TX]][2], 
                    ymax=Coordinate.Normo[[TX]][2] + SideLength, 
                    color="#f03b20", fill=NA, linewidth=3) +
           theme(plot.title = element_blank(),
                 plot.subtitle = element_blank(),
                 plot.caption = element_blank(),
                 axis.text = element_blank(),
                 axis.ticks = element_blank(),
                 axis.line = element_blank(),
                 plot.margin = margin(t=0, r=0, b=0, l=0, unit="cm"),
                 plot.background = element_rect(fill = "transparent", color = NA),
                 panel.background = element_rect(fill="transparent", color=NA),
                 panel.border = element_rect(fill = "transparent", color = "black")),
         file=paste0(DirInteg,"XeniumView_Sub/ForPaper/HypoxicNormoxic_HeatMap_TX",TX,
                     "_Hypo_x",paste(Coordinate.Hypo[[TX]], collapse="y"),
                     "_Normo_x",paste(Coordinate.Normo[[TX]], collapse="y"),
                     ".png"),
         width=SaveWidth, height=SaveHeight, dpi=500, bg="transparent")
}
Func.HypoxiaHeatmap.Rect(TX="01")
Func.HypoxiaHeatmap.Rect(TX="02")
Func.HypoxiaHeatmap.Rect(TX="11")
Func.HypoxiaHeatmap.Rect(TX="15")
Func.HypoxiaHeatmap.Rect(TX="16")
Func.HypoxiaHeatmap.Rect(TX="19")

Func.CellGroupHypoNormo <- function(TX){
  DataList = readRDS(
    file=paste0(DirInteg,"[MetaDataList][InitialClust][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".rds"))
  DF.MetaData_0 <- 
    DataList[["MetaData"]] %>% 
    dplyr::filter(sample == TX)
  DF.MetaData_1 <- 
    DF.MetaData_0 %>% 
    dplyr::mutate(
      area = 
        case_when(CoordX > Coordinate.Hypo[[TX]][1] &
                    CoordX < Coordinate.Hypo[[TX]][1] + SideLength &
                    CoordY > Coordinate.Hypo[[TX]][2] &
                    CoordY < Coordinate.Hypo[[TX]][2] + SideLength ~ "Hypoxia",
                  CoordX > Coordinate.Normo[[TX]][1] &
                    CoordX < Coordinate.Normo[[TX]][1] + SideLength &
                    CoordY > Coordinate.Normo[[TX]][2] &
                    CoordY < Coordinate.Normo[[TX]][2] + SideLength ~ "Normoxia",
                  TRUE ~ "Residual")) %>% 
    dplyr::filter(area %in% c("Hypoxia", "Normoxia"))
  DF.MetaData_2 <- 
    DF.MetaData_1 %>% 
    dplyr::select(c(cell_id, area)) %>% 
    dplyr::mutate(cell_id = str_remove(cell_id, pattern = "TX5K.._")) %>% 
    set_colnames(c("cell_id", "group"))
  write_csv(
    DF.MetaData_2,
    paste0(DirInteg,"XeniumView_Sub/ForPaper/[CellGroupTable][HypxicAndNormoxicLesion_TX5K_",TX,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
           "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".csv"))
}
Func.CellGroupHypoNormo(TX="01")
Func.CellGroupHypoNormo(TX="02")
Func.CellGroupHypoNormo(TX="11")
Func.CellGroupHypoNormo(TX="15")
Func.CellGroupHypoNormo(TX="16")
Func.CellGroupHypoNormo(TX="19")

Func.GenerateGreenAndSubclustcolored <- function(TX){
  XenView.Sub.FullLayer = readRDS(
    file=paste0(DirInteg,"XeniumView_Sub/[FigureRDS][SubClust_XeniumView_TX5K_",TX,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".rds"))
  XenView.Sub.AllGreenCAF = XenView.Sub.FullLayer
  XenView.Sub.AllGreenCAF@layers$geom_polygon...3 = NULL
  XenView.Sub.AllGreenCAF@layers$geom_polygon...2 = NULL
  saveRDS(
    XenView.Sub.AllGreenCAF,
    file=paste0(DirInteg,"XeniumView_Sub/ForPaper/[FigureRDS][SubClust_XeniumView_GreenCAF_TX5K_",TX,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".rds"))
  XenView.Subclust <- XenView.Sub.FullLayer
  XenView.Subclust@layers$geom_polygon...3 = NULL
  XenView.Subclust@layers$geom_polygon = NULL
  saveRDS(
    XenView.Subclust + 
      scale_fill_manual(values = setNames(cols2, paste0("CAF-",names(cols2)))),
    file=paste0(DirInteg,"XeniumView_Sub/ForPaper/[FigureRDS][SubClust_XeniumView_SubclustcoloredCAF_TX5K_",TX,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".rds"))
}
Func.GenerateGreenAndSubclustcolored(TX="01")
Func.GenerateGreenAndSubclustcolored(TX="02")
Func.GenerateGreenAndSubclustcolored(TX="11")
Func.GenerateGreenAndSubclustcolored(TX="15")
Func.GenerateGreenAndSubclustcolored(TX="16")
Func.GenerateGreenAndSubclustcolored(TX="19")

Func.HypoxiaGreenCAF <- function(TX){
  RDS <- readRDS(
    file=paste0(DirInteg,"XeniumView_Sub/ForPaper/[FigureRDS][SubClust_XeniumView_GreenCAF_TX5K_",TX,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".rds"))
  Length.OldX = diff(range(RDS@data$vertex_x))
  Length.OldY = diff(range(RDS@data$vertex_y))
  AspectRatio = max(Length.OldX, Length.OldY) / min(Length.OldX, Length.OldY)
  if(Length.OldX > Length.OldY){
    Coord.Hypo.Final = rev(Coordinate.Hypo[[TX]])   # rotated
    Coord.Normo.Final = rev(Coordinate.Normo[[TX]])   # rotated
    CoordY.Bar.Mag.Normo = Coord.Normo.Final[2]+SideLength*0.050
    CoordY.Bar.Mag.Hypo = Coord.Hypo.Final[2]+SideLength*0.050
    Coords.BarEnd.NotMag = c("x"=min(RDS@data$vertex_y) + Length.OldY*0.95,
                             "y"=min(RDS@data$vertex_x) + Length.OldX*0.05)
  } else {
    Coord.Normo.Final = Coordinate.Normo[[TX]]   # not rotated
    Coord.Hypo.Final = Coordinate.Hypo[[TX]]   # not rotated
    CoordY.Bar.Mag.Normo = Coord.Normo.Final[2]+SideLength*0.950
    CoordY.Bar.Mag.Hypo = Coord.Hypo.Final[2]+SideLength*0.950
    Coords.BarEnd.NotMag = c("x"=min(RDS@data$vertex_x) + Length.OldX*0.95,
                             "y"=min(RDS@data$vertex_y) + Length.OldY*0.95)
  }
  CommonTheme = 
    theme(plot.title = element_blank(),
          plot.subtitle = element_blank(),
          plot.caption = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          plot.margin = margin(t=0, r=0, b=0, l=0, unit="cm"),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill="black", color=NA),
          panel.border = element_rect(fill = "transparent", color = NA))
  # 1. Weakly magnified Green CAF
  Plot.greenCAF = 
    RDS + 
    annotate("rect", 
             xmin=Coord.Hypo.Final[1], xmax=(Coord.Hypo.Final[1])+SideLength,
             ymin=Coord.Hypo.Final[2], ymax=(Coord.Hypo.Final[2])+SideLength, 
             color="black", fill=NA, linewidth=2.8) +
    annotate("rect", 
             xmin=Coord.Hypo.Final[1], xmax=(Coord.Hypo.Final[1])+SideLength,
             ymin=Coord.Hypo.Final[2], ymax=(Coord.Hypo.Final[2])+SideLength, 
             color="skyblue", fill=NA, linewidth=2.1) +
    annotate("rect", 
             xmin=Coord.Normo.Final[1], xmax=(Coord.Normo.Final[1])+SideLength,
             ymin=Coord.Normo.Final[2], ymax=(Coord.Normo.Final[2])+SideLength, 
             color="black", fill=NA, linewidth=2.8) +
    annotate("rect", 
             xmin=Coord.Normo.Final[1], xmax=(Coord.Normo.Final[1])+SideLength,
             ymin=Coord.Normo.Final[2], ymax=(Coord.Normo.Final[2])+SideLength, 
             color="#f03b20", fill=NA, linewidth=2.1) +
    annotate("segment", color="white", linewidth=7.0,
             x=Coords.BarEnd.NotMag[["x"]]-1020, xend=Coords.BarEnd.NotMag[["x"]]+20, 
             y=Coords.BarEnd.NotMag[["y"]], yend=Coords.BarEnd.NotMag[["y"]]) +
    annotate("segment", color="black", linewidth=4.0,
             x=Coords.BarEnd.NotMag[["x"]]-1000, xend=Coords.BarEnd.NotMag[["x"]], 
             y=Coords.BarEnd.NotMag[["y"]], yend=Coords.BarEnd.NotMag[["y"]]) +
    CommonTheme
  ggsave(Plot.greenCAF,
         file=paste0(DirInteg,"XeniumView_Sub/ForPaper/HypoxicNormoxic_TX",TX,"_GreenCAF_WeakMag.png"), 
         width=10, height=10*AspectRatio, 
         dpi=500, bg="transparent", limitsize=FALSE)
  # 2. Strongly magnified Hypo/Normo area
  ggsave(plot = 
           RDS +
           annotate("segment", color="white", linewidth=6.0,
                    x=Coord.Normo.Final[1]+SideLength*0.698, xend=Coord.Normo.Final[1]+SideLength*0.952, 
                    y=CoordY.Bar.Mag.Normo, yend=CoordY.Bar.Mag.Normo) +
           annotate("segment", color="black", linewidth=4.0,
                    x=Coord.Normo.Final[1]+SideLength*0.700, xend=Coord.Normo.Final[1]+SideLength*0.950, 
                    y=CoordY.Bar.Mag.Normo, yend=CoordY.Bar.Mag.Normo) +
           coord_cartesian(
             xlim=c(Coord.Normo.Final[1], Coord.Normo.Final[1]+SideLength),
             ylim=c(Coord.Normo.Final[2], Coord.Normo.Final[2]+SideLength),
             ratio=1) +
           CommonTheme,
         width=10, height=10, dpi=500, bg="transparent", 
         file=paste0(DirInteg,"XeniumView_Sub/ForPaper/HypoxicNormoxic_TX",TX,"_GreenCAF_HighMag_Normoxic_x",paste(Coordinate.Normo[[TX]], collapse="y"),".png"))
  ggsave(plot = 
           RDS +
           annotate("segment", color="white", linewidth=6.0,
                    x=Coord.Hypo.Final[1]+SideLength*0.698, xend=Coord.Hypo.Final[1]+SideLength*0.952, 
                    y=CoordY.Bar.Mag.Hypo, yend=CoordY.Bar.Mag.Hypo) +
           annotate("segment", color="black", linewidth=4.0,
                    x=Coord.Hypo.Final[1]+SideLength*0.700, xend=Coord.Hypo.Final[1]+SideLength*0.950, 
                    y=CoordY.Bar.Mag.Hypo, yend=CoordY.Bar.Mag.Hypo) +
           coord_cartesian(
             xlim=c(Coord.Hypo.Final[1], Coord.Hypo.Final[1]+SideLength),
             ylim=c(Coord.Hypo.Final[2], Coord.Hypo.Final[2]+SideLength),
             ratio=1) +
           CommonTheme,
         width=10, height=10, dpi=500, bg="transparent", 
         file=paste0(DirInteg,"XeniumView_Sub/ForPaper/HypoxicNormoxic_TX",TX,"_GreenCAF_HighMag_Hypoxic_x",paste(Coordinate.Hypo[[TX]], collapse="y"),".png"))
}
Func.HypoxiaGreenCAF(TX="01")
Func.HypoxiaGreenCAF(TX="02")
Func.HypoxiaGreenCAF(TX="11")
Func.HypoxiaGreenCAF(TX="15")
Func.HypoxiaGreenCAF(TX="16")
Func.HypoxiaGreenCAF(TX="19")

Func.HypoxiaSubclustcoloredCAF <- function(TX){
  RDS <- readRDS(
    file=paste0(DirInteg,"XeniumView_Sub/ForPaper/[FigureRDS][SubClust_XeniumView_SubclustcoloredCAF_TX5K_",TX,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".rds"))
  Length.OldX = diff(range(RDS@data$vertex_x))
  Length.OldY = diff(range(RDS@data$vertex_y))
  AspectRatio = max(Length.OldX, Length.OldY)/min(Length.OldX, Length.OldY)
  if(Length.OldX > Length.OldY){
    Coord.Hypo.Final = rev(Coordinate.Hypo[[TX]])   # rotated
    Coord.Normo.Final = rev(Coordinate.Normo[[TX]])   # rotated
    CoordY.Bar.Mag.Normo = Coord.Normo.Final[2]+SideLength*0.050
    CoordY.Bar.Mag.Hypo = Coord.Hypo.Final[2]+SideLength*0.050
    Coords.BarEnd.NotMag = c("x"=min(RDS@data$vertex_y) + Length.OldY*0.95,
                             "y"=min(RDS@data$vertex_x) + Length.OldX*0.05)
  } else {
    Coord.Normo.Final = Coordinate.Normo[[TX]]   # not rotated
    Coord.Hypo.Final = Coordinate.Hypo[[TX]]   # not rotated
    CoordY.Bar.Mag.Normo = Coord.Normo.Final[2]+SideLength*0.950
    CoordY.Bar.Mag.Hypo = Coord.Hypo.Final[2]+SideLength*0.950
    Coords.BarEnd.NotMag = c("x"=min(RDS@data$vertex_x) + Length.OldX*0.95,
                             "y"=min(RDS@data$vertex_y) + Length.OldY*0.95)
  }
  CommonTheme = 
    theme(plot.title = element_blank(),
          plot.subtitle = element_blank(),
          plot.caption = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          plot.margin = margin(t=0, r=0, b=0, l=0, unit="cm"),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill="black", color=NA),
          panel.border = element_rect(fill = "transparent", color = NA))
  # 1. Weakly magnified Green CAF
  Plot.subclustcoloredCAF = 
    RDS + 
    annotate("rect", 
             xmin=Coord.Hypo.Final[1], xmax=(Coord.Hypo.Final[1])+SideLength,
             ymin=Coord.Hypo.Final[2], ymax=(Coord.Hypo.Final[2])+SideLength, 
             color="black", fill=NA, linewidth=2.8) +
    annotate("rect", 
             xmin=Coord.Hypo.Final[1], xmax=(Coord.Hypo.Final[1])+SideLength,
             ymin=Coord.Hypo.Final[2], ymax=(Coord.Hypo.Final[2])+SideLength, 
             color="skyblue", fill=NA, linewidth=2.1) +
    annotate("rect", 
             xmin=Coord.Normo.Final[1], xmax=(Coord.Normo.Final[1])+SideLength,
             ymin=Coord.Normo.Final[2], ymax=(Coord.Normo.Final[2])+SideLength, 
             color="black", fill=NA, linewidth=2.8) +
    annotate("rect", 
             xmin=Coord.Normo.Final[1], xmax=(Coord.Normo.Final[1])+SideLength,
             ymin=Coord.Normo.Final[2], ymax=(Coord.Normo.Final[2])+SideLength, 
             color="#f03b20", fill=NA, linewidth=2.1) +
    annotate("segment", color="white", linewidth=7.0,
             x=Coords.BarEnd.NotMag[["x"]]-1020, xend=Coords.BarEnd.NotMag[["x"]]+20, 
             y=Coords.BarEnd.NotMag[["y"]], yend=Coords.BarEnd.NotMag[["y"]]) +
    annotate("segment", color="black", linewidth=4.0,
             x=Coords.BarEnd.NotMag[["x"]]-1000, xend=Coords.BarEnd.NotMag[["x"]], 
             y=Coords.BarEnd.NotMag[["y"]], yend=Coords.BarEnd.NotMag[["y"]]) +
    CommonTheme
  ggsave(Plot.subclustcoloredCAF,
         file=paste0(DirInteg,"XeniumView_Sub/ForPaper/HypoxicNormoxic_TX",TX,"_SubclustcoloredCAF_WeakMag.png"), 
         width=10, height=10*AspectRatio, 
         dpi=500, bg="transparent", limitsize=FALSE)
  # 2. Strongly magnified Hypo/Normo area
  ggsave(plot = 
           RDS +
           annotate("segment", color="white", linewidth=6.0,
                    x=Coord.Normo.Final[1]+SideLength*0.698, xend=Coord.Normo.Final[1]+SideLength*0.952, 
                    y=CoordY.Bar.Mag.Normo, yend=CoordY.Bar.Mag.Normo) +
           annotate("segment", color="black", linewidth=4.0,
                    x=Coord.Normo.Final[1]+SideLength*0.700, xend=Coord.Normo.Final[1]+SideLength*0.950, 
                    y=CoordY.Bar.Mag.Normo, yend=CoordY.Bar.Mag.Normo) +
           coord_cartesian(
             xlim=c(Coord.Normo.Final[1], Coord.Normo.Final[1]+SideLength),
             ylim=c(Coord.Normo.Final[2], Coord.Normo.Final[2]+SideLength),
             ratio=1) +
           CommonTheme,
         width=10, height=10, dpi=500, bg="transparent", 
         file=paste0(DirInteg,"XeniumView_Sub/ForPaper/HypoxicNormoxic_TX",TX,"_SubclustcoloredCAF_HighMag_Normoxic_x",paste(Coordinate.Normo[[TX]], collapse="y"),".png"))
  ggsave(plot = 
           RDS +
           annotate("segment", color="white", linewidth=6.0,
                    x=Coord.Hypo.Final[1]+SideLength*0.698, xend=Coord.Hypo.Final[1]+SideLength*0.952, 
                    y=CoordY.Bar.Mag.Hypo, yend=CoordY.Bar.Mag.Hypo) +
           annotate("segment", color="black", linewidth=4.0,
                    x=Coord.Hypo.Final[1]+SideLength*0.700, xend=Coord.Hypo.Final[1]+SideLength*0.950, 
                    y=CoordY.Bar.Mag.Hypo, yend=CoordY.Bar.Mag.Hypo) +
           coord_cartesian(
             xlim=c(Coord.Hypo.Final[1], Coord.Hypo.Final[1]+SideLength),
             ylim=c(Coord.Hypo.Final[2], Coord.Hypo.Final[2]+SideLength),
             ratio=1) +
           CommonTheme,
         width=10, height=10, dpi=500, bg="transparent", 
         file=paste0(DirInteg,"XeniumView_Sub/ForPaper/HypoxicNormoxic_TX",TX,"_SubclustcoloredCAF_HighMag_Hypoxic_x",paste(Coordinate.Normo[[TX]], collapse="y"),".png"))
}
Func.HypoxiaSubclustcoloredCAF(TX="01")
Func.HypoxiaSubclustcoloredCAF(TX="02")
Func.HypoxiaSubclustcoloredCAF(TX="11")
Func.HypoxiaSubclustcoloredCAF(TX="15")
Func.HypoxiaSubclustcoloredCAF(TX="16")
Func.HypoxiaSubclustcoloredCAF(TX="19")


##  3. Hypoxic / Normoxic area (main figure) ####
Func.GenerateGreenAndSubclustcolored.main <- function(TX){
  XenView.Sub.FullLayer = readRDS(
    file=paste0(DirInteg,"XeniumView_Sub/[FigureRDS][SubClust_XeniumView_TX5K_",TX,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".rds"))
  XenView.Sub.DAPI_and_Green = XenView.Sub.FullLayer
  XenView.Sub.DAPI_and_Green@layers$geom_polygon...2 = NULL
  saveRDS(
    XenView.Sub.DAPI_and_Green,
    file=paste0(DirInteg,"XeniumView_Sub/ForPaper/[FigureRDS][SubClust_XeniumView_DAPI_and_Green_TX5K_",TX,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".rds"))
  XenView.Sub.DAPI_and_Subclust = XenView.Sub.FullLayer
  XenView.Sub.DAPI_and_Subclust@layers$geom_polygon = NULL
  saveRDS(
    XenView.Sub.DAPI_and_Subclust,
    file=paste0(DirInteg,"XeniumView_Sub/ForPaper/[FigureRDS][SubClust_XeniumView_DAPI_and_Subclust_TX5K_",TX,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".rds"))
}
Func.GenerateGreenAndSubclustcolored.main(TX = "16")

Func.HypoxiaAndCAF.Main <- function(TX){
  RDS.DAPI_and_Green <- readRDS(
    file=paste0(DirInteg,"XeniumView_Sub/ForPaper/[FigureRDS][SubClust_XeniumView_DAPI_and_Green_TX5K_",TX,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".rds"))
  Length.OldX = diff(range(RDS.DAPI_and_Green@data$vertex_x))
  Length.OldY = diff(range(RDS.DAPI_and_Green@data$vertex_y))
  AspectRatio = max(Length.OldX, Length.OldY) / min(Length.OldX, Length.OldY)
  if(Length.OldX > Length.OldY){
    Coord.Hypo.Final = rev(Coordinate.Hypo[[TX]])   # rotated
    Coord.Normo.Final = rev(Coordinate.Normo[[TX]])   # rotated
    CoordY.Bar.Mag.Normo = Coord.Normo.Final[2]+SideLength*0.050
    CoordY.Bar.Mag.Hypo = Coord.Hypo.Final[2]+SideLength*0.050
    Coords.BarEnd.NotMag = c("x"=min(RDS.DAPI_and_Green@data$vertex_y) + Length.OldY*0.95,
                             "y"=min(RDS.DAPI_and_Green@data$vertex_x) + Length.OldX*0.05)
  } else {
    Coord.Normo.Final = Coordinate.Normo[[TX]]   # not rotated
    Coord.Hypo.Final = Coordinate.Hypo[[TX]]   # not rotated
    CoordY.Bar.Mag.Normo = Coord.Normo.Final[2]+SideLength*0.950
    CoordY.Bar.Mag.Hypo = Coord.Hypo.Final[2]+SideLength*0.950
    Coords.BarEnd.NotMag = c("x"=min(RDS.DAPI_and_Green@data$vertex_x) + Length.OldX*0.95,
                             "y"=min(RDS.DAPI_and_Green@data$vertex_y) + Length.OldY*0.95)
  }
  Common.Cartesian.Hypo <- coord_cartesian(
    xlim=c(Coord.Hypo.Final[1], 
           Coord.Hypo.Final[1] + SideLength),
    ylim=c(Coord.Hypo.Final[2], 
           Coord.Hypo.Final[2] + SideLength),
    ratio=1)
  Common.Cartesian.Normo <- coord_cartesian(
    xlim=c(Coord.Normo.Final[1], 
           Coord.Normo.Final[1] + SideLength),
    ylim=c(Coord.Normo.Final[2], 
           Coord.Normo.Final[2] + SideLength),
    ratio=1)
  ScaleBar.White <- annotate(
    "segment", color="white", linewidth=6.0,
    x=Coord.Hypo.Final[1]+SideLength*0.698, xend=Coord.Hypo.Final[1]+SideLength*0.952, 
    y=CoordY.Bar.Mag.Hypo, yend=CoordY.Bar.Mag.Hypo)
  ScaleBar.Black <- annotate(
    "segment", color="black", linewidth=4.0,
    x=Coord.Hypo.Final[1]+SideLength*0.700, xend=Coord.Hypo.Final[1]+SideLength*0.950, 
    y=CoordY.Bar.Mag.Hypo, yend=CoordY.Bar.Mag.Hypo)
  CommonTheme = 
    theme(plot.title = element_blank(),
          plot.subtitle = element_blank(),
          plot.caption = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          plot.margin = margin(t=0, r=0, b=0, l=0, unit="cm"),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill="black", color=NA),
          panel.border = element_rect(fill = "transparent", color = NA))
  # 1. DAPI and Green
  ggsave(plot = 
           RDS.DAPI_and_Green +
           #ScaleBar.White + ScaleBar.Black + 
           Common.Cartesian.Normo +
           CommonTheme,
         width=10, height=10, dpi=500, bg="transparent", 
         file=paste0(DirInteg,"XeniumView_Sub/ForPaper/HypoxicNormoxic(main)_TX",TX,
                     "_DAPI_and_Green_HighMag_Normoxic_x",paste(Coordinate.Normo[[TX]], collapse="y"),".png"))
  ggsave(plot = 
           RDS.DAPI_and_Green + 
           #ScaleBar.White + ScaleBar.Black +
           Common.Cartesian.Hypo +
           CommonTheme,
         width=10, height=10, dpi=500, bg="transparent", 
         file=paste0(DirInteg,"XeniumView_Sub/ForPaper/HypoxicNormoxic(main)_TX",TX,
                     "_DAPI_and_Green_HighMag_Hypoxic_x",paste(Coordinate.Hypo[[TX]], collapse="y"),".png"))
  # 2. DAPI and CAF-8
  RDS.DAPI_and_Subclust <- readRDS(
    file=paste0(DirInteg,"XeniumView_Sub/ForPaper/[FigureRDS][SubClust_XeniumView_DAPI_and_Subclust_TX5K_",TX,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".rds"))
  ggsave(plot = 
           RDS.DAPI_and_Subclust +
           scale_fill_manual(
             values = setNames(c(NA, NA, NA, NA, NA, NA, NA, NA, "#b80058"),
                               paste0("CAF-", 0:8))) +
           #ScaleBar.White + ScaleBar.Black +
           Common.Cartesian.Normo +
           CommonTheme,
         width=10, height=10, dpi=500, bg="transparent", 
         file=paste0(DirInteg,"XeniumView_Sub/ForPaper/HypoxicNormoxic(main)_TX",TX,
                     "_DAPI_and_CAF8_HighMag_Normoxic_x",paste(Coordinate.Normo[[TX]], collapse="y"),".png"))
  ggsave(plot = 
           RDS.DAPI_and_Subclust +
           scale_fill_manual(
             values = setNames(c(NA, NA, NA, NA, NA, NA, NA, NA, "#b80058"),
                               paste0("CAF-", 0:8))) +
           #ScaleBar.White + ScaleBar.Black +
           Common.Cartesian.Hypo +
           CommonTheme,
         width=10, height=10, dpi=500, bg="transparent", 
         file=paste0(DirInteg,"XeniumView_Sub/ForPaper/HypoxicNormoxic(main)_TX",TX,
                     "_DAPI_and_CAF8_HighMag_Hypoxic_x",paste(Coordinate.Hypo[[TX]], collapse="y"),".png"))
  # 3. DAPI and Merge
  ggsave(plot = 
           RDS.DAPI_and_Subclust +
           scale_fill_manual(
             values = setNames(c("green", "green", "green", "green", "green", 
                                 "green", "green", "green", "#b80058"),
                               paste0("CAF-", 0:8))) +
           #ScaleBar.White + ScaleBar.Black +
           Common.Cartesian.Normo +
           CommonTheme,
         width=10, height=10, dpi=500, bg="transparent", 
         file=paste0(DirInteg,"XeniumView_Sub/ForPaper/HypoxicNormoxic(main)_TX",TX,
                     "_DAPI_and_Merged_HighMag_Normoxic_x",paste(Coordinate.Normo[[TX]], collapse="y"),".png"))
  ggsave(plot = 
           RDS.DAPI_and_Subclust +
           scale_fill_manual(
             values = setNames(c("green", "green", "green", "green", "green", 
                                 "green", "green", "green", "#b80058"),
                               paste0("CAF-", 0:8))) +
           #ScaleBar.White + ScaleBar.Black +
           Common.Cartesian.Hypo +
           CommonTheme,
         width=10, height=10, dpi=500, bg="transparent", 
         file=paste0(DirInteg,"XeniumView_Sub/ForPaper/HypoxicNormoxic(main)_TX",TX,
                     "_DAPI_and_Merged_HighMag_Hyopoxic_x",paste(Coordinate.Hypo[[TX]], collapse="y"),".png"))
}
Func.HypoxiaAndCAF.Main(TX = "16")

##  4. Niche                                 ####
Dir.XenViewNiche = paste0(DirInteg,"/NicheAnalysisData/XeniumView_Sub_Final/Final_",QCInfo.FileName,
                          "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,
                          "/[",Unclassified,"Unclassified_Labeled",Labels,"_Neighbor",n.neighbors,
                          "_npc",pc_use,"_k",k.niche,"]")
Xen.Niche.Point.01 = readRDS(file=paste0(Dir.XenViewNiche,"/[PointView_RDS_Final]_TX5K_01.rds"))
Xen.Niche.Point.02 = readRDS(file=paste0(Dir.XenViewNiche,"/[PointView_RDS_Final]_TX5K_02.rds")) + theme(legend.position = "none")
Xen.Niche.Point.19 = readRDS(file=paste0(Dir.XenViewNiche,"/[PointView_RDS_Final]_TX5K_19.rds")) + theme(legend.position = "none")
Xen.Niche.Point.11 = readRDS(file=paste0(Dir.XenViewNiche,"/[PointView_RDS_Final]_TX5K_11.rds")) + theme(legend.position = "none")
Xen.Niche.Point.16 = readRDS(file=paste0(Dir.XenViewNiche,"/[PointView_RDS_Final]_TX5K_16.rds")) + theme(legend.position = "none")
Xen.Niche.Point.15 = readRDS(file=paste0(Dir.XenViewNiche,"/[PointView_RDS_Final]_TX5K_15.rds")) + theme(legend.position = "none")

Func.Niche.Point.ForFig <- function(Plot, TXnum){
  pt.size <- 1
  xlim <- Plot@scales$scales[[1]]$limits
  ylim <- (Plot@scales$scales[[2]]$limits)*(-1)
  range.x <- max(xlim) - min(xlim)
  range.y <- max(ylim) - min(ylim)
  ScaleBarRt.x <- quantile(xlim, 0.93, names = F)
  ScaleBarRt.y <- quantile(ylim, 0.97, names = F)
  Plot@layers$geom_point$aes_params$size <- pt.size
  Plot@layers$geom_point$aes_params$size <- pt.size
  Plot@layers$geom_point...2$aes_params$stroke <- 2
  ggsave(
    Plot + 
      annotate(geom = "segment",
               x = ScaleBarRt.x-1000, 
               xend = ScaleBarRt.x, 
               y = ScaleBarRt.y, 
               yend = ScaleBarRt.y,
               color = "black",
               linewidth = 2) + 
      theme(legend.position = "none",
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank()),
    width = range.x/800, height = range.y/800,
    file = paste0("[NichePlusCAF8]_TX",TXnum,".png"), dpi=500)
}
Func.Niche.Point.ForFig(Plot = Xen.Niche.Point.01, TXnum = "01")
Func.Niche.Point.ForFig(Plot = Xen.Niche.Point.02, TXnum = "02")
Func.Niche.Point.ForFig(Plot = Xen.Niche.Point.19, TXnum = "19")
Func.Niche.Point.ForFig(Plot = Xen.Niche.Point.11, TXnum = "11")
Func.Niche.Point.ForFig(Plot = Xen.Niche.Point.16, TXnum = "16")
Func.Niche.Point.ForFig(Plot = Xen.Niche.Point.15, TXnum = "15")

##  5. Hypoxic / Normoxic point              ####

DF.MetaData <- readRDS(
  file=paste0(DirInteg,"[MetaDataList][InitialClust][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
              "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".rds"))$MetaData
SeuObj.CAF_0 = readRDS(
  paste0(DirInteg,"Objects/[SeuObj_ExtractedCAF][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
         QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2,".rds") )
DF.Coord.AllCase <- 
  data.frame(
    cell_id = colnames(SeuObj.CAF_0),
    sample = str_replace(SeuObj.CAF_0$orig.ident, pattern = "TX5K_", replacement = "TX"),
    CoordX = SeuObj.CAF_0$X,
    CoordY = SeuObj.CAF_0$Y,
    Subclusters = SeuObj.CAF_0$seurat_clusters)
NumOfSubclusts <- DF.Coord.AllCase$subclusters %>% unique() %>% length()

Func.Point.GenerateGreenAndSubclustcolored <- function(TX){
  # Individual coordinate data
  DF.Coord.Indiv_0 <- 
    dplyr::filter(DF.Coord.AllCase, sample == paste0("TX",TX))
  RangeX <- c(min(DF.Coord.Indiv_0$CoordX), max(DF.Coord.Indiv_0$CoordX))
  Length_X <- diff(RangeX)
  RangeY <- c(min(DF.Coord.Indiv_0$CoordY), max(DF.Coord.Indiv_0$CoordY))
  Length_Y <- diff(RangeY)
  WideOrLong <- case_when(
    Length_Y/Length_X >= 1 ~ c("WideOrLong"="Long", "ncol"=ceiling(NumOfSubclusts/2), "nrow"=2, "width"=ceiling(NumOfSubclusts/2)*4+1, "height"=2*8),
    Length_Y/Length_X < 1 ~ c("WideOrLong"="Wide", "ncol"=2, "nrow"=ceiling(NumOfSubclusts/2), "width"=2*8, "height"=ceiling(NumOfSubclusts/2)*4+1),
    TRUE ~ c("WideOrLong"="Square", "ncol"=ceiling(NumOfSubclusts/2), "nrow"=2, "width"=ceiling(NumOfSubclusts/2)*4+1, "height"=2*8 ) )
  AspectRatio = max(Length_X, Length_Y) / min(Length_X, Length_Y)
  
  # base plot
  if(WideOrLong[["WideOrLong"]] == "Long") {
    XenPoint.Sub_0 =
      ggplot(DF.Coord.Indiv_0, 
             aes(x = CoordX, y = CoordY, color = Subclusters)) +
      #geom_point() +
      scale_x_continuous(expand=expansion(mult=c(0, 0)),
                         breaks=seq(from=ceiling(RangeX[1]/1000) * 1000,
                                    to=floor(RangeX[2]/1000) * 1000,
                                    by=1000),
                         limits=RangeX) +
      scale_y_reverse(expand=expansion(mult=c(0, 0)),
                      breaks=seq(from=ceiling(RangeY[1]/1000) * 1000,
                                 to=floor(RangeY[2]/1000) * 1000,
                                 by=1000),
                      limits=RangeY)
  } else {
    XenPoint.Sub_0 =
      ggplot(DF.Coord.Indiv_0, 
             aes(x=CoordY, y=CoordX, color = Subclusters)) +
      #geom_point() +
      scale_x_continuous(expand=expansion(mult=c(0, 0)),
                         breaks=seq(from=ceiling(RangeY[1]/1000) * 1000,
                                    to=floor(RangeY[2]/1000) * 1000,
                                    by=1000),
                         limits=RangeY) +
      scale_y_continuous(expand=expansion(mult=c(0, 0)),
                         breaks=seq(from=ceiling(RangeX[1]/1000) * 1000,
                                    to=floor(RangeX[2]/1000) * 1000,
                                    by=1000),
                         limits=RangeX)
  }
  XenPoint.Sub_1 <- 
    XenPoint.Sub_0 +
    labs(title=paste0("Xenium point view, TX5K_",TX),
         x=NULL, y=NULL) +
    coord_fixed(ratio=1) + 
    theme(text = element_text(color="gray50", face="bold"),
          legend.position = "none",
          axis.line = element_blank(),
          plot.background = element_rect(fill = "transparent", colour = NA),
          panel.background = element_rect(fill = "black", colour = NA),
          panel.border = element_rect(fill = "transparent", colour = "black"),
          panel.grid = element_blank())
  
  fs::dir_create(path=paste0(DirInteg,"XeniumView_Point_Sub/ForPaper"))
  
  # Colored as green
  XenPoint.Sub_green <- 
    XenPoint.Sub_1 +
    scale_color_manual(values = rep("green",times=NumOfSubclusts))
  saveRDS(
    XenPoint.Sub_green,
    file=paste0(DirInteg,"XeniumView_Point_Sub/ForPaper/[FigureRDS][SubClust_XeniumView_GreenCAF_TX5K_",TX,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".rds"))

  # Subcluster color code
  XenPoint.Sub_multi <- 
    XenPoint.Sub_1 +
    scale_color_manual(values = cols2)
  saveRDS(
    XenPoint.Sub_multi,
    file=paste0(DirInteg,"XeniumView_Point_Sub/ForPaper/[FigureRDS][SubClust_XeniumView_SubclustcoloredCAF_TX5K_",TX,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".rds"))
}
Func.Point.GenerateGreenAndSubclustcolored(TX="01")
Func.Point.GenerateGreenAndSubclustcolored(TX="02")
Func.Point.GenerateGreenAndSubclustcolored(TX="11")
Func.Point.GenerateGreenAndSubclustcolored(TX="15")
Func.Point.GenerateGreenAndSubclustcolored(TX="16")
Func.Point.GenerateGreenAndSubclustcolored(TX="19")

Func.WeakAndStrongMag <- function(TX,Coloring){
  RDS <- readRDS(
    file=paste0(DirInteg,"XeniumView_Point_Sub/ForPaper/[FigureRDS][SubClust_XeniumView_",
                Coloring,"CAF_TX5K_",TX,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".rds"))
  Length.OldX = diff(range(RDS@data$CoordX))
  Length.OldY = diff(range(RDS@data$CoordY))
  AspectRatio = max(Length.OldX, Length.OldY) / min(Length.OldX, Length.OldY)
  if(Length.OldX > Length.OldY){
    Coord.Hypo.Final = rev(Coordinate.Hypo[[TX]])   # rotated
    Coord.Normo.Final = rev(Coordinate.Normo[[TX]])   # rotated
    CoordY.Bar.Mag.Normo = Coord.Normo.Final[2]+SideLength*0.050
    CoordY.Bar.Mag.Hypo = Coord.Hypo.Final[2]+SideLength*0.050
    Coords.BarEnd.NotMag = c("x"=min(RDS@data$CoordY) + Length.OldY*0.95,
                             "y"=min(RDS@data$CoordX) + Length.OldX*0.05)
  } else {
    Coord.Normo.Final = Coordinate.Normo[[TX]]   # not rotated
    Coord.Hypo.Final = Coordinate.Hypo[[TX]]   # not rotated
    CoordY.Bar.Mag.Normo = Coord.Normo.Final[2]+SideLength*0.950
    CoordY.Bar.Mag.Hypo = Coord.Hypo.Final[2]+SideLength*0.950
    Coords.BarEnd.NotMag = c("x"=min(RDS@data$CoordX) + Length.OldX*0.95,
                             "y"=min(RDS@data$CoordY) + Length.OldY*0.95)
  }
  CommonTheme = 
    theme(plot.title = element_blank(),
          plot.subtitle = element_blank(),
          plot.caption = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          plot.margin = margin(t=0, r=0, b=0, l=0, unit="cm"),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill="black", color=NA),
          panel.border = element_rect(fill = "transparent", color = NA))
  # 1. Weakly magnified Green CAF
  Plot.greenCAF = 
    RDS + 
    geom_point(size = 2) +
    annotate("rect", 
             xmin=Coord.Hypo.Final[1], xmax=(Coord.Hypo.Final[1])+SideLength,
             ymin=Coord.Hypo.Final[2], ymax=(Coord.Hypo.Final[2])+SideLength, 
             color="black", fill=NA, linewidth=2.8) +
    annotate("rect", 
             xmin=Coord.Hypo.Final[1], xmax=(Coord.Hypo.Final[1])+SideLength,
             ymin=Coord.Hypo.Final[2], ymax=(Coord.Hypo.Final[2])+SideLength, 
             color="skyblue", fill=NA, linewidth=2.1) +
    annotate("rect", 
             xmin=Coord.Normo.Final[1], xmax=(Coord.Normo.Final[1])+SideLength,
             ymin=Coord.Normo.Final[2], ymax=(Coord.Normo.Final[2])+SideLength, 
             color="black", fill=NA, linewidth=2.8) +
    annotate("rect", 
             xmin=Coord.Normo.Final[1], xmax=(Coord.Normo.Final[1])+SideLength,
             ymin=Coord.Normo.Final[2], ymax=(Coord.Normo.Final[2])+SideLength, 
             color="#f03b20", fill=NA, linewidth=2.1) +
    annotate("segment", color="white", linewidth=7.0,
             x=Coords.BarEnd.NotMag[["x"]]-1020, xend=Coords.BarEnd.NotMag[["x"]]+20, 
             y=Coords.BarEnd.NotMag[["y"]], yend=Coords.BarEnd.NotMag[["y"]]) +
    annotate("segment", color="black", linewidth=4.0,
             x=Coords.BarEnd.NotMag[["x"]]-1000, xend=Coords.BarEnd.NotMag[["x"]], 
             y=Coords.BarEnd.NotMag[["y"]], yend=Coords.BarEnd.NotMag[["y"]]) +
    CommonTheme
  ggsave(Plot.greenCAF,
         file=paste0(DirInteg,"XeniumView_Point_Sub/ForPaper/HypoxicNormoxic_TX",TX,"_",Coloring,"CAF_WeakMag.png"), 
         width=10, height=10*AspectRatio, 
         dpi=500, bg="transparent", limitsize=FALSE)
  
  # 2. Strongly magnified Hypo/Normo area
  ggsave(plot = 
           RDS + 
           geom_point(size = 10) +
           annotate("segment", color="white", linewidth=6.0,
                    x=Coord.Normo.Final[1]+SideLength*0.698, xend=Coord.Normo.Final[1]+SideLength*0.952, 
                    y=CoordY.Bar.Mag.Normo, yend=CoordY.Bar.Mag.Normo) +
           annotate("segment", color="black", linewidth=4.0,
                    x=Coord.Normo.Final[1]+SideLength*0.700, xend=Coord.Normo.Final[1]+SideLength*0.950, 
                    y=CoordY.Bar.Mag.Normo, yend=CoordY.Bar.Mag.Normo) +
           coord_cartesian(
             xlim=c(Coord.Normo.Final[1], Coord.Normo.Final[1]+SideLength),
             ylim=c(Coord.Normo.Final[2], Coord.Normo.Final[2]+SideLength),
             ratio=1) +
           CommonTheme,
         width=10, height=10, dpi=500, bg="transparent", 
         file=paste0(DirInteg,"XeniumView_Point_Sub/ForPaper/HypoxicNormoxic_TX",TX,"_",Coloring,"CAF_HighMag_Normoxic_x",paste(Coordinate.Normo[[TX]], collapse="y"),".png"))
  ggsave(plot = 
           RDS + 
           geom_point(size = 10) +
           annotate("segment", color="white", linewidth=6.0,
                    x=Coord.Hypo.Final[1]+SideLength*0.698, xend=Coord.Hypo.Final[1]+SideLength*0.952, 
                    y=CoordY.Bar.Mag.Hypo, yend=CoordY.Bar.Mag.Hypo) +
           annotate("segment", color="black", linewidth=4.0,
                    x=Coord.Hypo.Final[1]+SideLength*0.700, xend=Coord.Hypo.Final[1]+SideLength*0.950, 
                    y=CoordY.Bar.Mag.Hypo, yend=CoordY.Bar.Mag.Hypo) +
           coord_cartesian(
             xlim=c(Coord.Hypo.Final[1], Coord.Hypo.Final[1]+SideLength),
             ylim=c(Coord.Hypo.Final[2], Coord.Hypo.Final[2]+SideLength),
             ratio=1) +
           CommonTheme,
         width=10, height=10, dpi=500, bg="transparent", 
         file=paste0(DirInteg,"XeniumView_Point_Sub/ForPaper/HypoxicNormoxic_TX",TX,"_",Coloring,"CAF_HighMag_Hypoxic_x",paste(Coordinate.Hypo[[TX]], collapse="y"),".png"))
}
Func.WeakAndStrongMag(TX="01", "Green")
Func.WeakAndStrongMag(TX="01", "Subclustcolored")
Func.WeakAndStrongMag(TX="02", "Green")
Func.WeakAndStrongMag(TX="02", "Subclustcolored")
Func.WeakAndStrongMag(TX="11", "Green")
Func.WeakAndStrongMag(TX="11", "Subclustcolored")
Func.WeakAndStrongMag(TX="15", "Green")
Func.WeakAndStrongMag(TX="15", "Subclustcolored")
Func.WeakAndStrongMag(TX="16", "Green")
Func.WeakAndStrongMag(TX="16", "Subclustcolored")
Func.WeakAndStrongMag(TX="19", "Green")
Func.WeakAndStrongMag(TX="19", "Subclustcolored")











#################################################
##
##            Xenium mapping (Old)
#####
  TX="01"
  TX="02"
  TX="11"
  TX="15"
  TX="16"
  TX="19"
  #XenView.Ini.OnlyCAF = readRDS(
  #  file=paste0(DirInteg,"/XeniumView_Initial/",QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,
  #              "/[RDS][IniClustOnlyCAF_XeniumView_TX5K_",TX,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
  #              "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".rds"))
  XenView.Sub.FullLayer = readRDS(
    file=paste0(DirInteg,"XeniumView_Sub/[FigureRDS][SubClust_XeniumView_TX5K_",TX,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".rds"))
  XenView.Sub.AllGreenCAF = XenView.Sub.FullLayer
    XenView.Sub.AllGreenCAF@layers$geom_polygon...3 = NULL
    XenView.Sub.AllGreenCAF@layers$geom_polygon...2 = NULL
  XenView.Subclust <- XenView.Sub.FullLayer
    XenView.Subclust@layers$geom_polygon...3 = NULL
    XenView.Subclust@layers$geom_polygon...1 = NULL
  XenView.Sub.OnlyProlifCAF = XenView.Sub.FullLayer
    XenView.Sub.FullLayer@layers[[1]]$aes_params$alpha =0.8
    XenView.Sub.FullLayer@layers[[2]]$aes_params$alpha =0.8
    XenView.Sub.OnlyProlifCAF@layers[[1]] = NULL
    NumOfSubclust = length(unique(XenView.Sub.FullLayer@data$group))
    Vec.Colors.CAFsubclust = setNames(rep(NA, times=9), # gray20
                                      paste0("CAF-", 0:8))
    Vec.Colors.CAFsubclust["CAF-8"]="#b80058"
    XenView.Sub.OnlyProlifCAF = 
      XenView.Sub.OnlyProlifCAF + scale_fill_manual(values=Vec.Colors.CAFsubclust)
  Dir.XenViewNiche = paste0(DirInteg,"/NicheAnalysisData/XeniumView_Sub_Final/Final_",QCInfo.FileName,
                            "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,
                            "/[",Unclassified,"Unclassified_Labeled",Labels,"_Neighbor",n.neighbors,
                            "_npc",pc_use,"_k",k.niche,"]")
  Xen.Niche.Polygon = readRDS(file=paste0(Dir.XenViewNiche,"/[PolygonView_RDS]_TX5K_",TX,".rds") )
  
  SideLength = 400
  Coordinate.Area1 = list(
    "01" = c(2500, 2900), "02" = c(1300, 9500), "11" = c(2100, 3350), 
    "16" = c(4600, 5200), "15" = c(1250, 10550), "19" = c(7770, 6670))
  Coordinate.Area2 <- list(
    "01" = c(6200, 3100), "02" = c(3800, 2500), "11" = c(2650, 3000),
    "16" = c(4800, 3400), "15" = c(5400, 7450), "19" = c(8000, 750))
  Coordinate.Hypo <- list(
    "01" = c(9400, 1200), "02" = c(1400, 9400), "11" = c(3600, 3600),
    "16" = c(5200, 3600), "15" = c(800, 10400), "19" = c(7908, 8000))
  Coordinate.Normo <- list(
    "01" = #c(7600, 4400), 
      c(3800, 2000),
    "02" = c(2200, 5000), "11" = c(2200, 3200),
    "16" = c(4600, 4200), "15" = c(4600, 9200), "19" = c(7000, 4600))
  Coordinate.Cancer <- list("16" = c(3750, 4800))
  Coordinate.Vasc <- list("16" = c(1100, 4150))
  Coordinate.Plasma <- list("16" = c(4600, 4100))
  
  

  Func.Sub.LongVersion = function(TX){
    RDS = readRDS(
      file=paste0(DirInteg,"XeniumView_Sub/[FigureRDS][SubClust_XeniumView_TX5K_",TX,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                  "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".rds"))
    Length.OldX = diff(range(RDS@data$vertex_x))
    Length.OldY = diff(range(RDS@data$vertex_y))
    AspectRatio = max(Length.OldX, Length.OldY)/min(Length.OldX, Length.OldY)
    if(Length.OldX > Length.OldY){
      Coord.Hypo.Final = rev(Coord.Hypo)   # rotated
      Coord.Normo.Final = rev(Coord.Normo)   # rotated
      CoordY.Bar.Mag.Normo = Coord.Normo.Final[2]+SideLength*0.050
      CoordY.Bar.Mag.Hypo = Coord.Hypo.Final[2]+SideLength*0.050
      Coords.BarEnd.NotMag = c("x"=quantile(RDS@data$vertex_y, 0.95, names=F),
                               "y"=quantile(RDS@data$vertex_x, 0.02, names=F))
    } else {
      Coord.Normo.Final = Coord.Normo   # not rotated
      Coord.Hypo.Final = Coord.Hypo   # not rotated
      CoordY.Bar.Mag.Normo = Coord.Normo.Final[2]+SideLength*0.950
      CoordY.Bar.Mag.Hypo = Coord.Hypo.Final[2]+SideLength*0.950
      Coords.BarEnd.NotMag = c("x"=quantile(RDS@data$vertex_x, 0.95, names=F),
                               "y"=quantile(RDS@data$vertex_y, 0.98, names=F))
    }
    CommonTheme = 
      theme(plot.title = element_blank(),
            plot.subtitle = element_blank(),
            plot.caption = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            plot.margin = margin(t=0, r=0, b=0, l=0, unit="cm"),
            plot.background = element_rect(fill = "transparent", color = NA),
            panel.background = element_rect(fill="black", color=NA),
            panel.border = element_rect(fill = "transparent", color = NA))
    # 1. Weakly magnified Green CAF
    RDS.greenCAF = RDS
    RDS.greenCAF@layers$geom_polygon...3 = NULL # Nuclei
    RDS.greenCAF@layers$geom_polygon...2 = NULL # all CAF filled by multicolor
    Plot.greenCAF = 
      RDS.greenCAF + 
      annotate("rect", 
               xmin=Coord.Hypo.Final[1], xmax=(Coord.Hypo.Final[1])+SideLength,
               ymin=Coord.Hypo.Final[2], ymax=(Coord.Hypo.Final[2])+SideLength, 
               color="black", fill=NA, linewidth=2.8) +
      annotate("rect", 
               xmin=Coord.Hypo.Final[1], xmax=(Coord.Hypo.Final[1])+SideLength,
               ymin=Coord.Hypo.Final[2], ymax=(Coord.Hypo.Final[2])+SideLength, 
               color="skyblue", fill=NA, linewidth=2.1) +
      annotate("rect", 
               xmin=Coord.Normo.Final[1], xmax=(Coord.Normo.Final[1])+SideLength,
               ymin=Coord.Normo.Final[2], ymax=(Coord.Normo.Final[2])+SideLength, 
               color="black", fill=NA, linewidth=2.8) +
      annotate("rect", 
               xmin=Coord.Normo.Final[1], xmax=(Coord.Normo.Final[1])+SideLength,
               ymin=Coord.Normo.Final[2], ymax=(Coord.Normo.Final[2])+SideLength, 
               color="#f03b20", fill=NA, linewidth=2.1) +
      annotate("segment", color="white", linewidth=7.0,
               x=Coords.BarEnd.NotMag[["x"]]-1020, xend=Coords.BarEnd.NotMag[["x"]]+20, 
               y=Coords.BarEnd.NotMag[["y"]], yend=Coords.BarEnd.NotMag[["y"]]) +
      annotate("segment", color="black", linewidth=4.0,
               x=Coords.BarEnd.NotMag[["x"]]-1000, xend=Coords.BarEnd.NotMag[["x"]], 
               y=Coords.BarEnd.NotMag[["y"]], yend=Coords.BarEnd.NotMag[["y"]]) +
      CommonTheme
    ggsave(Plot.greenCAF,
           file=paste0("SubClusts_WeakMag_Green_TX",TX,".png"), 
           width=10, height=10*AspectRatio, 
           dpi=500, bg="transparent", limitsize=FALSE)
    # 2. Strongly magnified Hypo/Normo area
    ggsave(plot = 
             RDS.greenCAF +
             annotate("segment", color="white", linewidth=6.0,
                      x=Coord.Normo.Final[1]+SideLength*0.698, xend=Coord.Normo.Final[1]+SideLength*0.952, 
                      y=CoordY.Bar.Mag.Normo, yend=CoordY.Bar.Mag.Normo) +
             annotate("segment", color="black", linewidth=4.0,
                      x=Coord.Normo.Final[1]+SideLength*0.700, xend=Coord.Normo.Final[1]+SideLength*0.950, 
                      y=CoordY.Bar.Mag.Normo, yend=CoordY.Bar.Mag.Normo) +
             coord_cartesian(
               xlim=c(Coord.Normo.Final[1], Coord.Normo.Final[1]+SideLength),
               ylim=c(Coord.Normo.Final[2], Coord.Normo.Final[2]+SideLength),
               ratio=1) +
             CommonTheme,
          width=10, height=10, dpi=500, bg="transparent", 
          file=paste0("/Volumes/PortableSSD/[4]myR/TX5K",TX,"_x",paste(Coord.Normo, collapse="y"),"_Normo_OnlyCAF.png"))
    ggsave(plot = 
             RDS.greenCAF +
             annotate("segment", color="white", linewidth=6.0,
                      x=Coord.Hypo.Final[1]+SideLength*0.698, xend=Coord.Hypo.Final[1]+SideLength*0.952, 
                      y=CoordY.Bar.Mag.Hypo, yend=CoordY.Bar.Mag.Hypo) +
             annotate("segment", color="black", linewidth=4.0,
                      x=Coord.Hypo.Final[1]+SideLength*0.700, xend=Coord.Hypo.Final[1]+SideLength*0.950, 
                      y=CoordY.Bar.Mag.Hypo, yend=CoordY.Bar.Mag.Hypo) +
             coord_cartesian(
               xlim=c(Coord.Hypo.Final[1], Coord.Hypo.Final[1]+SideLength),
               ylim=c(Coord.Hypo.Final[2], Coord.Hypo.Final[2]+SideLength),
               ratio=1) +
             CommonTheme,
           width=10, height=10, dpi=500, bg="transparent", 
           file=paste0("/Volumes/PortableSSD/[4]myR/TX5K",TX,"_x",paste(Coord.Hypo, collapse="y"),"_Hypo_OnlyCAF.png"))
    # 3. multi colored weakly magnified
    RDS.multicoloredCAF = RDS
    RDS.multicoloredCAF@layers$geom_polygon...3 = NULL # Nuclei
    RDS.multicoloredCAF@layers$geom_polygon =NULL # all CAF filled green
    Plot.multicoloredCAF = 
      RDS.multicoloredCAF + 
      scale_fill_manual(values=c(setNames(cols2, paste0("CAF-",0:8)))) +
      annotate("rect", 
               xmin=Coord.Hypo.Final[1], xmax=(Coord.Hypo.Final[1])+SideLength,
               ymin=Coord.Hypo.Final[2], ymax=(Coord.Hypo.Final[2])+SideLength, 
               color="black", fill=NA, linewidth=2.8) +
      annotate("rect", 
               xmin=Coord.Hypo.Final[1], xmax=(Coord.Hypo.Final[1])+SideLength,
               ymin=Coord.Hypo.Final[2], ymax=(Coord.Hypo.Final[2])+SideLength, 
               color="skyblue", fill=NA, linewidth=2.1) +
      annotate("rect", 
               xmin=Coord.Normo.Final[1], xmax=(Coord.Normo.Final[1])+SideLength,
               ymin=Coord.Normo.Final[2], ymax=(Coord.Normo.Final[2])+SideLength, 
               color="black", fill=NA, linewidth=2.8) +
      annotate("rect", 
               xmin=Coord.Normo.Final[1], xmax=(Coord.Normo.Final[1])+SideLength,
               ymin=Coord.Normo.Final[2], ymax=(Coord.Normo.Final[2])+SideLength, 
               color="#f03b20", fill=NA, linewidth=2.1) +
      annotate("segment", color="white", linewidth=7.0,
               x=Coords.BarEnd.NotMag[["x"]]-1020, xend=Coords.BarEnd.NotMag[["x"]]+20, 
               y=Coords.BarEnd.NotMag[["y"]], yend=Coords.BarEnd.NotMag[["y"]]) +
      annotate("segment", color="black", linewidth=4.0,
               x=Coords.BarEnd.NotMag[["x"]]-1000, xend=Coords.BarEnd.NotMag[["x"]], 
               y=Coords.BarEnd.NotMag[["y"]], yend=Coords.BarEnd.NotMag[["y"]]) +
      CommonTheme
    ggsave(Plot.multicoloredCAF,
           file=paste0("SubClusts_WeakMag_MultiColor_TX",TX,".png"), 
           width=10, height=10*AspectRatio, 
           dpi=500, bg="transparent", limitsize=FALSE)
    # 4. 
    ggsave(plot = 
             RDS.multicoloredCAF + 
             annotate("segment", color="white", linewidth=6.0,
                      x=Coord.Normo.Final[1]+SideLength*0.698, xend=Coord.Normo.Final[1]+SideLength*0.952, 
                      y=CoordY.Bar.Mag.Normo, yend=CoordY.Bar.Mag.Normo) +
             annotate("segment", color="black", linewidth=4.0,
                      x=Coord.Normo.Final[1]+SideLength*0.700, xend=Coord.Normo.Final[1]+SideLength*0.950, 
                      y=CoordY.Bar.Mag.Normo, yend=CoordY.Bar.Mag.Normo) +
             scale_fill_manual(values=c(setNames(cols2, paste0("CAF-",0:8)))) +
             coord_cartesian(
               xlim=c(Coord.Normo.Final[1], Coord.Normo.Final[1]+SideLength),
               ylim=c(Coord.Normo.Final[2], Coord.Normo.Final[2]+SideLength),
               ratio=1) +
             CommonTheme,
           width=10, height=10, dpi=500, bg="transparent", 
           file=paste0("/Volumes/PortableSSD/[4]myR/TX5K",TX,"_x",paste(Coord.Normo.Final, collapse="y"),"_Normo_CAFsubclust.png"))
    ggsave(plot = 
             RDS.multicoloredCAF + 
             annotate("segment", color="white", linewidth=6.0,
                      x=Coord.Hypo.Final[1]+SideLength*0.698, xend=Coord.Hypo.Final[1]+SideLength*0.952, 
                      y=CoordY.Bar.Mag.Hypo, yend=CoordY.Bar.Mag.Hypo) +
             annotate("segment", color="black", linewidth=4.0,
                      x=Coord.Hypo.Final[1]+SideLength*0.700, xend=Coord.Hypo.Final[1]+SideLength*0.950, 
                      y=CoordY.Bar.Mag.Hypo, yend=CoordY.Bar.Mag.Hypo) +
             scale_fill_manual(values=c(setNames(cols2, paste0("CAF-",0:8)))) +
             coord_cartesian(
               xlim=c(Coord.Hypo.Final[1], Coord.Hypo.Final[1]+SideLength),
               ylim=c(Coord.Hypo.Final[2], Coord.Hypo.Final[2]+SideLength),
               ratio=1) +
             CommonTheme,
           width=10, height=10, dpi=500, bg="transparent", 
           file=paste0("/Volumes/PortableSSD/[4]myR/TX5K",TX,"_x",paste(Coord.Hypo.Final, collapse="y"),"_Hypo_CAFsubclusts.png"))
  }
  Func.GreenCAF.FromSub = function(RDS){
    Plot.GreenCAF_0 = RDS
    LengthX = diff(range(Plot.GreenCAF_0@data$vertex_x))
    LengthY = diff(range(Plot.GreenCAF_0@data$vertex_y))
    WideOrLong = case_when( LengthY/LengthX > 3/2 ~ c("WideOrLong"="Long", "width"=10, "height"=10*LengthY/LengthX),
                            LengthX/LengthY < 2/3 ~ c("WideOrLong"="Wide", "width"=10*LengthX/LengthY, "height"=10),
                            TRUE ~ c("WideOrLong"="Square", "width"=10, "height"=10*LengthY/LengthX) )
    Plot.GreenCAF_0@layers$geom_polygon...3 = NULL # Nuclei
    Plot.GreenCAF_0@layers$geom_polygon...2 = NULL # all CAF filled by multicolor
    Plot.GreenCAF_1 = 
      RDS.greenCAF_0 + 
      labs(title = NULL, subtitle = NULL, caption = NULL) +
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            plot.margin = margin(t=0, r=0, b=0, l=0, unit="cm"))
    return(Plot.GreenCAF_1)
  }
  Func.Cartesian = function(RDS, Coord, FileName){
    Length.OldX = diff(RDS@scales$get_scales("x")$limits)
    Length.OldY = diff(RDS@scales$get_scales("y")$limits)
    if(Length.OldX < Length.OldY){ # Not rotated
      InSeg.X = Coord[1]+SideLength*0.70
      InSeg.Xend = Coord[1]+SideLength*0.95
      InSeg.Y =Coord[2]+SideLength*0.05
      InSeg.Yend = Coord[2]+SideLength*0.05
      OutSeg.X = Coord[1]+SideLength*0.698
      OutSeg.Xend = Coord[1]+SideLength*0.952
      OutSeg.Y =Coord[2]+SideLength*0.05
      OutSeg.Yend = Coord[2]+SideLength*0.05
    } else {                      # Rotated
      InSeg.X = Coord[1]+SideLength*0.05 
      InSeg.Xend = Coord[1]+SideLength*0.05
      InSeg.Y =Coord[2]+SideLength*0.70
      InSeg.Yend = Coord[2]+SideLength*0.95
      OutSeg.X = Coord[1]+SideLength*0.05 
      OutSeg.Xend = Coord[1]+SideLength*0.05
      OutSeg.Y =Coord[2]+SideLength*0.698
      OutSeg.Yend = Coord[2]+SideLength*0.952
    }
    plot = 
      RDS +
      annotate("segment", color="white", linewidth=6.0,
               x = OutSeg.X, 
               xend = OutSeg.Xend, 
               y = OutSeg.Y, 
               yend = OutSeg.Yend) +
      annotate("segment", color="black", linewidth=4.0,
               x = InSeg.X, 
               xend = InSeg.Xend, 
               y = InSeg.Y, 
               yend = InSeg.Yend) +
      coord_cartesian(
        xlim=c(Coord[1], Coord[1]+SideLength),
        ylim=c(Coord[2], Coord[2]+SideLength),
        ratio=1) +
      theme(plot.margin = unit(c(0,0,0,0), "mm"),
            plot.title = element_blank(),
            plot.subtitle = element_blank(),
            plot.caption = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            plot.background = element_rect(fill = "transparent", color = NA),
            panel.background = element_rect(fill="black", color=NA),
            panel.border = element_rect(fill = "transparent", color = NA))
    ggsave(plot,
      width=10, height=10, dpi=500, bg="transparent", 
      file=paste0("/Volumes/PortableSSD/[4]myR/TX5K",TX,"_x",paste(Coord, collapse="y"),"_",FileName,".png"))
  }
  Func.NichePoint.Rect = function(Coord1, Coord2){
    Dir.XenViewNiche = paste0(DirInteg,"NicheAnalysisData/XeniumView_Sub_Final/Final_",QCInfo.FileName,
                              "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,
                              "/[",Unclassified,"Unclassified_Labeled",Labels,"_Neighbor",n.neighbors,
                              "_npc",pc_use,"_k",k.niche,"]")
    Xen.Niche.Point = 
      readRDS(file=paste0(Dir.XenViewNiche,"/[PointView_RDS_Final]_TX5K_",TX,".rds"))
    Range.OldX = range(Xen.Niche.Point@data$CoordX)
    Range.OldY = range(Xen.Niche.Point@data$CoordY)
    Length.OldX = diff(Range.OldX)
    Length.OldY = diff(Range.OldY)
    AspectRatio = max(Length.OldX, Length.OldY)/min(Length.OldX, Length.OldY)
    if(Length.OldX>Length.OldY){
      Coord1.NewXmin = Coord1[2]
      Coord1.NewXmax = Coord1.NewXmin + SideLength
      Coord1.NewYmin = Range.OldX[2] - Coord1[1]
      Coord1.NewYmax = Coord1.NewYmin - SideLength
      Coord2.NewXmin = Coord2[2]
      Coord2.NewXmax = Coord2.NewXmin + SideLength
      Coord2.NewYmin = Range.OldX[2] - Coord2[1]
      Coord2.NewYmax = Coord2.NewYmin - SideLength
    } else {
      Coord1.NewXmin = Coord1[1]
      Coord1.NewXmax = Coord1.NewXmin + SideLength
      Coord1.NewYmin = Coord1[2]
      Coord1.NewYmax = Coord1.NewYmin + SideLength
      Coord2.NewXmin = Coord2[1]
      Coord2.NewXmax = Coord2.NewXmin + SideLength
      Coord2.NewYmin = Coord2[2]
      Coord2.NewYmax = Coord2.NewYmin + SideLength
    }
    #Xen.Niche.Point@layers[[2]]$aes_params$size=2
    #Xen.Niche.Point@layers[[2]]$aes_params$stroke=1.2
    Plot.Niche.Point.Rect = 
      Xen.Niche.Point +
      annotate("rect", 
               xmin=Coord1.NewXmin, xmax=Coord1.NewXmax,
               ymin=Coord1.NewYmin, ymax=Coord1.NewYmax, 
               color="black", fill=NA, linewidth=2) +
      annotate("rect", 
               xmin=Coord2.NewXmin, xmax=Coord2.NewXmax,
               ymin=Coord2.NewYmin, ymax=Coord2.NewYmax, 
               color="black", fill=NA, linewidth=2) +
      theme(legend.position = "none",
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            plot.margin = margin(t=0, r=0, b=0, l=0))
    ggsave(
      plot = Plot.Niche.Point.Rect,
      file=paste0("/Volumes/PortableSSD/[4]myR/TX5K",TX,"_NichePointMap_Final",
                  "_Area1x",paste(Coord1, collapse="y"),
                  "_Area2x",paste(Coord2, collapse="y"),
                  ".png"),
      width=5, height=5*AspectRatio, dpi=400)
    }
  Func.HypoxiaHeatmap.Rect = function(Coord.Hypo, Coord.Normo){
    Xen.HypoxiaHeatmap = readRDS(
      file=paste0(DirInteg,"/XeniumView_Sub/",QCInfo.FileName,
                  "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,
                  "/ssGSVA/[Figure][ExtractedCAF_XenView_GSVA_HypoxicAreaAndProlifCAF_TX5K_",TX,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                  "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".rds"))
    Xen.HypoxiaHeatmap@layers$geom_point$aes_params$size = 5
    Xen.HypoxiaHeatmap@layers$geom_point$aes_params$stroke = 3
    LengthX = diff(range(Xen.HypoxiaHeatmap@layers$geom_point$data$X))
    LengthY = diff(range(Xen.HypoxiaHeatmap@layers$geom_point$data$Y))
    if(LengthY > LengthX){ # long
      SaveWidth=9
      SaveHeight=9/LengthX*LengthY
    } else {  # Wide
      SaveWidth=9/LengthY*LengthX
      SaveHeight=9
    }
    ggsave(plot = 
      Xen.HypoxiaHeatmap +
      annotate("rect", 
               xmin=Coord.Hypo[1], xmax=(Coord.Hypo[1])+SideLength,
               ymin=Coord.Hypo[2], ymax=(Coord.Hypo[2])+SideLength, 
               color="black", fill=NA, linewidth=4) +
      annotate("rect", 
               xmin=Coord.Hypo[1], xmax=(Coord.Hypo[1])+SideLength,
               ymin=Coord.Hypo[2], ymax=(Coord.Hypo[2])+SideLength, 
               color="#253494", fill=NA, linewidth=3) +
      annotate("rect", 
               xmin=Coord.Normo[1], xmax=(Coord.Normo[1])+SideLength,
               ymin=Coord.Normo[2], ymax=(Coord.Normo[2])+SideLength, 
               color="black", fill=NA, linewidth=4) +
      annotate("rect", 
               xmin=Coord.Normo[1], xmax=(Coord.Normo[1])+SideLength,
               ymin=Coord.Normo[2], ymax=(Coord.Normo[2])+SideLength,  
               color="#f03b20", fill=NA, linewidth=3) +
      theme(plot.title = element_blank(),
            plot.subtitle = element_blank(),
            plot.caption = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            plot.margin = margin(t=0, r=0, b=0, l=0, unit="cm"),
            plot.background = element_rect(fill = "transparent", color = NA),
            panel.background = element_rect(fill="transparent", color=NA),
            panel.border = element_rect(fill = "transparent", color = "black")),
      file=paste0("/Volumes/PortableSSD/[4]myR/TX5K",TX,"_HeatMap_",
                  "_Area1x",paste(Coord.Hypo, collapse="y"),
                  "_Area2x",paste(Coord.Normo, collapse="y"),
                  ".png"),
      width=SaveWidth, height=SaveHeight, dpi=500, bg="transparent")
  }
  # Weak magnified view Subclusts
  Func.Sub.LongVersion(TX)
  # Relative Hypoxic / Normoxic
    #Func.Cartesian(XenView.Ini.OnlyCAF, Coordinate.Hypo[[TX]], "Hypo_OnlyCAF")
    #Func.Cartesian(XenView.Ini.OnlyCAF, Coordinate.Normo[[TX]], "Normo_OnlyCAF")
  Func.Cartesian(XenView.Sub.AllGreenCAF, Coordinate.Hypo[[TX]], "Hypo_OnlyCAF")
  Func.Cartesian(XenView.Sub.AllGreenCAF, Coordinate.Normo[[TX]], "Normo_OnlyCAF")
  Func.Cartesian(XenView.Subclust, Coordinate.Hypo[[TX]], "Hypo_Subclust")
  Func.Cartesian(XenView.Subclust, Coordinate.Normo[[TX]], "Normo_Subclust")
  Func.Cartesian(XenView.Sub.OnlyProlifCAF, Coordinate.Hypo[[TX]], "Hypo_OnlyProlifCAF")
  Func.Cartesian(XenView.Sub.OnlyProlifCAF, Coordinate.Normo[[TX]], "Normo_OnlyProlifCAF")
  Func.Cartesian(XenView.Sub.FullLayer + scale_fill_manual(values=Vec.Colors.CAFsubclust), Coordinate.Hypo[[TX]], "Hypo_Merged")
  Func.Cartesian(XenView.Sub.FullLayer + scale_fill_manual(values=Vec.Colors.CAFsubclust), Coordinate.Normo[[TX]], "Normo_Merged")
  Func.Cartesian(Xen.Niche.Polygon, Coordinate.Hypo[[TX]], "Hypo_OnlyCAF_NicheView")
  Func.Cartesian(Xen.Niche.Polygon, Coordinate.Normo[[TX]], "Normo_OnlyCAF_NicheView")
  for(i in 1:6){
    TX=TXnumInteg[i]
    Func.HypoxiaHeatmap.Rect(Coordinate.Hypo[[TX]], Coordinate.Normo[[TX]])
    }
  for(i in 1:6){
    TX=TXnumInteg[i]
    Func.NichePoint.Rect(Coordinate.Hypo[[TX]], Coordinate.Normo[[TX]])}
  # Niche cancer
  Func.Cartesian(Xen.Niche.Polygon, Coordinate.Cancer[[TX]], "CancerNiche_niche")
  Func.Cartesian(XenView.Ini.OnlyCAF, Coordinate.Cancer[[TX]], "CancerNiche_OnlyCAF")
  Func.Cartesian(XenView.Sub.OnlyProlifCAF, Coordinate.Cancer[[TX]], "CancerNiche_OnlyProlifCAF")
  Func.Cartesian(XenView.Sub.FullLayer + scale_fill_manual(values=Vec.Colors.CAFsubclust), Coordinate.Cancer[[TX]], "CancerNiche_Merged")
  # Niche Vascular
  Func.Cartesian(Xen.Niche.Polygon, Coordinate.Vasc[[TX]], "VascNiche_niche")
  Func.Cartesian(XenView.Ini.OnlyCAF, Coordinate.Vasc[[TX]], "VascNiche_OnlyCAF")
  Func.Cartesian(XenView.Sub.OnlyProlifCAF, Coordinate.Vasc[[TX]], "VascNiche_OnlyProlifCAF")
  Func.Cartesian(XenView.Sub.FullLayer + scale_fill_manual(values=Vec.Colors.CAFsubclust), Coordinate.Vasc[[TX]], "VascNiche_Merged")
  # Niche T+Myelo
  Func.Cartesian(Xen.Niche.Polygon, Coordinate.TM[[TX]], "TMNiche")
  Func.Cartesian(XenView.Sub.FullLayer, Coordinate.TM[[TX]], "TMNiche2")
  Func.Cartesian(XenView.Ini, Coordinate.TM[[TX]], "TMNiche3")
  # Niche not cancer
  Func.Cartesian(Xen.Niche.Polygon, Coordinate.Plasma[[TX]], "PlasmaNiche_niche")
  Func.Cartesian(XenView.Sub.FullLayer, Coordinate.Plasma[[TX]], "PlasmaNiche_caf")
  Func.Cartesian(XenView.Ini, Coordinate.Plasma[[TX]], "PlasmaNiche_celltype")

  
  
  #XenView.Ini.FullLayer = readRDS(
  #  file=paste0(DirInteg,"/XeniumView_Initial/",QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,
  #              "/[RDS][IniClust_XeniumView_TX5K_",TX,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
  #              "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".rds"))
  #XenView.Ini.FullLayer@layers[[2]]$aes_params$alpha = 0.1
  #XenView.Ini.FullLayer@layers[[2]]$aes_params$alpha = 0.2
  # Area
  Func.Ini.OnlyCP.Rect = function(TX){
    Coord.Area1 = Coordinate.Area1[[TX]]
    Coord.Area2 = Coordinate.Area2[[TX]]
    XenView.Ini.OnlyCP = readRDS(
      file=paste0(DirInteg,"/XeniumView_Initial/",QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,
                  "/[RDS][IniClustOnlyCP_XeniumView_TX5K_",TX,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                  "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".rds"))
    List.MetaData = readRDS(file=paste0(DirInteg,"[MetaDataList][InitialClust][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                                      "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".rds"))
    XenPlot.Ini.ColCodeForFinal = 
      c(List.MetaData[["OrderedClustnumToColor"]], "Excluded_byQC"=NA)
    XenPlot.Ini.ColCodeForFinal[c("22","25","28","35","36","Excluded_byQC")] = "#3A3A3A"
    LengthX = diff(range(XenView.Ini.OnlyCP@data$CoordX, na.rm=TRUE))
    LengthY = diff(range(XenView.Ini.OnlyCP@data$CoordY, na.rm=TRUE))
    if(LengthY > LengthX){
        Coords.Bar.X = c("start" = min(XenView.Ini.OnlyCP@data$CoordX, na.rm=TRUE) + LengthX*0.94,
                         "end" = (min(XenView.Ini.OnlyCP@data$CoordX, na.rm=TRUE) + LengthX*0.94) - 1000)
        Coords.Bar.Y = c("start" = min(XenView.Ini.OnlyCP@data$CoordY, na.rm=TRUE) + LengthY*0.96,
                         "end" = min(XenView.Ini.OnlyCP@data$CoordY, na.rm=TRUE) + LengthY*0.96)
        ImageHeight = 10*LengthY/LengthX
        ImageWidth = 10
    } else {
        Coords.Bar.X = c("start" = min(XenView.Ini.OnlyCP@data$CoordX, na.rm=TRUE) + LengthX*0.04,
                         "end" = min(XenView.Ini.OnlyCP@data$CoordX, na.rm=TRUE) + LengthX*0.04)
        Coords.Bar.Y = c("start" = (min(XenView.Ini.OnlyCP@data$CoordY, na.rm=TRUE) + LengthY*0.94) - 1000,
                         "end" = min(XenView.Ini.OnlyCP@data$CoordY, na.rm=TRUE) + LengthY*0.94)
        ImageHeight = 10
        ImageWidth = 10*LengthX/LengthY
  }
  Plot.Ini.WeakMag = 
    XenView.Ini.OnlyCP +
    annotate("rect", 
             xmin=Coord.Area1[1], xmax=Coord.Area1[1]+SideLength,
             ymin=Coord.Area1[2], ymax=Coord.Area1[2]+SideLength, 
             color="white", fill=NA, linewidth=3) +
    annotate("rect", 
             xmin=Coord.Area1[1], xmax=Coord.Area1[1]+SideLength,
             ymin=Coord.Area1[2], ymax=Coord.Area1[2]+SideLength, 
             color="black", fill=NA, linewidth=2) +
    annotate("rect", 
             xmin=Coord.Area2[1], xmax=Coord.Area2[1]+SideLength,
             ymin=Coord.Area2[2], ymax=Coord.Area2[2]+SideLength, 
             color="white", fill=NA, linewidth=3) +
    annotate("rect", 
             xmin=Coord.Area2[1], xmax=Coord.Area2[1]+SideLength,
             ymin=Coord.Area2[2], ymax=Coord.Area2[2]+SideLength, 
             color="black", fill=NA, linewidth=2) +
    annotate("segment", 
             x=Coords.Bar.X[["start"]], xend=Coords.Bar.X[["end"]],
             y=Coords.Bar.Y[["start"]], yend=Coords.Bar.Y[["end"]],
             color="black", linewidth=5) +
    annotate("segment", 
             x=Coords.Bar.X[["start"]], xend=Coords.Bar.X[["end"]],
             y=Coords.Bar.Y[["start"]], yend=Coords.Bar.Y[["end"]],
             color="white", linewidth=2.5) +
    scale_fill_manual(values=XenPlot.Ini.ColCodeForFinal) +
    theme(
      plot.margin = unit(c(0,0,0,0), "mm"),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      plot.background = element_rect(fill="transparent", color=NA),
      panel.background = element_rect(fill="black", color=NA),
      panel.border = element_rect(fill="transparent", color="black"))
    ggsave(plot = Plot.Ini.WeakMag,
           width=ImageWidth, height=ImageHeight, dpi=500, bg="transparent",
           file=paste0("/Volumes/PortableSSD/[4]myR/TX5K",TX,"_NoZoomUp_",
                       "_Area1x",paste(Coord.Area1, collapse="y"),
                       "_Area2x",paste(Coord.Area2, collapse="y"),
                       "2.png"))
  } 
  for(i in 1:6){
    TX = TXnumInteg[i]
    Func.Ini.OnlyCP.Rect(TX)
  }
  for(i in 1:6){
    TX = TXnumInteg[i]
    XenView.Ini.OnlyCP = readRDS(
      file=paste0(DirInteg,"/XeniumView_Initial/",QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,
                  "/[RDS][IniClustOnlyCP_XeniumView_TX5K_",TX,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                  "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".rds"))
    Func.Cartesian(XenView.Ini.OnlyCP, Coordinate.Area1[[TX]], "Area1")
    Func.Cartesian(XenView.Ini.OnlyCP, Coordinate.Area2[[TX]], "Area2")
  }

  
  
  
  
  
  

