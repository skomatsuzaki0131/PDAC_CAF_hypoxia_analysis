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
Vec.colorcode = c("Acinar"="#ff00ff", "Ductal"="#ddbcff", "Islet"="#ffff00",
                  "CAF"="#00ff00","Mural"="#990000","Endothelial"="#ff0000",
                  "Lymphoid"="#00bfff", "Lymph_T"="#00bfff","Lymph_B"="#00bfff","Plasma"="#00bfff",
                  "Myeloid"="#ff8800", "Mast"="#ff8800", "Nerve"="#0000ff", 
                  "Proliferating_Immune"="#999999", "Unclassified"="#999999")
Vec.colorcode.Final = 
  c("Acinar"="#ff00ff", "Ductal_like_acinar"="#ddbcff", "Ductal"="#ddbcff",
   "Normal_ductal"="#FFFFCC","PanIN1"="#FFFF00","PanIN2"="#FD8D3C","PDAC"="#ff0000",
   "Islet"="#BF812D",
   "CAF"="#00ff00","Mural"="#990000","Endothelial"="#990000",
   "Lymphoid"="#00bfff", "Lymph_T"="#00bfff","Lymph_B"="#00bfff","Plasma"="#00bfff",
   "Myeloid"="#8300FF", "Mast"="#8300FF", "Nerve"="#0000ff", 
   "Proliferating_Immune"="#999999", "Unclassified"="#999999")
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
DirInteg = paste0('/Volumes/Extreme SSD/Analysis/Data/IntegAnalysis/',
                  "[",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                  "_nFeatRNA:",paste(nFeatRNA,collapse="~"),"_nCountRNA:",paste(nCountRNA,collapse="~"),"/")


###############################################################
###    Integration
### 0. NormalizeData() ...already have done ####
Func.ReadSeuObjForInteg = function(TX, Mag, nFeatRNA, nCountRNA, Dim1, Res1){
  TXdl=if(TX %in% c("27","28")){ "27_28" }else{ TX }
  Directory = paste0("/Volumes/Extreme SSD/Analysis/Data/TX5K_",TXdl,"/")
  Dir0 = "/Volumes/Extreme SSD/Analysis/Data/"
  SeuratObj = readRDS(paste0(Directory,"Objects/[SeuObj][TX5K_",TX,"]_Countable_Mag",Mag,
                             "_nFeat",nFeatRNA[1],"-",nFeatRNA[2],
                             "_nCount",nCountRNA[1],"-",nCountRNA[2],".rds") )
  colnames(SeuratObj) = paste(str_remove(SeuratObj$orig.ident, pattern="_"), 
                              colnames(SeuratObj), sep="_")
  return(SeuratObj)}
ObjList = list(Func.ReadSeuObjForInteg(TX="01", Mag=1, nFeatRNA=c(100, 900), nCountRNA=c(100, 1800)), # Tis 0
               Func.ReadSeuObjForInteg(TX="02", Mag=1, nFeatRNA=c(100, 900), nCountRNA=c(100, 1800)), # T1b IA
               Func.ReadSeuObjForInteg(TX="19", Mag=1, nFeatRNA=c(100, 900), nCountRNA=c(100, 1800)), # T1c IA
               Func.ReadSeuObjForInteg(TX="11", Mag=1, nFeatRNA=c(100, 900), nCountRNA=c(100, 1800)), # T2b IA
               Func.ReadSeuObjForInteg(TX="16", Mag=1, nFeatRNA=c(100, 900), nCountRNA=c(100, 1800)), # T1b IA
               Func.ReadSeuObjForInteg(TX="15", Mag=1, nFeatRNA=c(100, 900), nCountRNA=c(100, 1800)), # T3 IIA
               #Func.ReadSeuObjForInteg(TX="04", Mag=1, nFeatRNA=c(100, 900), nCountRNA=c(100, 1800)), # T1b IA
               Func.ReadSeuObjForInteg(TX="14", Mag=1, nFeatRNA=c(100, 900), nCountRNA=c(100, 1800)) # T1 IA
               )
CombinedObj_1 = Reduce(function(m,n) merge(m, y=n, add.cell.ids=NULL), ObjList)
rm(ObjList)
Fuc.MinorQC = function(nFeatRNA, nCountRNA){
  CombinedObj_1_QC = subset(CombinedObj_1, 
                          subset=nFeature_RNA>nFeatRNA[1] & nFeature_RNA<nFeatRNA[2] &
                                 nCount_RNA>nCountRNA[1] & nCount_RNA<nCountRNA[2])
  n_before = ncol(CombinedObj_1) %>% format(big.mark=",")
  Breaks.nFeat = seq(0, max(CombinedObj_1$nFeature_RNA), by=10)
  Count.nFeat = table(cut(CombinedObj_1$nFeature_RNA, breaks=Breaks.nFeat, include.lowest=T, right=T) )
  Vec.Frequent.nFeat = c(Breaks.nFeat[which.max(Count.nFeat)], Breaks.nFeat[which.max(Count.nFeat)+1] )
  Breaks.nCount = seq(0, max(CombinedObj_1$nCount_RNA), by=10)
  Count.nCount = table(cut(CombinedObj_1$nCount_RNA, breaks=Breaks.nCount, include.lowest=T, right=T) )
  Vec.Frequent.nCount = c(Breaks.nCount[which.max(Count.nCount)], Breaks.nCount[which.max(Count.nCount)+1] )
  n_after = ncol(CombinedObj_1_QC) %>% format(big.mark=",")
  pct_qc = formatC((100*ncol(CombinedObj_1_QC)/ncol(CombinedObj_1)),digits=3)
  p = ggarrange(VlnPlot(CombinedObj_1, features = c("nFeature_RNA"), pt.size = 0, layer = "counts", raster=T) & 
                labs(subtitle=paste0("BeforeQC\nn=",n_before,"cells\nMost frequent : ",Vec.Frequent.nFeat[1],"~",Vec.Frequent.nFeat[2]), x=NULL) &
                geom_hline(yintercept=Vec.Frequent.nFeat) &
                geom_hline(yintercept=nFeatRNA, linetype="dashed") &
                theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="none", 
                      panel.grid.major.y=element_line(color="gray90"), panel.grid.minor.y=element_line(color="gray90")),
              VlnPlot(CombinedObj_1, features = c("nCount_RNA"), pt.size = 0, layer = "counts", raster=T) & 
                labs(subtitle=paste0("BeforeQC\nn=",n_before,"cells\nMost frequent : ",Vec.Frequent.nCount[1],"~",Vec.Frequent.nCount[2]), x=NULL) &
                geom_hline(yintercept=Vec.Frequent.nCount) &
                geom_hline(yintercept=nCountRNA, linetype="dashed") &
                theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="none", 
                      panel.grid.major.y=element_line(color="gray90"), panel.grid.minor.y=element_line(color="gray90")),
              FeatureScatter(CombinedObj_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster=T) &
                geom_hline(yintercept=nFeatRNA, linetype="dashed") &
                geom_vline(xintercept=nCountRNA, linetype="dashed") &
                theme(legend.position="none"),
              VlnPlot(CombinedObj_1_QC, features = c("nFeature_RNA"), pt.size = 0, layer = "counts", raster=T) & 
                labs(subtitle=paste0("AfterQC\n",
                                     " nFeature_RNA : ",nFeatRNA[1],"~",nFeatRNA[2],"\n",
                                     " nCount_RNA : ",nCountRNA[1],"~",nCountRNA[2],"\nn=",n_after,"cells (",pct_qc,"%)"), x=NULL) &
                theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="none", 
                      panel.grid.major.y=element_line(color="gray90"), panel.grid.minor.y=element_line(color="gray90")),
              VlnPlot(CombinedObj_1_QC, features = c("nCount_RNA"), pt.size = 0, layer = "counts", raster=T) & 
                labs(subtitle=paste0("AfterQC\n",
                                     " nFeature_RNA : ",nFeatRNA[1],"~",nFeatRNA[2],"\n",
                                     " nCount_RNA : ",nCountRNA[1],"~",nCountRNA[2],"\nn=",n_after,"cells (",pct_qc,"%)"), x=NULL) &
                theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="none", 
                      panel.grid.major.y=element_line(color="gray90"), panel.grid.minor.y=element_line(color="gray90")),
              FeatureScatter(CombinedObj_1_QC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster=T) &
                theme(legend.position="none"),
              ncol=3, nrow=2, widths=c(1,1,1))
  fs::dir_create(path=paste0(DirInteg))
  ggsave(plot=p,
         file=paste0(DirInteg,"[Figure][QC_VlnPlot][",NumOfSamples,"case(",paste(TXnumInteg,collapse=","),")]_",
                     QCInfo.FileName,".png"),
         width=10, height=10, dpi=100, bg="white")
  return(CombinedObj_1_QC)
}
CombinedObj_2 = Fuc.MinorQC(nFeatRNA=c(100,900), nCountRNA=c(100,1800))
rm(CombinedObj_1)
beep(expr=NULL, sound=3)
### 1. FindVariableFeatures() ####
### 2. ScaleData() ####
### 3. RunPCA() ####
CombinedObj_2 <- CombinedObj_2 %>% 
                 FindVariableFeatures(selection.method="vst", nfeatures=2000) 
CombinedObj_2 <- ScaleData(CombinedObj_2, features=VariableFeatures(CombinedObj_2))
set.seed(123)
CombinedObj_2 <- RunPCA(CombinedObj_2, features=VariableFeatures(CombinedObj_2), npcs=50, verbose=FALSE)
### 4. RunHarmony() ... Harmony uses PCA reduction ####
library(harmony)
set.seed(123)
CombinedObj_3 = RunHarmony(CombinedObj_2,
                           group.by.vars="orig.ident",
                           reduction.use="pca",
                           dims.use=1:50,
                           assay.use="RNA")
beep(expr=NULL, sound=3)
fs::dir_create(path=paste0(DirInteg,"/Objects"))
saveRDS(CombinedObj_3, file=paste0(DirInteg,"Objects/[SeuObj][",NumOfSamples,"case(",paste(TXnumInteg,collapse=","),")]_",
                                   QCInfo.FileName,"_1stPCA50vf2000.rds"))
rm(CombinedObj_2)
### 5. Clustering ... FindNeighbors(), FindClusters() ####
ggsave(plot = ElbowPlot(object=CombinedObj_3, ndims=50) &
              labs(subtitle=paste0(NumOfSamples,"case Integration(",IntegArgo,", TX5K_",paste(TXnumInteg,collapse=","),")\n",
                                   QCInfo,"\n1stPCA50vf2000")) &
              scale_x_continuous(breaks=c(5,10,15,20,25,30,35,40,45)) &
              theme(panel.grid.major.x=element_line(color="gray80")),
       file = paste0(DirInteg,"[Figure][InitialClust_PCA][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg,collapse=","),")]","_1stPCA50vf2000.png"), 
       dpi = 200, width = 6, height = 6)
Func.IniClust = function(Dim1, Res1){
    CombinedObj_3 = readRDS(file=paste0(DirInteg,"Objects/[SeuObj][",NumOfSamples,"case(",paste(TXnumInteg,collapse=","),")]_",
                                        QCInfo.FileName,"_1stPCA50vf2000.rds"))
    old_fgm <- getOption("future.globals.maxSize")
    options(future.globals.maxSize = 128 * 1024^3)
    set.seed(123)
    CombinedObj_4 = RunUMAP(CombinedObj_3, reduction="harmony", dims=1:Dim1) %>% 
                    FindNeighbors(reduction="harmony", dims = 1:Dim1)
    options(future.globals.maxSize = old_fgm)
    set.seed(123)
    CombinedObj_4 = FindClusters(CombinedObj_4, dims = 1:Dim1, resolution = Res1)
    saveRDS(CombinedObj_4, paste0(DirInteg,"Objects/[SeuObj][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                                  QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".rds") )
}
Func.IniClust(Dim1=20, Res1=1) # 6 Case
Func.IniClust(Dim1=20, Res1=1) # 7 Case
Func.IniClust(Dim1=25, Res1=1) # 7 Case
beep(expr=NULL, sound=3)



###
###############################################################
###    Lower analysis after initial clustering  
### 1. General lower analysis ####
Func.GeneralLower.Ini = function(Dim1,Res1){
  ####################
  # 0. read seurat obj
  SeuObjIni_0 = readRDS(paste0(DirInteg,"Objects/[SeuObj][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                               QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".rds") )
  CommonTitle.Ini = paste0(NumOfSamples,"case Integration(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")\n",
                           QCInfo,"\n1stClust:PCA50vf2000_ClustDim",Dim1,"Res",formatC(Res1,digits=2,format="f"))
  NumOfClust = length(unique(SeuObjIni_0$seurat_clusters))
  CommonCaption.Ini = paste0(format(ncol(SeuObjIni_0), big.mark=",", scientific=F), " cells(AfterQC), ",NumOfClust, " clusters.")
  Ncol=ceiling(NumOfClust^0.5)
  ####################
  # 1. Dimplot: Fusion & By cluster
  DF.UMAP_0 = Embeddings(SeuObjIni_0, "umap") %>% as.data.frame() %>%
              rownames_to_column(var = "cell_id") %>% 
              dplyr::mutate("seurat_clusters" = SeuObjIni_0$seurat_clusters,
                            "sample" = SeuObjIni_0$orig.ident)
  set.seed(123)
  DF.UMAP_1 = dplyr::slice_sample(DF.UMAP_0, prop=1)
  DF.UMAPcenter = group_by(DF.UMAP_1, seurat_clusters) %>% 
                  dplyr::summarize(umap_1 = median(umap_1), umap_2 = median(umap_2))
  DimPlot_Base = 
    ggplot(DF.UMAP_1, aes(x=umap_1, y=umap_2, color=seurat_clusters)) +
      geom_point(size=0.3, alpha=0.3) +
      scale_x_continuous(breaks=c(-10,-5,0,5,10)) +
      scale_y_continuous(breaks=c(-10,-5,0,5,10)) +
      coord_fixed(ratio=1) + 
      theme(legend.position="none",
            axis.text = element_text(face = "bold"),
            axis.title = element_text(face = "bold"),
            panel.background = element_blank(),
            axis.line = element_line(color="black")) # + facet_wrap(. ~ seurat_clusters)
  ggsave(
    ggarrange(
      DimPlot_Base + geom_text_repel(
                      data=DF.UMAPcenter, 
                      aes(label=seurat_clusters),
                      color="black", box.padding=0.3, point.padding=0.2, segment.alpha=0.5),
      DimPlot_Base + facet_wrap(. ~ seurat_clusters),
      ncol=2, widths=c(3, 4)) %>% 
    annotate_figure(top=text_grob(CommonTitle.Ini), bottom=CommonCaption.Ini),
    file=paste0(DirInteg,"[Figure][InitialClust_DimPlot2][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".png"), 
    dpi=200, width=16, height=8, bg = "white" )
  ####################
  # 2. FeaturePlot　nFeatRNA and nCountRNA
  FeatPlot_ScaleX = scale_x_continuous(limits = range(as.data.frame(SeuObjIni_0@reductions[["umap"]]@cell.embeddings)$umap_1), breaks=c(-10,0,10))
  FeatPlot_ScaleY = scale_y_continuous(limits = range(as.data.frame(SeuObjIni_0@reductions[["umap"]]@cell.embeddings)$umap_2), breaks=c(-10,0,10))
  FeatPlot_RangeFeature = range(SeuObjIni_0$nFeature_RNA)
  FeatPlot_RangeCount = range(SeuObjIni_0$nCount_RNA)
  Plot.GeneralLower2 = 
    ggarrange(FeaturePlot(SeuObjIni_0, reduction="umap", label=TRUE, split.by="seurat_clusters", raster=FALSE,
                          features="nFeature_RNA") + FeatPlot_ScaleX + FeatPlot_ScaleY + 
                plot_layout(ncol=Ncol, nrow=ceiling(NumOfClust/Ncol), guides="collect") & 
                scale_color_viridis(option="magma", limits=FeatPlot_RangeFeature) &
                coord_fixed(ratio=1) & theme(legend.position="right") & DimPlot_Theme &
                patchwork::plot_annotation(title="nFeature_RNA"),
              FeaturePlot(SeuObjIni_0, reduction="umap", label=TRUE, split.by="seurat_clusters", raster=FALSE,
                          features="nCount_RNA") + FeatPlot_ScaleX + FeatPlot_ScaleY +
                plot_layout(ncol=Ncol, nrow=ceiling(NumOfClust/Ncol), guides="collect") & 
                scale_color_viridis(option="magma", limits=FeatPlot_RangeCount) &
                coord_fixed(ratio=1) & theme(legend.position="right") & DimPlot_Theme & 
              patchwork::plot_annotation(title="nCount_RNA"),
              ncol=2, widths=c(1, 1)) +
              patchwork::plot_annotation(title=CommonTitle.Ini)
  ggsave(Plot.GeneralLower2,
         file=paste0(DirInteg,"[Figure][InitialClust_DimPlot2_nFandC][",
                     NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                     "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".png"), 
        dpi=100, width=32, height=12, bg = "white" )
  # 3. FeaturePlot　Cell type marker
#  SeuObjIni_1 = JoinLayers(SeuObjIni_0)
#  MT.NormalizedExp = GetAssayData(SeuObjIni_1, assay="RNA", layer="data")
#  Func = function(Target){
#    df = cbind(as.data.frame(SeuObjIni_1@reductions[["umap"]]@cell.embeddings),
#               "seurat_clusters" = SeuObjIni_0$seurat_clusters,
#               "TargetGene"=MT.NormalizedExp[Target, ]) %>% 
#         dplyr::arrange(TargetGene)
#    ggplot(df, aes(x=umap_1, y=umap_2, color=TargetGene)) + 
#      geom_point(alpha=0.5, size=0.5) +
#      labs(color=NULL, title=Target) +
#      scale_color_gradientn(colors = c("gray90","gray90","gray90","gray90","gray90","#fee5d9","#fcbba1","#ef3b2c","#99000d")) +
#      coord_fixed(ratio=1) +
#      theme(axis.text = element_text(face = "bold"),
#            axis.title = element_text(face = "bold"),
#            legend.text = element_text(face = "bold"),
#            legend.title = element_text(face = "bold"),
#            panel.background = element_blank(),
#            axis.line = element_line(color = "black"))
#    }
#  ggsave(ggarrange(Func("AMY2A"),Func("EPCAM"),Func("INS"),
#                   Func("COL1A1"),Func("PECAM1"),Func("PTPRC"),
#                   nrow=2, ncol=3, align="hv") %>% 
#         annotate_figure(top=text_grob(CommonTitle.Ini), bottom=CommonCaption.Ini),
#         file=paste0(DirInteg,"[Figure][InitialClust_DimPlot2_CellTypeMarkers][",
#                     NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
#                     "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".png"), 
#         dpi=100, width=24, height=12, bg = "white" )
#  
  ###############################################
  # 3. Dot plot : Cell type markers & CAF markers
  DotPlotTheme = theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                       axis.text = element_text(face = "bold", size = 19),
                       legend.text = element_text(face = "bold", size = 15),
                       legend.title = element_text(face = "bold", size = 15),
                       axis.title = element_text(face = "bold", size = 17))
  DotPlotLgdLab = guides(size = guide_legend(title = "Percent\nExpressed"),
                         color = guide_colorbar(title = "Scaled\nAverage\nExpression",
                                                direction="vertical",
                                                frame.colour="black",
                                                ticks.colour="black"))
  DotPlotColScale = scale_color_gradient2(high="#f03b20",mid="gray90",low="#2c7fb8",midpoint=0)
  ggsave(ggarrange(
          DotPlot(SeuObjIni_0, cluster.idents=FALSE,
                  dot.scale=10,
                  features=Vec.MarkerGenes_CellType) &
          labs(y="Cluster", x=NULL) & DotPlotLgdLab & RotatedAxis() & DotPlotTheme & DotPlotColScale,
          DotPlot(SeuObjIni_0, cluster.idents=FALSE, 
                  dot.scale = 10, 
                  features=Vec.CAFmarkers) &
          labs(y="Cluster", x=NULL) & DotPlotLgdLab & RotatedAxis() & DotPlotTheme & DotPlotColScale,
          ncol=1, nrow=2) %>% 
          annotate_figure(top=text_grob(CommonTitle.Ini), bottom=CommonCaption.Ini),
          file=paste0(DirInteg,"[Figure][InitialClust_DotPlot_CellTypeAndCAFmarkers][",
                      NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                      "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".png"), 
          dpi=300, width=20, height=25, bg="white")
  ##############################
  # 4. Dot plot : Variable genes 
  # High SD genes in average expression
  DF.average = AverageExpression(SeuObjIni_0, group.by = "seurat_clusters", assay = "RNA", layer = "data") %>% data.frame()
  DF.average_SD = mutate(DF.average, StaDev = apply(DF.average, 1, sd)) %>% arrange(-StaDev)
  Vec.VarioubleGenes = rownames(DF.average_SD[1:40, ])
  # hierarchical clust
  rho = stats::cor(t(DF.average[Vec.VarioubleGenes, ]), method = "spearman")   # spearmanの相関係数ρ(rho)
  dist = stats::as.dist(1 - rho)
  #これひつよう？ -> attributes(dist)$class == "dist"
  tree = hclust(dist, method = "ward.D")
  #plot(tree)
  Pdot_variable = 
    DotPlot(SeuObjIni_0, dot.scale = 10, 
            features = Vec.VarioubleGenes[tree$order],
            cluster.idents = TRUE) &
    guides(size = guide_legend(title = "Percent\nExpressed"),
           color = guide_colorbar(title = "Scaled\nAverage\nExpression")) &
    RotatedAxis() &
    theme(text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.text = element_text(size = 15),
          axis.title = element_blank()) &
    scale_color_gradientn(colors = c("gray90","#fee5d9","#fcbba1","#ef3b2c","#99000d")) &
    coord_fixed(ratio=1) &
    ggtitle(paste0("Top variable 40 genes\n",CommonTitle.Ini))
  ggsave(Pdot_variable, 
         file=paste0(DirInteg,"[Figure][InitialClust_DotPlot_DotPlot_VariableGenes][",
                     NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                     "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".png"), 
        dpi = 300, width=10, height = 10, bg = "white") 
  # 5. CellGroupTable
  DF.CellGroupTable_0 = data.frame(cell_id=colnames(SeuObjIni_0),
                                   sample=SeuObjIni_0$orig.ident,
                                   group=SeuObjIni_0$seurat_clusters,
                                   row.names=NULL)
  DF.CellGroupTable_1 = DF.CellGroupTable_0 %>% 
    dplyr::arrange(sample, group) %>% 
    dplyr::mutate(group=paste0("Clust.",group))
  DF.CellGroupTable_1$cell_id = str_remove(DF.CellGroupTable_1$cell_id, pattern="TX5K.._")
  DF.CellGroupTable_1$sample = str_remove(DF.CellGroupTable_1$sample, pattern="TX5K_")
  Func.CellGroupTableSubsetAndWrite = function(TX){
    DF_0 = subset(DF.CellGroupTable_1, subset=sample==TX)
    DF_1 = select(DF_0, -sample)
    library("fs")
    fs::dir_create(path=paste0(DirInteg,"/CellGroupTable"))
    write.csv(DF_1, 
              file=paste0(DirInteg,"/CellGroupTable/[InitialClust_CellGroupTable_TX5K_",TX,"][",
                          NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                          "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".csv"),
              row.names=FALSE)
  }
  for(i in 1:NumOfSamples){Func.CellGroupTableSubsetAndWrite(TX=TXnumInteg[i])}
  ##############################
  # 6. DEG 
  SeuObjIni_0 = JoinLayers(SeuObjIni_0)
  DF.DEG_0 = FindAllMarkers(SeuObjIni_0, 
                            only.pos=T, logfc.threshold=0.8,
                            min.pct=0.2, min.cell.feature=0.01*(ncol(SeuObjIni_0)), 
                            min.cells.group=3, 
                            return.thresh=0.01)
  DF.DEG_1 = DF.DEG_0 %>% 
    dplyr::mutate(PctRatio = pct.1/pct.2) %>% 
    dplyr::select(-p_val) %>% 
    relocate(gene, .before=avg_log2FC)
  write.csv(DF.DEG_1, 
            file=paste0(DirInteg,"DEG_table/[InitialClust_log2FC-0.8_adjP-0.01_pct1-0.2][",
                        NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                        "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".csv"),
            row.names=FALSE)
  }
Func.GeneralLower.Ini(Dim1=20, Res1=1.0)
Func.GeneralLower.Ini(Dim1=25, Res1=1.0)

### 2. Cluster overlapping with 6 cases integration ####
Func.Overlap.Bar = function(Dim1,Res1){
  i=1
  TX=TXnumInteg[i]
  DirInteg_6cases = paste0('/Volumes/Extreme SSD/Analysis/Data/IntegAnalysis/',
                           "[6case(",paste(TXnumInteg[1:6],collapse=","),")]",
                           "_nFeatRNA:",paste(nFeatRNA,collapse="~"),"_nCountRNA:",paste(nCountRNA,collapse="~"),"/")
  DF.CellGroup_6Cases = read.csv(
    file=paste0(DirInteg_6cases,"/CellGroupTable/[InitialClust_CellGroupTable_TX5K_",TX,
                "][6case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:6],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim20Res1.csv") ) %>% 
                        set_colnames(c("cell_id", "Integ6cases") )
  DF.CellGroup_New = read.csv(file=paste0(DirInteg,"/CellGroupTable/[InitialClust_CellGroupTable_TX5K_",TX,"][",
                         NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                         "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".csv") ) %>% 
                     set_colnames(c("cell_id", "IntegNew") )
  DF.CellGroupJoined = left_join(DF.CellGroup_6Cases,
                                 DF.CellGroup_New,
                                 by="cell_id")
  DF.Count_0 = dplyr::count(DF.CellGroupJoined, Integ6cases, IntegNew, name="n") %>% 
               group_by(Integ6cases) %>% 
               dplyr::mutate(prop=n/sum(n)) %>% ungroup()
  DF.Count_0$Integ6cases = factor(DF.Count_0$Integ6cases, levels=unique(DF.CellGroup_6Cases$Integ6cases))
  DF.Count_1 = dplyr::arrange(DF.Count_0, Integ6cases, -prop)
  DF.Count_2 = DF.Count_1 %>% 
               group_by(Integ6cases) %>%
               dplyr::slice(1)
  Vec.NewClustOrder_0 = DF.Count_2$IntegNew
  Vec.NewClustOrder_1 = str_remove(Vec.NewClustOrder_0, pattern="Clust.")
  Vec.NewClustOrder_2 = unique(Vec.NewClustOrder_1) %>% as.numeric()
  Vec.AllClustNew = as.numeric(0:(length(unique(DF.CellGroup_New$IntegNew))-1))
  Vec.NewClustOrder_3 = c(Vec.NewClustOrder_2,
                          setdiff(Vec.AllClustNew, Vec.NewClustOrder_2))
  Str.NewClustOrder = paste(Vec.NewClustOrder_3, collapse=",")
  write.csv(Str.NewClustOrder, 
          file=paste0(DirInteg,"/[InitialClust_ClustOrder][",
                        NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                        "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".csv") )
}
Func.Overlap.Bar(Dim1=25, Res1=1.0)

### 3. General lower analysis final ####
Func.GeneralLower.Ini.Final = 
  function(Dim1,Res1,Order.Ini,
           Acinar,Ductal,Islet,CAF,Mural,Endothelial,Lymph_T,Lymph_B,
           Plasma,Myeloid,Mast,Nerve,Unclassified){
  ###################
  # 1. Reorder clusts
  SeuObjIni_0 = readRDS(paste0(DirInteg,"Objects/[SeuObj][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                               QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".rds") )
  CommonTitle.Ini = paste0(NumOfSamples,"case Integration(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")\n",
                           QCInfo,"\n1stClust:PCA50vf2000_ClustDim",Dim1,"Res",formatC(Res1,digits=2,format="f"))
  NumOfClust = length(unique(SeuObjIni_0$seurat_clusters))
  CommonCaption.Ini = paste0(format(ncol(SeuObjIni_0), big.mark=",", scientific=F), " cells(AfterQC), ",NumOfClust, " clusters.")
  DF.MetaData_0 = 
    Embeddings(SeuObjIni_0, "umap") %>% as.data.frame() %>%
    rownames_to_column(var = "cell_id") %>% 
    dplyr::mutate("CoordX"=SeuObjIni_0$X,
                  "CoordY"=SeuObjIni_0$Y,
                  "seurat_clusters" = SeuObjIni_0$seurat_clusters,
                  "sample" = str_remove(SeuObjIni_0$orig.ident, pattern="TX5K_"))
  Vec.OrderedCnumToCTLAB = 
    case_when(Order.Ini %in% List.Annot[["Acinar"]] ~ "Acinar cell",
              Order.Ini %in% List.Annot[["Ductal_like_acinar"]] ~ "Ductal-like acinar cell",
              Order.Ini %in% List.Annot[["Normal_ductal"]] ~ "Ductal cell (normal)",
              Order.Ini %in% List.Annot[["PanIN"]] ~ "Ductal cell (PanIN)",
              Order.Ini %in% List.Annot[["PDAC"]] ~ "Ductal cell (PDAC)",
              Order.Ini %in% List.Annot[["Islet"]] ~ "Islet cell",
              Order.Ini %in% List.Annot[["CAF"]] ~ "CAF",
              Order.Ini %in% List.Annot[["Mural"]] ~ "Mural cell",
              Order.Ini %in% List.Annot[["Endothelial"]] ~ "Endothelial cell/Pericyte",
              Order.Ini %in% List.Annot[["Lymph_T"]] ~ "T Lymphocyte",
              Order.Ini %in% List.Annot[["Lymph_B"]] ~ "B Lymphocyte",
              Order.Ini %in% List.Annot[["Plasma"]] ~ "Plasma celll",
              Order.Ini %in% List.Annot[["Myeloid"]] ~ "Myeloid cell",
              Order.Ini %in% List.Annot[["Mast"]] ~ "Mast cell",
              Order.Ini %in% List.Annot[["Nerve"]] ~ "Neuronal cell",
              Order.Ini %in% List.Annot[["Unclassified"]] ~ "Unclassified",
              TRUE ~ "Others" )
  names(Vec.OrderedCnumToCTLAB) = Order.Ini
  Vec.OrderedCnumToCT = 
    unlist(lapply(names(List.Annot), function(ct){
      setNames(rep(ct, length(List.Annot[[ct]])),
               List.Annot[[ct]]) 
      }) )
  Vec.OrderedCnumToLab = paste0("C",names(Vec.OrderedCnumToCTLAB),".",unname(Vec.OrderedCnumToCTLAB))
  names(Vec.OrderedCnumToLab) = names(Vec.OrderedCnumToCTLAB)
  DF.MetaData_1 = DF.MetaData_0 %>% 
                  dplyr::mutate(Label = Vec.OrderedCnumToLab[as.character(DF.MetaData_0$seurat_clusters)] )
  Vec.OrderedCnumToColor = Vec.colorcode.Final[unname(Vec.OrderedCnumToCT)]
  names(Vec.OrderedCnumToColor) = Order.Ini
  SeuObjIni_0$Label = Vec.OrderedCnumToLab[as.character(SeuObjIni_0$seurat_clusters)] %>% unname()
  ###################
  # 2. DotPlot
  SeuObjIni_0$Label = factor(SeuObjIni_0$Label, 
                             levels=Vec.OrderedCnumToLab[as.character(Order.Ini)])
  Idents(SeuObjIni_0) = "Label"
  DotPlotTheme = theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                       axis.text = element_text(face = "bold", size = 19),
                       legend.text = element_text(face = "bold", size = 15),
                       legend.title = element_text(face = "bold", size = 15),
                       axis.title = element_text(face = "bold", size = 17))
  DotPlotLgdLab = guides(size = guide_legend(title = "Percent\nExpressed"),
                         color = guide_colorbar(title = "Scaled\nAverage\nExpression",
                                                direction="vertical",
                                                frame.colour="black",
                                                ticks.colour="black"))
  DotPlotColScale = scale_color_gradient2(high="#f03b20",mid="gray90",low="#2c7fb8",midpoint=0)
  ggsave(
      DotPlot(SeuObjIni_0,
              dot.scale=10,
              features=c("AMY2A","CPA1",
                         "C6","HNF1B","FGFR3","SOX9","EPCAM","CDH1",
                         "MUC5B","TFF1","MUC1","MUC5AC","CLDN18","KRT19","CTSE","TNS4","SLC2A1","LAMC2",
                         "COL17A1",
                         "CDH3",
                         "FAM83A",
                         "MACC1",
                         "DKK1",#"CEACAM5",
                         "INS","SCG2",
                         "COL5A1","COL1A1","FAP","PDGFRB","ACTA2","TAGLN","RGS5","NOTCH3","MYH11",
                         "PECAM1","CD34","IL3RA","FLT4",
                         "PTPRC","CD3E","CD8A","GZMA","CD4","CD19","CD79A",
                         "CD27","CD38","TNFRSF17","MZB1",
                         "CSF1R","CD14","ITGAM","CD68","FCGR3A","ITGAX","LY75","FCER1A",
                         "KIT","CMA1",
                         "SCN7A","MPZ","NCAM1")) &
      labs(y=NULL, x=NULL, subtitle=CommonTitle.Ini, caption=CommonCaption.Ini) & 
      DotPlotLgdLab & RotatedAxis() & DotPlotTheme & DotPlotColScale,
  file=paste0(DirInteg,"[Figure][2InitialClust_DotPlot_CellTypeAndCAFmarkers_Final][",
              NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
              "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".png"), 
  dpi=300, width=20, height=12.5, bg="white")
  # 3. DimPlot
  set.seed(123)
  DF.UMAP_1 = dplyr::slice_sample(DF.MetaData_1, prop=1)
  DF.UMAPcenter = group_by(DF.UMAP_1, seurat_clusters) %>% 
    dplyr::summarize(umap_1 = median(umap_1), umap_2 = median(umap_2)) %>% 
    dplyr::mutate(Label = paste0("C",seurat_clusters,".",
                                 Vec.OrderedCnumToCT[as.character(seurat_clusters)] ) )
  #DimPlot_Base = 
    ggplot(DF.UMAP_1, aes(x=umap_1, y=umap_2, color=seurat_clusters)) +
    geom_point(size=0.3, alpha=0.2) +
    #scale_x_continuous(breaks=c(-10,-5,0,5,10)) +
    #scale_y_continuous(breaks=c(-10,-5,0,5,10)) +
    scale_color_manual(values=Vec.OrderedCnumToColor) +
    coord_fixed(ratio=1) + 
    theme(legend.position="none",
          axis.text = element_text(face = "bold"),
          axis.title = element_text(face = "bold"),
          panel.background = element_blank(),
          axis.line = element_line(color="black")) # + facet_wrap(. ~ seurat_clusters)
  ggsave(
    ggarrange(
      DimPlot_Base + geom_text_repel(
        data=DF.UMAPcenter, 
        aes(label=Annot),
        color="black", box.padding=0.3, point.padding=0.2, segment.alpha=0.5),
      DimPlot_Base + facet_wrap(. ~ seurat_clusters),
      ncol=2, widths=c(3, 4)) %>% 
      annotate_figure(top=text_grob(CommonTitle.Ini), bottom=CommonCaption.Ini),
    file=paste0(DirInteg,"[Figure][InitialClust_DimPlot2_Final][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".png"), 
    dpi=200, width=16, height=8, bg = "white" )
  ###################
  # 4. CellGroupTable
  DF.CellGroupTable_0 = data.frame(cell_id=colnames(SeuObjIni_0),
                                   sample=SeuObjIni_0$orig.ident,
                                   group=SeuObjIni_0$Label,
                                   row.names=NULL)
  DF.CellGroupTable_1 = DF.CellGroupTable_0 %>% 
                        dplyr::arrange(sample, group)
  DF.CellGroupTable_1$cell_id = str_remove(DF.CellGroupTable_1$cell_id, pattern="TX5K.._")
  DF.CellGroupTable_1$sample = str_remove(DF.CellGroupTable_1$sample, pattern="TX5K_")
  Func.CellGroupTableSubsetAndWrite = function(TX){
    DF_0 = subset(DF.CellGroupTable_1, subset=sample==TX)
    DF_1 = select(DF_0, -sample)
    library("fs")
    fs::dir_create(path=paste0(DirInteg,"/CellGroupTable"))
    write.csv(DF_1, 
              file=paste0(DirInteg,"/CellGroupTable/[InitialClust_CellGroupTable_TX5K_",TX,"_Final][",
                          NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                          "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".csv"),
              row.names=FALSE)
  }
  for(i in 1:NumOfSamples){Func.CellGroupTableSubsetAndWrite(TX=TXnumInteg[i])}
  ###################
  # 5. Cluster proportions of each sample
  CommonTheme_bar = theme(axis.text = element_text(face = "bold"),
                          axis.title = element_text(face = "bold"),
                          legend.text = element_text(face = "bold"),
                          legend.title = element_text(face = "bold"),
                          legend.key = element_blank(),
                          panel.background = element_blank(),
                          panel.grid = element_blank(),
                          axis.line = element_line(color = "gray20"))
  DF.CellTable_0 = table(DF.MetaData_1$seurat_clusters, DF.MetaData_1$sample) %>% as.data.frame() %>% 
    set_colnames(c("cluster","sample","cell"))
  DF.CellsBySample_Labs_0 = DF.CellTable_0 %>% 
    group_by(sample) %>% 
    dplyr::mutate(total=sum(cell)) %>% 
    ungroup() %>% rowwise() %>% 
    dplyr::mutate(prop=cell/total) %>% 
    dplyr::mutate(label=paste0(formatC(cell, big.mark=","),"\n",
                               "(",formatC(prop*100, format="f", digits=1),"%)")) %>% 
    ungroup() %>% group_by(sample) %>% 
    dplyr::mutate(cumsum=cumsum(prop)) %>% 
    dplyr::mutate(coord_y=1-cumsum+prop/2)
  DF.CellsByClust_Labs_0 = dplyr::arrange(DF.CellTable_0, cluster) %>% group_by(cluster) %>% 
    dplyr::mutate(total=sum(cell)) %>% 
    ungroup() %>% rowwise() %>% 
    dplyr::mutate(prop=cell/total) %>% 
    dplyr::mutate(label=paste0(formatC(cell, big.mark=","),"\n",
                               "(",formatC(prop*100, format="f", digits=1),"%)")) %>% 
    ungroup() %>% group_by(cluster) %>% 
    dplyr::mutate(cumsum=cumsum(prop)) %>% 
    dplyr::mutate(coord_y=1-cumsum+prop/2)
  DF.CellTable_0$cluster = factor(DF.CellTable_0$cluster,
                                  levels=Order.Ini)
  #BarPlot = 
    ggarrange(
    ggplot(DF.CellTable_0, aes(x=sample, y=cell, fill=cluster)) +
      geom_bar(stat="identity", position="fill", alpha=0.8, color="gray50") +
      #geom_text(data=DF.CellsBySample_Labs_0, 
      #          aes(y=coord_y, label=case_when(prop>0.05 ~label,
      #                                         prop>0.02 ~str_replace(label, pattern="\n", replace="  "),
      #                                         TRUE ~ ""))) +
      geom_text(data=dplyr::slice(group_by(DF.CellsByClust_Labs_0, sample), 1), 
                aes(label=formatC(total, big.mark=",")), y=1.03) +
      labs(y="Proportion", x="Sample", fill="Initial cluster") +
      scale_y_continuous(labels=scales::percent_format(),
                         expand=expansion(mult=c(0, 0.05))) +
      scale_fill_manual(values=Vec.OrderedCnumToColor) +
      CommonTheme_bar,
    ggplot(DF.CellTable_0, aes(x=cluster, y=cell, fill=sample)) +
      geom_bar(stat="identity", position="fill", color="gray50") +
      geom_text(data=DF.CellsByClust_Labs_0, 
                aes(y=coord_y, label=case_when(prop>0.10 ~label,
                                               prop>0.02 ~str_replace(label, pattern="\n", replace="  "),
                                               TRUE ~ ""))) +
      geom_text(data=dplyr::slice(group_by(DF.CellsByClust_Labs_0, cluster), 1), 
                aes(label=formatC(total, big.mark=",")), y=1.02) +
      labs(y="Proportion", x="Initial cluster", fill="Sample") +
      scale_y_continuous(labels=scales::percent_format(),
                         expand=expansion(mult=c(0, 0.05))) +
      CommonTheme_bar,
    ncol=1, nrow=2, align="hv") %>% 
    annotate_figure(top=text_grob(CommonTitle.Ini), bottom=CommonCaption.Ini)
  ggsave(BarPlot, 
         file=paste0(DirInteg,"[Figure][SubClust_BarPlot_CellProportion][",
                     NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                     "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"),
         dpi = 300, width = 13, height = 13, bg = "white") 
  ###################
  # 6. Save DataList
  DataList = list("OrderedClustnumToColor" = Vec.OrderedCnumToColor,
                  "OrderedClustnumToLab" = Vec.OrderedCnumToLab,
                  "OrderedClustnumToCelltype" = Vec.OrderedCnumToCT,
                  "MetaData"=DF.MetaData_1)
  saveRDS(DataList,
    file=paste0(DirInteg,"[MetaDataList][InitialClust][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".rds"))
  }
#Func.GeneralLower.Ini.Final(Dim1=20, Res1=1.0, Order.Ini=c(1,3,24,27,19,7,0,33,29,22,26,9,30,23,25,
#                                                           11,2,8,16,12,37,17,10,35,13,14,5,21,15,18,38,4,6,
#                                                           20,36,28,34,32,31,39))
Func.GeneralLower.Ini.Final(Dim1=20, Res1=1.0, 
                            Order.Ini=c(0,13,24,33,12,1,17,30,23,10,27,20,
                                        14,9,6,7,11,34,3,32,2,8,16,15,18,37,5,4,21,
                                        19,31,26,29,25,22,28,36,35),
                            List.Annot = list(Acinar=c(0,13,24), Ductal_like_acinar=c(33,1,12,17,30),
                                              Normal_ductal=23, PanIN=c(10,27),
                                              PDAC=c(20), Islet=c(14),
                                              CAF=c(9,6,7), Mural=c(11),
                                              Endothelial=c(34,3),Lymph_T=c(32,2,8,16),
                                              Lymph_B=c(15),Plasma=c(18,37),
                                              Myeloid=c(5,4,21,19,31),Mast=c(26),
                                              Nerve=c(29),Unclassified=c(36,35,25,22,28) ) )   # 6 Case
Func.GeneralLower.Ini.Final(Dim1=20, Res1=1.0, Order.Ini=c(0,1,21,26,33,2,3,24,7,28,20,15,4,5,14,12,36,38,13,16,32,
                                                           9,8,19,18,17,37,22,23,10,6,11,31,27,30,34,39,25,35,29))   # 7 Case 7個目＃４
Func.GeneralLower.Ini.Final(Dim1=20, Res1=1.0, Order.Ini=c(0,16,25,26,20,1,7,14,27,23,29,10,22,13,15,5,2,
                                                           11,35,3,32,4,12,17,19,36,6,9,21,8,24,28,30,34,31,18,33,37))   # 7 Case 7個目＃14 Dim1:20
Func.GeneralLower.Ini.Final(Dim1=25, Res1=1.0, Order.Ini=c(1,13,24,27,22,15,0,29,16,23,31,11,21,14,4,3,9,12,38,5,26,35,2,10,
                                                           32,19,17,18,39,6,8,7,25,34,30,28,40,36,37,20,33))   # 7 Case 7個目＃14 Dim1:25
Func.GeneralLower.Ini.Final(Dim1=20, Res1=1.0, Order.Ini=c(1,2,23,27,24,34,11,0,18,15,20,22,37,14,
                                                           7,3,6,10,36,12,16,32,5,26,21,38,19,17,9,8,33,28,29,4,35,39,
                                                           13,31,30,25))   # 8 Case
beep(expr=NULL, sound=3)


### 4. Cluster overlapping alluvial ####
Func.Overlap = function(Dim1,Res1,NewIntegCases,NewIntegCasesAlphabetical,NewIntegOrder,NewIntegCAFclust){
  library(ggalluvial)
  Order.6Cases=paste0("SixCases_C",
                      c(0,13,24,25,33,1,12,17,30,23,10,27,20,
                        14,9,6,7,11,34,3,32,2,8,16,15,18,37,5,38,4,21,
                        19,31,26,29,39,36,35,22,28) )
  Order.New=paste0(NewIntegCasesAlphabetical,"Cases_C",NewIntegOrder)
  DirInteg_6cases = paste0('/Volumes/Extreme SSD/Analysis/Data/IntegAnalysis/',
                    "[6case(",paste(TXnumInteg[1:6],collapse=","),")]",
                    "_nFeatRNA:",paste(nFeatRNA,collapse="~"),"_nCountRNA:",paste(nCountRNA,collapse="~"),"/")
  DirInteg_New = paste0('/Volumes/Extreme SSD/Analysis/Data/IntegAnalysis/',
                           "[7case(",paste(TXnumInteg[1:NewIntegCases],collapse=","),")]",
                           "_nFeatRNA:",paste(nFeatRNA,collapse="~"),"_nCountRNA:",paste(nCountRNA,collapse="~"),"/")
  Func = function(SampleNum){
  TX=TXnumInteg[[SampleNum]]
  DF.CellGroupTable_6cases_0 = read.csv(
              file=paste0(DirInteg_6cases,"/CellGroupTable/[InitialClust_CellGroupTable_TX5K_",TX,
                          "][6case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:6],collapse=","),")]",
                          "_1stPCA50vf2000_ClustDim20Res1.csv") )
  DF.CellGroupTable_New_0 = read.csv(
              file=paste0(DirInteg_New,"/CellGroupTable/[InitialClust_CellGroupTable_TX5K_",TX,
                          "][",NewIntegCases,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NewIntegCases],collapse=","),")]",
                          "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".csv") )
  Vec.CellID_CAF_6cases = subset(DF.CellGroupTable_6cases_0, subset=group%in%c("Clust.9","Clust.6","Clust.7"))$cell_id
  Vec.CellID_CAF_New = subset(DF.CellGroupTable_New_0, subset=group%in%paste("Clust",NewIntegCAFclust,sep="."))$cell_id
  DF.CellGroupTable_6cases_1 = dplyr::mutate(DF.CellGroupTable_6cases_0,
                                             group=str_replace(DF.CellGroupTable_6cases_0$group,
                                                               pattern="Clust.", replace="SixCases_C"))
  DF.CellGroupTable_New_1 = dplyr::mutate(DF.CellGroupTable_New_0,
                                             group=str_replace(DF.CellGroupTable_New_0$group,
                                                               pattern="Clust.", replace="SevCases_C"))
  DF.CellGroupTable_Joined_0 = inner_join(DF.CellGroupTable_6cases_1, DF.CellGroupTable_New_1, by="cell_id")
  DF.Alluvial_0 = count(DF.CellGroupTable_Joined_0, group.x, group.y, name="Freq")
  DF.Alluvial_0$group.x = factor(DF.Alluvial_0$group.x, levels=Order.6Cases)
  DF.Alluvial_0$group.y = factor(DF.Alluvial_0$group.y, levels=Order.New)
  DF.CellGroupTable_Joined_1 = subset(DF.CellGroupTable_Joined_0,
                                      subset=cell_id%in%c(Vec.CellID_CAF_6cases,Vec.CellID_CAF_New))
  DF.Alluvial_1 = count(DF.CellGroupTable_Joined_1, group.x, group.y, name="Freq")
  DF.Alluvial_1$group.x = factor(DF.Alluvial_1$group.x, levels=Order.6Cases)
  DF.Alluvial_1$group.y = factor(DF.Alluvial_1$group.y, levels=Order.New)
  ggplot(DF.Alluvial_1,
         aes(axis1=group.x, axis2=group.y, y=Freq)) +
    geom_alluvium(aes(fill=group.x), width=1/12, alpha=0.8, discern=TRUE) +
    geom_stratum(width=1/7, fill="grey90", color="grey40", discern=TRUE) +
    geom_text(stat="stratum", size=3, discern=TRUE,
              aes(label = case_when(after_stat(count) > 1500 ~ 
                                    paste0(str_remove(as.character(after_stat(stratum)),pattern=".*_"),"\n",
                                           "(",after_stat(count),")"),
                                    after_stat(count) > 200 ~ 
                                    paste0(str_remove(as.character(after_stat(stratum)),pattern=".*_"),
                                          "(",after_stat(count),")"),
                                    TRUE ~ ""))) +
    scale_x_discrete(limits=c(paste0("6 cases integ.\n(",length(Vec.CellID_CAF_6cases)," CAFs)"),
                              paste0("7 cases integ.\n(",length(Vec.CellID_CAF_New)," CAFs)")),
                     expand=c(0, 0)) +
    scale_y_continuous(expand=c(0, 0)) +
    labs(y=NULL, x="", fill="6 case integration", title=paste0("TX5K_",TX)) +
    theme(axis.text = element_text(face = "bold"),
          axis.title = element_text(face = "bold"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold"),
          legend.key = element_blank(),
          strip.text = element_text(face = "bold"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          plot.margin=unit(c(0.5,1,0,1),"lines"))
  }
  Plot.Alluvials =
    ggarrange(Func(1),Func(2),Func(3),Func(4),Func(5),Func(6),
              common.legend=T, legend="right") %>% 
    annotate_figure(bottom=text_grob(paste0("6 case integration: TX5K_",paste(TXnumInteg[1:6],collapse=","),
                                            "; CAF clusters are C9, 6, 7.\n",
                                            "7 case integration: TX5K_",paste(TXnumInteg[1:NewIntegCases],collapse=","),
                                            "; CAF clusters are C.",paste(NewIntegCAFclust,collapse=", "),"." ) ) )
  ggsave(Plot.Alluvials,
    file=paste0(DirInteg,"[Figure][InitialClust_ClusterOverlapping_Alluvial][",
                NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".png"), 
    dpi=300, width=20, height=12, bg="white")
}
Func.Overlap(Dim1=20, Res1=1.0, 
             NewIntegCases=7, NewIntegCasesAlphabetical="Sev", NewIntegCAFclust=c(15,5,2),
             NewIntegOrder=c(0,16,25,26,20,1,7,14,27,23,29,10,22,13,15,5,2,
                             11,35,3,32,4,12,17,19,36,6,9,21,8,24,28,30,34,31,18,33,37) )
Func.Overlap(Dim1=20, Res1=1.0, 
             NewIntegCases=8, NewIntegCasesAlphabetical="Eig", NewIntegCAFclust=c(7,3,6),
             NewIntegOrder=c(1,2,23,27,24,34,11,18,0,15,20,22,37,14,7,3,6,4,10,36,12,16,32,5,
                             26,21,38,19,17,9,8,33,28,29,35,39,13,31,30,25) )
Func.Overlap(Dim1=25, Res1=1.0, 
             NewIntegCases=7, NewIntegCasesAlphabetical="Sev", NewIntegCAFclust=c(4,3,9),
             NewIntegOrder=c(1,13,24,27,22,15,0,29,16,23,31,11,21,14,4,3,9,12,38,5,26,35,2,10,
                             32,19,17,18,39,6,8,7,25,34,30,28,40,36,37,20,33))

### 5. Xenium view of initial clusters ####
Dim1=20
Res1=1.0
Func.XenViewIniclust = function(Dim1, Res1){
  library(arrow)
  for(i in 1:NumOfSamples){
    TX = TXnumInteg[i]
    DF.Boundaries_0 = read_parquet(paste0("/Volumes/Extreme SSD/Data/TX5K_",
                                              case_when(TX %in% c("15","16") ~ "15_16",
                                                        TX %in% c("18","19","20","22") ~ "18_19_20_22",
                                                        TX %in% c("27","28") ~ "27_28",
                                                        TRUE ~ TX),
                                              "/cell_boundaries.parquet") )
    DF.NucBoundaries_0 = read_parquet(paste0("/Volumes/Extreme SSD/Data/TX5K_",
                                          case_when(TX %in% c("15","16") ~ "15_16",
                                                    TX %in% c("18","19","20","22") ~ "18_19_20_22",
                                                    TX %in% c("27","28") ~ "27_28",
                                                    TRUE ~ TX),
                                          "/nucleus_boundaries.parquet") )
    List.MetaData = readRDS(file=paste0(DirInteg,"[MetaDataList][InitialClust][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                                        "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".rds"))
    CAFclusts = names(List.MetaData[["OrderedClustnumToCelltype"]])[unname(List.MetaData[["OrderedClustnumToCelltype"]])=="CAF"]
    DF.CellGroup_0 = List.MetaData[["MetaData"]] %>% 
                     subset(subset=sample==TX)
    DF.CellGroup_1 = DF.CellGroup_0 %>% 
      dplyr::mutate(cell_id=str_remove(DF.CellGroup_0$cell_id, pattern="TX5K.*_"))
    DF.Boundaries_2 = DF.Boundaries_0 %>% 
      left_join(DF.CellGroup_1, by="cell_id") %>% 
      dplyr::mutate(seurat_clusters = case_when(is.na(seurat_clusters) ~ "Excluded_byQC",
                                                TRUE ~ seurat_clusters))
    DF.NucBoundaries_2 = DF.NucBoundaries_0 %>% 
      left_join(DF.CellGroup_1, by="cell_id") %>% 
      dplyr::mutate(seurat_clusters = case_when(is.na(seurat_clusters) ~ "Excluded_byQC",
                                                TRUE ~ seurat_clusters))
    CommonTitle.Ini = paste0(NumOfSamples,"case Integration(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")\n",
                             QCInfo,"\n1stClust:PCA50vf2000_ClustDim",Dim1,"Res",formatC(Res1,digits=2,format="f"))
    NumOfClust = length(unique(DF.CellGroup_1$seurat_clusters))
    CommonCaption.Ini = paste0(format(nrow(DF.CellGroup_1), big.mark=",", scientific=F),
                               " cells(AfterQC), ",NumOfClust, " clusters.")
    RangeX = c(min(DF.CellGroup_0$CoordX), max(DF.CellGroup_0$CoordX))
    RangeY = c(min(DF.CellGroup_0$CoordY), max(DF.CellGroup_0$CoordY))
    WideOrLong = case_when( diff(RangeY)/diff(RangeX) > 3/2 ~ c("WideOrLong"="Long", "width"=10, "height"=21),
                            diff(RangeY)/diff(RangeX) < 2/3 ~ c("WideOrLong"="Wide", "width"=21, "height"=10),
                            TRUE ~ c("WideOrLong"="Square", "width"=10, "height"=11 ) )
    XenView.Ini.FullLayer = 
      ggplot(DF.Boundaries_2, 
             aes(x=vertex_x, y=vertex_y, group=label_id)) +
      geom_polygon(aes(fill=seurat_clusters), alpha=1.0, color=NA) +
      #geom_polygon(aes(color=seurat_clusters), alpha=1.0, linewidth=0.3, fill=NA) +
      geom_polygon(data=DF.NucBoundaries_2, fill="blue", alpha=0.3, color=NA) +
      #geom_polygon(data=DF.NucBoundaries_2, color="blue", alpha=0.8, linewidth=0.1, fill=NA) +
      labs(title=paste0("Xenium view of initial clusters, TX5K_",TX),
           x=NULL, y=NULL,
           subtitle=CommonTitle.Ini, caption=CommonCaption.Ini) +
      scale_x_continuous(expand=expansion(mult=c(0, 0)),
                         breaks=seq(from=ceiling(RangeX[1]/1000) * 1000,
                                    to=floor(RangeX[2]/1000) * 1000,
                                    by=1000),
                         limits=c(RangeX[1],RangeX[2])) +
      scale_y_reverse(expand=expansion(mult=c(0, 0)),
                      breaks=seq(from=ceiling(RangeY[1]/1000) * 1000,
                                 to=floor(RangeY[2]/1000) * 1000,
                                 by=1000),
                      limits=c(RangeY[2],RangeY[1])) +
      scale_fill_manual(values=
                          #colorspace::lighten(c(List.MetaData[["OrderedClustnumToColor"]],
                          #                      "Excluded_byQC"="gray20"), amount=0.35)) +
                          c(List.MetaData[["OrderedClustnumToColor"]],
                            "Excluded_byQC"="gray20")) +
      #scale_color_manual(values=c(List.MetaData[["OrderedClustnumToColor"]],
      #                           "Excluded_byQC"="gray20")) +
      coord_fixed(ratio=1) + 
      #coord_cartesian(xlim=c(2500,2900), ylim=c(2500, 2900), ratio=1) + 　　#微調整用！
      theme(text = element_text(color="gray50", face="bold"),
            plot.caption = element_text(size = 5),
            legend.position = "none",
            axis.ticks = element_line(color="gray50"),
            axis.line = element_blank(),
            plot.background = element_rect(fill = "transparent", colour = NA),
            panel.background = element_rect(fill = "black", colour = NA),
            panel.border = element_rect(fill = "transparent", colour = "gray50"),
            panel.grid = element_blank())
    fs::dir_create(path=paste0(DirInteg,"/XeniumView_Initial/",QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1))
    saveRDS(XenView.Ini.FullLayer,
           file=paste0(DirInteg,"/XeniumView_Initial/",QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,
                       "/[RDS][IniClust_XeniumView_TX5K_",TX,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                       "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".rds"))
    # Only cytoplasm without Nucleus
    XenView.Ini.OnlyCP = XenView.Ini.FullLayer
    #XenView.Ini.OnlyCP@layers[[1]]$aes_params$alpha <- 1.0 # cell type fill
    XenView.Ini.OnlyCP@layers[[2]] <- NULL # Nucleus fill
    saveRDS(XenView.Ini.OnlyCP,
            file=paste0(DirInteg,"/XeniumView_Initial/",QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,
                        "/[RDS][IniClustOnlyCP_XeniumView_TX5K_",TX,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                        "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".rds"))
    ggsave(XenView.Ini.OnlyCP,
           file=paste0(DirInteg,"/XeniumView_Initial/",QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,
                       "/[Figure][IniClust_XeniumView_TX5K_",TX,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                       "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".png"), 
           width=as.numeric(WideOrLong[["width"]]), height=as.numeric(WideOrLong[["height"]]), 
           dpi=500, bg="transparent", limitsize=FALSE)
    # Only CAF cytoplasm and Nucreus
    ColorsForCAF = setNames(rep(NA, times=NumOfClust), 
                            as.character(0:(NumOfClust-1)))
    ColorsForCAF[CAFclusts] = "#00ff00"
    XenView.Ini.OnlyCAF = XenView.Ini.FullLayer + 
      scale_fill_manual(values = c(ColorsForCAF, "Excluded_byQC"=NA)) #+ 
      #scale_color_manual(values = c(ColorsForCAF, "Excluded_byQC"=NA))
    saveRDS(XenView.Ini.OnlyCAF,
            file=paste0(DirInteg,"/XeniumView_Initial/",QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,
                        "/[RDS][IniClustOnlyCAF_XeniumView_TX5K_",TX,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                        "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".rds"))
  }  
}
Func.XenViewIniclust(Dim1=20, Res1=1.0)

###
###############################################################
###
###    Extract CAFs
### 1. Pre-process ####
Func.ExtractCAF = function(Dim1, Res1, FibroClusts){
  library(harmony)
  # 0. Read Seurat object
  SeuObj_all = readRDS(paste0(DirInteg,"Objects/[SeuObj][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                              QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".rds") )
  SeuObj_CAF_0 = subset(SeuObj_all, seurat_clusters %in% FibroClusts)
  saveRDS(SeuObj_CAF_0,
          file=paste0(paste0(DirInteg,"Objects/[SeuObj_ExtractedCAF][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                             QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".rds")  ) )
  rm(SeuObj_all)
  # 1. NormalizeData, FindvariableFeatures
  SeuObj_CAF_1 = NormalizeData(SeuObj_CAF_0, verbose=FALSE) %>% 
                 FindVariableFeatures(selection.method="vst", nfeatures=2000)
  # 2. ScaleData, RunPCA
  SeuObj_CAF_2 = ScaleData(SeuObj_CAF_1, features = VariableFeatures(SeuObj_CAF_1))
  set.seed(123)
  SeuObj_CAF_2 = RunPCA(SeuObj_CAF_2, npcs = 50, # total number of PCs to compute and store (50 by default)
                        verbose = FALSE, # print the top genes associated with high/low loadings for the PC
                        features = VariableFeatures(SeuObj_CAF_1)) 
  # 3. RunHarmony()
  set.seed(123)
  SeuObj_CAF_3 <- RunHarmony(
    SeuObj_CAF_2,
    group.by.vars = "orig.ident",
    reduction.use = "pca",
    dims.use      = 1:50,
    assay.use     = "RNA")
  saveRDS(SeuObj_CAF_2,
          file=paste0(paste0(DirInteg,"Objects/[SeuObj_ExtractedCAF][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                             QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000.rds")  ) )
  CommonTitle.CAF = paste0(NumOfSamples,"case Integration(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")\n",
                           QCInfo,"\n1stClust:PCA50vf2000_ClustDim",Dim1,"Res",formatC(Res1,digits=2,format="f"))
  ggsave(plot = ElbowPlot(object = SeuObj_CAF_2, ndims=50) &
           labs(subtitle=paste0(CommonTitle.CAF,"\n2ndClust:PCA50vf2000"),
                caption=paste0("All CAFs, ",format(ncol(SeuObj_CAF_2), big.mark=",", scientific=F)," cells.")) &
           scale_x_continuous(breaks=c(5,10,15,20,25,30,35,40,45)) &
           theme(panel.grid.major.x=element_line(color="gray80")),
         file = paste0(DirInteg,"[Figure][SubClust_ElbowPlot][",
                       NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                       "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000.png"), 
         dpi = 100, width = 6, height = 6)
}
# Func.ExtractCAF(Dim1=20, Res1=1.0, FibroClusts=c(11,2,8,16))
Func.ExtractCAF(Dim1=20, Res1=1.0, FibroClusts=c(9,6,7)) # 6 cases
Func.ExtractCAF(Dim1=20, Res1=1.0, FibroClusts=c(15,5,2)) # 7 cases(last14)
Func.ExtractCAF(Dim1=25, Res1=1.0, FibroClusts=c(4,3,9)) # 7 cases(last14)
Func.ExtractCAF(Dim1=20, Res1=1.0, FibroClusts=c(7,3,6)) # 8 cases


###
###############################################################
###
###    Scoring CAFs
### 1. Module score ####
#  1-1. Calculate module score, scoring and ranking
Func.ModScore.Ini = function(Dim1, Res1, Dim2){
  #--------------------
  #  1. Calculate
  #####################
  SeuObj_CAF_0 = readRDS(paste0(paste0(DirInteg,"Objects/[SeuObj_ExtractedCAF][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                                       QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".rds")  ))
  CommonTitle.CAF = paste0(NumOfSamples,"case Integration(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")\n",
                           QCInfo,"\n1stClust:PCA50vf2000_ClustDim",Dim1,"Res",formatC(Res1,digits=2,format="f"))
  CommonCaption.CAF = paste0("Total ",format(ncol(SeuObj_CAF_0), big.mark=",", scientific=F), " CAFs.")
  SeuObj_CAF_1 = AddModuleScore(SeuObj_CAF_0,
                                features = list(DF.Markers$Gene[DF.Markers$Classification1=="myCAF"],
                                                DF.Markers$Gene[DF.Markers$Classification1=="iCAF"],
                                                c("CD74","HLA-DRA"),
                                                DF.Markers$Gene[DF.Markers$Classification1=="Buffa"&
                                                                DF.Markers$Panel5k=="Included"],
                                                DF.Markers$Gene[DF.Markers$Classification1=="WinterCore"&
                                                                  DF.Markers$Panel5k=="Included"],
                                                DF.Markers$Gene[DF.Markers$GeneType=="Mix3"&
                                                                DF.Markers$Classification1=="NormoCAF"&
                                                                DF.Markers$Panel5k=="Included"],
                                                DF.Markers$Gene[DF.Markers$GeneType=="Mix3"&
                                                                DF.Markers$Classification1=="HypoCAF"&
                                                                DF.Markers$Panel5k=="Included"],
                                                DF.Markers$Gene[DF.Markers$Classification1=="Proliferation"&
                                                                DF.Markers$Panel5k=="Included"],
                                                DF.Markers$Gene[DF.Markers$Classification1=="KEGGcellularsenescence"&
                                                                DF.Markers$Panel5k=="Included"],
                                                DF.Markers$Gene[DF.Markers$Classification1=="SASP"&
                                                                DF.Markers$Panel5k=="Included"],
                                                DF.Markers$Gene[DF.Markers$Classification1=="P53pathway"&
                                                                DF.Markers$Panel5k=="Included"]),
                                seed = 1234,
                                name = c("myCAF","iCAF","apCAF","BuffaOrig","WinterOrig","NormoCAF","HypoCAF",
                                         "Proliferation","Senescence","SASP","P53"))
  #saveRDS(SeuObj_CAF_1,
  #        file = paste0(DirInteg,"Objects/[SeuObj_ExtractedCAF_ModuleScore][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
  #                      QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".rds")  )
  DF.ModScores_0 = data.frame("cell_id"=colnames(SeuObj_CAF_1),
                              "myCAF_score"=SeuObj_CAF_1$myCAF1,
                              "iCAF_score"=SeuObj_CAF_1$iCAF2,
                              "apCAF_score"=SeuObj_CAF_1$apCAF3,
                              "BuffaOrig"=SeuObj_CAF_1$BuffaOrig4,
                              "WinterOrig"=SeuObj_CAF_1$WinterOrig5,
                              "NormoCAF"=SeuObj_CAF_1$NormoCAF6,
                              "HypoCAF"=SeuObj_CAF_1$HypoCAF7,
                              "Proliferation"=SeuObj_CAF_1$Proliferation8,
                              "Senescence"=SeuObj_CAF_1$Senescence9,
                              "SASP"=SeuObj_CAF_1$SASP10,
                              "P53"=SeuObj_CAF_1$P5311)
  Func.Ranking = function(ColName){
    CutOff = quantile(DF.ModScores_0[,ColName], probs=c(0.95, 0.9, 0.8, 0.7, 0.6), names=F)
    case_when(DF.ModScores_0[,ColName] > CutOff[1] ~ "R1_05pct",
              DF.ModScores_0[,ColName] > CutOff[2] ~ "R2_10pct",
              DF.ModScores_0[,ColName] > CutOff[3] ~ "R3_20pct",
              DF.ModScores_0[,ColName] > CutOff[4] ~ "R4_30pct",
              DF.ModScores_0[,ColName] > CutOff[5] ~ "R5_40pct",
              TRUE ~ "Residual")
  }
  DF.ModScores_1 = DF.ModScores_0 %>% 
    dplyr::mutate(myCAF_score_rank = Func.Ranking("myCAF_score")) %>% 
    dplyr::mutate(iCAF_score_rank = Func.Ranking("iCAF_score")) %>% 
    dplyr::mutate(apCAF_score_rank = Func.Ranking("apCAF_score")) %>% 
    dplyr::mutate(BuffaOrig_rank = Func.Ranking("BuffaOrig")) %>% 
    dplyr::mutate(WinterOrig_rank = Func.Ranking("WinterOrig")) %>% 
    dplyr::mutate(NormoCAF_rank = Func.Ranking("NormoCAF")) %>% 
    dplyr::mutate(HypoCAF_rank = Func.Ranking("HypoCAF")) %>% 
    dplyr::mutate(Proliferation_rank = Func.Ranking("Proliferation")) %>% 
    dplyr::mutate(Senescence_rank = Func.Ranking("Senescence")) %>% 
    dplyr::mutate(SASP_rank = Func.Ranking("SASP")) %>% 
    dplyr::mutate(P53_rank = Func.Ranking("P53"))
  write.csv(DF.ModScores_1,
            file=paste0(DirInteg,"/[DataTable_ExtractedCAF_ModuleScore][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                        QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".csv"),
            row.names=F)
  DF.ModScores_2 = DF.ModScores_1 %>% 
                   dplyr::mutate(sample=str_remove_all(DF.ModScores_0$cell_id, pattern="TX5K|_.*")) %>% 
                   dplyr::mutate(cell_id=str_remove_all(DF.ModScores_0$cell_id, pattern="TX5K.*_"))
  Vec.SignatsCellGroupTable = c("myCAF_score_rank","iCAF_score_rank","apCAF_score_rank",
                                "BuffaOrig_rank","WinterOrig_rank","Proliferation_rank")
  for(m in 1:NumOfSamples){
    TX = TXnumInteg[m]
    DF.ModScores_3 = subset(DF.ModScores_2, sample==TX)
    for(n in 1:length(Vec.SignatsCellGroupTable)){
      TargetSignat = Vec.SignatsCellGroupTable[n]
      DF.ModScores_4 = DF.ModScores_3[,c("cell_id",TargetSignat)] %>% 
                       set_colnames(c("cell_id","group"))
      DF.ModScores_5 = DF.ModScores_4 %>% 
                       dplyr::mutate(group=paste(str_remove_all(TargetSignat, pattern="_.*"),
                                                 DF.ModScores_4$group, sep="_") )
      DF.ALL = read.csv(file=paste0(DirInteg,"CellGroupTable/[InitialClust_CellGroupTable_TX5K_",TX,"_Final][",
                                    NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                                    "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".csv"))
      DF_NotCAF = subset(DF.ALL, subset = !(cell_id %in% DF.ModScores_5$cell_id) )
      DF.ModScores_6 = rbind(DF.ModScores_5, DF_NotCAF)
      write.csv(DF.ModScores_6, 
                file=paste0(DirInteg,"/CellGroupTable/[InitialClust_CellGroupTable_SigRank_TX5K",TX,"_",str_remove_all(TargetSignat, pattern="_.*"),"][",
                            NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                            "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".csv"),
                row.names=FALSE)
    }
  } 
  #--------------------------------------------
  # 2.  visualize by feature plot and vlnplot
  #############################################
  set.seed(123)
  SeuObj_CAF_2 = RunUMAP(SeuObj_CAF_1, reduction = "harmony", dims = 1:50)
  Vec.Features = c("myCAF1","iCAF2","apCAF3","BuffaOrig4","WinterOrig5","NormoCAF6","HypoCAF7",
                   "Proliferation8","Senescence9","SASP10","P5311")
  ggsave(plot=FeaturePlot(SeuObj_CAF_2, reduction="umap", label=F, raster=FALSE, order=T, features=Vec.Features) &
         scale_color_viridis(option="turbo"),
         file=paste0(DirInteg,"/[Figure][ExtractedCAF_SignatureScores_FeaturePlot]_countable_mag",Mag,
                     "_nFeat",paste(nFeatRNA,collapse="-"),"_nCount",paste(nCountRNA,collapse="-"),
                     "_1stPCA50vf2000ClustDim",Dim1,"Res",Res1,"]_2ndPCA50vf2000ClustDim",Dim2,"2.png"),
         dpi=300, width=20, height=15, bg="white")
  FeaturePlot_ScaleX = scale_x_continuous(breaks=c(-10,0,10), limits=c(-10,10))
  FeaturePlot_ScaleY = scale_y_continuous(breaks=c(-10,0,10), limits=c(-10,10))
  FeaturePlot_Theme = theme(panel.grid.major.x=element_line(color="gray95"),
                            panel.grid.major.y=element_line(color="gray95"))
  p.combine.single = FeaturePlot(SeuObj_CAF_2, reduction="umap", label=TRUE, 
                                 raster=FALSE, order=T, ncol=1, combine=FALSE,
                                 features=Vec.Features)
  p.combine.list = lapply(p.combine.single, function(x){
      title_text = x$labels$title
      x = x + ggtitle(NULL)
      x + coord_fixed(ratio=1) + labs(y=paste0(title_text,"\n\numap_2")) +
      FeaturePlot_ScaleX + FeaturePlot_ScaleY + FeaturePlot_Theme + scale_color_viridis(option = "turbo")})
  p.combine = patchwork::wrap_plots(p.combine.list, ncol=1)
  splits   <- sort(unique(SeuObj_CAF_2$seurat_clusters))
  n_split  <- length(splits)
  feat_lims <- lapply(Vec.Features, function(f) {
    v <- FetchData(SeuObj_CAF_2, vars = f)[, 1]
    range(v, na.rm = TRUE)  })
  names(feat_lims) <- Vec.Features
  p.byclust.single <- FeaturePlot(SeuObj_CAF_2, reduction="umap", label=TRUE, split.by="seurat_clusters",
                                  raster=FALSE, order=TRUE, keep.scale=NULL, features=Vec.Features)
  p.byclust.list <- lapply(seq_along(p.byclust.single), function(i){
    feat_index <- ((i - 1) %/% n_split) + 1
    feat       <- Vec.Features[feat_index]
    lims       <- feat_lims[[feat]]
    p.byclust.single[[i]] + 
      coord_fixed(ratio = 1) +
      FeaturePlot_ScaleX +
      FeaturePlot_ScaleY +
      FeaturePlot_Theme +
      scale_color_viridis_c(option="turbo", limits = lims)
    } )
  p.byclust.wrap = patchwork::wrap_plots(p.byclust.list, ncol=length(splits))
  ggsave(ggarrange(p.combine, p.byclust.wrap,
                   ncol=2, widths=c(1.5, length(unique(SeuObj_CAF_2$seurat_clusters))), align="h") %>% 
         annotate_figure(top=text_grob(CommonTitle.CAF)),
         file=paste0(DirInteg,"[Figure][ExtractedCAF_SignatureScores_FeaturePlot2][",
                     NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                     QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"2.png"), 
        dpi=200, width=48, height=36, bg = "white" )
  ggsave(VlnPlot(SeuObj_CAF_2, feature=Vec.Features, pt.size=0, ncol=1) +
         plot_annotation(title=CommonTitle.CAF),
         file=paste0(DirInteg,"[Figure][ExtractedCAF_SignatureScores_VlnPlot][",
                     NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                     QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"2.png"),
         dpi=200, width=0.5+1*(length(unique(SeuObj_CAF_2$seurat_clusters))), height=20, bg = "white" )
  #--------------------------------------------
  # 3. dot plot of each genes
  #############################################
  DotPlotTheme = theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                       axis.text = element_text(face = "bold", size = 19),
                       legend.text = element_text(face = "bold", size = 15),
                       legend.title = element_text(face = "bold", size = 15),
                       axis.title = element_text(face = "bold", size = 17))
  DotPlotLgdLab = guides(size = guide_legend(title = "Percent\nExpressed"),
                         color = guide_colorbar(title = "Scaled\nAverage\nExpression",
                                                direction="vertical",
                                                frame.colour="black",
                                                ticks.colour="black"))
  DotPlotColScale = scale_color_gradient2(high="#f03b20",mid="gray90",low="#2c7fb8",midpoint=0)
  ggsave(ggarrange(
    DotPlot(SeuObj_CAF_1, cluster.idents=FALSE,
            dot.scale=10,
            features=c(Vec.Mix300.NormoCAF,"MAPK12","TRPV4","TAP1","FUCA1") ) &
    geom_vline(xintercept=c(8.5, 82.5)) & annotate(geom="text",
                                                   label="Senescence genes\nwhen logFCcutoff=0.8",
                                                   x=Inf, y=Inf, hjust=1.1, vjust=1.1) &
    labs(y="Cluster", x=NULL) & DotPlotLgdLab & RotatedAxis() & DotPlotTheme & DotPlotColScale,
    DotPlot(SeuObj_CAF_1, cluster.idents=FALSE, 
            dot.scale = 10, 
            features=c(Vec.Mix300.HypoCAF,"WEE1","RACGAP1","CDK2","TPX2") ) &
    geom_vline(xintercept=c(2.5, 35.5)) & annotate(geom="text",
                                                   label="proliferation genes\nwhen logFCcutoff=0.8",
                                                   x=Inf, y=Inf, hjust=1.1, vjust=1.1) &
    labs(y="Cluster", x=NULL) & DotPlotLgdLab & RotatedAxis() & DotPlotTheme & DotPlotColScale,
    ncol=1, nrow=2) %>% 
    annotate_figure(top=text_grob(CommonTitle.CAF)),
    file=paste0(DirInteg,"[Figure][ExtractedCAF_SignatureScores__Dotplot]_countable_mag",Mag,
                "_nFeat",paste(nFeatRNA,collapse="-"),"_nCount",paste(nCountRNA,collapse="-"),
                "_1stPCA50vf2000ClustDim",Dim1,"Res",Res1,"]_2ndPCA50vf2000ClustDim",Dim2,".png"),
    dpi=300, width=20, height=15, bg="white")
  }
Func.ModScore.Ini(Dim1=20, Res1=1.0, Dim2=40)

# 1-2. Bar plot visualizing proportion
Func.ModScore.Ini.RanksAndSamples = function(Dim1, Res1){
  Vec.Signatures.rank.list = c("myCAF_score_rank","iCAF_score_rank","apCAF_score_rank","BuffaOrig_rank","WinterOrig_rank",
                               "NormoCAF_rank","HypoCAF_rank","Proliferation_rank","Senescence_rank")
  DF.Scores_0 = read.csv(file=paste0(DirInteg,"/[DataTable_ExtractedCAF_ModuleScore][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                                     QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".csv"), row.names=1)
  DF.Scores_1 = data.frame("cell_id"=str_remove(rownames(DF.Scores_0), pattern="TX5K.._"),
                           "sample"=str_remove_all(rownames(DF.Scores_0), pattern="TX5K|_.*"),
                           DF.Scores_0)
  for(i in 1:length(Vec.Signatures.rank.list)){
    DF.ScoreTarget_0 = DF.Scores_1[,c("sample",Vec.Signatures.rank.list[i],str_remove(Vec.Signatures.rank.list[i],"_rank"))] %>% 
                       set_colnames(c("sample", "rank", "score"))
    DF.CellTable_0 = table(DF.ScoreTarget_0$rank, DF.Scores_1$sample) %>% as.data.frame() %>% 
                     set_colnames(c("rank","sample","cell"))
    DF.CellsBySample_Labs_0 = DF.CellTable_0 %>% 
      group_by(sample) %>% 
      dplyr::mutate(total=sum(cell)) %>% 
      ungroup() %>% rowwise() %>% 
      dplyr::mutate(prop=cell/total) %>% 
      dplyr::mutate(label=paste0(formatC(cell, big.mark=","),"\n",
                                 "(",formatC(prop*100, format="f", digits=1),"%)")) %>% 
      ungroup() %>% group_by(sample) %>% 
      dplyr::mutate(cumsum=cumsum(prop)) %>% 
      dplyr::mutate(coord_y=1-cumsum+prop/2)
    DF.CellsByClust_Labs_0 = dplyr::arrange(DF.CellTable_0, rank) %>% group_by(rank) %>% 
      dplyr::mutate(total=sum(cell)) %>% 
      ungroup() %>% rowwise() %>% 
      dplyr::mutate(prop=cell/total) %>% 
      dplyr::mutate(label=paste0(formatC(cell, big.mark=","),"\n",
                                 "(",formatC(prop*100, format="f", digits=1),"%)")) %>% 
      ungroup() %>% group_by(rank) %>% 
      dplyr::mutate(cumsum=cumsum(prop)) %>% 
      dplyr::mutate(coord_y=1-cumsum+prop/2)
    CommonTheme_bar = theme(axis.text = element_text(face = "bold"),
                            axis.title = element_text(face = "bold"),
                            legend.text = element_text(face = "bold"),
                            legend.title = element_text(face = "bold"),
                            legend.key = element_blank(),
                            panel.background = element_blank(),
                            panel.grid = element_blank(),
                            axis.line = element_line(color = "gray20"))
    Vec.LabChange = c("R1_05pct"="Rank1 (5%)",
                      "R2_10pct"="Rank2 (10%)",
                      "R3_20pct"="Rank3 (20%)",
                      "R4_30pct"="Rank4 (30%)",
                      "R5_40pct"="Rank5 (40%)",
                      "Residual"="Residual")
    CommonTitle.CAF = paste0(NumOfSamples,"case Integration(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")\n",
                             QCInfo,"\n1stClust:PCA50vf2000_ClustDim",Dim1,"Res",formatC(Res1,digits=2,format="f"))
    CommonCaption.CAF = paste0("Total ",format(nrow(DF.Scores_0), big.mark=",", scientific=F), " CAFs.")
    BarPlot = 
      ggarrange(
      ggplot(DF.CellTable_0, aes(x=sample, y=cell, fill=rank)) +
        geom_bar(stat="identity", position="fill", color="gray50") +
        geom_text(data=DF.CellsBySample_Labs_0, 
                  aes(y=coord_y, label=case_when(prop>0.08 ~label,
                                                 prop>0.02 ~str_replace(label, pattern="\n", replace="  "),
                                                TRUE ~ ""))) +
        geom_text(data=dplyr::slice(group_by(DF.CellsByClust_Labs_0, sample), 1), 
                  aes(label=formatC(total, big.mark=",")), y=1.03) +
        labs(y="Proportion", x="Sample", fill=Vec.Signatures.rank.list[i]) +
        scale_fill_manual(values=rev(c("#002051","#3c4d6e","#7f7c75","#bbaf71", "#fdea45","#FFFFFF")),
                          labels=Vec.LabChange) +
        scale_y_continuous(labels=scales::percent_format(),
                           expand=expansion(mult=c(0, 0.05))) +
        CommonTheme_bar,
      ggplot(DF.CellTable_0, aes(x=rank, y=cell, fill=sample)) +
        geom_bar(stat="identity", position="fill", color="gray50") +
        geom_text(data=DF.CellsByClust_Labs_0, 
                  aes(y=coord_y, label=case_when(prop>0.10 ~label,
                                                 prop>0.02 ~str_replace(label, pattern="\n", replace="  "),
                                                 TRUE ~ ""))) +
        geom_text(data=dplyr::slice(group_by(DF.CellsByClust_Labs_0, rank), 1), 
                  aes(label=formatC(total, big.mark=",")), y=1.02) +
        labs(y="Proportion", x=Vec.Signatures.rank.list[i], fill="Sample") +
        scale_y_continuous(labels=scales::percent_format(),
                           expand=expansion(mult=c(0, 0.05))) +
        scale_x_discrete(labels=Vec.LabChange) +
        CommonTheme_bar,
      ncol=1, nrow=2, align="hv") %>% 
      annotate_figure(top=text_grob(CommonTitle.CAF), bottom=CommonCaption.CAF)
    ggsave(BarPlot, 
           file=paste0(DirInteg,"[Figure][InitialClust_ModuleScores_Barplot_CellProportion_",Vec.Signatures.rank.list[i],"][",
                       NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                       "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"),
           dpi = 300, width = 13, height = 13, bg = "white") 
    }
}
Func.ModScore.Ini.RanksAndSamples(Dim1=20, Res1=1.0)

# 1-3. Xenium view of ranking
Func.XenViewModScores = function(Dim1, Res1){
  DF.Scores_0 = read.csv(file=paste0(DirInteg,"/[DataTable_ExtractedCAF_ModuleScore][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                                     QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".csv"), row.names=1)
  DF.Scores_1 = data.frame("cell_id"=str_remove(rownames(DF.Scores_0), pattern="TX5K.._"),
                           "sample"=str_remove_all(rownames(DF.Scores_0), pattern="TX5K|_.*"),
                           DF.Scores_0)
  library(arrow)
  for(i in 1:NumOfSamples){
  TX = TXnumInteg[i]
  DF.Boundaries.ALL = read_parquet(paste0("/Volumes/Extreme SSD/Data/TX5K_",
                                        case_when(TX %in% c("15","16") ~ "15_16",
                                                  TX %in% c("18","19","20","22") ~ "18_19_20_22",
                                                  TX %in% c("27","28") ~ "27_28",
                                                  TRUE ~ TX),
                                        "/cell_boundaries.parquet") )
  DF.Scores_2 = subset(DF.Scores_1, subset=sample==TX)
  CommonSubtitle.CAF = paste0(NumOfSamples,"case Integration(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")\n",
                           QCInfo,"\n1stClust:PCA50vf2000_ClustDim",Dim1,"Res",formatC(Res1,digits=2,format="f"))
  CommonCaption.CAF.InSample = paste0("Total ",format(nrow(DF.Scores_2), big.mark=",", scientific=F), " CAFs in TX5K_",TX," (afterQC).")
  DF.Boundaries.CAF_0 = subset(DF.Boundaries.ALL, subset = cell_id %in% DF.Scores_2$cell_id)
  DF.Boundaries.CAF_1 = left_join(DF.Boundaries.CAF_0,
                                DF.Scores_2, by="cell_id")
  DF.Boundaries.NotCAF = subset(DF.Boundaries.ALL, subset = !(cell_id %in% DF.Scores_2$cell_id) )
  RangeX = c(min(DF.Boundaries.CAF_1$vertex_x), max(DF.Boundaries.CAF_1$vertex_x))
  RangeY = c(min(DF.Boundaries.CAF_1$vertex_y), max(DF.Boundaries.CAF_1$vertex_y))
  WideOrLong = if( diff(RangeY)>diff(RangeX) ){c("WideOrLong"="Long", "ncol"=2, "nrow"=1)} else {c("WideOrLong"="Wide", "ncol"=1, "nrow"=2)}
  Vec.Colnames.Gradient = c("myCAF_score","iCAF_score","apCAF_score","Hypoxia","NormoCAF","HypoCAF","Proliferation","Senescence")
 for(PlotNumber in 1:length(Vec.Colnames.Gradient)){
  DF.SubsetForPlot.Gradient = DF.Boundaries.CAF_1[,c("cell_id","vertex_x","vertex_y","label_id",Vec.Colnames.Gradient[PlotNumber])] %>% 
                              magrittr::set_colnames(c("cell_id","vertex_x","vertex_y","label_id","Target_signature"))
  Vec.ScoreRange = c(min(DF.SubsetForPlot.Gradient$Target_signature), max(DF.SubsetForPlot.Gradient$Target_signature))
  Mid = case_when(Vec.Colnames.Gradient[PlotNumber]=="Proliferation" ~ 0.6,
                  TRUE ~ 0.6)
  Plot.XenView.Gradient = 
    ggplot(DF.SubsetForPlot.Gradient, 
           aes(x=vertex_x, y=vertex_y, group=label_id)) +
    geom_polygon(aes(fill=Target_signature), color=NA, linewidth=0.1) +
    geom_polygon(data=DF.Boundaries.NotCAF, color=NA, fill="gray20",linewidth=0.1) +
    #scale_fill_identity(name="cell_type", breaks=unname(Vec.colorcode), 
    #                    labels=names(Vec.colorcode), guide="legend") +
    labs(x=NULL, y=NULL, subtitle="Gradient", fill=Vec.Colnames.Gradient[PlotNumber]) +
    scale_x_continuous(expand=expansion(mult=c(0, 0)),
                       breaks=seq(from=ceiling(RangeX[1]/1000) * 1000,
                                  to=floor(RangeX[2]/1000) * 1000,
                                  by=1000),
                       limits=c(RangeX[1],RangeX[2])) +
    scale_y_reverse(expand=expansion(mult=c(0, 0)),
                    breaks=seq(from=ceiling(RangeY[1]/1000) * 1000,
                               to=floor(RangeY[2]/1000) * 1000,
                               by=1000),
                    limits=c(RangeY[2],RangeY[1])) +
    coord_fixed(ratio=1) + 
    theme(text = element_text(face = "bold"),
          plot.caption = element_text(size = 5),
          #legend.key = element_blank(),
          #legend.position = "none",
          axis.ticks = element_line(color="gray10"),
          axis.line = element_line(color="gray10"),
          panel.background = element_rect(fill = "black", colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA),
          panel.grid = element_blank()) +
      scale_fill_gradientn(
        colors = viridis::viridis(256, option = "inferno"),
        values = scales::rescale(c(Vec.ScoreRange[1], Vec.ScoreRange[2]*Mid ,Vec.ScoreRange[2])),
        limits = c(Vec.ScoreRange[1], Vec.ScoreRange[2]))
    
  Vec.Colnames.Rank = paste0(Vec.Colnames.Gradient, "_rank")
  DF.SubsetForPlot.Rank = DF.Boundaries.CAF_1[,c("cell_id","vertex_x","vertex_y","label_id",Vec.Colnames.Rank[PlotNumber])] %>% 
    magrittr::set_colnames(c("cell_id","vertex_x","vertex_y","label_id","Target_signature"))
  DF.SubsetForPlot.Rank$Target_signature = factor(DF.SubsetForPlot.Rank$Target_signature,
                                                  levels=rev(c("R1_05pct","R2_10pct","R3_20pct","R4_30pct","R5_40pct","Residual")))
  Plot.XenView.Rank = 
    ggplot(DF.SubsetForPlot.Rank, 
           aes(x=vertex_x, y=vertex_y, group=label_id)) +
    geom_polygon(aes(fill=Target_signature), color=NA, linewidth=0.1) +
    geom_polygon(data=DF.Boundaries.NotCAF, color=NA, fill="gray20",linewidth=0.1) +
    labs(x=NULL, y=NULL, subtitle="Rank", fill=Vec.Colnames.Rank[PlotNumber]) +
    scale_x_continuous(expand=expansion(mult=c(0.02, 0.02)),
                       breaks=seq(from=ceiling(RangeX[1]/1000) * 1000,
                                  to=floor(RangeX[2]/1000) * 1000,
                                  by=1000),
                       limits=c(RangeX[1],RangeX[2])) +
    scale_y_reverse(expand=expansion(mult=c(0.02, 0.02)),
                    breaks=seq(from=ceiling(RangeY[1]/1000) * 1000,
                               to=floor(RangeY[2]/1000) * 1000,
                               by=1000),
                    limits=c(RangeY[2],RangeY[1])) +
    coord_fixed(ratio=1) + theme(text = element_text(face = "bold"),
                                 plot.caption = element_text(size = 5),
                                 #legend.key = element_blank(),
                                 #legend.position = "none",
                                 axis.ticks = element_line(color="gray10"),
                                 axis.line = element_line(color="gray10"),
                                 panel.background = element_rect(fill = "black", colour = NA),
                                 plot.background = element_rect(fill = "transparent", colour = NA),
                                 panel.grid = element_blank()) +
    scale_fill_manual(values=c("#002051","#3c4d6e","#7f7c75","#bbaf71", "#fdea45","#fdea45"))
  library(fs)
  fs::dir_create(path=paste0(DirInteg,"/XeniumView_Initial/","[",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                                QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1))
  ggsave(ggarrange(Plot.XenView.Gradient, Plot.XenView.Rank,
                   ncol=as.numeric(WideOrLong[["ncol"]]), 
                   nrow=as.numeric(WideOrLong[["nrow"]]), 
                   align="hv") %>% 
          annotate_figure(top=paste0("TX5K_",TX,"\n",CommonSubtitle.CAF,"\n",Vec.Colnames.Gradient[PlotNumber]),
                          bottom=CommonCaption.CAF.InSample),
         file=paste0(DirInteg,"/XeniumView_Initial/","[",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                     QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"/[Figure][InitialClust_SignatureScores_XenView_",Vec.Colnames.Gradient[PlotNumber],"-TX5K",TX,"]_countable_mag",Mag,
                     "_nFeat",paste(nFeatRNA,collapse="-"),"_nCount",paste(nCountRNA,collapse="-"),
                     "_PCA50vf2000_ClustDim",Dim1,"Res",Res1,".png"),,
         dpi=300, width=20, height=20)
 }
  }
}
Func.XenViewModScores(Dim1=20, Res1=1.0)

### 2. GSVA ####
# 2-1. Calculate GSVA score, scoring and ranking, save
Func.GSVA.Initial = function(Dim1, Res1){
  library(msigdbr)
  library(GSVA)
  # 1. Get gene sets
  CategoryH = msigdbr(species = "Homo sapiens", category = "H")
  CategoryH_list = CategoryH %>% split(.$gs_name) %>%
                   lapply(function(x) unique(x$gene_symbol))
  CategoryC2 = msigdbr(species="Homo sapiens",  category="C2")
  AddPathways = list(
    #NormoCAF = Vec.Mix300.NormoCAF,
    #HypoCAF = Vec.Mix300.HypoCAF,
    BuffaOrig = Vec.BuffaOrig,
    WinterOrig = Vec.WinterOrig,
    myCAF = DF.Markers$Gene[DF.Markers$Classification1=="myCAF"],
    iCAF = DF.Markers$Gene[DF.Markers$Classification1=="iCAF"],
    apCAF = c("CD74","HLA-DRA"),
    Proliferation = DF.Markers$Gene[DF.Markers$Classification1=="Proliferation"&
                                    DF.Markers$Panel5k=="Included"],
    KEGG_Cell_cycle = CategoryC2 %>% filter(gs_name == "KEGG_CELL_CYCLE") %>% pull(gene_symbol) %>% unique(),
    REACTOME_Cellular_senescence = CategoryC2 %>% filter(gs_name == "REACTOME_CELLULAR_SENESCENCE") %>% pull(gene_symbol) %>% unique() )
  Pathways = c(AddPathways, CategoryH_list)
  str(CategoryH_list)
  Pathways.FilterUnobservedGenes = 
    lapply(Pathways, function(gs) {
    intersect(gs, Vec.AllGenes)
    })
  # 2.Calculate GSVA score
  set.seed(1234)
  SeuObj_CAF_0 = readRDS(file=paste0(DirInteg,"Objects/[SeuObj_ExtractedCAF][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                                     QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000.rds")  )
  SeuObj_CAF_0[["RNA"]] = JoinLayers(SeuObj_CAF_0[["RNA"]], layers="data") 
  expr = GetAssayData(SeuObj_CAF_0, assay="RNA", layer="data") %>% 
         as.matrix()
  param = ssgseaParam(exprData=expr, 
                      geneSets=Pathways.FilterUnobservedGenes,
                      minSize=1, maxSize=Inf, normalize=TRUE)  # default
  # 3.save with other meta data
  MT.gsva_res = gsva(param, verbose=TRUE)
  MT.gsva_res_trans = t(MT.gsva_res)
  colnames(MT.gsva_res_trans) = paste0(colnames(MT.gsva_res_trans), "_ssGSVA")
  DF.gsva_res_trans_0 = as.data.frame(MT.gsva_res_trans) %>% 
                        rownames_to_column(var="full_id")
  DF.gsva_res_trans_1 = cbind("sample"=str_remove(DF.gsva_res_trans_0$full_id, pattern="_.*"),
                              "cell_id"=str_remove(DF.gsva_res_trans_0$full_id, pattern="TX5K.*_"),
                              "X"=SeuObj_CAF_0$X,
                              "Y"=SeuObj_CAF_0$Y) %>% 
                        cbind(DF.gsva_res_trans_0) %>% 
                        dplyr::select(-full_id)
  write.csv(DF.gsva_res_trans_1,
            file=paste0(DirInteg,"/[DataTable_ExtractedCAF_ssGSVA][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                        QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".csv"),
            row.names=FALSE)
}
Func.GSVA.Initial(Dim1=20, Res1=1.0)

# 2-2. visualize via Xenium plot
Func.GSVA.Initial.XenView = function(Dim1, Res1){
  library(arrow)
  for(i in seq_along(TXnumInteg)){
  TX = TXnumInteg[i]
  DF.Boundaries.ALL = read_parquet(paste0("/Volumes/Extreme SSD/Data/TX5K_",
                                          case_when(TX %in% c("15","16") ~ "15_16",
                                                    TX %in% c("18","19","20","22") ~ "18_19_20_22",
                                                    TX %in% c("27","28") ~ "27_28",
                                                    TRUE ~ TX),
                                          "/cell_boundaries.parquet") )
  Vec.CellIDsAfterQC = read.csv(file=paste0(DirInteg,"CellGroupTable/[InitialClust_CellGroupTable_TX5K_",TX,"_Final][",
                                                    NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                                                    "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".csv"))$cell_id
  DF.Boundaries.ALL = subset(DF.Boundaries.ALL, subset=cell_id%in%Vec.CellIDsAfterQC)
  DF.Boundaries.CAF_0 = subset(DF.Boundaries.ALL, subset=cell_id%in%DF.gsva_res_trans_1$cell_id)
  DF.Boundaries.NotCAF_0 = subset(DF.Boundaries.ALL, subset= !(cell_id%in%DF.gsva_res_trans_1$cell_id))
  Func.GSVAScoreRanking = function(ScoreColmn){case_when(ScoreColmn < quantile(ScoreColmn, probs=seq(0,1,0.2), names=F)[2] ~ "Rank1",
                                                        ScoreColmn < quantile(ScoreColmn, probs=seq(0,1,0.2), names=F)[3] ~ "Rank2",
                                                        ScoreColmn < quantile(ScoreColmn, probs=seq(0,1,0.2), names=F)[4] ~ "Rank3",
                                                        ScoreColmn < quantile(ScoreColmn, probs=seq(0,1,0.2), names=F)[5] ~ "Rank4",
                                                        TRUE ~ "Rank5")}
  DF.gsva_res_trans_1 = 
    read.csv(file=paste0(DirInteg,"/[DataTable_ExtractedCAF_ssGSVA][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                         QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".csv"))
  DF.GSVAscore_1 = dplyr::mutate(DF.gsva_res_trans_1,
                                 "BuffaOrig_rank"=Func.GSVAScoreRanking(DF.gsva_res_trans_1$BuffaOrig_ssGSVA),
                                 "WinterOrig_rank"=Func.GSVAScoreRanking(DF.gsva_res_trans_1$WinterOrig_ssGSVA),
                                 "myCAF_rank"=Func.GSVAScoreRanking(DF.gsva_res_trans_1$myCAF_ssGSVA),
                                 "iCAF_rank"=Func.GSVAScoreRanking(DF.gsva_res_trans_1$iCAF_ssGSVA),
                                 "apCAF_rank"=Func.GSVAScoreRanking(DF.gsva_res_trans_1$apCAF_ssGSVA),
                                 "Proliferation_rank"=Func.GSVAScoreRanking(DF.gsva_res_trans_1$Proliferation_ssGSVA),
                                 "G2M_rank"=Func.GSVAScoreRanking(DF.gsva_res_trans_1$HALLMARK_G2M_CHECKPOINT_ssGSVA),
                                 "E2F_rank"=Func.GSVAScoreRanking(DF.gsva_res_trans_1$HALLMARK_E2F_TARGETS_ssGSVA),
                                 "MITOTIC_SPINDLE_rank"=Func.GSVAScoreRanking(DF.gsva_res_trans_1$HALLMARK_MITOTIC_SPINDLE_ssGSVA),
                                 "MYC_TARGETS_V1_rank"=Func.GSVAScoreRanking(DF.gsva_res_trans_1$HALLMARK_MYC_TARGETS_V1_ssGSVA),
                                 "MYC_TARGETS_V2_rank"=Func.GSVAScoreRanking(DF.gsva_res_trans_1$HALLMARK_MYC_TARGETS_V2_ssGSVA),
                                 "P53_rank"=Func.GSVAScoreRanking(DF.gsva_res_trans_1$HALLMARK_P53_PATHWAY_ssGSVA),
                                 "IFNalpha_rank"=Func.GSVAScoreRanking(DF.gsva_res_trans_1$HALLMARK_INTERFERON_ALPHA_RESPONSE_ssGSVA),
                                 "IFNgamma_rank"=Func.GSVAScoreRanking(DF.gsva_res_trans_1$HALLMARK_INTERFERON_GAMMA_RESPONSE_ssGSVA),
                                 "KEGG_Cell_cycle_rank"=Func.GSVAScoreRanking(DF.gsva_res_trans_1$KEGG_CellCycle_ssGSVA),
                                 "REACTOME_Cellular_senesence_rank"=Func.GSVAScoreRanking(DF.gsva_res_trans_1$REACTOME_Cellular_Senesence_ssGSVA))
  DF.Boundaries.CAF_GSVA = left_join(DF.Boundaries.CAF_0,
                                     DF.GSVAscore_1, by="cell_id")
  RangeX = c(min(DF.Boundaries.ALL$vertex_x), max(DF.Boundaries.ALL$vertex_x))
  RangeY = c(min(DF.Boundaries.ALL$vertex_y), max(DF.Boundaries.ALL$vertex_y))
  WideOrLong = case_when( diff(RangeY)/diff(RangeX) > 3/2 ~ c("WideOrLong"="Long", "ncol"=length(pathways), "nrow"=1, "width"=length(pathways)*10, "height"=20),
                          diff(RangeY)/diff(RangeX) < 2/3 ~ c("WideOrLong"="Wide", "ncol"=1, "nrow"=length(pathways), "width"=20, "height"=length(pathways)*10),
                          TRUE ~ c("WideOrLong"="Square", "ncol"=length(pathways), "nrow"=1, "width"=length(pathways)*10, "height"=20) )
  CommonTheme.XenView.GSVA = 
    theme(text = element_text(face = "bold"),
          plot.caption = element_text(size = 5),
          legend.position = case_when(WideOrLong[["WideOrLong"]]%in%c("Long","Square") ~ "bottom",
                                      WideOrLong[["WideOrLong"]]=="Wide" ~ "right"),
          axis.ticks = element_line(color="gray10"),
          axis.line = element_line(color="gray10"),
          panel.background = element_rect(fill = "black", colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA),
          panel.grid = element_blank())
  CommonScaleX = scale_x_continuous(
                     expand=expansion(mult=c(0.02, 0.02)),
                     breaks=seq(from=ceiling(RangeX[1]/1000) * 1000,
                                to=floor(RangeX[2]/1000) * 1000,
                                by=1000),
                     limits=c(RangeX[1],RangeX[2]))
  CommonScaleY = scale_y_reverse(
                    expand=expansion(mult=c(0.02, 0.02)),
                    breaks=seq(from=ceiling(RangeY[1]/1000) * 1000,
                               to=floor(RangeY[2]/1000) * 1000,
                               by=1000),
                    limits=c(RangeY[2],RangeY[1]))
  PlotList.Gradient = list()
  PlotList.Rank = list()
  for(PlotNumber in 1:length(pathways)){
    Vec.SignatureNames = c("BuffaOrig","WinterOrig","myCAF","iCAF","apCAF","Proliferation",
                           "G2M","E2F","MITOTIC_SPINDLE","MYC_TARGETS_V1","MYC_TARGETS_V2",
                           "P53","IFNalpha","IFNgamma","KEGG_CellCycle","KEGG_CellularSenes")
    DF.SubsetForPlot.Gradient = DF.Boundaries.CAF_GSVA[,c("cell_id",
                                                          "vertex_x",
                                                          "vertex_y",
                                                          "label_id",
                                                          paste0(Vec.SignatureNames[PlotNumber],"_ssGSVA") )] %>% 
      magrittr::set_colnames(c("cell_id","vertex_x","vertex_y","label_id","Target_signature"))
    Plot.XenView.Gradient = 
      ggplot(DF.SubsetForPlot.Gradient, 
             aes(x=vertex_x, y=vertex_y, group=label_id)) +
      geom_polygon(aes(fill=Target_signature), color=NA, linewidth=0.1) +
      geom_polygon(data=DF.Boundaries.NotCAF_0, color=NA, fill="gray15") +
      labs(x=NULL, y=NULL, title=paste0("ssGSVA_",Vec.SignatureNames[PlotNumber],"_gradient"), 
           fill=paste0("ssGSVA\n",Vec.SignatureNames[PlotNumber]) ) +
      coord_fixed(ratio=1) + CommonScaleX + CommonScaleY + CommonTheme.XenView.GSVA +
      scale_fill_viridis(option="inferno")
    DF.SubsetForPlot.Rank = DF.Boundaries.CAF_GSVA[,c("cell_id",
                                                      "vertex_x",
                                                      "vertex_y",
                                                      "label_id",
                                                      paste0(Vec.SignatureNames[PlotNumber],"_rank") )] %>% 
      magrittr::set_colnames(c("cell_id","vertex_x","vertex_y","label_id","Target_signature"))
    PlotList.Gradient[[PlotNumber]] = Plot.XenView.Gradient
    Plot.XenView.Rank = 
      ggplot(DF.SubsetForPlot.Rank, 
             aes(x=vertex_x, y=vertex_y, group=label_id)) +
      geom_polygon(aes(fill=Target_signature), color=NA, linewidth=0.1) +
      geom_polygon(data=DF.Boundaries.NotCAF_0, color=NA, fill="gray15") +
      labs(x=NULL, y=NULL, title=paste0("ssGSVA_",Vec.SignatureNames[PlotNumber],"_rank"),
           fill=paste0("ssGSVA\n",paste0(Vec.SignatureNames[PlotNumber],"_rank") ) ) +
      coord_fixed(ratio=1) + CommonScaleX + CommonScaleY + CommonTheme.XenView.GSVA +
      scale_fill_manual(values=c("#002051","#3c4d6e","#7f7c75","#bbaf71", "#fdea45"))
    PlotList.Rank[[PlotNumber]] = Plot.XenView.Rank
  }
  Plots.Gradient = wrap_plots(PlotList.Gradient, 
                              ncol=as.numeric(WideOrLong[["ncol"]]), 
                              nrow=as.numeric(WideOrLong[["nrow"]]))
  Plots.Rank = wrap_plots(PlotList.Rank, 
                              ncol=as.numeric(WideOrLong[["ncol"]]), 
                              nrow=as.numeric(WideOrLong[["nrow"]]))
  fs::dir_create(path=paste0(DirInteg,"/XeniumView_Initial/",QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"/ssGSVA"))
  ggsave(plot=Plots.Gradient,
         file=paste0(DirInteg,"/XeniumView_Initial/",QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,
                     "/ssGSVA/[Figure][ExtractedCAF_XeniumView_GSVA(Gradient)_TX5K_",TX,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                     "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".png"), 
         width=as.numeric(WideOrLong[["width"]]), height=as.numeric(WideOrLong[["height"]]), 
         dpi=100, limitsize=FALSE)
  ggsave(plot=Plots.Rank,
         file=paste0(DirInteg,"/XeniumView_Initial/",QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,
                     "/ssGSVA/[Figure][ExtractedCAF_XeniumView_GSVA(Rank)_TX5K_",TX,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                     "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".png"), 
         width=as.numeric(WideOrLong[["width"]]), height=as.numeric(WideOrLong[["height"]]), 
         dpi=100, limitsize=FALSE)
  }
}
Func.GSVA.Initial.XenView(Dim1=20, Res1=1.0)


###
###############################################################
###
###    CAFs subclustering
### 1. Subclustering ####
Func.Subclust = function(Dim1, Res1, Dim2, Res2){
  # 1. RunUMAP(), FindNeighbors() : Reduction="harmony"
  SeuObj.CAF_0 = readRDS(paste0(DirInteg,"Objects/[SeuObj_ExtractedCAF][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                                QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000.rds")  )
  set.seed(123)
  SeuObj.CAF_1 = RunUMAP(SeuObj.CAF_0, reduction="harmony", dims=1:Dim2) %>% 
                 FindNeighbors(reduction="harmony", dims=1:Dim2)
  # 2. FindClusters() using SNN graph generated by Findneighbors()
  SeuObj.CAF_2 = FindClusters(SeuObj.CAF_1, dims=1:Dim2, resolution=Res2, random.seed=123)
  saveRDS(SeuObj.CAF_2, 
          paste0(DirInteg,"Objects/[SeuObj_ExtractedCAF][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                 QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2,".rds") ) 
}
#Func.Subclust(Dim1=20, Res1=1.0, Dim2=15, Res2=1.0) ####   <- 6cases BAD
#Func.Subclust(Dim1=20, Res1=1.0, Dim2=15, Res2=0.5) ####   <- 6cases BAD
Func.Subclust(Dim1=20, Res1=1.0, Dim2=20, Res2=0.5)
Func.Subclust(Dim1=20, Res1=1.0, Dim2=30, Res2=0.5)
Func.Subclust(Dim1=20, Res1=1.0, Dim2=40, Res2=0.5)
#Func.Subclust(Dim1=20, Res1=1.0, Dim2=50, Res2=0.5) ####   <- 6cases BAD
Func.Subclust(Dim1=25, Res1=1.0, Dim2=15, Res2=0.5)
Func.Subclust(Dim1=25, Res1=1.0, Dim2=25, Res2=0.5)
Func.Subclust(Dim1=25, Res1=1.0, Dim2=35, Res2=0.5)
beep(expr=NULL, sound=8)







###
###
###
###############################################################
###
###    Lower analysis sub-clustered CAFs
### 1. General analysis ####
Func.Subclust.GeneralLower = function(Dim1, Res1, Dim2, Res2){
  ##############################
  # 0. read seurat obj
  SeuObj.CAF_0 = readRDS(paste0(DirInteg,"Objects/[SeuObj_ExtractedCAF][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                                QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2,".rds") )
  CommonTitle.Sub = paste0(NumOfSamples,"case Integration(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")\n",
                           QCInfo,"\n1stClust:PCA50vf2000_ClustDim",Dim1,"Res",formatC(Res1,digits=2,format="f"),
                           "\n2ndClust:PCA50vf2000_ClustDim",Dim2,"Res",formatC(Res2,digits=2,format="f"))
  NumOfSubclust = length(unique(SeuObj.CAF_0$seurat_clusters))
  CommonCaption.Sub = paste0("Total ",format(ncol(SeuObj.CAF_0), big.mark=",", scientific=F), " CAFs, ",NumOfSubclust, " clusters.")
  CommonCAFLab = paste0("CAF-", 0:(NumOfSubclust-1))
  names(CommonCAFLab) = as.character(0:(NumOfSubclust-1))
  set.seed(123)
  ##############################
  # 1. Dimplot: Fusion & By cluster
  DF.UMAP_0 = Embeddings(SeuObj.CAF_0, "umap") %>% as.data.frame() %>%
    rownames_to_column(var = "full_id") %>% 
    dplyr::mutate("cell_id"=str_remove(colnames(SeuObj.CAF_0), pattern="TX5K.*_"),
                  "seurat_clusters" = SeuObj.CAF_0$seurat_clusters,
                  "sample" = str_remove(SeuObj.CAF_0$orig.ident, pattern="TX5K_"))
  DF.UMAP_1 = dplyr::slice_sample(DF.UMAP_0, prop=1)
  DF.UMAPcenter = group_by(DF.UMAP_1, seurat_clusters) %>% 
    dplyr::summarize(umap_1 = median(umap_1), umap_2 = median(umap_2)) %>% 
    dplyr::mutate(label = CommonCAFLab[seurat_clusters])
  DimPlot_ScaleX = scale_x_continuous(limits = range(DF.UMAP_0$umap_1), breaks=c(-10,-5,0,5,10))
  DimPlot_ScaleY = scale_y_continuous(limits = range(DF.UMAP_0$umap_2), breaks=c(-10,-5,0,5,10))
  DimPlot = 
    ggplot(DF.UMAP_1, aes(x=umap_1, y=umap_2, color=seurat_clusters)) +
    geom_point(size=0.2, alpha=0.2) +
    #geom_text(data=DF.UMAPcenter, 
    #          aes(label=label),
    #          color="black") +
    DimPlot_ScaleX + DimPlot_ScaleY +
    coord_fixed(ratio=1) + 
    theme(legend.position="none",
          axis.text = element_text(face = "bold"),
          axis.title = element_text(face = "bold"),
          panel.background = element_blank(),
          axis.line = element_line(color="black")) # + facet_wrap(. ~ seurat_clusters)
  DimPlot.Split.list = list()
  ColPalette = hue_pal()(NumOfSubclust)
  for(i in 0:(NumOfSubclust-1)){
    Plot.Split = 
      ggplot(subset(DF.UMAP_1, subset=seurat_clusters==i), 
             aes(x=umap_1, y=umap_2)) +
      geom_point(data=DF.UMAP_1, color="gray90", size=0.2, alpha=0.2) +
      geom_point(color=ColPalette[i+1], size=0.2, alpha=0.2) +
      DimPlot_ScaleX + DimPlot_ScaleY +
      coord_fixed(ratio=1) + 
      theme(legend.position="none",
            axis.text = element_text(face = "bold"),
            axis.title = element_text(face = "bold"),
            panel.background = element_blank(),
            axis.line = element_line(color="black")) #+ facet_wrap(. ~ seurat_clusters)
    DimPlot.Split.list[[paste0("CAF_",i)]] = Plot.Split
  }
  ggsave(
    plot=ggarrange(DimPlot, 
                   wrap_plots(DimPlot.Split.list),
                   ncol=2) %>% 
         annotate_figure(top=text_grob(CommonTitle.Sub), bottom=CommonCaption.Sub),
    file=paste0(DirInteg,"[Figure][SubClust_DimPlot2][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,
                "_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"), 
    dpi=200, width=16, height=8, bg = "white" )
  ##############################
  # 2. FeaturePlot
  FeatPlot_RangeFeature = range(SeuObj.CAF_0$nFeature_RNA)
  FeatPlot_RangeCount = range(SeuObj.CAF_0$nCount_RNA)
  FeatPlot_RangeFeature = range(SeuObj.CAF_0$nFeature_RNA)
  FeatPlot_RangeCount = range(SeuObj.CAF_0$nCount_RNA)
  Ncol=ceiling(NumOfSubclust^0.5)
  Plot.GeneralLower2 = 
    ggarrange(FeaturePlot(SeuObj.CAF_0, reduction="umap", label=TRUE, split.by="seurat_clusters", raster=FALSE,
                          features="nFeature_RNA") + DimPlot_ScaleX + DimPlot_ScaleY + 
                plot_layout(ncol=Ncol, nrow=ceiling(NumOfSubclust/Ncol), guides="collect") & 
                scale_color_viridis(option="magma", limits=FeatPlot_RangeFeature) &
                coord_fixed(ratio=1) & theme(legend.position="right") & #DimPlot_Theme &
                patchwork::plot_annotation(title="nFeature_RNA"),
              FeaturePlot(SeuObj.CAF_0, reduction="umap", label=TRUE, split.by="seurat_clusters", raster=FALSE,
                          features="nCount_RNA") + DimPlot_ScaleX + DimPlot_ScaleY +
                plot_layout(ncol=Ncol, nrow=ceiling(NumOfSubclust/Ncol), guides="collect") & 
                scale_color_viridis(option="magma", limits=FeatPlot_RangeCount) &
                coord_fixed(ratio=1) & theme(legend.position="right") & #DimPlot_Theme & 
                patchwork::plot_annotation(title="nCount_RNA"),
              ncol=2, widths=c(1, 1)) +
    patchwork::plot_annotation(title=CommonTitle.Sub)
  ggsave(Plot.GeneralLower2,
         file=paste0(DirInteg,"[Figure][SubClust_DimPlot2_nFandC][",
                     NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                     "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"), 
         dpi=100, width=32, height=12, bg = "white" )
  ##############################
  # 3. Dot plot : Cell type markers & CAF markers
  DotPlotTheme = theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                       axis.text = element_text(face = "bold", size = 19),
                       legend.text = element_text(face = "bold", size = 15),
                       legend.title = element_text(face = "bold", size = 15),
                       axis.title = element_text(face = "bold", size = 17))
  DotPlotLgdLab = guides(size = guide_legend(title = "Percent\nExpressed"),
                         color = guide_colorbar(title = "Scaled\nAverage\nExpression",
                                                direction="vertical",
                                                frame.colour="black",
                                                ticks.colour="black"))
  DotPlotColScale = scale_color_gradient2(high="#f03b20",mid="gray90",low="#2c7fb8",midpoint=0)
  ggsave(ggarrange(
    DotPlot(SeuObj.CAF_0, cluster.idents=FALSE,
            dot.scale=10,
            features=Vec.MarkerGenes_CellType) &
      labs(y="Cluster", x=NULL) & DotPlotLgdLab & RotatedAxis() & DotPlotTheme & DotPlotColScale,
    DotPlot(SeuObj.CAF_0, cluster.idents=FALSE, 
            dot.scale = 10, 
            features=Vec.CAFmarkers) &
      labs(y="Cluster", x=NULL) & DotPlotLgdLab & RotatedAxis() & DotPlotTheme & DotPlotColScale,
    ncol=1, nrow=2) %>% 
      annotate_figure(top=text_grob(CommonTitle.Sub), bottom=CommonCaption.Sub),
    file=paste0(DirInteg,"[Figure][SubClust_DotPlot_CellTypeAndCAFmarkers][",
                NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"), 
    dpi=300, width=20, height=25, bg="white")
  ##############################
  # 4. Dot plot : Senescence & Proliferation markers
  ggsave(ggarrange(
    ggarrange(
        DotPlot(SeuObj.CAF_0, cluster.idents=FALSE, dot.scale=10,
                features=Vec.BuffaOrig) &
        labs(y="Cluster", x=NULL, subtitle="Buffa's hypoxia") & DotPlotLgdLab & RotatedAxis() & DotPlotTheme & DotPlotColScale,
        DotPlot(SeuObj.CAF_0, cluster.idents=FALSE, dot.scale=10,
                features=Vec.WinterOrig ) &
        labs(y="Cluster", x=NULL, subtitle="Winter's hypoxia") & DotPlotLgdLab & RotatedAxis() & DotPlotTheme & DotPlotColScale,
        ncol=2, nrow=1, align="hv"),
    DotPlot(SeuObj.CAF_0, cluster.idents=FALSE,
            dot.scale=10,
            features=unique(subset(DF.Markers, subset=GeneType=="Senescence"&
                                     Panel5k=="Included")$Gene) ) &
      labs(y="Cluster", x=NULL, subtitle="Senescence signatures (KEGG cellular senescence, SASP, HumanHallmark P53 pathway)") & DotPlotLgdLab & RotatedAxis() & DotPlotTheme & DotPlotColScale,
    DotPlot(SeuObj.CAF_0, cluster.idents=FALSE, 
            dot.scale = 10, 
            features=subset(DF.Markers, subset=GeneType=="CAFpaper"&
                              Classification1=="Proliferation"&
                              Panel5k=="Included")$Gene ) &
      labs(y="Cluster", x=NULL, subtitle="Proliferation genes") & DotPlotLgdLab & RotatedAxis() & DotPlotTheme & DotPlotColScale,
    ncol=1, nrow=3) %>% 
      annotate_figure(top=text_grob(CommonTitle.Sub)),
    file=paste0(DirInteg,"[Figure][SubClust_DotPlot_SenescenceAndProliferation][",
                NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"), 
    dpi=300, width=48, height=35, bg="white")
  ##############################
  # 5. Dot plot : Variable genes 
  # High SD genes in average expression
  DF.average = AverageExpression(SeuObj.CAF_0, group.by="seurat_clusters", assay="RNA", layer="data") %>% data.frame()
  DF.average_SD = mutate(DF.average, StaDev = apply(DF.average, 1, sd)) %>% arrange(-StaDev)
  Vec.VarioubleGenes = rownames(DF.average_SD[1:40, ])
  # hierarchical clust
  rho = stats::cor(t(DF.average[Vec.VarioubleGenes, ]), method = "spearman")   # spearmanの相関係数ρ(rho)
  dist = stats::as.dist(1 - rho)
  #これひつよう？ -> attributes(dist)$class == "dist"
  tree = hclust(dist, method = "ward.D")
  #plot(tree)
  Pdot_variable = DotPlot(SeuObj.CAF_0, dot.scale = 10, 
                          features = Vec.VarioubleGenes[tree$order],
                          cluster.idents = TRUE) &
    labs(x=NULL, y ="Cluster", caption=paste0(CommonTitle.Sub, "Clusters and genes are hierarchicaly clustered.") ) &
    guides(size = guide_legend(title = "Percent\nExpressed"),
           color = guide_colorbar(title = "Scaled\nAverage\nExpression")) &
    RotatedAxis() &
    theme(text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 15)) &
    scale_color_gradientn(colors = c("gray90","#fee5d9","#fcbba1","#ef3b2c","#99000d")) &
    coord_fixed(ratio=1) &
    ggtitle(paste0("Top variable 40 genes\n",CommonTitle.Sub))
  ggsave(Pdot_variable, 
         file=paste0(DirInteg,"[Figure][SubClust_DotPlot_VariableGenes][",
                     NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                     "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"),
         dpi = 300, width = 13, height = 10, bg = "white") 
  ##############################
  # 6. Sample proportion by custer
  CommonTheme_bar = theme(axis.text = element_text(face = "bold"),
                      axis.title = element_text(face = "bold"),
                      legend.text = element_text(face = "bold"),
                      legend.title = element_text(face = "bold"),
                      legend.key = element_blank(),
                      panel.background = element_blank(),
                      panel.grid = element_blank(),
                      axis.line = element_line(color = "gray20"))
  DF.CellTable_0 = table(SeuObj.CAF_0$seurat_clusters, SeuObj.CAF_0$orig.ident) %>% as.data.frame() %>% 
                      set_colnames(c("cluster","sample","cell")) %>% 
                      dplyr::mutate(cluster = paste0("CAF-",cluster))
  DF.CellsBySample_Labs_0 = DF.CellTable_0 %>% 
                           group_by(sample) %>% 
                           dplyr::mutate(total=sum(cell)) %>% 
                           ungroup() %>% rowwise() %>% 
                           dplyr::mutate(prop=cell/total) %>% 
                           dplyr::mutate(label=paste0(formatC(cell, big.mark=","),"\n",
                                                      "(",formatC(prop*100, format="f", digits=1),"%)")) %>% 
                           ungroup() %>% group_by(sample) %>% 
                           dplyr::mutate(cumsum=cumsum(prop)) %>% 
                           dplyr::mutate(coord_y=1-cumsum+prop/2)
  DF.CellsByClust_Labs_0 = dplyr::arrange(DF.CellTable_0, cluster) %>% group_by(cluster) %>% 
                           dplyr::mutate(total=sum(cell)) %>% 
                           ungroup() %>% rowwise() %>% 
    dplyr::mutate(prop=cell/total) %>% 
    dplyr::mutate(label=paste0(formatC(cell, big.mark=","),"\n",
                               "(",formatC(prop*100, format="f", digits=1),"%)")) %>% 
    ungroup() %>% group_by(cluster) %>% 
    dplyr::mutate(cumsum=cumsum(prop)) %>% 
    dplyr::mutate(coord_y=1-cumsum+prop/2)
  BarPlot = ggarrange(
    ggplot(DF.CellTable_0, aes(x=sample, y=cell, fill=cluster)) +
      geom_bar(stat="identity", position="fill", color="gray50") +
      geom_text(data=DF.CellsBySample_Labs_0, 
                aes(y=coord_y, label=case_when(prop>0.05 ~label,
                                               prop>0.02 ~str_replace(label, pattern="\n", replace="  "),
                                               TRUE ~ ""))) +
      geom_text(data=dplyr::slice(group_by(DF.CellsByClust_Labs_0, sample), 1), 
                aes(label=formatC(total, big.mark=",")), y=1.03) +
      labs(y="Proportion", x="Sample", fill="Sub-cluster") +
      scale_y_continuous(labels=scales::percent_format(),
                         expand=expansion(mult=c(0, 0.05))) +
      CommonTheme_bar,
    ggplot(DF.CellTable_0, aes(x=cluster, y=cell, fill=sample)) +
      geom_bar(stat="identity", position="fill", color="gray50") +
      geom_text(data=DF.CellsByClust_Labs_0, 
                aes(y=coord_y, label=case_when(prop>0.10 ~label,
                                               prop>0.02 ~str_replace(label, pattern="\n", replace="  "),
                                               TRUE ~ ""))) +
      geom_text(data=dplyr::slice(group_by(DF.CellsByClust_Labs_0, cluster), 1), 
                aes(label=formatC(total, big.mark=",")), y=1.02) +
      labs(y="Proportion", x="Sub-cluster", fill="Sample") +
      scale_y_continuous(labels=scales::percent_format(),
                         expand=expansion(mult=c(0, 0.05))) +
      CommonTheme_bar,
    ncol=1, nrow=2, align="hv") %>% 
      annotate_figure(top=text_grob(CommonTitle.Sub), bottom=CommonCaption.Sub)
  ggsave(BarPlot, 
         file=paste0(DirInteg,"[Figure][SubClust_BarPlot_CellProportion][",
                     NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                     "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"),
         dpi = 300, width = 13, height = 13, bg = "white") 
  }
#Func.Subclust.GeneralLower(Dim1=20, Res1=1.0, Dim2=15, Res2=1.0) ####   <- BAD
#Func.Subclust.GeneralLower(Dim1=20, Res1=1.0, Dim2=15, Res2=0.5) ####   <- BAD
Func.Subclust.GeneralLower(Dim1=20, Res1=1.0, Dim2=20, Res2=0.5)
Func.Subclust.GeneralLower(Dim1=20, Res1=1.0, Dim2=30, Res2=0.5)
Func.Subclust.GeneralLower(Dim1=20, Res1=1.0, Dim2=40, Res2=0.5)
#Func.Subclust.GeneralLower(Dim1=20, Res1=1.0, Dim2=50, Res2=0.5) ####   <- BAD
Func.Subclust.GeneralLower(Dim1=25, Res1=1.0, Dim2=15, Res2=0.5)
Func.Subclust.GeneralLower(Dim1=25, Res1=1.0, Dim2=25, Res2=0.5)
Func.Subclust.GeneralLower(Dim1=25, Res1=1.0, Dim2=35, Res2=0.5)
beep(expr=NULL, sound=2)


### 2. Module score (insert) ####
Func.ModuleScore = function(Dim1, Res1, Dim2, Res2, MidPoint.Proliferation, MidPoint.BuffaOrig, MidPoint.WinterOrig, MidPoint.myCAF_score, MidPoint.iCAF_score, MidPoint.apCAF_score){
  #-------------------------------------------------------------------
  # 1.Scatter plot: Correlations between RNA quality and proliferation
  SeuObj.CAF_0 = readRDS(paste0(DirInteg,"Objects/[SeuObj_ExtractedCAF][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                                QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2,".rds") )
  NumObSubclust = length(unique(SeuObj.CAF_0$seurat_clusters))
  CommonTitle.Sub = paste0(NumOfSamples,"case Integration(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")\n",
                           QCInfo,"\n1stClust:PCA50vf2000_ClustDim",Dim1,"Res",formatC(Res1,digits=2,format="f"),
                           "\n2ndClust:PCA50vf2000_ClustDim",Dim2,"Res",formatC(Res2,digits=2,format="f"))
  CommonCaption.Sub = paste0("Total ",format(ncol(SeuObj.CAF_0), big.mark=",", scientific=F), " CAFs, ",NumObSubclust, " clusters.")
  CommonTheme.Vln = theme(legend.position="none",
                          axis.text = element_text(face = "bold"),
                          axis.title = element_text(face = "bold"),
                          legend.text = element_text(face = "bold"),
                          legend.title = element_text(face = "bold"),
                          panel.background = element_blank(),
                          panel.grid = element_blank(),
                          axis.line = element_line(color = "black"))
  DF.ModScore_0 = read.csv(paste0(DirInteg,"/[DataTable_ExtractedCAF_ModuleScore][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                                  QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".csv"), row.names=1 )
  DF.ModScore_1 = DF.ModScore_0[, 1:(ncol(DF.ModScore_0)/2)]
  MT.Correlation_0 = cor(DF.ModScore_1)
  if ( identical(colnames(SeuObj.CAF_0), rownames(DF.ModScore_0)) ){ 
    SeuObj.CAF_0@meta.data = cbind(SeuObj.CAF_0@meta.data, DF.ModScore_0) 
    DF.ModScore_1 = data.frame(SeuObj.CAF_0@reductions[["umap"]]@cell.embeddings[,c("umap_1","umap_2")],
                               "nCount_RNA"=SeuObj.CAF_0$nCount_RNA,
                               "nFeature_RNA"=SeuObj.CAF_0$nFeature_RNA,
                               "subclust"=SeuObj.CAF_0$seurat_clusters,
                               DF.ModScore_0)
  }
  Func.Scatter = function(TargetGene, ProlifCAFclust){
    DF.ModScore_2 = cbind(DF.ModScore_1,
                          "TargetGene"=GetAssayData(JoinLayers(SeuObj.CAF_0), assay="RNA", layer="count")[TargetGene, ])
    DF.ModScore_3 = subset(DF.ModScore_2, subset=subclust==ProlifCAFclust)
    ggarrange(
      ggplot(DF.ModScore_3, aes(x=Proliferation, y=nFeature_RNA, color=TargetGene)) + 
        geom_point() + labs(color=TargetGene) + scale_color_viridis(option="turbo") +
        geom_text(label=paste0("CAF-",ProlifCAFclust), x=0, y=900) +
        scale_y_continuous(limits=nFeatRNA, breaks=seq(nFeatRNA[1], nFeatRNA[2], by=100)),
      ggplot(DF.ModScore_3, aes(x=Proliferation, y=nCount_RNA, color=TargetGene)) + 
        geom_point() + labs(color=TargetGene) + scale_color_viridis(option="turbo") +
        geom_text(label=paste0("CAF-",ProlifCAFclust), x=0, y=1800) +
        scale_y_continuous(limits=nCountRNA, breaks=seq(nCountRNA[1], nCountRNA[2], by=100)),
      ncol=2, nrow=1)
    }
  Plot.Scatter = ggarrange(
      Func.Scatter("TOP2A", ProlifCAFclust=8),
      Func.Scatter("BIRC5", ProlifCAFclust=8),
      Func.Scatter("MKI67", ProlifCAFclust=8),
      Func.Scatter("PGK1", ProlifCAFclust=8),
      ncol=1) %>% 
    annotate_figure(top=CommonTitle.Sub)
  ggsave(Plot.Scatter,
         file=paste0(DirInteg,"[Figure][SubClust_ProlifCAFandQuality][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                     "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,
                     "_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"), 
         dpi=200, width=12, height=16, bg = "white" )
  #-------------------------------------------------------------------
  # 2. Violin plot : Module scores of sub-clusters
  Func.Violin = function(ScoreName){
    DF = DF.ModScore_1[,c("subclust",ScoreName)] %>% set_colnames(c("subclust","Target"))
    NamedVec.colors = hue_pal()(NumObSubclust)
    names(NamedVec.colors) = 0:(NumObSubclust-1)
    DF.Average = summarise(group_by(DF, subclust),
                           Mean = mean(Target)) %>% 
                 dplyr::arrange(Mean)
    DF$subclust = factor(DF$subclust, levels=as.character(DF.Average$subclust))
    ggplot(DF, aes(x=subclust, y=Target, fill=subclust)) + 
      geom_violin() + 
      geom_boxplot(width=0.2, fill="white", outliers=F) +
      labs(x="CAF sub-cluster", y=ScoreName) +
      scale_fill_manual(values=NamedVec.colors) +
      CommonTheme.Vln
  }
  VlnPlot.MSbySubclust = 
    ggarrange(
      Func.Violin(ScoreName="Proliferation"),
      Func.Violin(ScoreName="BuffaOrig"),
      Func.Violin(ScoreName="WinterOrig"),
      Func.Violin(ScoreName="myCAF_score"),
      Func.Violin(ScoreName="iCAF_score"), 
      Func.Violin(ScoreName="apCAF_score"), 
    ncol=1, nrow=6, align="hv") %>% 
    annotate_figure(top=CommonTitle.Sub, bottom=paste0(CommonCaption.Sub, "Subclusters are sorted by mean value.") )
  ggsave(VlnPlot.MSbySubclust,
         file=paste0(DirInteg,"[Figure][SubClust_ModuleScoresBySubclust][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                     "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"), 
         dpi=300, width=5, height=20, bg = "white" )
  
  #-------------------------------------------------------------------
  # 3. Feature plot 
  CommonTheme.Feat = theme(axis.text = element_text(face = "bold"),
                           axis.title = element_text(face = "bold"),
                           legend.text = element_text(face = "bold"),
                           legend.title = element_text(face = "bold"),
                           panel.background = element_blank(),
                           panel.grid = element_blank(),
                           axis.line = element_line(color = "black"))
  Vec.ColorScale = c("#253494", "#2c7fb8","#41b6c4", "#a1dab4", "#ffffcc", "#fecc5c", "#fd8d3c","#f03b20", "#bd0026")
  Vec.Range.Proliferation = c(min(DF.ModScore_1$Proliferation), max(DF.ModScore_1$Proliferation))
  Vec.Range.BuffaOrig = c(min(DF.ModScore_1$BuffaOrig), max(DF.ModScore_1$BuffaOrig))
  Vec.Range.WinterOrig = c(min(DF.ModScore_1$WinterOrig), max(DF.ModScore_1$WinterOrig))
  Vec.Range.myCAF_score = c(min(DF.ModScore_1$myCAF_score), max(DF.ModScore_1$myCAF_score))
  Vec.Range.iCAF_score = c(min(DF.ModScore_1$iCAF_score), max(DF.ModScore_1$iCAF_score))
  Vec.Range.apCAF_score = c(min(DF.ModScore_1$apCAF_score), max(DF.ModScore_1$apCAF_score))
  Probs.Proliferation = c(quantile(c(0, MidPoint.Proliferation),probs=c(0, 0.25, 0.5, 0.75)),
                          quantile(c(MidPoint.Proliferation, 1),probs=c(0, 0.25, 0.5, 0.75, 1.0)))
  Probs.BuffaOrig = c(quantile(c(0, MidPoint.BuffaOrig),probs=c(0, 0.25, 0.5, 0.75)),
                          quantile(c(MidPoint.BuffaOrig, 1),probs=c(0, 0.25, 0.5, 0.75, 1.0)))
  Probs.WinterOrig = c(quantile(c(0, MidPoint.WinterOrig),probs=c(0, 0.25, 0.5, 0.75)),
                    quantile(c(MidPoint.WinterOrig, 1),probs=c(0, 0.25, 0.5, 0.75, 1.0)))
  Probs.myCAF_score = c(quantile(c(0, MidPoint.myCAF_score),probs=c(0, 0.25, 0.5, 0.75)),
                          quantile(c(MidPoint.myCAF_score, 1),probs=c(0, 0.25, 0.5, 0.75, 1.0)))
  Probs.iCAF_score = c(quantile(c(0, MidPoint.iCAF_score),probs=c(0, 0.25, 0.5, 0.75)),
                          quantile(c(MidPoint.iCAF_score, 1),probs=c(0, 0.25, 0.5, 0.75, 1.0)))
  Probs.apCAF_score = c(quantile(c(0, MidPoint.apCAF_score),probs=c(0, 0.25, 0.5, 0.75)),
                          quantile(c(MidPoint.apCAF_score, 1),probs=c(0, 0.25, 0.5, 0.75, 1.0)))
  FeatPlot.ModScore = 
    ggarrange(
      ggplot(dplyr::arrange(DF.ModScore_1, Proliferation), aes(x=umap_1, y=umap_2, color=Proliferation)) + 
        geom_point(size=0.1) + coord_fixed(ratio=1) + CommonTheme.Feat + 
        scale_color_gradientn(colors=Vec.ColorScale, 
                              values=scales::rescale(quantile(Vec.Range.Proliferation, probs=Probs.Proliferation)),
                              limits=Vec.Range.Proliferation,
                              guide=guide_colorbar(frame.colour="black", ticks.colour="black")),
      ggplot(dplyr::arrange(DF.ModScore_1, BuffaOrig), aes(x=umap_1, y=umap_2, color=BuffaOrig)) + 
        geom_point(size=0.1) + coord_fixed(ratio=1) + CommonTheme.Feat + 
        scale_color_gradientn(colors=Vec.ColorScale, 
                                values=scales::rescale(quantile(Vec.Range.BuffaOrig, probs=Probs.BuffaOrig)),
                              limits=Vec.Range.BuffaOrig,
                              guide=guide_colorbar(frame.colour="black", ticks.colour="black")),
      ggplot(dplyr::arrange(DF.ModScore_1, WinterOrig), aes(x=umap_1, y=umap_2, color=WinterOrig)) + 
        geom_point(size=0.1) + coord_fixed(ratio=1) + CommonTheme.Feat + 
        scale_color_gradientn(colors=Vec.ColorScale, 
                              values=scales::rescale(quantile(Vec.Range.WinterOrig, probs=Probs.WinterOrig)),
                              limits=Vec.Range.WinterOrig,
                              guide=guide_colorbar(frame.colour="black", ticks.colour="black")),
      ggplot(dplyr::arrange(DF.ModScore_1, myCAF_score), aes(x=umap_1, y=umap_2, color=myCAF_score)) + 
        geom_point(size=0.1) + coord_fixed(ratio=1) + CommonTheme.Feat + labs(color="myCAF") +
        scale_color_gradientn(colors=Vec.ColorScale, 
                              values=scales::rescale(quantile(Vec.Range.myCAF_score, probs=Probs.myCAF_score)),
                              limits=Vec.Range.myCAF_score,
                              guide=guide_colorbar(frame.colour="black", ticks.colour="black")),
      ggplot(dplyr::arrange(DF.ModScore_1, iCAF_score), aes(x=umap_1, y=umap_2, color=iCAF_score)) + 
        geom_point(size=0.1) + coord_fixed(ratio=1) + CommonTheme.Feat + labs(color="iCAF") +
        scale_color_gradientn(colors=Vec.ColorScale,
                              values=scales::rescale(quantile(Vec.Range.iCAF_score, probs=Probs.iCAF_score)),
                              limits=Vec.Range.iCAF_score,
                              guide=guide_colorbar(frame.colour="black", ticks.colour="black")),
      ggplot(dplyr::arrange(DF.ModScore_1, apCAF_score), aes(x=umap_1, y=umap_2, color=apCAF_score)) + 
        geom_point(size=0.1) + coord_fixed(ratio=1) + CommonTheme.Feat + labs(color="apCAF") +
        scale_color_gradientn(colors=Vec.ColorScale,
                              values=scales::rescale(quantile(Vec.Range.apCAF_score, probs=Probs.apCAF_score)),
                              limits=Vec.Range.apCAF_score,
                              guide=guide_colorbar(frame.colour="black", ticks.colour="black")),
    ncol=2, nrow=3, align="hv") %>% 
    annotate_figure(top=CommonTitle.Sub, 
                    bottom=paste0(CommonCaption.Sub, "\nMidpoints of scores are shifted; Proliferation:",MidPoint.Proliferation,
                                  ", BuffaOriginal:",MidPoint.BuffaOrig,", WinterOriginal:",MidPoint.WinterOrig,
                                  ", myCAF:",MidPoint.myCAF_score,", iCAF:",MidPoint.iCAF_score,", apCAF:",MidPoint.apCAF_score) )
  ggsave(plot=FeatPlot.ModScore,
         file=paste0(DirInteg,"[Figure][SubClust_ModuleScores_FeaturePlot][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                     "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"), 
         dpi=300, width=10, height=8, bg = "white" )
  FeatPlot.ModScore.Byclust = ggarrange(
    ggplot(dplyr::arrange(DF.ModScore_1, Proliferation), 
           aes(x=umap_1, y=umap_2, color=Proliferation)) + geom_point(size=0.1) + 
        scale_color_gradientn(colors=Vec.ColorScale, 
                              values=scales::rescale(quantile(Vec.Range.Proliferation, probs=Probs.Proliferation)),
                              limits=Vec.Range.Proliferation,
                              guide=guide_colorbar(frame.colour="black", ticks.colour="black") ) + 
        CommonTheme.Feat + coord_fixed(ratio=1) + facet_wrap(subclust ~ ., nrow=1),
    ggplot(dplyr::arrange(DF.ModScore_1, BuffaOrig), 
           aes(x=umap_1, y=umap_2, color=BuffaOrig)) + geom_point(size=0.1, alpha=0.3) + 
        scale_color_gradientn(colors=Vec.ColorScale, 
                              values=scales::rescale(quantile(Vec.Range.BuffaOrig, probs=Probs.BuffaOrig)),
                              limits=Vec.Range.BuffaOrig,
                              guide=guide_colorbar(frame.colour="black", ticks.colour="black") ) + 
        CommonTheme.Feat + coord_fixed(ratio=1) + facet_wrap(subclust ~ ., nrow=1),
    ggplot(dplyr::arrange(DF.ModScore_1, WinterOrig), 
           aes(x=umap_1, y=umap_2, color=WinterOrig)) + geom_point(size=0.1, alpha=0.3) + 
      scale_color_gradientn(colors=Vec.ColorScale, 
                            values=scales::rescale(quantile(Vec.Range.WinterOrig, probs=Probs.WinterOrig)),
                            limits=Vec.Range.WinterOrig,
                            guide=guide_colorbar(frame.colour="black", ticks.colour="black") ) + 
      CommonTheme.Feat + coord_fixed(ratio=1) + facet_wrap(subclust ~ ., nrow=1),
    ggplot(dplyr::arrange(DF.ModScore_1, myCAF_score), 
           aes(x=umap_1, y=umap_2, color=myCAF_score)) + geom_point(size=0.1, alpha=0.3) + labs(color="myCAF") +
        scale_color_gradientn(colors=Vec.ColorScale,
                              values=scales::rescale(quantile(Vec.Range.myCAF_score, probs=Probs.myCAF_score)),
                              limits=Vec.Range.myCAF_score,
                              guide=guide_colorbar(frame.colour="black", ticks.colour="black") ) + 
        CommonTheme.Feat + coord_fixed(ratio=1) + facet_wrap(subclust ~ ., nrow=1),
    ggplot(dplyr::arrange(DF.ModScore_1, iCAF_score), 
           aes(x=umap_1, y=umap_2, color=iCAF_score)) + geom_point(size=0.1, alpha=0.3) + labs(color="iCAF") +
        scale_color_gradientn(colors=Vec.ColorScale,
                              values=scales::rescale(quantile(Vec.Range.iCAF_score, probs=Probs.iCAF_score)),
                              limits=Vec.Range.iCAF_score,
                              guide=guide_colorbar(frame.colour="black", ticks.colour="black") ) + 
        CommonTheme.Feat + coord_fixed(ratio=1) + facet_wrap(subclust ~ ., nrow=1),
    ggplot(dplyr::arrange(DF.ModScore_1, apCAF_score), 
           aes(x=umap_1, y=umap_2, color=apCAF_score)) + geom_point(size=0.1, alpha=0.3) + labs(color="apCAF") +
        scale_color_gradientn(colors=Vec.ColorScale, 
                              values=scales::rescale(quantile(Vec.Range.apCAF_score, probs=Probs.apCAF_score)),
                              limits=Vec.Range.apCAF_score,
                              guide=guide_colorbar(frame.colour="black", ticks.colour="black") ) + 
        CommonTheme.Feat + coord_fixed(ratio=1) + facet_wrap(subclust ~ ., nrow=1),
    ncol=1, align="hv") %>% 
    annotate_figure(top=CommonTitle.Sub, 
                    bottom=paste0(CommonCaption.Sub, "\nMidpoints of scores are shifted; Proliferation:",MidPoint.Proliferation,
                                  ", BuffaOriginal:",MidPoint.BuffaOrig,", WinterOriginal:",MidPoint.WinterOrig,
                                  ", myCAF:",MidPoint.myCAF_score,", iCAF:",MidPoint.iCAF_score,", apCAF:",MidPoint.apCAF_score) )
  ggsave(plot=FeatPlot.ModScore.Byclust,
         file=paste0(DirInteg,"[Figure][SubClust_ModuleScores_FeaturePlot_ByClust][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                     "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"), 
         dpi=300, width=16, height=12, bg = "white" )
}
Func.ModuleScore(Dim1=20, Res1=1.0, Dim2=20, Res2=0.5, MidPoint.Proliferation=0.65, MidPoint.BuffaOrig=0.65, MidPoint.WinterOrig=0.65, MidPoint.myCAF_score=0.60, MidPoint.iCAF_score=0.5, MidPoint.apCAF_score=0.75)
Func.ModuleScore(Dim1=20, Res1=1.0, Dim2=30, Res2=0.5, MidPoint.Proliferation=0.65, MidPoint.BuffaOrig=0.65, MidPoint.WinterOrig=0.65, MidPoint.myCAF_score=0.60, MidPoint.iCAF_score=0.5, MidPoint.apCAF_score=0.75)
Func.ModuleScore(Dim1=20, Res1=1.0, Dim2=40, Res2=0.5, MidPoint.Proliferation=0.65, MidPoint.BuffaOrig=0.65, MidPoint.WinterOrig=0.65, MidPoint.myCAF_score=0.60, MidPoint.iCAF_score=0.5, MidPoint.apCAF_score=0.75)


### 3. ssGSVA score (insert) ####
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
MidPoint.iCAF=0.6
MidPoint.apCAF=0.999
Func.gsvaScore.Sub = function(Dim1, Res1, Dim2, Res2, method,
                              MidPoint.Proliferation, MidPoint.KEGG_CellCycle,
                              MidPoint.E2F, MidPoint.G2M, MidPoint.MITOTIC_SPINDLE,
                              MidPoint.MYC_TARGETS_V1, MidPoint.MYC_TARGETS_V2,
                              MidPoint.BuffaOrig, MidPoint.WinterOrig, MidPoint.HYPOXIA, 
                              MidPoint.myCAF, MidPoint.iCAF, MidPoint.apCAF){
  #-------------------------------------------------------------------
  # 0. Read data
  SeuObj.CAF_0 = readRDS(paste0(DirInteg,"Objects/[SeuObj_ExtractedCAF][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                                QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2,".rds") )
  NumObSubclust = length(unique(SeuObj.CAF_0$seurat_clusters))
  CommonTitle.Sub = paste0(NumOfSamples,"case Integration(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")\n",
                           QCInfo,"\n1stClust:PCA50vf2000_ClustDim",Dim1,"Res",formatC(Res1,digits=2,format="f"),
                           "\n2ndClust:PCA50vf2000_ClustDim",Dim2,"Res",formatC(Res2,digits=2,format="f"))
  CommonCaption.Sub = paste0("Total ",format(ncol(SeuObj.CAF_0), big.mark=",", scientific=F), " CAFs, ",NumObSubclust, " clusters.")
  DF.gsvaScore_0 = read.csv(file=paste0(DirInteg,"/[DataTable_ExtractedCAF_ssGSVA][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                                       QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".csv"))
  rownames(DF.gsvaScore_0) = paste(DF.gsvaScore_0$sample,
                                   DF.gsvaScore_0$cell_id, 
                                   sep="_")
  DF.gsvaScore_0 = dplyr::mutate(DF.gsvaScore_0,
                                 full_id = rownames(DF.gsvaScore_0))
  #-------------------------------------------------------------------
  # 1.Heatmap of all pathways
  DF.Meta_0 = data.frame("full_id"=colnames(SeuObj.CAF_0),
                         "seurat_clusters"=paste0("CAF-", SeuObj.CAF_0$seurat_clusters),
                         SeuObj.CAF_0@reductions[["umap"]]@cell.embeddings)
  if( identical(DF.Meta_0$full_id, DF.gsvaScore_0$full_id) ){
    DF.Meta_1 = inner_join(DF.Meta_0, DF.gsvaScore_0, by="full_id")
  }
  DF.ScoreMean_0 = DF.Meta_1 %>% 
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
  if(method=="SpearmanAverage"){
    d_pathway = stats::as.dist(1 - cor(DF.Zscore_0, method="spearman") )
    class(d_pathway) = "dist"
    Hclust_pathway = stats::hclust(d_pathway, method="average")
    d_clusters = stats::as.dist(1 - cor(t(DF.Zscore_0), method="spearman") )
    class(d_clusters) = "dist"
    Hclust_clusters = stats::hclust(d_clusters, method="average")
    MethodLabel = "Hierarchical clustering was performed using average linkage based on Spearman correlation distances."
  } else if (method=="EuclideanWardd2"){               # ward.D2なら距離はEuclideanで！！
    d_pathway = dist(t(DF.Zscore_0), method="euclidean")
    class(d_pathway) = "dist"
    Hclust_pathway <- hclust(d_pathway, method = "ward.D2")
    d_clusters = dist(DF.Zscore_0, method="euclidean")
    class(d_clusters) = "dist"
    Hclust_clusters <- hclust(d_clusters, method = "ward.D2")
    MethodLabel = "Hierarchical clustering was performed using Ward’s method (ward.D2) with Euclidean distances."
  }
  Vec.PathwayOrder = colnames(DF.Zscore_0)[Hclust_pathway$order]
  DF.Zscore_2$pathway = factor(DF.Zscore_2$pathway,
                               levels=Vec.PathwayOrder)
  Vec.ClustOrder = rownames(DF.Zscore_0)[Hclust_clusters$order]
  DF.Zscore_2$seurat_clusters = factor(DF.Zscore_2$seurat_clusters,
                                       levels=Vec.ClustOrder) 
  # dendrogram
  library(ggdendro)
  dend_pathway = as.dendrogram(Hclust_pathway)
  dend_data = dendro_data(dend_pathway)
  p_dend = 
    ggplot(dend_data$segments) +
    geom_segment(aes(x=y, y=x, xend=yend, yend=xend)) +
    scale_y_continuous(limits = c(0.5, length(Hclust_pathway$labels) + 0.5),
                       breaks = seq_along(Hclust_pathway$labels),
                       labels = Hclust_pathway$labels,
                       expand = expansion(mult=c(0,0))) +
    scale_x_reverse() +
    theme_void()
  p_heat = 
  ggplot(DF.Zscore_2, 
         aes(x=seurat_clusters, y=pathway, fill=zscore)) +
    geom_tile(color="gray50") +
    labs(subtitle=CommonTitle.Sub,
         x=NULL, y=NULL, fill="Scaled mean ssGSEA ES") +
    scale_x_discrete(expand=expansion(mult=c(0,0))) +
    scale_y_discrete(expand=expansion(mult=c(0,0)),
                     position="right",
                     label=c("BuffaOrig"="Buffa_original",
                             "WinterOrig"="Winter_original")) +
    scale_fill_gradientn(limits=c(max(abs(range(DF.Zscore_2$zscore))),
                                  -max(abs(range(DF.Zscore_2$zscore)))),
                         colors=c("#4dac26","#b8e186","#f7f7f7","#f1b6da","#d01c8b"),
                         guide=guide_colorbar(#direction="horizontal",
                                              title.position="top",
                                              title.hjust=0.5,
                                              frame.colour="black",
                                              ticks.colour="black")) +
    theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
          legend.position="bottom",
          axis.text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust=1),
          axis.title = element_text(face = "bold"),
          legend.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold"),
          panel.border = element_rect(fill = NA, color = "black"))
  ggsave(plot = p_dend + p_heat + 
                plot_layout(widths = c(1.2, 4)) +
                plot_annotation(caption = MethodLabel),
         file=paste0(DirInteg,"[Figure][SubClust_ssGSVAscores_Heatmap][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                     "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,
                     "_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"), 
         dpi=500, width=6, height=12, bg = "white" )
  #-------------------------------------------------------------------
  # 2.Scatter plot: Correlations between RNA quality and proliferation
  CommonTheme.Vln = theme(legend.position="none",
                          axis.text = element_text(face = "bold"),
                          axis.title = element_text(face = "bold"),
                          legend.text = element_text(face = "bold"),
                          legend.title = element_text(face = "bold"),
                          panel.background = element_blank(),
                          panel.grid = element_blank(),
                          axis.line = element_line(color = "black"))
  if ( identical(colnames(SeuObj.CAF_0), rownames(DF.gsvaScore_0)) ){ 
    SeuObj.CAF_0@meta.data = cbind(SeuObj.CAF_0@meta.data, DF.gsvaScore_0) 
    DF.gsvaScore_2 = data.frame(SeuObj.CAF_0@reductions[["umap"]]@cell.embeddings[,c("umap_1","umap_2")],
                               "nCount_RNA"=SeuObj.CAF_0$nCount_RNA,
                               "nFeature_RNA"=SeuObj.CAF_0$nFeature_RNA,
                               "subclust"=SeuObj.CAF_0$seurat_clusters,
                               DF.gsvaScore_0)
  }
  Func.Scatter = function(TargetGene, ProlifCAFclust){
    DF.gsvaScore_3 = cbind(DF.gsvaScore_2,
                          "TargetGene"=GetAssayData(JoinLayers(SeuObj.CAF_0), assay="RNA", layer="count")[TargetGene, ])
    DF.gsvaScore_4 = subset(DF.gsvaScore_3, subset=subclust==ProlifCAFclust)
    ggarrange(
      ggplot(DF.gsvaScore_4, aes(x=Proliferation_ssGSVA, y=nFeature_RNA, color=TargetGene)) + 
        geom_point() + labs(color=TargetGene) + scale_color_viridis(option="turbo") +
        geom_text(label=paste0("CAF-",ProlifCAFclust), x=0, y=900) +
        scale_y_continuous(limits=nFeatRNA, breaks=seq(nFeatRNA[1], nFeatRNA[2], by=100)),
      ggplot(DF.gsvaScore_4, aes(x=Proliferation_ssGSVA, y=nCount_RNA, color=TargetGene)) + 
        geom_point() + labs(color=TargetGene) + scale_color_viridis(option="turbo") +
        geom_text(label=paste0("CAF-",ProlifCAFclust), x=0, y=1800) +
        scale_y_continuous(limits=nCountRNA, breaks=seq(nCountRNA[1], nCountRNA[2], by=100)),
      ncol=2, nrow=1)
  }
  Plot.Scatter = ggarrange(
    Func.Scatter("TOP2A", ProlifCAFclust=8),
    Func.Scatter("BIRC5", ProlifCAFclust=8),
    Func.Scatter("MKI67", ProlifCAFclust=8),
    Func.Scatter("PGK1", ProlifCAFclust=8),
    ncol=1) %>% 
    annotate_figure(top=CommonTitle.Sub)
  ggsave(Plot.Scatter,
         file=paste0(DirInteg,"tiral[Figure][SubClust_ProlifCAFandQuality(ssGSVA)][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                     "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,
                     "_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"), 
         dpi=200, width=12, height=16, bg = "white" )
  #-------------------------------------------------------------------
  # 3. Violin plot : Module scores of sub-clusters
  Func.Violin = function(ScoreName){
    DF = 
      DF.gsvaScore_2[,c("subclust", ScoreName)] %>% 
      set_colnames(c("subclust", "Target"))
    NamedVec.colors = hue_pal()(NumObSubclust)
    names(NamedVec.colors) = 0:(NumObSubclust-1)
    DF.Average = summarise(group_by(DF, subclust),
                           Mean = mean(Target)) %>% 
      dplyr::arrange(Mean)
    DF$subclust = factor(DF$subclust, levels=as.character(DF.Average$subclust))
    ggplot(DF, aes(x=subclust, y=Target, fill=subclust)) + 
      geom_violin() + 
      geom_boxplot(width=0.10, fill="white", outliers=F) +
      labs(x="CAF sub-cluster", 
           y=str_remove(ScoreName, pattern="_ssGSVA") %>% 
             str_replace(pattern="HALLMARK_", replacement="HALLMARK\n") %>% 
             str_replace(pattern="KEGG_", replacement="KEGG\n") %>% 
             str_replace(pattern="REACTOME_", replacement="REACTOME\n") %>% 
             str_replace(pattern="Orig", replacement="_original") %>% 
             str_replace(pattern="MYC_TARGETS_V", replacement="MYC_v") %>% 
             str_replace(pattern="INTERFERON", replacement="IFN")) +
      scale_fill_manual(values=NamedVec.colors) +
      CommonTheme.Vln
  }
  VlnPlot.GSVAbySubclust = 
    ggarrange(
      Func.Violin(ScoreName=Vec.SignatureNames[1]),Func.Violin(ScoreName=Vec.SignatureNames[9]),
      Func.Violin(ScoreName=Vec.SignatureNames[2]),Func.Violin(ScoreName=Vec.SignatureNames[10]),
      Func.Violin(ScoreName=Vec.SignatureNames[3]),Func.Violin(ScoreName=Vec.SignatureNames[11]),
      Func.Violin(ScoreName=Vec.SignatureNames[4]),Func.Violin(ScoreName=Vec.SignatureNames[12]),
      Func.Violin(ScoreName=Vec.SignatureNames[5]),Func.Violin(ScoreName=Vec.SignatureNames[13]),
      Func.Violin(ScoreName=Vec.SignatureNames[6]),Func.Violin(ScoreName=Vec.SignatureNames[14]),
      Func.Violin(ScoreName=Vec.SignatureNames[7]),Func.Violin(ScoreName=Vec.SignatureNames[15]),
      Func.Violin(ScoreName=Vec.SignatureNames[8]),Func.Violin(ScoreName=Vec.SignatureNames[16]),
      Func.Violin(ScoreName=Vec.SignatureNames[17]),
      ncol=2, nrow=9, align="hv") %>% 
    annotate_figure(top=CommonTitle.Sub, bottom=paste0(CommonCaption.Sub, "Subclusters are sorted by mean value.") )
  ggsave(VlnPlot.GSVAbySubclust,
         file=paste0(DirInteg,"[Figure][SubClust_ssGSVAscores_Violin][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                     "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"), 
         dpi=300, width=8, height=21, bg = "white" )
  #-------------------------------------------------------------------
  # 4. Feature plot 
  CommonTheme.Feat = theme(axis.text = element_text(face = "bold"),
                           axis.title = element_text(face = "bold"),
                           legend.text = element_text(face = "bold"),
                           legend.title = element_text(face = "bold"),
                           panel.background = element_blank(),
                           panel.grid = element_blank(),
                           axis.line = element_line(color = "black"))
  Vec.ColorScale = c("#253494", "#2c7fb8","#41b6c4", "#a1dab4", "#ffffcc", "#fecc5c", "#fd8d3c","#f03b20", "#bd0026")
  Func.FeatPlotGSVA = function(ColumnName, MidPoint){
    DF.Score_0 = DF.gsvaScore_2[ , c("umap_1","umap_2",paste0(ColumnName,"_ssGSVA"))] %>% 
                 set_colnames(c("umap_1","umap_2","TargetSignature")) %>% 
                 dplyr::arrange(TargetSignature)
    Plot = 
      ggplot(DF.Score_0, 
           aes(x=umap_1, y=umap_2, color=TargetSignature)) + 
      geom_point(size=0.1) + 
      labs(color= ColumnName %>% 
                  str_replace(pattern="HALLMARK_", replacement="HALLMARK\n") %>% 
                  str_replace(pattern="KEGG_", replacement="KEGG\n") %>% 
                  str_replace(pattern="REACTOME_", replacement="REACTOME\n") %>% 
                  str_replace(pattern="Orig", replacement="_original") %>% 
                  str_replace(pattern="MYC_TARGETS_V", replacement="MYC_v") %>% 
                  str_replace(pattern="INTERFERON", replacement="IFN")) + 
      coord_fixed(ratio=1) + 
      CommonTheme.Feat + 
      scale_color_gradientn(
        colors=Vec.ColorScale, 
        values=scales::rescale(quantile(range(DF.Score_0$TargetSignature), 
                                        probs=c(quantile(c(0, MidPoint),probs=c(0, 0.25, 0.5, 0.75)),
                                                quantile(c(MidPoint, 1),probs=c(0, 0.25, 0.5, 0.75, 1.0))))),
        limits=range(DF.Score_0$TargetSignature),
        guide=guide_colorbar(frame.colour="black", ticks.colour="black"))
    return(Plot)
  }
  FeatPlot.gsvaScore = 
    ggarrange(
      Func.FeatPlotGSVA("Proliferation", MidPoint.Proliferation),
      Func.FeatPlotGSVA("KEGG_Cell_cycle", MidPoint.KEGG_CellCycle),
      Func.FeatPlotGSVA("HALLMARK_E2F_TARGETS", MidPoint.E2F),
      Func.FeatPlotGSVA("HALLMARK_G2M_CHECKPOINT", MidPoint.G2M),
      Func.FeatPlotGSVA("HALLMARK_MITOTIC_SPINDLE", MidPoint.MITOTIC_SPINDLE),
      Func.FeatPlotGSVA("HALLMARK_MYC_TARGETS_V1", MidPoint.MYC_TARGETS_V1),
      Func.FeatPlotGSVA("HALLMARK_MYC_TARGETS_V2", MidPoint.MYC_TARGETS_V2),
      Func.FeatPlotGSVA("BuffaOrig", MidPoint.BuffaOrig),
      Func.FeatPlotGSVA("WinterOrig", MidPoint.WinterOrig),
      Func.FeatPlotGSVA("HALLMARK_HYPOXIA", MidPoint.HYPOXIA),
      Func.FeatPlotGSVA("myCAF", MidPoint.myCAF),
      Func.FeatPlotGSVA("iCAF", MidPoint.iCAF),
      Func.FeatPlotGSVA("apCAF", MidPoint.apCAF),
      ncol=3, nrow=5, align="hv") %>% 
    annotate_figure(top=CommonTitle.Sub, 
                    bottom=paste0(CommonCaption.Sub, "\nMidpoints of scores are shifted; Proliferation:",MidPoint.Proliferation,
                                  ", BuffaOriginal:",MidPoint.BuffaOrig,", WinterOriginal:",MidPoint.WinterOrig,
                                  ", myCAF:",MidPoint.myCAF,", iCAF:",MidPoint.iCAF,", apCAF:",MidPoint.apCAF) )
  ggsave(plot=FeatPlot.gsvaScore,
         file=paste0(DirInteg,"[Figure][SubClust_ssGSVAscores_FeaturePlot][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                     "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"), 
         dpi=300, width=15, height=15, bg = "white" )
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
  
  #-------------------------------------------------------------------
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
  }
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


### 4. Cell cycle scoring ####
Func.CCscoring = function(Dim1, Res1, Dim2, Res2){
  # 4-1. Read seurat object
  SeuObj.CAF_0 = readRDS(paste0(DirInteg,"Objects/[SeuObj_ExtractedCAF][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                                QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2,".rds") )
  CommonTitle.Sub = paste0(NumOfSamples,"case Integration(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")\n",
                           QCInfo,"\n1stClust:PCA50vf2000_ClustDim",Dim1,"Res",formatC(Res1,digits=2,format="f"),
                           "\n2ndClust:PCA50vf2000_ClustDim",Dim2,"Res",formatC(Res2,digits=2,format="f"))
  NumOfClust = length(unique(SeuObj.CAF_0$seurat_clusters))
  CommonCaption.Sub = paste0("Total ",format(ncol(SeuObj.CAF_0), big.mark=",", scientific=F), " CAFs, ",NumOfClust, " clusters.")
  SeuObj.CAF_1 = SeuObj.CAF_0 %>% 
                 CellCycleScoring(s.features=cc.genes$s.genes,
                                  g2m.features=cc.genes$g2m.genes,
                                  set.ident=T)
  DF.CCscores_0 = data.frame("subclust"=SeuObj.CAF_1$seurat_clusters,
                             "S.Score"=SeuObj.CAF_1$S.Score,
                             "G2M.Score"=SeuObj.CAF_1$G2M.Score,
                             "phase"=SeuObj.CAF_1$Phase)
  DF.AveScore_0 = group_by(DF.CCscores_0, subclust) %>% 
                  summarise(Ave.S=mean(S.Score),
                            Ave.G2M=mean(G2M.Score))
  Num.MostHighSscore = DF.AveScore_0$subclust[DF.AveScore_0$Ave.S==max(DF.AveScore_0$Ave.S)]
  Num.MostHighG2Mscore = DF.AveScore_0$subclust[DF.AveScore_0$Ave.G2M==max(DF.AveScore_0$Ave.G2M)]
  DF.FreqOfPhaseSubclust = table(DF.CCscores_0$subclust, DF.CCscores_0$phase) %>% 
                           as.data.frame() %>% 
                           set_colnames(c("Subclust","Phase","Cell"))
  CommonTheme.Vln = theme(legend.position="none",
                          axis.text = element_text(face="bold", color="black", size=25),
                          axis.title = element_text(face="bold", color="black", size=10),
                          legend.text = element_text(face = "bold"),
                          legend.title = element_text(face = "bold"),
                          plot.background = element_rect(fill="transparent", color=NA),
                          panel.background = element_rect(fill="white", color=NA),
                          panel.grid = element_blank(),
                          panel.border = element_rect(fill=NA, color = "black"))
  CommonCAFLab = paste0("CAF-", 0:(NumOfClust-1))
  names(CommonCAFLab) = as.character(0:(NumOfClust-1))
  # 4-2. Cell cycle scores by subclust, violin plot
  Plot.CCscore.Vln = 
    ggarrange(
      ggplot(DF.CCscores_0, aes(x=subclust, y=S.Score, fill=subclust)) + 
        geom_violin(scale="width", linewidth=1.0) + 
        geom_boxplot(width=0.1, fill="white", outliers=F, linewidth=1.0) +
        scale_x_discrete(labels=CommonCAFLab) +
        scale_fill_manual(values=cols2) +
        labs(x=NULL) + CommonTheme.Vln + theme(axis.text.x = element_text(angle=45, hjust=1)),
      ggplot(DF.CCscores_0, aes(x=subclust, y=G2M.Score, fill=subclust)) + 
        geom_violin(scale="width", linewidth=1.0) + 
        geom_boxplot(width=0.1, fill="white", outliers=F, linewidth=1.0) +
        scale_x_discrete(labels=CommonCAFLab) +
        scale_fill_manual(values=cols2) +
        labs(x=NULL) + CommonTheme.Vln + theme(axis.text.x = element_text(angle=45, hjust=1)),
      ncol=1, align="hv") %>% 
      annotate_figure(top=CommonTitle.Sub,
                      bottom=CommonCaption.Sub)
  ggsave(plot=Plot.CCscore.Vln,
         file=paste0(DirInteg,"[Figure][SubClust_CellCycleScores_BySubclust_VlnPlot][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                     "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"), 
         dpi=500, width=8, height=10, bg = "transparent" )
  # 4-3. Cell cycle related genes by subclust, dot plot
  ggsave(plot=
           DotPlot(SeuObj.CAF_0, features=c(cc.genes$s.genes, cc.genes$g2m.genes)) &
           labs(x=NULL, y=NULL, 
                   title=CommonTitle.Sub,
                   caption="43 genes are not in Custom Xenium 5k Panel.") &
           guides(size = guide_legend(title = "Percent\nExpressed"),
                  color = guide_colorbar(title = "Scaled\nAverage\nExpression",
                                         frame.colour="black",
                                         ticks.colour="black")) &
           RotatedAxis() &
           geom_vline(xintercept=20.5) &
           theme(text = element_text(face = "bold"),
                 axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                 axis.text = element_text(size = 15),
                 axis.title = element_text(size = 15)) &
           scale_color_gradientn(colors = c("gray90","#fee5d9","#fcbba1","#ef3b2c","#99000d")) &
           coord_fixed(ratio=1),
         file=paste0(DirInteg,"[Figure][SubClust_CellCycleScores_BySubclust_DotPlot][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                     "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"), 
         dpi=300, width=12, height=8, bg = "white" )
  # 4-3. S.score & G2/M.score, scatter plot highlight top subclust
  CommonTheme.Scatter = theme(legend.position="none",
                              axis.text = element_text(face="bold", color="black", size=15),
                              axis.title = element_text(face="bold", color="black"),
                              plot.background = element_rect(fill="transparent", color=NA),
                              panel.background = element_rect(fill="white", color=NA),
                              panel.border = element_rect(fill=NA, color=NA),
                              panel.grid = element_blank(),
                              axis.line = element_line(color = "black"))
  Plot.Scores.Scatter =
    ggarrange(
      ggplot(DF.CCscores_0, aes(S.Score, G2M.Score)) +
        geom_point(color="grey80", size=0.3, alpha=0.5) +
        geom_point(data=subset(DF.CCscores_0, subclust==Num.MostHighSscore), 
                   shape=21, fill="red", color="black", 
                   size=1.2, stroke=0.1) +
        CommonTheme.Scatter +
        labs(title=paste0("Most high S",
                          if(Num.MostHighSscore==Num.MostHighG2Mscore){" & G2M"}else{""}," score subcluster is highlighted\n",
                          "CAF-",Num.MostHighSscore," : ",
                          format(nrow(subset(DF.CCscores_0, subclust==Num.MostHighSscore)), big.mark=",", scientific=F),
                          " cells." ) ),
      if(Num.MostHighSscore==Num.MostHighG2Mscore){ 
        NULL }else{
        ggplot(DF.CCscores_0, aes(S.Score, G2M.Score)) +
        geom_point(color="grey80", size=0.3, alpha=0.5) +
        geom_point(data=subset(DF.CCscores_0, subclust==Num.MostHighG2Mscore), color="red", size=0.6) +
        CommonTheme.Scatter +
        labs(title=paste0("Most high G2M score subcluster is highlighted\n",
                          "CAF-",Num.MostHighG2Mscore," : ",
                          format(nrow(subset(DF.CCscores_0, subclust==Num.MostHighG2Mscore)), big.mark=",", scientific=F),
                          " cells." ) )
        },
      ncol=2, nrow=1, align="hv") %>% 
    annotate_figure(top=paste0("Cell-cycle scatter plot\n",CommonTitle.Sub), bottom=CommonCaption.Sub)
  ggsave(plot=Plot.Scores.Scatter,
         file=paste0(DirInteg,"[Figure][SubClust_CellCycleScores_SscoreByG2Mscore_ScatterPlot][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                     "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"), 
         dpi=300, width=10, height=6, bg = "transparent" )
  # 4-4. Phase breakdown of subclust, barplot.
  DF.FreqOfPhaseSubclust_1 = ungroup(DF.FreqOfPhaseSubclust) %>% 
                             group_by(Subclust) %>% 
                             dplyr::mutate(Total = sum(Cell) ) %>% ungroup() %>% rowwise() %>% 
                             dplyr::mutate(Prop = Cell/Total) %>% 
                             dplyr::mutate(Label = paste0(Phase,
                                                        "\n",format(Cell, big.mark=",", scientific=F),
                                                        "\n(",format(Prop*100, digits=2),"%)" ) ) %>% 
                             group_by(Subclust) %>% 
                             dplyr::mutate(Cumsum= cumsum(Prop)) %>% 
                             dplyr::mutate(CoordY= 1-Cumsum+Prop/2) %>% 
                             dplyr::rename(subclust=Subclust, phase=Phase) %>% ungroup()
  DF.FreqOfPhaseSubclust_2 = group_by(DF.FreqOfPhaseSubclust_1, subclust) %>% slice(1) %>% ungroup()
  Plot.Phase.barplot = 
    ggplot(DF.CCscores_0, aes(x=subclust)) +
    geom_bar(stat="count", position="fill", color="gray50", aes(fill=phase)) +
    geom_text(data=DF.FreqOfPhaseSubclust_1,
              aes(label=case_when(Prop<0.1 ~ "", TRUE ~ Label), 
                  y=CoordY)) +
    geom_text(data=DF.FreqOfPhaseSubclust_2,
              aes(label=scales::comma(Total),
                  x=subclust), 
              y=1.01, vjust=0, hjust=0.5,
              inherit.aes = FALSE) +
    labs(x=NULL, y="Cell proportion",
         title=CommonTitle.Sub, caption=CommonCaption.Sub) +
    scale_y_continuous(limits=c(0, 1.05),
                       expand=expansion(mult=c(0, 0) ),
                       labels=scales::percent_format()) +
    scale_x_discrete(limits=as.character(0:(NumOfClust-1)),
                     labels=CommonCAFLab) +
    scale_fill_manual(values=c("G1"="#ffffb3","G2M"="#8dd3c7","S"="#bebada")) +
    theme(axis.text = element_text(face = "bold"),
          axis.title = element_text(face = "bold"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(color = "black"))
  ggsave(plot=Plot.Phase.barplot,
         file=paste0(DirInteg,"[Figure][SubClust_CellCycleScores_PhaseBreakDownOfSubclust_BarPlot][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                     "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"), 
         dpi=300, width=10, height=6, bg = "white" )
}
Func.CCscoring(Dim1=20, Res1=1.0, Dim2=20, Res2=0.5)
Func.CCscoring(Dim1=20, Res1=1.0, Dim2=30, Res2=0.5)
Func.CCscoring(Dim1=20, Res1=1.0, Dim2=40, Res2=0.5)

### 5. Xenium view ####
Func.XenViewSubclust = function(Dim1, Res1, Dim2, Res2){
  library(arrow)
  for(i in 1:NumOfSamples){
    TX = TXnumInteg[i]
    List.MetaData = readRDS(file=paste0(DirInteg,"[MetaDataList][InitialClust][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                                        "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".rds"))
    DF.CellGroup_0 = List.MetaData[["MetaData"]] %>% 
      subset(subset=sample==TX)
    RangeX = c(min(DF.CellGroup_0$CoordX), max(DF.CellGroup_0$CoordX))
    Length_X = diff(RangeX)
    RangeY = c(min(DF.CellGroup_0$CoordY), max(DF.CellGroup_0$CoordY))
    Length_Y = diff(RangeY)
    DF.Boundaries.ALL_0 = read_parquet(paste0("/Volumes/Extreme SSD/Data/TX5K_",
                                            case_when(TX %in% c("15","16") ~ "15_16",
                                                      TX %in% c("18","19","20","22") ~ "18_19_20_22",
                                                      TX %in% c("27","28") ~ "27_28",
                                                      TRUE ~ TX),
                                            "/cell_boundaries.parquet") )
    DF.CellGroup_0 = read.csv(file=paste0(DirInteg,"/CellGroupTable/[SubClust_CellGroupTable_TX5K_",TX,"][",
                                        NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                                        "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2,".csv"))
    DF.CellGroup_1 = DF.CellGroup_0 %>% 
                     dplyr::mutate(CellType= str_remove(DF.CellGroup_0$group, pattern="-.*") %>% 
                                             str_replace(pattern="Clust.*", replace="NotCAF"))
    NumOfSubclusts = subset(DF.CellGroup_1, subset=CellType=="CAF")[ ,"group"] %>% 
                     unique() %>% length()
    Vec.NumOfCells = table(subset(DF.CellGroup_1, subset=CellType=="CAF")$group) %>% 
                     as.numeric()
    Vec.Cell_id.CAF = subset(DF.CellGroup_1, subset=CellType=="CAF")$cell_id
    CommonTitle.Sub = paste0(NumOfSamples,"case Integration(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")\n",
                             QCInfo,"\n1stClust:PCA50vf2000_ClustDim",Dim1,"Res",formatC(Res1,digits=2,format="f"),
                             "\n2ndClust:PCA50vf2000_ClustDim",Dim2,"Res",formatC(Res2,digits=2,format="f"))
    CommonCaption.Sub = paste0("Total ",format(length(Vec.Cell_id.CAF), big.mark=",", scientific=F), " CAFs, ",NumOfSubclusts, " sub-clusters.")
    DF.Boundaries.ALL_1 = left_join(DF.Boundaries.ALL_0, DF.CellGroup_1, by="cell_id")
    DF.Boundaries.CAF_0 = subset(DF.Boundaries.ALL_1, 
                                 subset = cell_id %in% Vec.Cell_id.CAF)
    DF.Boundaries.NotCAF_0 = subset(DF.Boundaries.ALL_1, 
                                    subset = !(cell_id %in% Vec.Cell_id.CAF) )
    DF.NucBoundaries_0 = read_parquet(paste0("/Volumes/Extreme SSD/Data/TX5K_",
                                             case_when(TX %in% c("15","16") ~ "15_16",
                                                       TX %in% c("18","19","20","22") ~ "18_19_20_22",
                                                       TX %in% c("27","28") ~ "27_28",
                                                       TRUE ~ TX),
                                             "/nucleus_boundaries.parquet") )
    DF.NucBoundaries_2 = DF.NucBoundaries_0 %>% 
      left_join(DF.CellGroup_1, by="cell_id") %>% 
      dplyr::mutate(group = case_when(is.na(group) ~ "Excluded_byQC",
                                      TRUE ~ group))
    WideOrLong = case_when( Length_Y/Length_X >= 1 ~ c("WideOrLong"="Long", "ncol"=ceiling(NumOfSubclusts/2), "nrow"=2, "width"=ceiling(NumOfSubclusts/2)*4+1, "height"=2*8),
                            Length_Y/Length_X < 1 ~ c("WideOrLong"="Wide", "ncol"=2, "nrow"=ceiling(NumOfSubclusts/2), "width"=2*8, "height"=ceiling(NumOfSubclusts/2)*4+1),
                            TRUE ~ c("WideOrLong"="Square", "ncol"=ceiling(NumOfSubclusts/2), "nrow"=2, "width"=ceiling(NumOfSubclusts/2)*4+1, "height"=2*8 ) )
    AspectRatio = max(Length_X, Length_Y) / min(Length_X, Length_Y)
    # all sub-clusters, save as RDS file
    if(WideOrLong[["WideOrLong"]] == "Long") {
      XenView.Sub.FullLayer_0 =
        ggplot(DF.Boundaries.CAF_0, 
               aes(x=vertex_x, y=vertex_y, group=label_id)) +
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
      XenView.Sub.FullLayer_0 =
        ggplot(DF.Boundaries.CAF_0, 
             aes(x=vertex_y, y=vertex_x, group=label_id)) +
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
    XenView.Sub.FullLayer_1 =
      XenView.Sub.FullLayer_0 +
      geom_polygon(fill="#00ff00", color=NA, alpha=1.0) +
      geom_polygon(aes(fill=group), color=NA, alpha=1.0) +
      #geom_polygon(aes(color=group), fill=NA, alpha=1.0, linewidth=0.3) +
      #geom_polygon(data=DF.Boundaries.NotCAF_0, fill="gray20", color=NA, alpha=0.9) +
      #geom_polygon(data=DF.Boundaries.NotCAF_0, color="gray20", fill=NA, alpha=1.0, linewidth=0.3) +
      geom_polygon(data=DF.NucBoundaries_2, fill="blue", alpha=0.3, color=NA) +
      #geom_polygon(data=DF.NucBoundaries_2, color="blue", alpha=0.8, linewidth=0.1, fill=NA) +
      labs(title=paste0("Xenium view of initial clusters, TX5K_",TX),
           x=NULL, y=NULL,
           subtitle=CommonTitle.Sub, caption=CommonCaption.Sub) +
      coord_fixed(ratio=1) + 
      #coord_cartesian(xlim=c(2500,2900), ylim=c(2500, 2900), ratio=1) + 　　#微調整用！
      theme(text = element_text(color="gray50", face="bold"),
          plot.caption = element_text(size = 5),
          legend.position = "none",
          axis.ticks = element_line(color="gray50"),
          axis.line = element_blank(),
          plot.background = element_rect(fill = "transparent", colour = NA),
          panel.background = element_rect(fill = "black", colour = NA),
          panel.border = element_rect(fill = "transparent", colour = "gray50"),
          panel.grid = element_blank())
    saveRDS(XenView.Sub.FullLayer_1,
           file=paste0(DirInteg,"XeniumView_Sub/[FigureRDS][SubClust_XeniumView_TX5K_",TX,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                       "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".rds"))
    # Only cytoplasm without Nucleus
    XenView.Sub.OnlyCP = XenView.Sub.FullLayer_1
    XenView.Sub.OnlyCP@layers[[3]] <- NULL # Nucleus fill
    XenView.Sub.OnlyCP@layers[[1]] <- NULL # CAF green fill
    #XenView.Sub.OnlyCP@layers[[1]]$aes_params$alpha <- 1.0 # CAF fill
    #XenView.Sub.OnlyCP@layers[[3]]$aes_params$alpha <- 1.0 # Not CAF fill
    #XenView.Sub.OnlyCP@layers[[6]] <- NULL # Nucleus color
    #XenView.Sub.OnlyCP@layers[[5]] <- NULL # Nucleus fill
    #XenView.Sub.OnlyCP@layers[[4]] <- NULL # Not CAF color
    #XenView.Sub.OnlyCP@layers[[2]] <- NULL # CAF color
    saveRDS(XenView.Sub.OnlyCP,
            file=paste0(DirInteg,"/XeniumView_Sub/[FigureRDS][SubClustOnlyCP_XeniumView_TX5K_",TX,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                        "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".rds"))
    ggsave(plot=XenView.Sub.OnlyCP,
           filename=paste0(DirInteg,"/XeniumView_Sub/[Figure][SubClustOnlyCP_XeniumView_TX5K_",TX,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                       "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"),
           width=as.numeric(WideOrLong[["width"]]), height=as.numeric(WideOrLong[["height"]]), 
           dpi=500, bg="transparent", limitsize=FALSE)
    # By sub-cluster
    Plot.List = list()
    for(TargetSubclust in 0:(NumOfSubclusts-1)){
      NamedVec.FillCols = rep("gray90", times=NumOfSubclusts)
      NamedVec.FillCols[TargetSubclust+1] = "#ff0000"
      names(NamedVec.FillCols) = paste0("CAF-",0:(NumOfSubclusts-1))
      NamedVec.LgdLabs = paste0("CAF-",0:(NumOfSubclusts-1),
                                " (",format(Vec.NumOfCells, big.mark=","),")")
      names(NamedVec.LgdLabs) = paste0("CAF-",0:(NumOfSubclusts-1))
      Plot = 
      ggplot(DF.Boundaries.CAF_0, 
             aes(x=vertex_x, y=vertex_y, group=label_id)) +
        geom_polygon(aes(fill=group), color=NA, linewidth=0.1) +
        geom_polygon(data=DF.Boundaries.NotCAF_0, color=NA, fill="gray50",linewidth=0.1) +
        labs(x=NULL, y=NULL, title=NamedVec.LgdLabs[[TargetSubclust+1]]) +
        scale_x_continuous(expand=expansion(mult=c(0, 0)),
                           breaks=seq(from=ceiling(RangeX[1]/1000) * 1000,
                                      to=floor(RangeX[2]/1000) * 1000,
                                      by=1000),
                           limits=c(RangeX[1],RangeX[2])) +
        scale_y_reverse(expand=expansion(mult=c(0, 0)),
                        breaks=seq(from=ceiling(RangeY[1]/1000) * 1000,
                                   to=floor(RangeY[2]/1000) * 1000,
                                   by=1000),
                        limits=c(RangeY[2],RangeY[1])) +
        coord_fixed(ratio=1) + 
        theme(text = element_text(face = "bold"),
              plot.caption = element_text(size = 5),
              legend.position = "none",
              axis.ticks = element_line(color="gray10"),
              axis.line = element_blank(),
              panel.background = element_rect(fill = "transparent", colour = "black"),
              plot.background = element_rect(fill = "transparent", colour = NA),
              panel.grid = element_blank()) +
        scale_fill_manual(values=NamedVec.FillCols)
      Plot.List[[length(Plot.List)+1]] = Plot
    }
    Plot.XenViewSubclust = 
      wrap_plots(Plot.List, ncol=as.numeric(WideOrLong[["ncol"]]), nrow=as.numeric(WideOrLong[["nrow"]])) +
      plot_annotation(title=paste0("TX5K_",TX),
                      subtitle=CommonTitle.Sub,
                      caption=paste0(CommonCaption.Sub,"\n",paste(unname(NamedVec.LgdLabs), collapse=", "))) &
      theme(plot.background  = element_rect(fill="transparent", color=NA),
            panel.background = element_rect(fill="transparent", color="black"))
    fs::dir_create(path=paste0(DirInteg,"/XeniumView_Sub"))
    ggsave(Plot.XenViewSubclust,
           file=paste0(DirInteg,"XeniumView_Sub/[Figure][SubClust_XeniumView_TX5K_",TX,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                       "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"), 
           width=as.numeric(WideOrLong[["width"]]), height=as.numeric(WideOrLong[["height"]]), 
           dpi=400, bg="transparent", limitsize=FALSE)
  }
}
Func.XenViewSubclust(Dim1=20, Res1=1.0, Dim2=30, Res2=0.5)

### 6. Spatial distribution (kNN) ####
Func.kNN = function(Dim1, Res1, Dim2, Res2, Num.k, i, Eliminate){
  library(FNN)
  library(Matrix)
  TX=TXnumInteg[i]
  TXdl=if(TX %in% c("27","28")){ "27_28" }else{ TX }
  Directory = paste0("/Volumes/Extreme SSD/Analysis/Data/TX5K_",TXdl,"/")
  DirInteg = paste0('/Volumes/Extreme SSD/Analysis/Data/IntegAnalysis/',
                    "[",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                    "_nFeatRNA:",paste(nFeatRNA,collapse="~"),"_nCountRNA:",paste(nCountRNA,collapse="~"),"/")
  DF.Groups_0 = read.csv(paste0(DirInteg,"/CellGroupTable/[SubClust_CellGroupTable_TX5K_",TX,"][",
                                NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2,".csv")) %>% 
                subset(subset= !(group %in% Eliminate) )
  Vec.Annotation = str_remove(DF.Groups_0$group, pattern="-.*") %>% 
                   str_replace(pattern="Clust.*", replacement="NotCAF")
  DF.Groups_1 = dplyr::mutate(DF.Groups_0, CorN = Vec.Annotation)
  Seu = readRDS(paste0(Directory,"Objects/[SeuObj][TX5K_",TX,"]_Countable_Mag",Mag,
     　                "_nFeat",nFeatRNA[1],"-",nFeatRNA[2],
                       "_nCount",nCountRNA[1],"-",nCountRNA[2],".rds") ) %>% 
        subset(cells = DF.Groups_0$cell_id)
  DF.Coords = data.frame("cell_id"=colnames(Seu), 
                         "X"=Seu$X, "Y"=Seu$Y)
  DF.Meta = inner_join(DF.Groups_1, 
                       DF.Coords,    # DF.Coords is ordered according to cell_id
                       by="cell_id")
  MT.Coords = as.matrix(DF.Meta[ , c("X","Y"), drop=FALSE])  #transform to matrix !!
  Vec.ClustOfCell = DF.Meta[["group"]]
  Vec.Clusters = unique(Vec.ClustOfCell)
  # kNN indices (exclude self by default in numeric)
  # calculate distance All - ALL, and explore nearest k cells
  MT.NNindex = FNN::get.knn(MT.Coords, k=Num.k)$nn.index # row x column = cell number x k , values are cell index
  Vec.CAFclusts = subset(DF.Meta, subset=CorN=="CAF")[ , "group"] %>% unique()
  Vec.IndexCAFs = which(Vec.ClustOfCell %in% Vec.CAFclusts)
  MT.countNN_0 = matrix(0L, 
                    nrow=length(Vec.CAFclusts), 
                    ncol=length(Vec.Clusters),
                    dimnames = list(Vec.CAFclusts, Vec.Clusters))  # start matrix
  for(i in Vec.IndexCAFs){
    a = Vec.ClustOfCell[i]  # which cluster the CAF cell (index=i) is included ?
    b = Vec.ClustOfCell[MT.NNindex[i, ]] # which cluster Nearest k cells are included ?
    tab = table(b) # How many cells by cluster in near 30 cells
    MT.countNN_0[a, names(tab)] = MT.countNN_0[a, names(tab)] + as.integer(tab)
    # Cluster a に属する細胞たちのnearest k cellに、Cluster bの細胞が合計で何回現れたか
  }
  DF.countNN_log = log1p(MT.countNN_0) %>% 
                   as.data.frame() %>% 
                   rownames_to_column(var="CAF_subcluster") %>% 
                   pivot_longer(-CAF_subcluster, 
                                names_to="Neigbors",
                                values_to="Log1pCounts")
  DF.countNN_log$Neigbors = factor(DF.countNN_log$Neigbors,
                                   levels=Vec.Clusters)
  # row-normalized proportion
  row_sums = rowSums(MT.countNN_0) # 各クラスタに属する細胞の近傍細胞の合計 ≒  クラスター細胞数　x k値
  MT.countNN_rowNorm = sweep(MT.countNN_0, 1, row_sums, "/") # normalized to each cluster size
  DF.countNN_rowNorm = MT.countNN_rowNorm %>% 
    as.data.frame() %>% 
    rownames_to_column(var="CAF_subcluster") %>% 
    pivot_longer(-CAF_subcluster, 
                 names_to="Neigbors",
                 values_to="RowNormedProp")
  DF.countNN_rowNorm$Neigbors = factor(DF.countNN_rowNorm$Neigbors,
                                       levels=Vec.Clusters)
  CommonTitle.Sub = paste0(NumOfSamples,"case Integration(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")\n",
                           QCInfo,"\n1stClust:PCA50vf2000_ClustDim",Dim1,"Res",formatC(Res1,digits=2,format="f"),
                           "\n2ndClust:PCA50vf2000_ClustDim",Dim2,"Res",formatC(Res2,digits=2,format="f"))
  CommonTheme = theme(axis.text = element_text(face = "bold"),
                      axis.text.x=element_text(angle=45, hjust=1),
                      axis.title = element_text(face = "bold"),
                      legend.text = element_text(face = "bold"),
                      legend.title = element_text(face = "bold"),
                      axis.line = element_blank(),
                      panel.border = element_rect(fill = NA, color = "black"))
  Plot.Tile = 
  ggarrange(
    ggplot(DF.countNN_log, 
           aes(y=CAF_subcluster, x=Neigbors, fill=Log1pCounts)) +
      geom_tile(color="gray50") +
      labs(title=paste0("TX5K_",TX,", kNN nighbor counts (log1p)"),
           y="Focal (CAF sub-clusters)", 
           fill="Neighbor\ncell\ncounts\n(log1p)") +
      scale_x_discrete(expand=expansion(mult=c(0,0))) +
      scale_y_discrete(expand=expansion(mult=c(0,0))) +
      scale_fill_gradientn(colors=c("#253494", "#2c7fb8","#41b6c4", "#a1dab4", "#ffffcc",
                                    "#fecc5c", "#fd8d3c","#f03b20", "#bd0026"),
                           guide=guide_colorbar(direction="vertical",
                                                frame.colour="black",
                                                ticks.colour="black")) + CommonTheme,
    ggplot(DF.countNN_rowNorm, 
           aes(y=CAF_subcluster, x=Neigbors, fill=RowNormedProp)) +
      geom_tile(color="gray50") +
      labs(title=paste0("TX5K_",TX,", kNN neighbor composition of focal CAF)"),
           y="Focal (CAF sub-clusters)", 
           fill="Neighbor\nproportion\n(within\nkNN of\nfocal CAF)") +
      scale_x_discrete(expand=expansion(mult=c(0,0))) +
      scale_y_discrete(expand=expansion(mult=c(0,0))) +
      scale_fill_gradientn(colors=c("#253494", "#2c7fb8","#41b6c4", "#a1dab4", "#ffffcc",
                                    "#fecc5c", "#fd8d3c","#f03b20", "#bd0026"),
                           guide=guide_colorbar(direction="vertical",
                                                frame.colour="black",
                                                ticks.colour="black")) + CommonTheme,
    ncol=1, nrow=2, align="hv") %>% 
    annotate_figure(top=text_grob(CommonTitle.Sub),
                    bottom=text_grob(paste0("Heatmaps showing the total counts (log1p) of k-nearest neighbors and the row-normalized neighbor\n",
                                     "composition for each CAF subcluster. For the lower panel, values represent the proportion of\n",
                                     "neighboring cells within the kNN of each focal CAF subcluster.")))
  ggsave(Plot.Tile,
         file=paste0(DirInteg,"[Figure][SubClust_kNN_TX5K_",TX,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                     "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"), 
         width=12, height=8, dpi=400)
  
  }
Func.kNN(Dim1=20, Res1=1.0, Dim2=30, Res2=0.5, Num.k=30, i=1, Eliminate=c("Clust.4","Clust.5","Clust.8"))
Func.kNN(Dim1=20, Res1=1.0, Dim2=30, Res2=0.5, Num.k=30, i=2, Eliminate=c("Clust.4","Clust.5","Clust.8"))
Func.kNN(Dim1=20, Res1=1.0, Dim2=30, Res2=0.5, Num.k=30, i=3, Eliminate=c("Clust.4","Clust.5","Clust.8"))
Func.kNN(Dim1=20, Res1=1.0, Dim2=30, Res2=0.5, Num.k=30, i=4, Eliminate=c("Clust.4","Clust.5","Clust.8"))
Func.kNN(Dim1=20, Res1=1.0, Dim2=30, Res2=0.5, Num.k=30, i=5, Eliminate=c("Clust.4","Clust.5","Clust.8"))
Func.kNN(Dim1=20, Res1=1.0, Dim2=30, Res2=0.5, Num.k=30, i=6, Eliminate=c("Clust.4","Clust.5","Clust.8"))

### 8. Spatial distribution (seurat build niche assay) ####
### 8-1. Calculate kNN count and assign niche number ####
  Unclassified = "Retain" # or Unclassified = "Exclude"
  Labels = "seurat_clusters" # or Labels = "cell_type"
  n.neighbors = 30
Func.NicheAssay = function(Dim1, Res1, Dim2, Res2, n.neighbors, Unclassified, Labels){
  # Eliminate unclassified clusters (optional)
  List.InitialClustData = readRDS(file=paste0(DirInteg,"[MetaDataList][InitialClust][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                                              "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".rds"))
  Add = List.InitialClustData$OrderedClustnumToLab
  names(Add) = paste0("Clust.",names(Add))
  List.InitialClustData[["ForPlotLab"]] = Add
  EliminateClusters = names(List.InitialClustData$OrderedClustnumToCelltype)[as.character(List.InitialClustData$OrderedClustnumToCelltype) %in% "Unclassified"]
  DF.Niche = data.frame()
  DF.CellidAndNiche = data.frame()
  for(i in 1:NumOfSamples){
    TX=TXnumInteg[i]
    message("Processing TX = ", TX)
    DF.CGtable_Sub_0 = 
      read.csv(file=paste0(DirInteg,"/CellGroupTable/[SubClust_CellGroupTable_TX5K_",TX,"][",
                           NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                           "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2,".csv")) 
    DF.CGtable_Sub_1 = DF.CGtable_Sub_0 %>% 
      dplyr::mutate(cell_type = str_remove(group, pattern="Clust."))
    DF.CGtable_Sub_1$cell_type[!grepl("^CAF-", DF.CGtable_Sub_1$cell_type)] = 
      List.InitialClustData$OrderedClustnumToCelltype[DF.CGtable_Sub_1$cell_type[!grepl("^CAF-", DF.CGtable_Sub_1$cell_type)]]
    DF.CGtable_Sub_2 = DF.CGtable_Sub_1 %>% 
      dplyr::mutate(Exclusion = case_when(group %in% paste0("Clust.",EliminateClusters) ~ "Exclude",
                                          TRUE ~ "Retain") )
    DF.CGtable_ForAssay = DF.CGtable_Sub_2 %>% 
      dplyr::rename(seurat_clusters = group) %>% 
      subset(subset = Exclusion %in% if(Unclassified=="Exclude"){
                                         "Retain" } else { c("Exclude", "Retain") } )
    Xen_0 = readRDS(paste0("/Volumes/Extreme SSD/Analysis/TEP/Objects/[XenObj]_TX5K_",
                           if(TX %in% c("27","28")){ "27_28" }else{ TX },
                           "_countable_coords.rds") ) %>% 
            subset(cells=DF.CGtable_ForAssay$cell_id)
    DF.CGtable_ForAssay$cell_id = 
      factor(DF.CGtable_ForAssay$cell_id, levels=colnames(Xen_0))
    DF.CGtable_ForAssay = dplyr::arrange(DF.CGtable_ForAssay, cell_id)
    if(identical(colnames(Xen_0), 
                 as.character(DF.CGtable_ForAssay$cell_id) )){
    Xen_0$LabelsForNicheAssay = DF.CGtable_ForAssay[[Labels]]
  }
    set.seed(123)
    Xen_1 = BuildNicheAssay(
      object = Xen_0,
      fov = "fov",
      group.by = "LabelsForNicheAssay",
      neighbors = n.neighbors)
    # save niche label data
    DF.CellidAndNiche_byCase = data.frame(
      full_id = paste0("TX5K",TX,"_",colnames(Xen_1)),
      cell_id = colnames(Xen_1),
      niches = Xen_1$niches,
      row.names = NULL)
    DF.CellidAndNiche = rbind(DF.CellidAndNiche, DF.CellidAndNiche_byCase)
    # save kNN count data
    library(Matrix)
    MT.kNNcounts_0 = Xen_1[["niche"]]@counts
    MT.kNNcounts_1 = t(MT.kNNcounts_0)
    rownames(MT.kNNcounts_1) = paste0("TX5K",TX,"_",rownames(MT.kNNcounts_1))
    if(i == 1){
      MT.kNNcounts_All = Matrix::Matrix(
        0,
        nrow = 0,
        ncol = length(colnames(MT.kNNcounts_1)),
        sparse = TRUE,
        dimnames = list(NULL, colnames(MT.kNNcounts_1))
      )
    }
    MT.kNNcounts_All = rbind(MT.kNNcounts_All,
                             MT.kNNcounts_1)
    #ImageDimPlot(Xen_1, 
    #             group.by="niches", 
    #             crop=TRUE, 
    #             alpha=1, 
    #             fov="fov")
    DF.count_0 = table(niche = Xen_1$niches,
                       LabelsForNicheAssay = Xen_1$LabelsForNicheAssay) %>% as.data.frame() %>% 
                 dplyr::rename(n_cell = Freq) %>% 
                 dplyr::arrange(niche) %>% 
                 dplyr::mutate(niche = paste0("Niche_",TX,"_",niche))
    DF.count_1 = DF.count_0 %>% group_by(niche) %>% 
      dplyr::mutate(total_byniche = sum(n_cell)) %>% 
      dplyr::mutate(prop_byniche = n_cell/total_byniche) %>% 
      ungroup() %>% group_by(LabelsForNicheAssay) %>% 
      dplyr::mutate(total_byclust = sum(n_cell)) %>% 
      dplyr::mutate(prop_byclust = n_cell/total_byclust)
    DF.count_2 = cbind(sample=paste0("TX5K_",TX),
                       DF.count_1)
    DF.Niche = rbind(DF.Niche, DF.count_2)
    }
  fs::dir_create(path=paste0(DirInteg,"/NicheAnalysisData"))
  saveRDS(MT.kNNcounts_All,
          file=paste0(DirInteg,"/NicheAnalysisData/[CountDataRDS][NicheAnalysis_",Unclassified,"Unclassified_Labeled",Labels,"_Neighbor",n.neighbors,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                      "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".rds"))
  write.csv(DF.CellidAndNiche,
            file=paste0(DirInteg,"/NicheAnalysisData/[DataTable][NicheAnalysis_CellidAndNiche_",Unclassified,"Unclassified_Labeled",Labels,"_Neighbor",n.neighbors,"_TX5K",TX,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                        "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".csv"),
            row.names=F)
  write.csv(DF.Niche,
         file=paste0(DirInteg,"NicheAnalysisData/[DataTable][NicheAnalysis_CellCount_",Unclassified,"Unclassified_Labeled",Labels,"_Neighbor",n.neighbors,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                     "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".csv"),
         row.names=F)
}
Func.NicheAssay(Dim1=20, Res1=1.0, Dim2=30, Res2=0.5, n.neighbors=20, Unclassified="Retain", Labels="seurat_clusters")
Func.NicheAssay(Dim1=20, Res1=1.0, Dim2=30, Res2=0.5, n.neighbors=30, Unclassified="Retain", Labels="seurat_clusters")
Func.NicheAssay(Dim1=20, Res1=1.0, Dim2=30, Res2=0.5, n.neighbors=35, Unclassified="Retain", Labels="seurat_clusters")
Func.NicheAssay(Dim1=20, Res1=1.0, Dim2=30, Res2=0.5, n.neighbors=40, Unclassified="Retain", Labels="seurat_clusters")
Func.NicheAssay(Dim1=20, Res1=1.0, Dim2=30, Res2=0.5, n.neighbors=40, Unclassified="Exclude", Labels="seurat_clusters")


### 8-2a. k-means clustering method ####
Dim1=20
Res1=1.0
Dim2=30
Res2=0.5
n.neighbors=40
Func.Kmeans1 = function(Dim1, Res1, Dim2, Res2, n.neighbors, Unclassified, Labels){
  # load count data
  CommonTitle.Sub = paste0(NumOfSamples,"case Integration(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")\n",
                           QCInfo,"\n1stClust:PCA50vf2000_ClustDim",Dim1,"Res",formatC(Res1,digits=2,format="f"),
                           "\n2ndClust:PCA50vf2000_ClustDim",Dim2,"Res",formatC(Res2,digits=2,format="f"))
  MT.kNNcounts_0 = readRDS(
    file=paste0(DirInteg,"/NicheAnalysisData/[CountDataRDS][NicheAnalysis_",Unclassified,"Unclassified_Labeled",Labels,"_Neighbor",n.neighbors,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".rds"))
  MT.kNNcounts_1 = MT.kNNcounts_0 / rowSums(MT.kNNcounts_0)
  MT.kNNcounts_2 = scale(MT.kNNcounts_1)
  set.seed(123)
  pca = prcomp(MT.kNNcounts_2, rank.=30)
  var_ratio <- pca$sdev^2 / sum(pca$sdev^2)
  cum_var   <- cumsum(var_ratio)
  LeftAxisMax = 0.25
  DF.pca <- data.frame(
    PC = seq_along(var_ratio),
    Variance = var_ratio,
    Cumulative = cum_var) %>% 
    dplyr::mutate(Cumulative_scale = Cumulative*LeftAxisMax)
  Plot = 
  ggplot(DF.pca, aes(x = PC)) +
    geom_vline(xintercept = c(5, 10, 15, 20, 25, 30, 35, 40), color = "gray90") +
    annotate(geom="rect", 
             xmin=0, xmax=45, ymin=0.6*LeftAxisMax, ymax=0.8*LeftAxisMax,
             alpha=0.2) +
    annotate(geom="rect", 
             xmin=0, xmax=45, ymin=0.65*LeftAxisMax, ymax=0.75*LeftAxisMax,
             alpha=0.2) +
    geom_col(aes(y = Variance), fill = "grey70", width = 0.8) +
    geom_line(aes(y = Cumulative_scale), color = "black", linewidth = 1) +
    geom_point(aes(y = Cumulative_scale), color = "black", size = 2) +
    labs(title="PCs and Explained variance ratio (decision of pc_use)",
         subtitle=paste0(CommonTitle.Sub,"\n",
                        "UnclassifiedClusters:",Unclassified,"_Label:",Labels,"_Neighbor:",n.neighbors),
         caption="Explained variance ratio : SD^2 / sum(SD^2)",
         x="Principal component (PC)") +
    scale_x_continuous(breaks=seq(5, 45, by=5),
                       expand = expansion(mult=c(0,0))) +
    scale_y_continuous(
      expand = expansion(mult=c(0,0)),
      name = "Explained variance ratio",
      sec.axis = sec_axis(~ ./LeftAxisMax,
                          name = "Cumulative explained variance",
                          breaks = seq(0, 1, by=0.1))) +
    theme_classic()
  ggsave(plot=Plot,
         file=paste0(DirInteg,"NicheAnalysisData/[Figure][KmeansClust_PCandExplain_",Unclassified,"Unclassified_Labeled",Labels,"_Neighbor",n.neighbors,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                     "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"),
         width=5, height=5, dpi=200)
  saveRDS(pca,
         file=paste0(DirInteg,"NicheAnalysisData/[RDS][KmeansClust_PCAres][",Unclassified,"Unclassified_Labeled",Labels,"_Neighbor",n.neighbors,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                     "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".rds"))
}
Func.Kmeans1(Dim1=20, Res1=1.0, Dim2=30, Res2=0.5, n.neighbors=20, Unclassified="Retain", Labels="seurat_clusters")
Func.Kmeans1(Dim1=20, Res1=1.0, Dim2=30, Res2=0.5, n.neighbors=30, Unclassified="Retain", Labels="seurat_clusters")
Func.Kmeans1(Dim1=20, Res1=1.0, Dim2=30, Res2=0.5, n.neighbors=35, Unclassified="Retain", Labels="seurat_clusters")
Func.Kmeans1(Dim1=20, Res1=1.0, Dim2=30, Res2=0.5, n.neighbors=40, Unclassified="Retain", Labels="seurat_clusters")
Func.Kmeans1(Dim1=20, Res1=1.0, Dim2=30, Res2=0.5, n.neighbors=40, Unclassified="Exclude", Labels="seurat_clusters")

Func.Kmeans2 = function(Dim1, Res1, Dim2, Res2, n.neighbors, Unclassified, Labels, pc_use){
  set.seed(123)
  pca = readRDS(file=paste0(DirInteg,"NicheAnalysisData/[RDS][KmeansClust_PCAres][",Unclassified,"Unclassified_Labeled",Labels,"_Neighbor",n.neighbors,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                            "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".rds"))
  Ncell = nrow(pca$x)
  NcellSampled = 50000
  SampledIndex = sample.int(Ncell, size=min(NcellSampled, Ncell))
  pcaSampled = pca$x[SampledIndex, 1:pc_use]
  ks = 2:24
  wss = 
    sapply(ks, function(k){
          kmeans(
            pcaSampled,
            centers=k,
            nstart=20,
            iter.max=50)$tot.withinss
    } )
  png(paste0(DirInteg,"/NicheAnalysisData/[Figure][KmeansClust_ElbowPlot_npc",pc_use,"][",Unclassified,"Unclassified_Labeled",Labels,"_Neighbor",n.neighbors,
             "][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
             "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"),
      width=1800, height=1500, res=300)
  plot(ks, wss/wss[1], type="b",
       xlab = "k\n(Number of spatial niches)",
       ylab = "Relative within-cluster sum of squares (WSS) ",
       main = "Elbow plot for decision of k (number of niches)",
       sub = paste0(Unclassified,"Unclassified_Labeled",Labels,"_Neighbor",n.neighbors,
                    "_Sampling:",NcellSampled) )
  abline(v = c(5, 10, 15, 20), col="gray80")
  dev.off()
}
Func.Kmeans2(pc_use=15, Dim1=20, Res1=1.0, Dim2=30, Res2=0.5, n.neighbors=40, Unclassified="Retain", Labels="seurat_clusters")
Func.Kmeans2(pc_use=15, Dim1=20, Res1=1.0, Dim2=30, Res2=0.5, n.neighbors=40, Unclassified="Exclude", Labels="seurat_clusters")

Func.Kmeans3 = function(Dim1, Res1, Dim2, Res2, n.neighbors, Unclassified, Labels, pc_use, k.niche){
  set.seed(123)
  pca = readRDS(file=paste0(
                  DirInteg,"NicheAnalysisData/[RDS][KmeansClust_PCAres][",Unclassified,"Unclassified_Labeled",Labels,"_Neighbor",n.neighbors,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                  "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".rds"))
  km = kmeans(pca$x[, 1:pc_use], 
              centers = k.niche, 
              nstart = 50, 
              iter.max = 100)
  DF.KmeansRes = data.frame(
                  cell_id = names(km$cluster),
                  niche = as.character(km$cluster))
  message("**MyNote** iter = ", km$iter, " <100 ならok")
  message("**MyNote** ifault = ", km$ifault, " 0が望ましい")
  #message(table(km$cluster) # 偏りないように
  write.csv(DF.KmeansRes, file=paste0(
        DirInteg,"/NicheAnalysisData/[KmeansClust][",Unclassified,"Unclassified_Labeled",Labels,"_Neighbor",n.neighbors,
        "_npc",pc_use,"_k",k.niche,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
        "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".csv"),
        row.names = F)
}
Func.Kmeans3(pc_use=10, k.niche=10, Dim1=20, Res1=1.0, Dim2=30, Res2=0.5, n.neighbors=40, Unclassified="Retain", Labels="seurat_clusters")
Func.Kmeans3(pc_use=15, k.niche=10, Dim1=20, Res1=1.0, Dim2=30, Res2=0.5, n.neighbors=40, Unclassified="Retain", Labels="seurat_clusters")
Func.Kmeans3(pc_use=15, k.niche=10, Dim1=20, Res1=1.0, Dim2=30, Res2=0.5, n.neighbors=40, Unclassified="Exclude", Labels="seurat_clusters")
Func.Kmeans3(pc_use=15, k.niche=11, Dim1=20, Res1=1.0, Dim2=30, Res2=0.5, n.neighbors=40, Unclassified="Exclude", Labels="seurat_clusters")

Func.Kmeans.Visualize = function(Dim1, Res1, Dim2, Res2, n.neighbors, Unclassified, Labels, pc_use, k.niche){
k.niche = 11
pc_use = 15
Dim1=20 
Res1=1.0 
Dim2=30 
Res2=0.5
n.neighbors=40
Unclassified="Exclude"
Labels="seurat_clusters"
method="SpearmanAverage"
CAForder = c(4,8,0,1,2,5,3,7,6)
  # 0-1. load data
  List.InitialClustData = readRDS(file=paste0(DirInteg,"[MetaDataList][InitialClust][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                                              "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".rds"))
  SeuObj.CAF_0 = readRDS(paste0(DirInteg,"Objects/[SeuObj_ExtractedCAF][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                                QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2,".rds") )
  NumOfCAFsubclusts = length(unique(SeuObj.CAF_0$seurat_clusters))
  DF.KmeansRes = read.csv(file=paste0(
      DirInteg,"/NicheAnalysisData/[KmeansClust][",Unclassified,"Unclassified_Labeled",Labels,"_Neighbor",n.neighbors,
      "_npc",pc_use,"_k",k.niche,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
      "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".csv"),
      row.names=1)
  CommonTitle.Sub = paste0(NumOfSamples,"case Integration(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")\n",
                           QCInfo,"\n1stClust:PCA50vf2000_ClustDim",Dim1,"Res",formatC(Res1,digits=2,format="f"),
                           "\n2ndClust:PCA50vf2000_ClustDim",Dim2,"Res",formatC(Res2,digits=2,format="f"))
  # 0-2. data molding
  DF.CellGroup_CAF = data.frame(
    cell_id = colnames(SeuObj.CAF_0),
    seurat_subclusters = paste0("CAF-", SeuObj.CAF_0$seurat_clusters) )
  DF.MetaData = List.InitialClustData$MetaData %>% 
    dplyr::mutate(cell_type = List.InitialClustData[["OrderedClustnumToCelltype"]][as.character(seurat_clusters)]) %>% 
    inner_join(rownames_to_column(DF.KmeansRes, var="cell_id"), by="cell_id") %>% 
    left_join(DF.CellGroup_CAF, by="cell_id") %>% 
    dplyr::mutate(seurat_subclusters = 
                    ifelse( !is.na(seurat_subclusters), seurat_subclusters, Label)) %>% 
    dplyr::mutate(niche = paste0("niche_",niche))
  ########################
  ### 1. Proportion heatmap
  DF.CellCount_0 = 
    table(DF.MetaData$seurat_subclusters, DF.MetaData$niche) %>% as.data.frame() %>% 
    set_colnames(c("seurat_subclusters", "niche", "ncell")) %>% 
    group_by(seurat_subclusters) %>% 
    dplyr::mutate(total_inSubclust = sum(ncell)) %>% 
    dplyr::mutate(prop_inSubclust = ncell/total_inSubclust) %>% 
    ungroup() %>% group_by(niche) %>% 
    dplyr::mutate(total_inNiche = sum(ncell)) %>% 
    dplyr::mutate(prop_inNiche = ncell/total_inNiche)
  # hclust
  DF.PropInNiche_0 = DF.CellCount_0[ , c("niche", "seurat_subclusters", "prop_inNiche")]
  DF.PropInNiche_1 = DF.PropInNiche_0 %>% 
    pivot_wider(names_from=seurat_subclusters,
                values_from=prop_inNiche,
                values_fill=0) %>% 
    column_to_rownames("niche")
  if(method=="SpearmanAverage"){
      cor_niche = cor(t(DF.PropInNiche_1), method="spearman", use="pairwise.complete.obs")
      d_niche = stats::as.dist(1 - cor_niche)
      Hclust_niche = stats::hclust(d_niche, method="average")
      MethodLabel = "Hierarchical clustering was performed using average linkage based on Spearman correlation distances."
  } else if (method=="EuclideanWardd2"){               # ward.D2なら距離はEuclideanで！！
      d_niche = dist(DF.PropInNiche_1, method="euclidean")
      Hclust_niche = stats::hclust(d_niche, method = "ward.D2")
      MethodLabel = "Hierarchical clustering was performed using Ward’s method (ward.D2) with Euclidean distances."
  }
  Vec.nicheorder = (Hclust_niche$labels)[Hclust_niche$order]
  Vec.LabelOrder=c(paste0("CAF-", CAForder), 
                   unname(List.InitialClustData$OrderedClustnumToLab))
  DF.CellCount_0$niche = factor(DF.CellCount_0$niche, levels=Vec.nicheorder)
  DF.CellCount_0$seurat_subclusters = factor(DF.CellCount_0$seurat_subclusters, levels=Vec.LabelOrder)
  # dendrogram
  library(ggdendro)
  dend_niche = as.dendrogram(Hclust_niche)
  dend_data = dendro_data(dend_niche)
  p_dend = 
    ggplot(dend_data$segments) +
    geom_segment(aes(x=y, y=x, xend=yend, yend=xend)) +
    scale_y_continuous(limits = c(0.5, length(Hclust_niche$labels) + 0.5),
                       breaks = seq_along(Hclust_niche$labels),
                       labels = Hclust_niche$labels,
                       expand = expansion(mult=c(0,0))) +
    scale_x_reverse() +
    theme_void()
  # Heatmap
  Vec.NicheAnnot = case_when(
    k.niche==10  ~ c("niche_1"="N1.Nerve Niche",        "niche_2"="N2.Atrophic Acinar Niche",
                     "niche_3"="N3.T cell-rich Niche",  "niche_4"="N4.Normal Acinar Niche",
                     "niche_5"="N5.Vascular-rich Niche","niche_6"="N6.Islet Niche",
                     "niche_7"="N7.Cancer Niche",       "niche_8"="N8.B cell-rich Niche",
                     "niche_9"="N9.Normal~PanIN Niche", "niche_10"="N10.Myeloid-rich Niche",
                     "niche_11"=""),
    k.niche==11  ~ c("niche_1"="N1.Vascular-rich Niche",  "niche_2"="N2.PanIN Niche",
                     "niche_3"="N3.T+Myelo Niche",        "niche_4"="N4.Islet Niche",
                     "niche_5"="N5.Atrophic Acinar Niche","niche_6"="N6.Cancer Niche",
                     "niche_7"="N7.Plasma cell Niche",    "niche_8"="N8.Nerve Niche",
                     "niche_9"="N9.T+B Niche",         "niche_10"="N10.Normal duct Niche",
                     "niche_11"="N11.Normal Acinar Niche"))
  Vec.NicheColor = case_when(
    k.niche==10  ~ c("niche_1"="#b2df8a",  "niche_2"="#a6cee3",
                     "niche_3"="#33a02c",  "niche_4"="#1f78b4",
                     "niche_5"="#fb9a99",  "niche_6"="#ffd000",
                     "niche_7"="#e31a1c",  "niche_8"="#ff7f00",
                     "niche_9"="#cab2d6", "niche_10"="#6a3d9a",
                     "niche_11"=""),
    k.niche==11  ~ c("niche_1"="#b2df8a",  "niche_2"="#6a3d9a",
                     "niche_3"="#fb9a99",  "niche_4"="#b15928",
                     "niche_5"="#a6cee3",  "niche_6"="#e31a1c",
                     "niche_7"="#ffd000",  "niche_8"="#ff7f00",
                     "niche_9"="#cab2d6", "niche_10"="#33a02c",
                     "niche_11"="#1f78b4"))
  CommonScaleX = scale_x_discrete(expand=expansion(mult=c(0,0)))
  CommonScaleY = scale_y_discrete(expand=expansion(mult=c(0,0)),
                                  labels=Vec.NicheAnnot,
                                  position="right")
  CommonScaleFill = scale_fill_viridis_c(option="magma",
                                         trans="sqrt",
                                         guide=guide_colorbar(
                                         #direction="horizontal",
                                         #title.position="top",
                                         #title.hjust=0.5,
                                         frame.colour="black",
                                         ticks.colour="black"))
  CommonTheme.HM = theme(plot.margin = margin(t=0, r=0.3, b=0, l=0, unit="cm"),
                      legend.position="bottom",
                      text = element_text(face="bold", color="black"),
                      axis.text.x = element_text(angle=45, hjust=1, face="bold",color="black"),
                      axis.text.y = element_text(face="bold",color="black"),
                      legend.text = element_text(angle=45, hjust=1, face="bold",color="black"),
                      panel.border = element_rect(fill = NA, color = "black"))
  p_heat1 =
    ggplot(DF.CellCount_0, aes(x=seurat_subclusters, y=niche, fill=prop_inNiche)) +
    geom_tile(color="gray20") + 
    labs(x=NULL, y=NULL, fill="Proportion of clusters\nwithin each niche (%)") +
    CommonScaleX + CommonScaleY +
    scale_fill_viridis_c(option="magma",
                         trans="sqrt",
                         labels=scales::percent_format(),
                         guide=guide_colorbar(
                         #direction="horizontal", title.position="top", 
                         title.vjust=1,
                         frame.colour="black",
                         ticks.colour="black")) + 
    CommonTheme.HM
  p_heat2 = 
    ggplot(DF.CellCount_0, aes(x=seurat_subclusters, y=niche, fill=prop_inSubclust)) +
    geom_tile(color="gray20") + 
    labs(x=NULL, y=NULL, fill="Proportion of niches\nwithin each cluster (%)") +
    CommonScaleX + CommonScaleY +
    scale_fill_viridis_c(option="mako",
                         trans="sqrt",
                         labels=scales::percent_format(),
                         guide=guide_colorbar(
                         #direction="horizontal", title.position="top",
                         title.vjust=1,
                         frame.colour="black",
                         ticks.colour="black") ) + 
    CommonTheme.HM
  ggsave(plot = p_dend + p_heat1 + p_heat2 +
         plot_layout(widths = c(0.5, 4, 4)) +
         plot_annotation(
           subtitle=paste0(CommonTitle.Sub,"\n",
                           "UnclassifiedClusts:",Unclassified,"_Label:",Labels,"_N.neighbors:",n.neighbors,"\n",
                           "npc",pc_use,"_k",k.niche),
           caption=paste0(MethodLabel) ),
         file=paste0(
            DirInteg,"/[Figure][Niche_Kmeans_ProportionHeatmap][",Unclassified,"Unclassified_Labeled",Labels,"_Neighbor",n.neighbors,
            "_npc",pc_use,"_k",k.niche,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
            "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"),
         height=6, width=20, dpi=400)
  ########################
  ###  2. Bar plot
  DF.MetaData$niche = factor(DF.MetaData$niche, 
                             levels=rev(Vec.nicheorder))
  DF.MetaData$seurat_subclusters = factor(DF.MetaData$seurat_subclusters,
                                          levels=Vec.LabelOrder)
  CommonTheme.Bar = theme(plot.margin = margin(t=0, r=0.3, b=0, l=0, unit="cm"),
                         text = element_text(face="bold", color="black"),
                         plot.subtitle = element_text(size=7),
                         axis.text.x = element_text(angle=45, hjust=1, face="bold",color="black"),
                         axis.text.y = element_text(face="bold",color="black"),
                         panel.border = element_rect(fill = NA, color = "black", linewidth=0.7))
  Plot.BarNiche = 
  ggplot(DF.MetaData, aes(y=seurat_subclusters, fill=niche)) + 
    geom_bar(stat="count", position="fill", color="gray50", linewidth=0.3,
             width=1) +
    labs(subtitle=paste0(CommonTitle.Sub,"\n",
                         "UnclassifiedClusts:",Unclassified,"_Label:",Labels,"_N.neighbors:",n.neighbors,"\n",
                         "npc",pc_use,"_k",k.niche),
         y=NULL, x="Niche proportion in cluster", 
         fill=paste0("Niche (k.niche=",k.niche,")")) +
    scale_x_continuous(expand=expansion(mult=c(0,0)), 
                       labels=scales::percent_format()) +
    scale_y_discrete(expand=expansion(mult=c(0,0))) +
    scale_fill_manual(#values=c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3',
                      #         '#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd' ),
                      values=Vec.NicheColor,
                      labels=Vec.NicheAnnot) +
    geom_hline(yintercept=NumOfCAFsubclusts+0.5) +
    CommonTheme.Bar
  ggsave(plot = Plot.BarNiche,
         file=paste0(
           DirInteg,"/[Niche_Kmeans_NichePropInClust][",Unclassified,"Unclassified_Labeled",Labels,"_Neighbor",n.neighbors,
           "_npc",pc_use,"_k",k.niche,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
           "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"),
         height=7, width=5, dpi=400)
  ggplot(DF.MetaData, aes(y=niche, fill=cell_type)) + 
    geom_bar(stat="count", position="fill", color="gray50") +
    theme(axis.text.x = element_text(angle=45, hjust=1))
  #ggplot(DF.MetaData, aes(x=niche, fill=seurat_subclusters)) + 
  #  geom_bar(stat="count", position="fill", color="gray50") +
  #  theme(axis.text.x = element_text(angle=45, hjust=1))

  ########################
  ###  3. Xenium view (point)
  for(i in 1:NumOfSamples){
  TX = TXnumInteg[i]
  DF.MetaData_subset_0 = subset(DF.MetaData, subset=sample==TX)
  RangeX = c(min(DF.MetaData_subset_0$CoordX), max(DF.MetaData_subset_0$CoordX))
  RangeY = c(min(DF.MetaData_subset_0$CoordY), max(DF.MetaData_subset_0$CoordY))
  WideOrLong = case_when( diff(RangeY)/diff(RangeX) > 3/2 ~ c("WideOrLong"="Long", "width"=5, "height"=10.5),
                          diff(RangeY)/diff(RangeX) < 2/3 ~ c("WideOrLong"="Wide", "width"=10.5, "height"=5),
                          TRUE ~ c("WideOrLong"="Square", "width"=5, "height"=5.5 ) )
  Xen.Niche.Point = 
  ggplot(subset(DF.MetaData, subset=sample==TX),
         aes(x=CoordX, y=CoordY, color=niche)) + 
    geom_point(size=0.2) + 
    labs(title=paste0("TX5K_",TX),
         subtitle=paste0(CommonTitle.Sub,"\n",
                         "UnclassifiedClusts:",Unclassified,"_Label:",Labels,"_N.neighbors:",n.neighbors,"\n",
                         "npc",pc_use,"_k",k.niche),
         x=NULL, y=NULL, color=paste0("Niche (k.niche=",k.niche,")") ) +
    scale_x_continuous(expand=expansion(mult=c(0.00, 0.00)),
                       breaks=seq(from=ceiling(RangeX[1]/1000) * 1000,
                                  to=floor(RangeX[2]/1000) * 1000,
                                  by=1000),
                       limits=c(RangeX[1],RangeX[2])) +
    scale_y_reverse(expand=expansion(mult=c(0.00, 0.00)),
                    breaks=seq(from=ceiling(RangeY[1]/1000) * 1000,
                               to=floor(RangeY[2]/1000) * 1000,
                               by=1000),
                    limits=c(RangeY[2],RangeY[1])) +
    scale_color_manual(#values=c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3',
                       #         '#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd' ),
                       values= Vec.NicheColor,
                       label = Vec.NicheAnnot) +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    coord_fixed(ratio=1) +
    theme(text = element_text(face = "bold"),
          axis.text = element_text(color="black"),
          plot.subtitle = element_text(size=5),
          plot.caption = element_text(size=5),
          legend.position = "right",
          legend.key = element_blank(),
          axis.ticks = element_line(color="black"),
          axis.line = element_blank(),
          plot.background = element_rect(fill = "transparent", colour = NA),
          panel.background = element_rect(fill = "transparent", colour = NA),
          panel.border = element_rect(fill = "transparent", colour = "black"),
          panel.grid = element_blank())
  Dir.XenViewNiche = paste0(DirInteg,"/NicheAnalysisData/XeniumView_Sub/",QCInfo.FileName,
                            "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,
                            "/[",Unclassified,"Unclassified_Labeled",Labels,"_Neighbor",n.neighbors,
                            "_npc",pc_use,"_k",k.niche,"]")
  fs::dir_create(path=Dir.XenViewNiche)
  ggsave(plot=Xen.Niche.Point,
         file=paste0(Dir.XenViewNiche,"/[PointView]_TX5K_",TX,".png"), 
         width=as.numeric(WideOrLong[["width"]]), height=as.numeric(WideOrLong[["height"]]), 
         dpi=400, bg="transparent")
  saveRDS(Xen.Niche.Point +
          geom_point(data=subset(DF.MetaData, 
                                 sample==TX & seurat_subclusters=="CAF-8"),
                     color="black", 
                     shape=21),
         file=paste0(Dir.XenViewNiche,"/[PointView_RDS]_TX5K_",TX,".rds"))
  }

  ########################
  ###  4. Xenium view (polygon)
  Func.Xen = function(Dim1, Res1){
    library(arrow)
    for(i in 1:NumOfSamples){
      TX = TXnumInteg[i]
      DF.Boundaries_0 = read_parquet(paste0("/Volumes/Extreme SSD/Data/TX5K_",
                                            case_when(TX %in% c("15","16") ~ "15_16",
                                                      TX %in% c("18","19","20","22") ~ "18_19_20_22",
                                                      TX %in% c("27","28") ~ "27_28",
                                                      TRUE ~ TX),
                                            "/cell_boundaries.parquet") )
      DF.CellNiche_0 = DF.KmeansRes %>% 
        rownames_to_column(var="cell_id") %>% 
        subset(subset = grepl(paste0("^TX5K",TX,"_"), cell_id)) %>% 
        dplyr::mutate(cell_id = str_remove(cell_id, pattern=paste0("TX5K",TX,"_"))) %>% 
        dplyr::mutate(niche = paste0("niche_",niche))
      DF.Boundaries_1 = DF.Boundaries_0 %>% 
        left_join(DF.CellNiche_0, by="cell_id") %>% 
        dplyr::mutate(niche2 = case_when(is.na(niche) ~ "Excluded",
                                        TRUE ~ niche))
      RangeX = c(min(subset(DF.Boundaries_1, !niche2=="Excluded")$vertex_x), 
                 max(subset(DF.Boundaries_1, !niche2=="Excluded")$vertex_x))
      RangeY = c(min(subset(DF.Boundaries_1, !niche2=="Excluded")$vertex_y), 
                 max(subset(DF.Boundaries_1, !niche2=="Excluded")$vertex_y))
      WideOrLong = case_when( diff(RangeY)/diff(RangeX) > 3/2 ~ c("WideOrLong"="Long", "width"=10, "height"=21),
                              diff(RangeY)/diff(RangeX) < 2/3 ~ c("WideOrLong"="Wide", "width"=21, "height"=10),
                              TRUE ~ c("WideOrLong"="Square", "width"=10, "height"=11 ) )
      DF.Boundaries_2 = DF.Boundaries_1　%>% 
        subset(vertex_x>=RangeX[1] & vertex_x<=RangeX[2] &
               vertex_y>=RangeY[1] & vertex_y<=RangeY[2])
      Xen.Niche.Polygon = 
      ggplot(DF.Boundaries_2, 
             aes(x=vertex_x, y=vertex_y, group=label_id)) +
        geom_polygon(aes(fill=niche2), color=NA) +
        labs(title=paste0("Xenium view of niches, TX5K_",TX),
             x=NULL, y=NULL,
             subtitle=CommonTitle.Sub) +
        scale_x_continuous(expand=expansion(mult=c(0, 0)),
                           breaks=seq(from=ceiling(RangeX[1]/1000) * 1000,
                                      to=floor(RangeX[2]/1000) * 1000,
                                      by=1000),
                           limits=c(RangeX[1],RangeX[2])) +
        scale_y_reverse(expand=expansion(mult=c(0, 0)),
                        breaks=seq(from=ceiling(RangeY[1]/1000) * 1000,
                                   to=floor(RangeY[2]/1000) * 1000,
                                   by=1000),
                        limits=c(RangeY[2],RangeY[1])) +
        scale_fill_manual(values=c(Vec.NicheColor,"Excluded"="gray40"),
                          label = c(Vec.NicheAnnot,"Excluded"="Excluded")) +
        coord_fixed(ratio=1) +
        theme(text = element_text(face = "bold"),
              plot.caption = element_text(size = 5),
              legend.position = "none",
              axis.ticks = element_line(color="gray50"),
              axis.line = element_blank(),
              panel.background = element_rect(fill = "transparent", colour = "gray50"),
              plot.background = element_rect(fill = "transparent", colour = NA),
              panel.grid = element_blank())
      Dir.XenViewNiche = paste0(DirInteg,"/NicheAnalysisData/XeniumView_Sub/",QCInfo.FileName,
                                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,
                                "/[",Unclassified,"Unclassified_Labeled",Labels,"_Neighbor",n.neighbors,
                                "_npc",pc_use,"_k",k.niche,"]")
      fs::dir_create(path=Dir.XenViewNiche)
      ggsave(plot=Xen.Niche.Polygon,
             file=paste0(Dir.XenViewNiche,"/[PolygonView]_TX5K_",TX,".png"), 
             width=as.numeric(WideOrLong[["width"]]), height=as.numeric(WideOrLong[["height"]]), dpi=700, bg="transparent")
      saveRDS(Xen.Niche.Polygon,
              file=paste0(Dir.XenViewNiche,"/[PolygonView_RDS]_TX5K_",TX,".rds") )
    }
  }

}

### 8-2b. counted in each sample ( old ... cannot use) ####
method="SpearmanAverage"  # or method="EuclideanWardd2"
Func.VisualizeNiche = function(Dim1, Res1, Dim2, Res2, n.neighbors, Unclassified, Labels, method){
  library(ggdendro)
  CommonTitle.Sub = paste0(NumOfSamples,"case Integration(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")\n",
                           QCInfo,"\n1stClust:PCA50vf2000_ClustDim",Dim1,"Res",formatC(Res1,digits=2,format="f"),
                           "\n2ndClust:PCA50vf2000_ClustDim",Dim2,"Res",formatC(Res2,digits=2,format="f"))
  # read and shape label data
  List.InitialClustData = readRDS(file=paste0(DirInteg,"[MetaDataList][InitialClust][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                                                "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".rds"))
  Add = List.InitialClustData$OrderedClustnumToLab
  names(Add) = paste0("Clust.",names(Add))
  List.InitialClustData[["ForPlotLab"]] = Add
  EliminateClusters = names(List.InitialClustData$OrderedClustnumToCelltype)[as.character(List.InitialClustData$OrderedClustnumToCelltype) %in% "Unclassified"]
  # load niche data
  DF.Niche_0 = read.csv(file=paste0(DirInteg,"NicheAnalysisData/[DataTable][NicheAnalysis_",Unclassified,"Unclassified_Labeled",Labels,"_Neighbor",n.neighbors,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                                    "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".csv"))# %>% 
               #subset(subset = !(LabelsForNicheAssay %in% paste0("Clust.", EliminateClusters)) )
  DF.ByNiche_0 = DF.Niche_0[ , c("niche", "LabelsForNicheAssay", "prop_byniche")]
  DF.ByNiche_1 = DF.ByNiche_0 %>% 
    pivot_wider(names_from=LabelsForNicheAssay,
                values_from=prop_byniche,
                values_fill=0) %>% 
    column_to_rownames("niche")
  if(method=="SpearmanAverage"){
    cor_niche = cor(t(DF.ByNiche_1), method="spearman", use="pairwise.complete.obs")
    d_niche = stats::as.dist(1 - cor_niche)
    Hclust_niche = stats::hclust(d_niche, method="average")
    MethodLabel = "Hierarchical clustering was performed using average linkage based on Spearman correlation distances."
  } else if (method=="EuclideanWardd2"){               # ward.D2なら距離はEuclideanで！！
    d_niche = dist(DF.ByNiche_1, method="euclidean")
    Hclust_niche = stats::hclust(d_niche, method = "ward.D2")
    MethodLabel = "Hierarchical clustering was performed using Ward’s method (ward.D2) with Euclidean distances."
  }
  nicheorder = (Hclust_niche$labels)[Hclust_niche$order]
  NumOfCAFsubclusts = grepl(pattern="^CAF-", unique(DF.Niche_0$LabelsForNicheAssay)) %>% 
                      sum() # How many clusters are named started with string "CAF"
  Vec.OrderedClustnum = c(paste0("CAF-",0:(NumOfCAFsubclusts-1)),
                          paste0("Clust.",names(List.InitialClustData$OrderedClustnumToLab) ) )
  DF.Niche_0$niche = factor(DF.Niche_0$niche, levels=nicheorder)
  DF.Niche_0$LabelsForNicheAssay = factor(DF.Niche_0$LabelsForNicheAssay, levels=Vec.OrderedClustnum)
  # dendrogram
  dend_niche = as.dendrogram(Hclust_niche)
  dend_data = dendro_data(dend_niche)
  p_dend = 
    ggplot(dend_data$segments) +
    geom_segment(aes(x=y, y=x, xend=yend, yend=xend)) +
    scale_y_continuous(limits = c(0.5, length(Hclust_niche$labels) + 0.5),
                       breaks = seq_along(Hclust_niche$labels),
                       labels = Hclust_niche$labels,
                       expand = expansion(mult=c(0,0))) +
    scale_x_reverse() +
    theme_void()
  # heatmap
  CommonScaleX = scale_x_discrete(expand=expansion(mult=c(0,0)),
                            labels=List.InitialClustData$ForPlotLab)
  CommonScaleY = scale_y_discrete(expand=expansion(mult=c(0,0)),
                            position="right")
  CommonScaleFill = scale_fill_viridis_c(option="magma",
                                         trans="sqrt",
                                         guide=guide_colorbar(
                                           #direction="horizontal",
                                           #title.position="top",
                                           #title.hjust=0.5,
                                           frame.colour="black",
                                           ticks.colour="black"))
  CommonTheme = theme(plot.margin = margin(t=0, r=0.3, b=0, l=0, unit="cm"),
                      legend.position="bottom",
                      axis.text = element_text(face = "bold"),
                      axis.text.x = element_text(angle=45, hjust=1),
                      axis.title = element_text(face = "bold"),
                      legend.text = element_text(face = "bold", angle=45, hjust=1),
                      legend.title = element_text(face = "bold"),
                      panel.border = element_rect(fill = NA, color = "black"))
  # for annotating niche
  p_heat1 =
    ggplot(DF.Niche_0, aes(x=LabelsForNicheAssay, y=niche, fill=prop_byniche)) +
    geom_tile(color="gray10") + 
    labs(x=NULL, y=NULL, fill="Cell proportion\nin niche",
         caption=paste0("Niches were hierarchically clustered based on the proportions of cell clusters composing each niche.")) +
    scale_fill_viridis_c(option="magma",
                         trans="sqrt",
                         guide=guide_colorbar(
                           #direction="horizontal",
                           #title.position="top",
                           #title.hjust=0.5,
                           frame.colour="black",
                           ticks.colour="black")) +
    CommonScaleX + CommonScaleY + CommonTheme
  # for annotating CAF
  p_heat2 =
  ggplot(DF.Niche_0, aes(x=LabelsForNicheAssay, y=niche, fill=prop_byclust)) +
    geom_tile(color="gray10") + 
    labs(x=NULL, y=NULL, fill="Cell proportion\nin cluster",
         caption=paste0("Heatmap shows the distribution of each cluster across niches,\n",
                        "normalized by the total number of cells in each cluster within each sample.")) +
    scale_fill_viridis_c(option="mako",
                         #trans="sqrt",
                         guide=guide_colorbar(
                           #direction="horizontal",
                           #title.position="top",
                           #title.hjust=0.5,
                           frame.colour="black",
                           ticks.colour="black")) +
    CommonScaleX + CommonScaleY + CommonTheme
  # save
  ggsave(plot = p_dend + p_heat1 + p_heat2 +
                plot_layout(widths = c(0.5, 4, 4)) +
                plot_annotation(
                    subtitle=paste0(CommonTitle.Sub,
                                    "\nN.neighbors = ",n.neighbors),
                    caption=paste0(MethodLabel,
                                   "\nUnclassified clusters (",paste(EliminateClusters, collapse=","), ") were eliminated.") ),
         file=paste0(DirInteg,"[Figure][NicheAnalysis_",Unclassified,"Unclassified_Labeled",Labels,"_Neighbor",n.neighbors,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                     "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"),
         height=8, width=20, dpi=400)

}
Func.VisualizeNiche(Dim1=20, Res1=1.0, Dim2=30, Res2=0.5, n.neighbors=40, method="SpearmanAverage")
Func.VisualizeNiche(Dim1=20, Res1=1.0, Dim2=30, Res2=0.5, n.neighbors=35, method="SpearmanAverage")
Func.VisualizeNiche(Dim1=20, Res1=1.0, Dim2=30, Res2=0.5, n.neighbors=30, Unclassified="Retain", Labels="seurat_clusters", method="SpearmanAverage")
Func.VisualizeNiche(Dim1=20, Res1=1.0, Dim2=30, Res2=0.5, n.neighbors=20, method="SpearmanAverage")

Func.XenViewNiche = function(Dim1, Res1, Dim2, Res2, n.neighbors, Unclassified, Labels){
List.InitialClustData = readRDS(file=paste0(DirInteg,"[MetaDataList][InitialClust][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                                            "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".rds"))
EliminateClusters = names(List.InitialClustData$OrderedClustnumToCelltype)[as.character(List.InitialClustData$OrderedClustnumToCelltype) %in% "Unclassified"]
# load niche data
DF.Niche_0 = read.csv(file=paste0(DirInteg,"NicheAnalysisData/[DataTable][NicheAnalysis_",Unclassified,"Unclassified_Labeled",Labels,"_Neighbor",n.neighbors,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                                  "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".csv"))
DF.Niche_1 = DF.Niche_0# %>% 
  #subset(subset = !(seurat_clusters %in% paste0("Clust.", EliminateClusters)) )
for(i in 1:NumOfSamples){
TX=TXnumInteg[i]
DF.NicheRes = read.csv(
  file=paste0(DirInteg,"NicheAnalysisData/[DataTable][NicheAnalysis_",Unclassified,"Unclassified_Labeled",Labels,"_Neighbor",n.neighbors,"_TX5K",TX,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
              "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".csv"))
DF.IDandNiche_0 = List.InitialClustData$MetaData %>% 
     rename(full_id = cell_id) %>% 
     dplyr::mutate(cell_id = str_remove(full_id, pattern="TX5K.*_")) %>% 
     subset(sample==TX)
DF.IDandNiche_1 = inner_join(DF.IDandNiche_0, 
                             DF.NicheRes, 
                             by="cell_id")
DF.IDandNiche_2 = DF.IDandNiche_1 %>% 
  subset(subset = !(seurat_clusters %in% paste0("Clust.", EliminateClusters)) ) %>% 
  dplyr::mutate(niches = as.character(niches))
RangeX = c(min(DF.IDandNiche_2$CoordX), max(DF.IDandNiche_2$CoordX))
RangeY = c(min(DF.IDandNiche_2$CoordY), max(DF.IDandNiche_2$CoordY))
WideOrLong = case_when( diff(RangeY)/diff(RangeX) > 3/2 ~ c("WideOrLong"="Long", "width"=5, "height"=10.5),
                        diff(RangeY)/diff(RangeX) < 2/3 ~ c("WideOrLong"="Wide", "width"=10.5, "height"=5),
                        TRUE ~ c("WideOrLong"="Square", "width"=5, "height"=5.5 ) )
Xen.Niche = 
ggplot(DF.IDandNiche_2, aes(x=CoordX, y=CoordY, color=niches)) +
  geom_point(size=0.05) +
  labs(subtitle=paste0("Distribution of niches, TX5K_",TX),
       x=NULL, y=NULL, 
       caption=paste0("N.neighbors=",n.neighbors,"cells") ) +
  scale_x_continuous(expand=expansion(mult=c(0, 0)),
                     breaks=seq(from=ceiling(RangeX[1]/1000) * 1000,
                                to=floor(RangeX[2]/1000) * 1000,
                                by=1000),
                     limits=c(RangeX[1],RangeX[2])) +
  scale_y_reverse(expand=expansion(mult=c(0, 0)),
                  breaks=seq(from=ceiling(RangeY[1]/1000) * 1000,
                             to=floor(RangeY[2]/1000) * 1000,
                             by=1000),
                  limits=c(RangeY[2],RangeY[1])) +
  coord_fixed(ratio=1) + 
  theme(text = element_text(face = "bold"),
        plot.caption = element_text(size = 5),
        legend.position = "right",
        legend.key = element_blank(),
        axis.ticks = element_line(color="gray50"),
        axis.line = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = "gray50"),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.grid = element_blank())
fs::dir_create(path=paste0(DirInteg,"/XeniumView_Sub/",QCInfo.FileName,
                           "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,
                           "/Niche"))
ggsave(plot=Xen.Niche,
       file=paste0(DirInteg,"/XeniumView_Sub/",QCInfo.FileName,
                   "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,
                   "/Niche/[Figure][XenView_Niche_TX5K_",TX,"_",Unclassified,"Unclassified_Labeled",Labels,"_neighbor",n.neighbors,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")].png"), 
       width=as.numeric(WideOrLong[["width"]]), height=as.numeric(WideOrLong[["height"]]), dpi=400, bg="transparent")
}
}
Func.XenViewNiche(Dim1=20, Res1=1.0, Dim2=30, Res2=0.5, n.neighbors=40)
Func.XenViewNiche(Dim1=20, Res1=1.0, Dim2=30, Res2=0.5, n.neighbors=35)
Func.XenViewNiche(Dim1=20, Res1=1.0, Dim2=30, Res2=0.5, n.neighbors=30)
Func.XenViewNiche(Dim1=20, Res1=1.0, Dim2=30, Res2=0.5, n.neighbors=20)

  
  
### 9. Spatial distribution だめ ####
function(){
library(FNN)
TX=TXnumInteg[4]
TX=TXnumInteg[5]
TX=TXnumInteg[6]
TX=TXnumInteg[1]
#SeuObj.CAF_0 = readRDS(paste0(DirInteg,"Objects/[SeuObj_ExtractedCAF][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
#                              QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2,".rds") )
DF.CellGroupTable_0 = read.csv(
  file=paste0(DirInteg,"/CellGroupTable/[SubClust_CellGroupTable_TX5K_",TX,"][",
              NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
              "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2,".csv") )
NumObSubclust = length(unique(DF.CellGroupTable_0$group))
CommonCaption.Sub = paste0("Total ",format(nrow(DF.CellGroupTable_0), big.mark=",", scientific=F), " CAFs, ",NumObSubclust, " clusters.")
CommonTitle.Sub = paste0(NumOfSamples,"case Integration(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")\n",
                         QCInfo,"\n1stClust:PCA50vf2000_ClustDim",Dim1,"Res",formatC(Res1,digits=2,format="f"),
                         "\n2ndClust:PCA50vf2000_ClustDim",Dim2,"Res",formatC(Res2,digits=2,format="f"))
Func.GSVAScoreRanking = function(ScoreColmn){case_when(ScoreColmn < quantile(ScoreColmn, probs=seq(0,1,0.2), names=F)[2] ~ "Rank1",
                                                       ScoreColmn < quantile(ScoreColmn, probs=seq(0,1,0.2), names=F)[3] ~ "Rank2",
                                                       ScoreColmn < quantile(ScoreColmn, probs=seq(0,1,0.2), names=F)[4] ~ "Rank3",
                                                       ScoreColmn < quantile(ScoreColmn, probs=seq(0,1,0.2), names=F)[5] ~ "Rank4",
                                                       TRUE ~ "Rank5")}
DF.ssGSVAres_0 = 
  read.csv(file=paste0(DirInteg,"/[DataTable_ExtractedCAF_ssGSVA][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                       QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".csv"))
DF.ssGSVAres_1 = dplyr::mutate(DF.ssGSVAres_0,
                               "BuffaOrig_rank"=Func.GSVAScoreRanking(DF.ssGSVAres_0$BuffaOrig_ssGSVA),
                               "WinterOrig_rank"=Func.GSVAScoreRanking(DF.ssGSVAres_0$WinterOrig_ssGSVA),
                               "HALLMARK_HYPOXIA_rank"=Func.GSVAScoreRanking(DF.ssGSVAres_0$HALLMARK_HYPOXIA_ssGSVA))
DF.ssGSVAres_2 = dplyr::select(DF.ssGSVAres_1,
                               c(sample, cell_id, X, Y, 
                                 BuffaOrig_ssGSVA, WinterOrig_ssGSVA, HALLMARK_HYPOXIA_ssGSVA,
                                 BuffaOrig_rank, WinterOrig_rank, HALLMARK_HYPOXIA_rank))
DF.ssGSVAres_3 = inner_join(DF.CellGroupTable_0,
                             DF.ssGSVAres_2,
                             by="cell_id")
DF.ssGSVAres_4 = subset(DF.ssGSVAres_3, subset=sample==paste0("TX5K",TX))
nn = FNN::get.knn(DF.ssGSVAres_4[,c("X","Y")], 
                  k=1000)      
isRank5 = DF.ssGSVAres_4$WinterOrig_rank %in% c("Rank5","Rank4")
Vec.Rank5_frac = numeric(nrow(DF.ssGSVAres_4)) 
for (i in seq_len(nrow(DF.ssGSVAres_4))) {
  index_0 = nn$nn.index[i, ]
  dist = nn$nn.dist[i, ]
  index_1 = index_0[dist <= 300]   # dist : smoothing scale (um)
  Vec.Rank5_frac[i] = 
    if (length(index_1)==0) { 0
  } else { mean(isRank5[index_1])}
}
summary(Vec.Rank5_frac)
hist(Vec.Rank5_frac, breaks=50) #その細胞の周囲 300µm 以内にいる細胞のうち、低酸素スコア上位 25%（Q4）が占める割合
DF.ssGSVAres_4$Rank5_frac_r300 = as.numeric(Vec.Rank5_frac)

Percentile_High = 0.85
Threshold = quantile(Vec.Rank5_frac, Percentile_High, na.rm=TRUE, names=FALSE)
DF.ssGSVAres_4$HypoxiaHigh = Vec.Rank5_frac >= Threshold
sum(DF.ssGSVAres_4$HypoxiaHigh)        # ← 必ず 0 以外になる
range(DF.ssGSVAres_4$Rank5_frac_r300)

coords <- as.matrix(DF.ssGSVAres_4[, c("X","Y")])
HypoxiaHigh_index <- which(DF.ssGSVAres_4$HypoxiaHigh)
coords_HypoxiaHigh <- coords[HypoxiaHigh_index, , drop = FALSE]

link_um <- 125      # ★50だと細切れになるのでまず75を推奨
k2 <- 30

DF.ssGSVAres_4$HypoxiaPatch <- NA_integer_

if (nrow(coords_HypoxiaHigh) >= 2) {
  
  nn2 <- FNN::get.knn(coords_HypoxiaHigh,
                      k = min(k2, nrow(coords_HypoxiaHigh)-1))
  
  edges <- do.call(rbind, lapply(seq_len(nrow(coords_HypoxiaHigh)), function(i) {
    j <- nn2$nn.index[i, ]
    d <- nn2$nn.dist[i, ]
    j <- j[d <= link_um]
    if (length(j) == 0) return(NULL)
    cbind(i, j)
  }))
  
  if (!is.null(edges) && nrow(edges) > 0) {
    g <- igraph::graph_from_edgelist(edges, directed = FALSE)
    comp <- igraph::components(g)$membership
    DF.ssGSVAres_4$HypoxiaPatch[HypoxiaHigh_index] <- comp
  }
}

table(DF.ssGSVAres_4$HypoxiaPatch, useNA = "ifany")

patch_sizes <- table(DF.ssGSVAres_4$HypoxiaPatch)

min_patch_size <- 50
keep_patches <- as.integer(names(patch_sizes[patch_sizes >= min_patch_size]))

DF.ssGSVAres_4$HypoxiaNiche <- ifelse(
  DF.ssGSVAres_4$HypoxiaPatch %in% keep_patches,
  DF.ssGSVAres_4$HypoxiaPatch,
  NA_integer_
)

table(DF.ssGSVAres_4$HypoxiaNiche, useNA="ifany")

isCAF8 <- DF.ssGSVAres_4$group == "CAF-8"
tab_niche <- table(CAF8 = isCAF8,
                   InNiche = !is.na(DF.ssGSVAres_4$HypoxiaNiche))
tab_niche
fisher.test(tab_niche)

ggplot(DF.ssGSVAres_4, aes(x = X, y = Y)) +
  geom_point(color = "gray90", size = 0.15) +
  geom_point(
    data = subset(DF.ssGSVAres_4, !is.na(HypoxiaNiche)),
    aes(color = factor(HypoxiaNiche)),
    size = 0.35
  ) +
  geom_point(
    data = subset(DF.ssGSVAres_4, group == "CAF-8"),
    shape = 21, stroke = 0.35, size = 0.8,
    fill = NA, color = "black"
  ) +
  scale_y_reverse() +
  coord_fixed() +
  theme_classic() +
  labs(color = "Hypoxia niche", x = "X (µm)", y = "Y (µm)", title=paste0("TX5K",TX))
}

### 10. Spatial prolif-CAF distribution merged onto tiled hypoxia heatmap ####
Dim1=20
Res1=1.0
Dim2=30
Res2=0.5
Num.GridLen = 200
Num.MinCells = 10
prolifCAFclust="8"
Fun.XenHypoxiaGrid = function(Dim1, Res1, Dim2, Res2, prolifCAFclust, Num.GridLen, Num.MinCells){
  Plot = function(TX){
    DF.CellGroupTable_0 = read.csv(
      file=paste0(DirInteg,"/CellGroupTable/[SubClust_CellGroupTable_TX5K_",TX,"][",
                  NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                  "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2,".csv") )
    NumOfSubclust = length(unique(DF.CellGroupTable_0$group))
    CommonCaption.Sub = paste0("Total ",format(nrow(DF.CellGroupTable_0), big.mark=",", scientific=F), " CAFs, ",NumOfSubclust, " clusters.")
    CommonTitle.Sub = paste0(NumOfSamples,"case Integration(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")\n",
                             QCInfo,"\n1stClust:PCA50vf2000_ClustDim",Dim1,"Res",formatC(Res1,digits=2,format="f"),
                             "\n2ndClust:PCA50vf2000_ClustDim",Dim2,"Res",formatC(Res2,digits=2,format="f"))
    DF.ssGSVAres_0 = 
      read.csv(file=paste0(DirInteg,"/[DataTable_ExtractedCAF_ssGSVA][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                           QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".csv"))
    DF.ssGSVAres_1 = DF.ssGSVAres_0 %>% 
      dplyr::mutate("BuffaOrig_z"=as.numeric(scale(BuffaOrig_ssGSVA)),
                    "WinterOrig_z"=as.numeric(scale(WinterOrig_ssGSVA)),
                    "HALLMARK_HYPOXIA_z"=as.numeric(scale(HALLMARK_HYPOXIA_ssGSVA)) )
    DF.ssGSVAres_2 = DF.ssGSVAres_1 %>% 
      dplyr::select(c(sample, cell_id, X, Y, 
                      BuffaOrig_z, BuffaOrig_ssGSVA, 
                      WinterOrig_z, WinterOrig_ssGSVA, 
                      HALLMARK_HYPOXIA_z, HALLMARK_HYPOXIA_ssGSVA) )
    DF.ssGSVAres_3 = inner_join(DF.CellGroupTable_0,
                                DF.ssGSVAres_2,
                                by="cell_id")
    DF.ssGSVAres_4 = subset(DF.ssGSVAres_3, 
                            subset=sample==paste0("TX5K",TX))
    DF.Grid = DF.ssGSVAres_4 %>% 
      dplyr::mutate(
        GridX.edge = floor(X / Num.GridLen),
        GridX.center = GridX.edge + 0.5,
        GridY.edge = floor(Y / Num.GridLen),
        GridY.center = GridY.edge + 0.5,
        Grid_id = paste(GridX.edge, GridY.edge, sep="_"))
    DF.Grid_sum_0 = 
      DF.Grid %>% 
      group_by(Grid_id, GridX.center, GridY.center) %>% 
      summarise(
        n_cells = n(),
        mean_local_hypoxia_z = mean(WinterOrig_z),
        mean_local_hypoxia = mean(WinterOrig_ssGSVA),
        med_local_hypoxia_z = median(WinterOrig_z),
        med_local_hypoxia = median(WinterOrig_ssGSVA),
        n_prolifCAF = sum(group == paste0("CAF-", prolifCAFclust), na.rm=TRUE),
        frac_caf8 = n_prolifCAF / n_cells,
        .groups = "drop") %>% 
      dplyr::mutate(
        mean_local_hypoxia_z = ifelse(n_cells < Num.MinCells, NA, mean_local_hypoxia_z),
        frac_caf8 = ifelse(n_cells < Num.MinCells, NA, frac_caf8) )
    RangeGridX = range(DF.Grid$GridX.edge)
    RangeGridY = range(DF.Grid$GridY.edge)
    Vec.AllGridId = as.vector(outer(seq(RangeGridX[1], RangeGridX[2]), 
                                    seq(RangeGridY[1], RangeGridY[2]), paste, sep="_"))
    Vec.DroppedGridIds = setdiff(Vec.AllGridId, unique(DF.Grid_sum_0$Grid_id))
    DF.naGrid = 
      DF.Grid_sum_0[rep(NA_integer_, length(Vec.DroppedGridIds)), ] %>% 
      dplyr::mutate(Grid_id = Vec.DroppedGridIds,
                    GridX.center = str_remove(Grid_id, pattern="_.*"),
                    GridX.center = as.numeric(GridX.center)+0.5,
                    GridY.center = str_remove(Grid_id, pattern=".*_"),
                    GridY.center = as.numeric(GridY.center)+0.5)
    DF.Grid_sum_1 = 
      rbind(DF.Grid_sum_0,DF.naGrid)
    DF.ssGSVAres_5 = DF.ssGSVAres_4 %>% 
      dplyr::mutate(GridX = X/Num.GridLen,
                    GridY = Y/Num.GridLen)
    RangeX = c(min(DF.ssGSVAres_4$X), max(DF.ssGSVAres_4$X))
    RangeY = c(min(DF.ssGSVAres_4$Y), max(DF.ssGSVAres_4$Y))
    WideOrLong = case_when(
      diff(RangeY)/diff(RangeX) > 4/3 ~ c("WorL"="Long", "width"=10, "height"=21, 
                                          "colbar"="vertical", "LgdTitlePos"="left", "LdgPos"="right",
                                          "LgdTitleAngle"=90, "LdgTextAngle"=0),
      diff(RangeY)/diff(RangeX) < 3/4 ~ c("WorL"="Wide", "width"=21, "height"=10,
                                          "colbar"="horizontal", "LgdTitlePos"="top", "LdgPos"="bottom",
                                          "LgdTitleAngle"=0, "LdgTextAngle"=45),
      TRUE ~ c("WorL"="Square", "width"=10, "height"=11, 
               "colbar"="vertical", "LgdTitlePos"="left", "LdgPos"="right",
               "LgdTitleAngle"=90, "LdgTextAngle"=0 ) )
  Xen.HypoxiaHeatmap = 
    ggplot() +
      geom_tile(
        data=DF.Grid_sum_1,
        aes(x=GridX.center, y=GridY.center, fill=mean_local_hypoxia_z),
        color = "gray60") + 
      geom_point(
        data=subset(DF.ssGSVAres_5, group == paste0("CAF-", prolifCAFclust)),
        aes(x=GridX, y=GridY), color="black",
        size=3.0, shape=21, stroke=1.5, fill=NA) +
      labs(title=paste0("TX5K_", TX),
        fill="Local mean of Winter's Hypoxia (ssGSVA)", x=NULL, y=NULL) +
      scale_x_continuous(expand=expansion(mult=c(0,0)),
                         labels = function(x) x*Num.GridLen) +
      scale_y_reverse(expand=expansion(mult=c(0,0)),
                      labels = function(y) y*Num.GridLen) + 
      coord_cartesian(xlim = RangeX/Num.GridLen, ylim = RangeY/Num.GridLen,
                      ratio=1) +
      scale_fill_gradientn(colors=c("#bd0026", "#f03b20", "#fd8d3c", "#fecc5c", "gray95",
                                    "#a1dab4", "#41b6c4", "#2c7fb8", "#253494"),
                           na.value="gray80",
                           limits=c(-max(abs(DF.Grid_sum_1$mean_local_hypoxia_z), na.rm=TRUE), 
                                    max(abs(DF.Grid_sum_1$mean_local_hypoxia_z), na.rm=TRUE)),
                           guide=guide_colorbar(
                             direction=WideOrLong[["colbar"]],
                             title.position=WideOrLong[["LgdTitlePos"]],
                             title.hjust=0.5,
                             title.vjust=0.5,
                             frame.colour="black",
                             ticks.colour="black")) +
      theme(legend.position = WideOrLong[["LdgPos"]],
            legend.background = element_rect(fill = "transparent", color = NA),
            text = element_text(size = 25, face = "bold"),
            legend.text = element_text(angle = as.numeric(WideOrLong[["LdgTextAngle"]]), hjust = 1),
            legend.title = element_text(angle = as.numeric(WideOrLong[["LgdTitleAngle"]])),
            plot.background = element_rect(fill = "transparent", color = NA),
            panel.background = element_rect(fill = "green", color = NA),
            panel.grid = element_blank(),
            panel.border = element_rect(fill = "transparent", color = "black") )
    fs::dir_create(path=paste0(DirInteg,"/XeniumView_Sub/",QCInfo.FileName,
                               "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,
                               "/ssGSVA"))
    ggsave(plot=Xen.HypoxiaHeatmap,
           file=paste0(DirInteg,"/XeniumView_Sub/",QCInfo.FileName,
                       "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,
                       "/ssGSVA/[Figure][ExtractedCAF_XenView_GSVA_HypoxicAreaAndProlifCAF_TX5K_",TX,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                       "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".png"), 
           width=as.numeric(WideOrLong[["width"]]), height=as.numeric(WideOrLong[["height"]]), 
           dpi=400, limitsize=FALSE, bg="transparent")
    saveRDS(Xen.HypoxiaHeatmap,
            file=paste0(DirInteg,"/XeniumView_Sub/",QCInfo.FileName,
                        "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,
                        "/ssGSVA/[Figure][ExtractedCAF_XenView_GSVA_HypoxicAreaAndProlifCAF_TX5K_",TX,"][",NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                        "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,".rds"))
    Vec.HypoxiaScoreGrid <- 
      setNames(DF.Grid_sum_1$mean_local_hypoxia_z, DF.Grid_sum_1$Grid_id)
    set.seed(1234)
    DF <- 
      DF.Grid %>% 
      dplyr::mutate(GridScore = Vec.HypoxiaScoreGrid[Grid_id],
                    GridScore_bin = cut(GridScore, breaks = 10)) %>% 
      dplyr::slice_sample(n = nrow(.)) %>% 
      dplyr::filter(!is.na(GridScore))
    ggplot(DF, aes(y = group)) +
      geom_bar(stat = "count", position = "fill",
               aes(fill = GridScore_bin)) +
      scale_fill_brewer(palette = "RdBu")
  ggplot(DF, aes(x=X, y=Y)) +
    geom_point(aes(color=WinterOrig_z)) +
    geom_point(
      data=subset(DF.ssGSVAres_5, group == paste0("CAF-", 7)),
      color="black",
      #size=3.0, 
      shape=21, stroke=1.5, fill=NA) +
    scale_y_reverse() +
    scale_color_gradientn(colors = c("red4","red", "orange","gray99", "skyblue", "blue", "blue4"),
                          limits = c(-max(abs(DF$WinterOrig_z)), max(abs(DF$WinterOrig_z)))) +
    coord_fixed(ratio = 1)
  Vec.CAForder = paste0("CAF-", c(4,8,0,1,2,5,3,7,6))
  DF.Barcode <- 
    DF %>% 
    dplyr::mutate(group = factor(group, levels = Vec.CAForder)) %>% 
    group_by(group) %>% 
    arrange(GridScore, .by_group = TRUE) %>% 
    dplyr::mutate(
      cell_order = row_number(),
      n_in_group = n(),
      xmin = (cell_order - 1) / n_in_group,
      xmax = cell_order / n_in_group,
      group_num = as.numeric(group),
      ymin = group_num - 0.45,
      ymax = group_num + 0.45,
      GridHypoxiaRank = case_when(
        GridScore > quantile(DF$GridScore, 0.9) ~ "R1",
        GridScore > quantile(DF$GridScore, 0.8) ~ "R2",
        GridScore > quantile(DF$GridScore, 0.7) ~ "R3",
        GridScore > quantile(DF$GridScore, 0.6) ~ "R4",
        GridScore > quantile(DF$GridScore, 0.5) ~ "R5",
        GridScore > quantile(DF$GridScore, 0.4) ~ "R6",
        GridScore > quantile(DF$GridScore, 0.3) ~ "R7",
        GridScore > quantile(DF$GridScore, 0.2) ~ "R8",
        GridScore > quantile(DF$GridScore, 0.1) ~ "R9",
        TRUE ~"R10"),
      GridHypoxiaRank = factor(GridHypoxiaRank,
                               levels = paste0("R",1:10))) %>% 
    ungroup()
  range(DF.Barcode$cell_order)
  ggplot(DF.Barcode) +
    geom_rect(aes(xmin = xmin, xmax = xmax,
                  ymin = ymin, ymax = ymax,
                  fill = GridScore)) +
    scale_fill_gradientn(
      colors=c("#bd0026", "#f03b20", "#fd8d3c", "#fecc5c", "gray95",
               "#a1dab4", "#41b6c4", "#2c7fb8", "#253494"),
      limits=c(-max(abs(DF.Barcode$GridScore), na.rm = TRUE), 
               max(abs(DF.Barcode$GridScore), na.rm = TRUE))) +
    scale_y_continuous(
      breaks = seq_along(Vec.CAForder),
      labels = Vec.CAForder) +
    scale_x_continuous(
      limits = c(0, 1),
      expand = expansion(mult = c(0, 0))) +
    theme(
      plot.background = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA),
      panel.border = element_rect(fill = "transparent", color = "black"),
      panel.grid = element_blank())+
  
  ggplot(DF.Barcode) +
    geom_rect(aes(xmin = xmin, xmax = xmax,
                  ymin = ymin, ymax = ymax,
                  fill = GridHypoxiaRank)) +
    scale_y_continuous(
      breaks = seq_along(Vec.CAForder),
      labels = Vec.CAForder) +
    scale_fill_manual(
      values=rev(c("#bd0026", "#f03b20", "#fd8d3c", "#fecc5c",
                   "gray95", "gray95",
                   "#a1dab4", "#41b6c4", "#2c7fb8", "#253494")))
    
  }
  for (i in 1:NumOfSamples){ Plot(TX=TXnumInteg[i]) }
}
Fun.XenHypoxiaGrid(Dim1=20, Res1=1.0, Dim2=30, Res2=0.5, prolifCAFclust="8", Num.GridLen=200, Num.MinCells=10)

### 11. Pseudo-time analysis  ####
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

DataSetName = "WholeCAF"
Seu.CAF_0 = readRDS(paste0(DirInteg,"Objects/[SeuObj_ExtractedCAF][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                           QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2,".rds") )
DataSetName = "rmCAF6"
Seu.CAF_rmCAF6 <- 
  Seu.CAF_0 %>% 
  subset(seurat_clusters != 6)

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







#Because a distinct normal fibroblast cluster was not identified, trajectory rooting was guided using spatial niche information. CAFs localized within the normal acinar niche were highly enriched in State 7, which was therefore selected as the root state for pseudotime ordering.





library("slingshot")
library("monocle3")
library("SeuratWrappers")
Seu.CAF_0 = readRDS(paste0(DirInteg,"Objects/[SeuObj_ExtractedCAF][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                           QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2,".rds") )
CommonTitle.Sub = paste0(NumOfSamples,"case Integration(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")\n",
                         QCInfo,"\n1stClust:PCA50vf2000_ClustDim",Dim1,"Res",formatC(Res1,digits=2,format="f"),
                         "\n2ndClust:PCA50vf2000_ClustDim",Dim2,"Res",formatC(Res2,digits=2,format="f"))
NumOfSublust = length(unique(Seu.CAF_0$seurat_clusters))
#CommonCaption.Sub = paste0("Total ",format(ncol(Seu.CAF_0), big.mark=",", scientific=F), " CAFs, ",NumOfSublust, " clusters.")
Seu.CAF_0$seurat_clusters = factor(Seu.CAF_0$seurat_clusters, levels=c(4,8,0,1,2,5,3,7,6))
Idents(Seu.CAF_0) = "seurat_clusters"

Seu.CAF_0 <- JoinLayers(Seu.CAF_0)
sceobj_0 <- as.SingleCellExperiment(Seu.CAF_0)
sceobj_1 <- slingshot(
  sceobj_0,
  clusterLabels = "seurat_clusters",
  reducedDim = "UMAP")
saveRDS(sceobj_1, "SceObj_wholeSeuObj.rds")
plot(reducedDims(sceobj_1)$UMAP,
     col = as.numeric(colData(sceobj_1)$seurat_clusters),
     pch = 16,
     asp = 1)
lines(SlingshotDataSet(sceobj_1),
      lwd = 3,
      col = "black")
pt <- slingPseudotime(sceobj_1)[, 1]
plot(reducedDims(sceobj_1)$UMAP,
     col = hcl.colors(100, "viridis")[cut(pt, 100)],
     pch = 16,
     asp = 1)
lines(SlingshotDataSet(sceobj_1),
      lwd = 3,
      col = "black")
# input embeddings direct
MT.UmapEmbeddings <- Embeddings(Seu.CAF_0, "umap")
Vec.ClustLabs <- Idents(Seu.CAF_0)
sceobj_modified <- slingshot(
  MT.UmapEmbeddings,
  clusterLabels = Vec.ClustLabs)
saveRDS(sceobj_modified, "SceObj_onlyEmbeddings.rds")
plot(MT.UmapEmbeddings,
     col = as.numeric(Vec.ClustLabs),
     pch = 16,
     asp = 1)
lines(SlingshotDataSet(sceobj_modified),
      lwd = 3,
      col = "black")
length(slingLineages(sceobj_modified))
slingLineages(sceobj_modified)
# define end clust ("8")
sceobj_modified_end8 <- slingshot(
  MT.UmapEmbeddings,
  clusterLabels = Vec.ClustLabs,
  end.clus = "8")
saveRDS(sceobj_modified_end8, "SceObj_onlyEmbeddings_endclust8.rds")
sceobj_modified_end8 <- readRDS("SceObj_onlyEmbeddings_endclust8.rds")
plot(MT.UmapEmbeddings,
     col = as.numeric(Vec.ClustLabs),
     pch = 16,
     asp = 1)
lines(SlingshotDataSet(sceobj_modified_end8),
      lwd = 3,
      col = "black")
length(slingLineages(sceobj_modified_end8))
slingLineages(sceobj_modified_end8)
pt <- slingPseudotime(sceobj_modified_end8)

summary(pt)
pt1 <- slingPseudotime(sceobj_modified_end8)[,1]

plot(MT.UmapEmbeddings,
     col = hcl.colors(100, "viridis")[cut(pt1, 100)],
     pch = 16,
     asp = 1)

lines(SlingshotDataSet(sceobj_modified_end8),
      lwd = 3)
# separately drawn
par(mfrow = c(1, 2),
    
    mar = c(4,4,2,1))

for(i in 1:2){
  
  pti <- pt[, i]
  
  cols <- hcl.colors(100, "viridis")
  
  pt_cut <- cut(pti, breaks = 100)
  
  plot(MT.UmapEmbeddings,
       
       col = ifelse(is.na(pti), "grey90", cols[pt_cut]),
       
       pch = 16,
       
       asp = 1,
       
       main = paste0("Lineage ", i))
  
  lines(SlingshotDataSet(sceobj_modified_end8),
        
        lwd = 3,
        
        col = "black")
  
}

par(mfrow = c(1,1))

DF.pt <- data.frame(
  
  cluster = as.character(Vec.ClustLabs),
  
  pt1 = pt[, 1],
  
  pt2 = pt[, 2]
  
)

aggregate(pt1 ~ cluster, DF.pt, median, na.rm = TRUE)

aggregate(pt2 ~ cluster, DF.pt, median, na.rm = TRUE)

slingLineages(sceobj_modified_end8)

# Reduction PCA
MT.PCAEmbeddings <- Embeddings(Seu.CAF_0, "pca")[,1:20]
sds_pca <- slingshot(
  MT.PCAEmbeddings,
  clusterLabels = Vec.ClustLabs,
  end.clus = "8")



###############################################################
###
###    DEG analysis
###  1. Check DEGs(FindAllMarkers) ####
Func.FindAllMarkers = function(Dim1, Res1, Dim2, Res2){
  library(presto)
  # read seurat obj
  SeuObj_0 = readRDS(paste0(DirInteg,"Objects/[SeuObj_ExtractedCAF][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                          QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2,".rds") )
  SeuObj_1 = JoinLayers(SeuObj_0)
  # check DEGs and save
  set.seed(123)
  df_0 = FindAllMarkers(SeuObj_1, slot="data", only.pos=FALSE,
                        min.pct = 0.00, min.cell.feature = 0, min.cells.group = 0, 
                        logfc.threshold = 0, return.thresh = 1.00,
                        test.use = "wilcox")
  library("fs")
  fs::dir_create(path=paste0(DirInteg,"/DEG_table"))
  write.csv(df_0, file = paste0(DirInteg, "DEG_table/[SubclustsOfCAFs_DEGs_All][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                                QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2,".csv"), row.names = T) 
  df_1 = subset(df_0, subset=p_val_adj<0.05 & avg_log2FC>=0.8 & pct.1>0.5) %>% 
         dplyr::mutate(pct1per2 = pct.1/pct.2)
  write.csv(df_1, file = paste0(DirInteg, "DEG_table/[SubclustsOfCAFs_DEGs_log2FC-0.8_adjP-0.05_pct1-0.5][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                                QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2,".csv"), row.names = T) 
  # volcano plot
  df_2 = df_0 %>% rowwise() %>% 
         dplyr::mutate(minuslog10AdjP = (-1)*log10(p_val_adj))
  df_2$minuslog10AdjP[is.infinite(df_2$minuslog10AdjP)] = 310
  df_2$colors_by_p = ifelse(df_2$p_val_adj<0.05 & (df_2$avg_log2FC>1 | df_2$avg_log2FC<(-1) ),
                            paste0("CAF-",df_2$cluster), "ns")
  P.Volc = 
    ggplot(df_2,aes(x=avg_log2FC, y=minuslog10AdjP)) +
    geom_point(data=subset(df_2, subset=colors_by_p=="ns"), color="gray70", size=0.3) +
    geom_point(data=subset(df_2, subset=colors_by_p!="ns"), aes(color=colors_by_p), size=0.3) +
    labs(x="log2(Fold change)", y="-log10(Adj.p-value)", title="DEGs of CAF sub-clusters",
         subtitle=paste0(NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")\n",
                         QCInfo.FileName,"\n1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2)) +
    scale_y_continuous(limits=c(0, 330))+
    scale_x_continuous(breaks=c(-10, -6, -3, -1, 0, 1, 3, 6, 10)) +
    #scale_color_manual(values=c("ns"="red")) +
    facet_grid( cluster ~ .) +
    geom_hline(yintercept=0, color="gray50") + geom_vline(xintercept=0, color="gray50") +
    theme(legend.position="none",
          axis.text = element_text(face = "bold"),
          axis.text.x = element_text(angle=45, hjust=1),
          axis.title = element_text(face = "bold"),
          legend.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"),
          strip.background = element_rect(fill="white", color="black"),
          panel.background = element_blank(),
          panel.grid = element_line(color = "gray90"),
          panel.border = element_rect(fill = NA, color = "black"))
  ggsave(P.Volc,
         file=paste0(DirInteg, "/[Figure][SubclustOfCAF_Volcano][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                     QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2,".png"), 
         width=8, height=8, dpi=300)
}
Func.FindAllMarkers(Dim1=20, Res1=1.0, Dim2=20, Res2=0.5)
Func.FindAllMarkers(Dim1=20, Res1=1.0, Dim2=30, Res2=0.5)
Func.FindAllMarkers(Dim1=20, Res1=1.0, Dim2=40, Res2=0.5)


###  2. pseudo-bulk GSEA ####
Func.GSEA.Subclust = function(Dim1, Res1, Dim2, Res2){
  library(org.Hs.eg.db)
  library(msigdbr)
  library(clusterProfiler)
  #SeuObj_0 = readRDS(paste0(DirInteg,"Objects/[SeuObj_ExtractedCAF][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
  #                          QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2,".rds") )
  H = clusterProfiler::read.gmt("/Volumes/PortableSSD/[4]myR/[1]dataset/h.all.v2024.1.Hs.entrez.gmt")
  #BW = clusterProfiler::read.gmt('/Volumes/PortableSSD/[4]myR/[1]dataset/c2.all.v2025.1.Hs.entrez.gmt.txt')
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
  DEG_0 = read.csv(file = paste0(DirInteg, "DEG_table/[SubclustsOfCAFs_DEGs_All][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                                 QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2,".csv"), 
                   row.names=1)
  DF.Sym_and_En_0 = clusterProfiler::bitr(DEG_0$gene, fromType="SYMBOL", toType="ENTREZID", 
                                          OrgDb=org.Hs.eg.db, drop=FALSE)
  #manually fill NA cell, failed to map
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
}
Func.GSEA.Subclust(Dim1=20, Res1=1.0, Dim2=20, Res2=0.5)
Func.GSEA.Subclust(Dim1=20, Res1=1.0, Dim2=30, Res2=0.5)
Func.GSEA.Subclust(Dim1=20, Res1=1.0, Dim2=40, Res2=0.5)

Func.visualizeGSEA.Subclust = function(Dim1, Res1, Dim2, Res2){}
  library(ggtext)
  DF.ResGSEA_0 = read.csv(
    file=paste0(DirInteg, "/[DataTable_ExtractedCAF_PseudoBulkGSEA][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2,".csv"),
    row.names=NULL) %>% 
    dplyr::mutate(Subclust = case_when(Sign=="UP" ~ paste0("CAF",Subclust,"_UP"),
                                       Sign=="DOWN" ~ paste0("CAF",Subclust,"_DOWN")))
  NumOfSubclust = (DF.ResGSEA_0$Subclust %>% unique() %>% length())/2
  ###### BarPlot
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
  #clean_gs_label <- function(x){
  #  x %>%
  #    str_remove("^(HALLMARK_|KEGG_)") %>%
  #    str_replace_all("_", " ") %>%
  #    str_squish() %>%
  #    str_to_lower() %>%
  #    str_to_sentence() %>%
  #    str_replace_all("\\bTnfa\\b", "TNF-α") %>%
  #    str_replace_all("\\bNfkb\\b", "NFkB") %>%
  #    str_replace_all("\\bIl6 jak stat3\\b", "IL-6/JAK/STAT3") %>%
  #    str_replace_all("\\bTgf beta\\b", "TGF-β") %>%
  #    str_replace_all("\\bMycaf\\b", "myCAF") %>%
  #    str_replace_all("\\bKras\\b", "KRAS") %>%
  #    str_replace_all("\\bIcaf\\b", "iCAF") %>%
  #    str_replace_all("\\bDna\\b", "DNA") %>%
  #    str_replace_all("\\bUv\\b", "UV") %>%
  #    str_replace_all("\\bMyc\\b", "MYC") %>%
  #    str_replace_all("\\bWnt beta catenin\\b", "Wnt/β-catenin") %>%
  #    str_replace_all("\\b dn\\b", " down") %>%
  #    str_replace_all("\\bnfkb\\b", "NFkB") %>%
  #    str_replace_all("\\bIl2 stat5\\b", "IL-2/STAT5") %>%
  #    str_replace_all("\\bG2m\\b", "G2/M") %>%
  #    str_replace_all("\\bE2f\\b", "E2F") %>%
  #    str_replace_all("\\bMtorc1\\b", "mTORC1")
  #}
  List.Plot = list()
  m = as.numeric(quantile(abs(DF.ResGSEA_1$NES), 0.975, na.rm=TRUE))
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
    scale_fill_distiller(type="div", palette="RdBu", direction=-1,
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
      legend.position = "right",
      plot.background = element_rect(fill="transparent", color=NA),
      panel.background = element_rect(fill="white", color=NA))
  ggsave(plot=Plot.All,
         width=18, height=10, dpi=500, filename = "subclust_cluster_level_gsea.png", bg="transparent")
  
  
  
  
  
  ###### Heatmap
  DF.ResGSEA_1 = DF.ResGSEA_0 %>% 
    dplyr::select(Subclust, ID, NES) %>% 
    pivot_wider(names_from = Subclust,
                values_from = NES)
  DF.ResGSEA_2 = column_to_rownames(DF.ResGSEA_1, var="ID")
  MT.cor.Subclust_0 = DF.ResGSEA_2 %>% 
    cor(method="spearman",
        use="pairwise.complete.obs")
  MT.cor.Pathway_0 = t(DF.ResGSEA_2) %>% 
    cor(method="spearman",
        use="pairwise.complete.obs")
  MT.Dist.pathway = as.dist(1-MT.cor.Pathway_0)
  Hclust.pathway = hclust(MT.Dist.pathway, method = "average")
  plot(Hclust.pathway)
  Vec.PathwayOrder = rownames(DF.ResGSEA_2)[Hclust.pathway$order]
  DF.ResGSEA_0$ID = factor(DF.ResGSEA_0$ID,
                                levels=Vec.PathwayOrder)
  #Vec.ClustOrder = rownames(DF.Zscore_0)[Hclust_clusters$order]
  #DF.Zscore_2$seurat_clusters = factor(DF.Zscore_2$seurat_clusters,
  #                                     levels=Vec.ClustOrder)
  GeomTileColor = "gray40"
  subset(DF.ResGSEA_0, subset=Sign=="UP") %>% 
    ggplot(aes(x=Subclust, y=ID)) +
    geom_tile(aes(fill=case_when(Significance=="ns" ~ NA,
                                 TRUE ~ NES)),
              color=GeomTileColor) +
    geom_text(aes(label=case_when(Significance == "ns" ~ NA,
                                  TRUE ~ formatC(NES, digit=2, width=2)))) +
    labs(x=NULL, y=NULL, fill="NES") +
    #scale_fill_viridis_c(option="turbo", trans=scales::pseudo_log_trans(base = 2)) +
    #scale_fill_gradientn(colors=c("#ffffcc","#fecc5c", "#fd8d3c","#f03b20", "#bd0026"),
    #                     trans=scales::pseudo_log_trans(base = 2)) +
    scale_fill_gradient(low="white", high="red",
                        trans=scales::pseudo_log_trans(base = 2),
                        na.value = "gray80") +
    scale_y_discrete(position="right") +
    theme(legend.position="bottom",
          plot.margin = unit(c(0,0,0,5), "mm"),
          axis.text.y = element_blank(),
          axis.text.x = element_text(angle=45, hjust=1)) +
  subset(DF.ResGSEA_0, subset=Sign=="DOWN") %>% 
    ggplot(aes(x=Subclust, y=ID)) +
    geom_tile(aes(fill=case_when(Significance=="ns" ~ NA,
                                 TRUE ~ NES)),
              color=GeomTileColor) +
    labs(x=NULL, y=NULL, fill="NES") +
    #scale_fill_viridis_c(option="turbo", trans=scales::pseudo_log_trans(base = 2)) +
    #scale_fill_gradientn(colors=c("#253494", "#2c7fb8","#41b6c4", "#a1dab4", "#ffffcc"),
    #                     trans=scales::pseudo_log_trans(base = 2)) +
    scale_fill_gradient(high="white", low="blue",
                        trans=scales::pseudo_log_trans(base = 2),
                        na.value = "gray80") +
    theme(legend.position="bottom",
          plot.margin = unit(c(0,5,0,1), "mm"),
          axis.text.y = element_text(hjust=0.5),
          axis.text.x = element_text(angle=45, hjust=1))

  
  
  
  



  BarPlot.list = lapply(0:(NumObClust-1), function(TargetClust){
    DF.ResGSEA_TargetClust = subset(DF.ResGSEA,
                                    subset=Subclust==TargetClust) %>% 
                             dplyr::arrange(NES)
    DF.ResGSEA_TargetClust$ID = factor(DF.ResGSEA_TargetClust$ID,
                                       levels=DF.ResGSEA_TargetClust$ID)
    #BarPlot.single = 
    ggplot(DF.ResGSEA_TargetClust, aes(x=NES, y=ID)) +
      geom_bar(stat="identity", aes(fill=NES), color="gray30") +
      labs(subtitle=paste0("CAF-",TargetClust), y=NULL) +
      scale_x_continuous(limits=c(-max(abs(DF.ResGSEA_TargetClust$NES)), 
                                  max(abs(DF.ResGSEA_TargetClust$NES)))) +
      scale_fill_gradient2(limits=c(-max(abs(DF.ResGSEA_TargetClust$NES)), 
                                    max(abs(DF.ResGSEA_TargetClust$NES))),
                           high="red",mid="white",low="steelblue",midpoint=0,
                           guide=guide_colorbar(#direction="vertical",
                             #title.position="top",
                             #title.hjust=0.5,
                             frame.colour="black",
                             ticks.colour="black")) +
      geom_vline(xintercept=0, color="gray10") +
      theme(#text = element_text(size = 20),
        axis.text = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        #axis.text.y = element_blank(),
        #axis.ticks.y = element_blank(),
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.key = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"))
  })
  p=ggarrange(plotlist=BarPlot.list, 
              ncol=ceiling(NumObClust/5), nrow=5,
              align="hv") %>% 
    annotate_figure(top=paste0("GSEA of Human hallmarks (MsigDB)\n",
                               NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")\n",
                               QCInfo.FileName,"\n1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2) )
  ggsave(plot=p, 
         file=paste0(DirInteg,"[Figure][",NumOfSamples,"case(",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]_",
                     QCInfo.FileName,"_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2,".png"),
         width=12, height=12, dpi=400, bg="white")



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
###
###############################################################
###
###    Ligand-Reseptor analysis
###    1. CellChat ####
###    2. 
###
###############################################################
###############################################################
Dim1=20; Res1=1.0; Dim2=30; Res2=0.5
library(CellChat)
TX=TXnumInteg[i]
DF.CGtable_Sub_0 = read.csv(file=paste0(DirInteg,"/CellGroupTable/[SubClust_CellGroupTable_TX5K_",TX,"][",
                                        NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                                        "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000ClustDim",Dim2,"Res",Res2,".csv"),)
Directory = paste0("/Volumes/Extreme SSD/Analysis/Data/TX5K_",
                   if(TX %in% c("27","28")){ "27_28/" }else{ paste0(TX,"/") })
SeuObj_0 = readRDS(paste0(Directory,"Objects/[SeuObj][TX5K_",TX,"]_Countable_Mag",Mag,
                           "_nFeat",nFeatRNA[1],"-",nFeatRNA[2],
                           "_nCount",nCountRNA[1],"-",nCountRNA[2],".rds") ) %>% 
           subset(cell=DF.CGtable_Sub_0$cell_id)
DF.CGtable_Sub_0$cell_id = factor(DF.CGtable_Sub_0$cell_id, levels=colnames(SeuObj_0))
DF.Meta = dplyr::arrange(DF.CGtable_Sub_0, cell_id) %>% 
          column_to_rownames(var="cell_id") %>% 
          dplyr::mutate(samples = as.factor(paste0("TX5K_",TX)))
SeuObj_1 = NormalizeData(SeuObj_0, verbose=FALSE)
rm(SeuObj_0)

DispClusts = c("CAF-0","CAF-1","CAF-2","CAF-3","CAF-4","CAF-5",#"CAF-6",
               "CAF-7","CAF-8",
               #"Clust.0","Clust.13","Clust.24","Clust.25","Clust.33","Clust.1","Clust.12",
               #"Clust.17","Clust.30",
               "Clust.23","Clust.27","Clust.10","Clust.20"#,"Clust.14"
               )
DF.Meta_EandC = subset(DF.Meta,
                       subset=group%in%DispClusts)
SeuObj_EandC_1 = subset(SeuObj_1,cells=rownames(DF.Meta_EandC) )
MT.ExpInput_EandC = GetAssayData(SeuObj_EandC_1, assay="RNA", layer="data")
CellChat_EandC_0 = createCellChat(object=MT.ExpInput_EandC, meta=DF.Meta_EandC, group.by="group")
CellChat_EandC_0@DB = CellChatDB.human
# pre-process
CellChat_EandC_1 = subsetData(CellChat_EandC_0)
CellChat_EandC_1 = identifyOverExpressedGenes(CellChat_EandC_1)
CellChat_EandC_1 = identifyOverExpressedInteractions(CellChat_EandC_1)
names(CellChat_EandC_1@LR)
coords_EandC <- data.frame("x"=SeuObj_EandC_1$X,
                           "y"=SeuObj_EandC_1$Y,
                           row.names=colnames(SeuObj_EandC_1))
CellChat_EandC_1@images$spatial = list(coords = coords_EandC)
# communication probability
set.seed(123)
CellChat_EandC_2 = computeCommunProb(CellChat_EandC_1, 
                                     raw.use=TRUE,
                                     distance.use=TRUE,
                                     interaction.range=200)
# filtering
CellChat_EandC_3 = filterCommunication(CellChat_EandC_2, min.cells = 10)
# compute communication probability pathway level
CellChat_EandC_3 = computeCommunProbPathway(CellChat_EandC_3)
# aggregate networks ( total Cell-cell networks )
CellChat_EandC_3 = aggregateNet(CellChat_EandC_3)
groupSize_EandC = as.numeric(table(CellChat_EandC_3@idents)[DispClusts])
# circle plot 1
MT.count = CellChat_EandC_3@net$count[DispClusts, DispClusts, drop=FALSE]
MT.weight = CellChat_EandC_3@net$weight[DispClusts, DispClusts, drop=FALSE]

netVisual_circle(MT.count,
                 vertex.weight = groupSize_EandC,
                 weight.scale = TRUE,
                 label.edge = FALSE)

netVisual_circle(MT.weight,
                 vertex.weight = groupSize_EandC,
                 weight.scale = TRUE,
                 label.edge = FALSE)
# circle plot 2 
caf <- paste0("CAF-", 0:8)
MT.W.CAFtoNON = MT.weight
MT.W.CAFtoNON[!rownames(MT.W.CAFtoNON)%in%caf, ] = 0
MT.W.CAFtoNON[ , colnames(MT.W.CAFtoNON)%in%caf] = 0
netVisual_circle(MT.W.CAFtoNON,
                 vertex.weight = groupSize_EandC,
                 weight.scale = TRUE,
                 label.edge = FALSE)
MT.W.NONtoCAF = MT.weight
MT.W.NONtoCAF[rownames(MT.W.NONtoCAF)%in%caf, ] = 0
MT.W.NONtoCAF[ , !colnames(MT.W.NONtoCAF)%in%caf] = 0
netVisual_circle(MT.W.NONtoCAF,
                 vertex.weight = groupSize_EandC,
                 weight.scale = TRUE,
                 label.edge = FALSE)
#CellChat_EandC_4 = CellChat_EandC_3
#CellChat_EandC_4@net$weight = MT.W.NONtoCAF
#netVisual_chord_cell(CellChat_EandC_4,
#                     net="weight",
#                     vertex.weight = groupSize_EandC,
#                     lab.cex = 0.6,
#                     small.gap = 1)

# 
CellChat_EandC_3@idents = factor(CellChat_EandC_3@idents, levels=DispClusts)
CellChat_EandC_3@meta$group = factor(CellChat_EandC_3@meta$group, levels=DispClusts)
CellChat_EandC_3 <- netAnalysis_computeCentrality(CellChat_EandC_3, slot.name = "netP")
netAnalysis_signalingRole_network(CellChat_EandC_3, slot.name = "netP",
                                  signaling="FGF")
names(CellChat_EandC_3@netP[["centr"]])    # pathway list

Func.SigRoleNetWork = function(pwy){
centr = CellChat_EandC_3@netP$centr[[pwy]]
clust_names = names(centr$outdeg)
names(centr$flowbet) = clust_names
names(centr$info) = clust_names
DF.Centr_0 = data.frame("cluster"    = DispClusts,
                        "Sender"     = centr$outdeg[DispClusts],
                        "Receiver"   = centr$indeg[DispClusts],
                        "Mediator"   = centr$flowbet[DispClusts],
                        "Influencer" = centr$info[DispClusts],
                        check.names = FALSE)
DF.Centr_0_long <- DF.Centr_0 %>%
  pivot_longer(-cluster, names_to="Role", values_to="value") %>%
  group_by(Role) %>%
  dplyr::mutate(value_scaled = {
                  r <- range(value, na.rm=TRUE)
                  if (diff(r) == 0) 0 else (value - r[1]) / diff(r)
                  }
                ) %>% ungroup()
DF.Centr_1 = DF.Centr_0_long %>%
             dplyr::mutate(cluster = factor(cluster, levels = DispClusts),
                           Role = factor(Role, levels = c("Influencer","Mediator","Receiver","Sender")))
ggplot(DF.Centr_1, aes(x=cluster, y=Role, fill=value_scaled)) +
  geom_tile(color="gray95") +
  scale_fill_gradientn(colors=c("#ffffff","#edf8fb","#ccece6","#99d8c9",
                                "#66c2a4","#41ae76","#238b45","#005824"),
                      #value,low="white", high="darkgreen",
                      limits=c(0,1), name="Importance\n(scaled)",
                      guide=guide_colorbar(direction="vertical",
                                           title.position="top",
                                           title.hjust=0.5,
                                           frame.colour="black",
                                           ticks.colour="black")) +
  labs(title=paste0("Signaling role heatmap: ", pwy), x=NULL, y=NULL) +
  scale_x_discrete(expand=expansion(mult=c(0,0))) +
  scale_y_discrete(expand=expansion(mult=c(0,0))) +
  theme_classic() +
  theme(axis.text = element_text(face = "bold"),
        axis.text.x = element_text(angle=45, hjust=1),
        axis.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.background = element_blank())
}
Func.SigRoleNetWork(pwy="COLLAGEN")
Func.SigRoleNetWork(pwy="FN1")
Func.SigRoleNetWork(pwy="LAMININ")
Func.SigRoleNetWork(pwy="TENASCIN")
Func.SigRoleNetWork(pwy="FGF")
Func.SigRoleNetWork(pwy="PDGF")
Func.SigRoleNetWork(pwy="TGFb")
Func.SigRoleNetWork(pwy="NOTCH")
Func.SigRoleNetWork(pwy="IGF")
Func.SigRoleNetWork(pwy="EGF")
Func.SigRoleNetWork(pwy="GAS")
Func.SigRoleNetWork(pwy="KLK")

  CommonTitle.Sub = paste0(NumOfSamples,"case Integration(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")\n",
                           QCInfo,"\n1stClust:PCA50vf2000_ClustDim",Dim1,"Res",formatC(Res1,digits=2,format="f"),
                           "\n2ndClust:PCA50vf2000_ClustDim",Dim2,"Res",formatC(Res2,digits=2,format="f"))
  NumOfSublust = length(unique(SeuObj.CAF_0$seurat_clusters))
  CommonCaption.Sub = paste0("Total ",format(ncol(Seu.CAF_0), big.mark=",", scientific=F), " CAFs, ",NumOfSublust, " clusters.")
  netAnalysis_signalingRole_network(CellChat_EandC_3, slot.name = "netP",
                                    signaling=pwy)
  Clusts.Epi = c("Clust.23","Clust.10","Clust.27","Clust.20")
  DF.LR.CtoE = subsetCommunication(
    object = CellChat_EandC_3,
    sources.use = c("CAF-0","CAF-1","CAF-2","CAF-3","CAF-4","CAF-5",#"CAF-6",
                    "CAF-7","CAF-8"),
    targets.use = Clusts.Epi,
    slot.name = "net" ) %>% 
    dplyr::rename("CAF"=source, "EPI"=target) %>% 
    dplyr::mutate(direction="CAF_to_EPI") %>% 
    dplyr::arrange(pathway_name, interaction_name_2)
  DF.LR.EtoC = subsetCommunication(
    object = CellChat_EandC_3,
    sources.use = Clusts.Epi,
    targets.use = c("CAF-0","CAF-1","CAF-2","CAF-3","CAF-4","CAF-5",#"CAF-6",
                    "CAF-7","CAF-8"),
    slot.name = "net") %>% 
    dplyr::rename("CAF"=target, "EPI"=source) %>% 
    dplyr::mutate(direction="EPI_to_CAF") %>% 
    dplyr::arrange(pathway_name, interaction_name_2)
  DF.LR.Both = rbind(DF.LR.CtoE, DF.LR.EtoC) %>% 
    dplyr::mutate(interaction_name_3 = 
                  str_replace(interaction_name_2,
                              pattern=" - ", replacement=" -> ") ) %>% 
    dplyr::mutate(EPI = factor(.$EPI, levels=Clusts.Epi) ) %>% 
    dplyr::mutate(interaction_name_3 = factor(.$interaction_name_3,
                                              levels=unique(.$interaction_name_3) ) )
  library("ggh4x")
  Plot.Dot = 
    ggplot(DF.LR.Both, 
         aes(x=CAF, y=interaction_name_3, color=pathway_name)) +
    geom_point(aes(size=prob)) + 
    geom_point(aes(size=prob), color="gray50", shape=21) + 
    labs(subtitle=CommonTitle.Sub,
         x=NULL, y=NULL,size="Communication\nprobability",
         caption=CommonCaption.Sub) +
    ggh4x::facet_grid2(direction ~ EPI, 
                       scales = "free_y",
                       space = "free_y") + 
    theme(axis.text = element_text(face = "bold"),
          axis.text.x = element_text(angle=45, hjust=1),
          axis.title = element_text(face = "bold"),
          legend.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold"),
          legend.key = element_blank(),
          strip.text = element_text(face = "bold"),
          strip.background = element_rect(color="black"),
          panel.background = element_blank(),
          panel.grid = element_line(color = "gray90"),
          panel.border = element_rect(fill = NA, color = "black"))
  ggsave(plot=Plot.Dot, 
         file=paste0(DirInteg,"[Figure][SubClust_LRanalysis_Dotplot][",
                     NumOfSamples,"case(",IntegArgo,", TX5K_",paste(TXnumInteg[1:NumOfSamples],collapse=","),")]",
                     "_1stPCA50vf2000_ClustDim",Dim1,"Res",Res1,"_2ndPCA50vf2000_ClustDim",Dim2,"Res",Res2,".png"), 
         width=12, height=15, dpi=400, bg="white")
  
  # 確認
  head(DF.LR)
VlnPlot(Seu.CAF_1, 
        feature=c("ITGAV","ITGB4","CD74"), 
        pt.size=0, raster=FALSE)
  
  
  P1 = ggplot(DF.LR.CtoE, 
         aes(x=CAF, y=interaction_name_3, color=pathway_name)) +
    geom_point(aes(size=prob)) + 
    geom_point(aes(size=prob), color="gray50", shape=21) + 
    labs(x=NULL, y=NULL) +
    facet_grid(direction ~ EPI) + CommonTheme + coord_fixed(ratio=1)
  P2 = ggplot(DF.LR.EtoC, 
         aes(x=CAF, y=interaction_name_3, color=pathway_name)) +
    geom_point(aes(size=prob)) + 
    geom_point(aes(size=prob), color="gray50", shape=21) + 
    labs(x=NULL, y=NULL) +
    facet_grid(direction ~ EPI) + CommonTheme + coord_fixed(ratio=1)
  P1.Lgd = get_legend(P1)
  P2.Lgd = get_legend(P2)
  
  Plot.Dot = ggarrange(P1 / P2) + plot_layout(widths = c(1))# +
  #annotate_figure(top=CommonTitle.Sub,
  #                bottom=CommonCaption.Sub)
  
  
  # Fig.5g に相当
  netVisual_bubble(CellChat_EandC_3, 
                   sources.use = c("CAF-0","CAF-1","CAF-2","CAF-3","CAF-4","CAF-5",#"CAF-6",
                                   "CAF-7","CAF-8"),
                   targets.use = c("Clust.23","Clust.27","Clust.10","Clust.20"))
  
  
  ###################################
  MT.ExpInput = GetAssayData(SeuObj_1, assay="RNA", layer="data")
  CellChat_0 = createCellChat(object=MT.ExpInput, meta=DF.Meta, group.by="group")
  CellChat_0@DB = CellChatDB.human
  
  # pre-process
  CellChat_1 = subsetData(CellChat_0)
  CellChat_1 = identifyOverExpressedGenes(CellChat_1)
  CellChat_1 = identifyOverExpressedInteractions(CellChat_1)
  names(CellChat_1@LR)
  # communication probability
  set.seed(123)
  CellChat_2 = computeCommunProb(CellChat_1, raw.use=TRUE)
  # filtering
  CellChat_3 <- filterCommunication(CellChat_2, min.cells = 10)
  # communication probability, pathway level
  CellChat_3 <- computeCommunProbPathway(CellChat_3)
  # aggregate networks (total Cell-cell networks)
  CellChat_3 <- aggregateNet(CellChat_3)
  groupSize <- as.numeric(table(CellChat_3@idents))
  netVisual_circle(CellChat_3@net$count,
                   vertex.weight = groupSize,
                   weight.scale = TRUE,
                   label.edge = FALSE)
  
  netVisual_circle(CellChat_3@net$weight,
                   vertex.weight = groupSize,
                   weight.scale = TRUE,
                   label.edge = FALSE)
  
  caf <- paste0("CAF-", 0:8)
  levels(CellChat_3@idents)  # 念のため存在確認
  
  netVisual_circle(CellChat_3@net$weight,
                   sources.use = caf,
                   vertex.weight = groupSize,
                   weight.scale = TRUE,
                   label.edge = FALSE)
  
  netVisual_circle(CellChat_3@net$weight,
                   targets.use = caf,
                   vertex.weight = groupSize,
                   weight.scale = TRUE,
                   label.edge = FALSE)
  
  CellChat_3 <- netAnalysis_computeCentrality(CellChat_3, slot.name = "netP")
  netAnalysis_signalingRole_network(CellChat_3, slot.name = "netP")
  
  ###############################################################