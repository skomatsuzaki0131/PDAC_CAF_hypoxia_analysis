library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggfortify)

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 2)
Dir.OutPCA = "/Volumes/PortableSSD/[3]Graduate school/[3]RNA-seq/analysis/PCA/"
ColorCode = c("#607DF0", "#91BFFA","#FEA0A0","#FF7F7F")
List.Colors = list(H=adjustcolor("#607DF0", alpha.f=0.6),
                   HN=adjustcolor("#7FBFFF", alpha.f=0.6),
                   NH=adjustcolor("#FFA3A3", alpha.f=0.6),#FFDDDCもいい
                   N=adjustcolor("red", alpha.f=0.6))
List.SampleName = list(N = expression(bold("N-CAF")),
                       NH = expression(bold("N-CAF in 1% "*O[2]*"")),
                       HN = expression(bold("H-CAF in 21% "*O[2]*"")),
                       H = expression(bold("H-CAF")))

# Data読み込み、成形
Data = read.csv("/Volumes/PortableSSD/[4]myR/[1]dataset/rnaseq_allTPM.csv", header = TRUE, sep = ",")
tData = as.data.frame(t(Data))
NumOfSamples = ncol(Data)-1
Vec.Patient = c(rep("23-0904",times=12), rep("23-1115",times=12), rep("23-1122",times = 12))
Vec.Oxygen = c(rep(c("H","H","H","HN","HN","HN","N","N","N","NH","NH","NH"), times = 3))
Vec.Sample = interaction(Vec.Patient, Vec.Oxygen, sep = "")

# 単位分散0の遺伝子を除去する
Vec.GeneNames = tData[1, ]
colnames(tData) = Vec.GeneNames
Vec.DeleteGenes = c()       #単位分散が0の遺伝子をこのvectorに入れていく
for (x in tData) {           #tDataの1列目から順に
  data5 = sapply(x[2:length(x)], as.numeric)     #1. xの2~末列までnumeric型にする
  data5_var = var(data5)         #2. var()で単位分散を求める
  if (data5_var == 0) {          #もしもdata5_varが0ならば
      Vec.DeleteGenes <- c(Vec.DeleteGenes, x[1])    #Vec.DeleteGenesに,GeneNameであるx[1]を付加する
  }
}
for (x in Vec.DeleteGenes) {tData = dplyr::select(tData, -x)}   #tDataからVec.DeleteGenesに含まれる遺伝子の列を削除

DF.TPM_2 = sapply(tData[2:(NumOfSamples+1), ], as.numeric)
DF.TPM_3 = t(DF.TPM_2) %>% as.data.frame()

Func.Filter.PCA = function(Coef){
DF.TPM_4 = DF.TPM_3[rowSums(DF.TPM_3 >= 1) >= (ncol(DF.TPM_3))*Coef, ] %>% 
           t() %>% as.data.frame()
PCAres = prcomp(DF.TPM_4, scale. = TRUE)
Vec.explained_variance <- PCAres$sdev^2 / sum(PCAres$sdev^2)
DF.TPM_5 = PCAres$x %>% as.data.frame()
Func.Plot = function(PCpair){
DF.TPM_6 = cbind(DF.TPM_5[ ,PCpair[1]],
                 DF.TPM_5[ ,PCpair[2]],
                 Patient = Vec.Patient,
                 Oxygen = Vec.Oxygen,
                 Sample = Vec.Sample) %>% as.data.frame()
DF.TPM_6$V1 = as.numeric(DF.TPM_6$V1)
DF.TPM_6$V2 = as.numeric(DF.TPM_6$V2)
DF.TPM_6$Oxygen = factor(DF.TPM_6$Oxygen, levels=c("H","HN","NH","N"))
DF.TPM_6$Patient = DF.TPM_6$Patient %>% 
                   str_replace(pattern="23-0904", replacement="Pair#2") %>% 
                   str_replace(pattern="23-1115", replacement="Pair#8") %>% 
                   str_replace(pattern="23-1122", replacement="Pair#9")
ggplot(DF.TPM_6, aes(x=V1, y=V2, color=Oxygen, shape=Patient)) +
  geom_point() +
  scale_color_manual(values=ColorCode,
                     labels=List.SampleName) +
  labs(x=paste0("PC", PCpair[1], " (",round(Vec.explained_variance[PCpair[1]],3)*100,"%)"),
       y=paste0("PC", PCpair[2], " (",round(Vec.explained_variance[PCpair[2]],3)*100,"%)"),
       shape="CAF pair") +
  theme(text = element_text(face = "bold"),
        legend.key = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"))
}
PlotPCA = 
  ggarrange(Func.Plot(PCpair=c(1,2)),
            Func.Plot(PCpair=c(1,3)),
            Func.Plot(PCpair=c(2,3)),
            ncol=3, nrow=1, common.legend=T, legend="right") %>% 
  annotate_figure(top=text_grob(paste0("ExpressingSample>=",Coef), face="bold"))
ggsave(plot=PlotPCA, 
       file=paste0(Dir.OutPCA,"[Figure][PCA]_GeneCutOff_ExpressingSampleRatio_",Coef,".png"),
       width=8.5, height=2.8, dpi=200, bg="white")
}
for(i in c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)[1:10]){
  Func.Filter.PCA(Coef=i)}

#######################################################

DF.TPM_4 = DF.TPM_3[rowSums(DF.TPM_3 >= 1) >= (ncol(DF.TPM_3))*1.0, ] %>% 
  t() %>% as.data.frame()
PCAres = prcomp(DF.TPM_4, scale. = TRUE)
Vec.explained_variance <- PCAres$sdev^2 / sum(PCAres$sdev^2)
DF.TPM_5 = PCAres$x %>% as.data.frame()
PCpair = c(1,3)
DF.TPM_6 = cbind(DF.TPM_5[ ,PCpair[1]],
                 DF.TPM_5[ ,PCpair[2]],
                 Patient = Vec.Patient,
                 Oxygen = Vec.Oxygen,
                 Sample = Vec.Sample) %>% as.data.frame()
DF.TPM_6$V1 = as.numeric(DF.TPM_6$V1)
DF.TPM_6$V2 = as.numeric(DF.TPM_6$V2)
DF.TPM_6$Oxygen = factor(DF.TPM_6$Oxygen, levels=c("H","HN","NH","N"))
DF.TPM_6$Patient = DF.TPM_6$Patient %>% 
  str_replace(pattern="23-0904", replacement="Pair#2") %>% 
  str_replace(pattern="23-1115", replacement="Pair#8") %>% 
  str_replace(pattern="23-1122", replacement="Pair#9")
Plot =
  ggplot(DF.TPM_6, aes(x=V1, y=V2, shape=Patient)) +
  geom_point(aes(color=Oxygen), size=4) +
  scale_color_manual(values=ColorCode,
                     labels=List.SampleName) +
  labs(x=paste0("PC1 (",round(Vec.explained_variance[PCpair[1]],3)*100,"%)"),
       y=paste0("PC2 (",round(Vec.explained_variance[PCpair[2]],3)*100,"%)"),
       shape="CAF pair", caption=paste0("GeneFiltering: ExpressingSample=AllOf36samples & var()>0.")) +
  theme(text = element_text(face="bold", size=15),
        plot.caption = element_text(hjust=0),
        legend.key = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"))
ggsave(plot=Plot, 
       file=paste0(Dir.OutPCA,"[Figure][PCA]_GeneFiltering_OnlyExpressingSample.png"),
       width=6, height=4, dpi=400, bg="white")



summary(pca_res)
##                            PC1     PC2      PC3
## Standard deviation     66.8197 51.8652 41.52666 : 「標準偏差」そのPCの標準偏差。固有値の平方根。
## Proportion of Variance  0.2496  0.1504  0.09641 : 「寄与率」各PCで何%説明できるか
## Cumulative Proportion   0.2496  0.4000  0.49644 : 「累積寄与率」PC1,2,...PCXまでの寄与率の合計
## 累積寄与率は70~80%を超える点までの主成分を選択すると良い

round(pca_res$rotation, 3)


round(pca_res$x, 3)
#主成分得点
##           PC1     PC2     PC3     PC4     PC5     PC6
## [1,]  -16.026  12.769  57.799  46.045 -25.629  15.895
## [2,]  -15.486  13.851  57.941  43.061 -32.286  13.627
## [3,]  -10.679  29.172  61.675  35.492 -69.376  13.188


##############################################
##############################################
##                                          ##
##                3 dimension               ##
##                                          ##
##############################################
##############################################

library(rgl)
library(scatterplot3d)
Vec.Color = List.Colors[[Vec.Oxygen[1]]]
for(i in 1:35){
  Colornext = List.Colors[[Vec.Oxygen[i+1]]]
  Vec.Color = c(Vec.Color, Colornext)}
DF.TPM_7 = cbind(DF.TPM_5[ ,1:3],
                 Patient = Vec.Patient,
                 Oxygen = Vec.Oxygen,
                 Color = Vec.Color,
                 Sample = Vec.Sample) %>% as.data.frame()
plot3d(x=DF.TPM_7$PC1,
       y=DF.TPM_7$PC2,
       z=DF.TPM_7$PC3,
       col=Vec.Color,
       shapes="x",
       size=1, 
       type=c("n"),
       #lwd=10,
       #radius=3, 
       add=F,
       xlab=paste0("PC1 (",round(Vec.explained_variance[1],3)*100,"%)"),
       ylab=paste0("PC2 (",round(Vec.explained_variance[2],3)*100,"%)"),
       zlab=paste0("PC3 (",round(Vec.explained_variance[3],3)*100,"%)"))
for (i in 1:36) {
  points3d(
    x = DF.TPM_7$PC1[i], 
    y = DF.TPM_7$PC2[i], 
    z = DF.TPM_7$PC3[i], 
    col = Vec.Color[i], 
    pch = 16, 
    size = 10
  )
}
bgplot3d({plot.new(); title("main")})
rgl.snapshot(filename="a.png")
snapshot3d(filename="b.png", webshot=F)

png("scatterplot3d_highres.png", width = 1800, height = 1600, res = 400)
scatterplot3d(x=DF.TPM_7$PC1,
       y=DF.TPM_7$PC3,
       z=DF.TPM_7$PC2,
       xlim=c(-max(abs(DF.TPM_7$PC1)),max(abs(DF.TPM_7$PC1))),
       ylim=c(-max(abs(DF.TPM_7$PC3)),max(abs(DF.TPM_7$PC3))),
       zlim=c(-max(abs(DF.TPM_7$PC2)),max(abs(DF.TPM_7$PC2))),
       bg=DF.TPM_7$Color, pch=21,
       cex.symbols=1.5,
       angle=105, scale.y=1.4,
       grid=T, box=T, 
       #lab = c(x=3, y=3), lab.z=3, cex.lab=1,
       mar=c(3,3,3,3),
       xlab=paste0("PC1 (",round(Vec.explained_variance[1],3)*100,"%)"),
       ylab=paste0("PC3 (",round(Vec.explained_variance[3],3)*100,"%)"),
       zlab=paste0("PC2 (",round(Vec.explained_variance[2],3)*100,"%)"))
dev.off()
