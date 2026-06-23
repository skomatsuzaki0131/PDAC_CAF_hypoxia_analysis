library(ggplot2)
library(dplyr)
library(magrittr)
library(cowplot)
library(tidyverse)
library(ggsignif)
library(legendry)
library(exactRankTests)
#library(psych)

DF.Morpho.0207Nor = read.csv('/Volumes/PortableSSD/[4]myR/[1]dataset/CAFmorphology_Results_CAF24-0207NorP2d2.csv') %>% 
                    dplyr::select(X, AR) %>% data.frame(Pair="Pair13", Oxygen="Normo")
DF.Morpho.0207Hyp = read.csv('/Volumes/PortableSSD/[4]myR/[1]dataset/CAFmorphology_Results_CAF24-0207HypP3d2.csv') %>% 
                    dplyr::select(X, AR) %>% data.frame(Pair="Pair13", Oxygen="Hypo")
DF.Morpho.0403Nor = read.csv('/Volumes/PortableSSD/[4]myR/[1]dataset/CAFmorphology_Results_CAF24-0403NorP3d2.csv') %>% 
                    dplyr::select(X, AR) %>% data.frame(Pair="Pair16", Oxygen="Normo")
DF.Morpho.0403Hyp = read.csv('/Volumes/PortableSSD/[4]myR/[1]dataset/CAFmorphology_Results_CAF24-0403HypP4d2.csv') %>% 
                    dplyr::select(X, AR) %>% data.frame(Pair="Pair16", Oxygen="Hypo")
DF.Morpho.0417Nor = read.csv('/Volumes/PortableSSD/[4]myR/[1]dataset/CAFmorphology_Results_CAF24-0417NorP2d2.csv') %>% 
                    dplyr::select(X, AR) %>% data.frame(Pair="Pair17", Oxygen="Normo")
DF.Morpho.0417Hyp = read.csv('/Volumes/PortableSSD/[4]myR/[1]dataset/CAFmorphology_Results_CAF24-0417HypP1d2.csv') %>% 
                    dplyr::select(X, AR) %>% data.frame(Pair="Pair17", Oxygen="Hypo")
DF.Morphology = bind_rows(DF.Morpho.0207Nor, DF.Morpho.0207Hyp, 
                          DF.Morpho.0403Nor, DF.Morpho.0403Hyp,
                          DF.Morpho.0417Nor, DF.Morpho.0417Hyp)
#psych::describeBy(Data.morphology, group = Data.morphology$Oxygen)


# Shapiro test
# if p<0.05  --> not normal distibution
# --> use non-parametric test (i.e. wilcoxon)
Shapiro_Morphology = shapiro.test(DF.Morphology$AR)
Shapiro_Morphology  # p = 2.041e-13  -> use wilcoxon

# non-parametiric test
Func.Wilcox.Morpho = function(PairNumber){
  wilcox.exact(DF.Morphology[DF.Morphology$Pair==PairNumber & DF.Morphology$Oxygen=="Normo","AR"],
               DF.Morphology[DF.Morphology$Pair==PairNumber & DF.Morphology$Oxygen=="Hypo","AR"], paired=F) }
Pvalue.0207 = Func.Wilcox.Morpho(PairNumber="Pair13")$p.value # *** p = 2.957699e-08
Pvalue.0403 = Func.Wilcox.Morpho(PairNumber="Pair16")$p.value # *** p = 3.986602e-05
Pvalue.0417 = Func.Wilcox.Morpho(PairNumber="Pair17")$p.value # *** p = 1.357134e-07

DF.Morphology_2 = DF.Morphology
DF.Morphology_2$Oxygen2 = factor(DF.Morphology_2$Oxygen, levels=c("Normo", "Hypo"))
DF.Morphology_2$Int = interaction(DF.Morphology_2$Oxygen2,
                                  DF.Morphology_2$Pair, sep=".") %>% 
                      factor(levels=c("Normo.Pair13", "Hypo.Pair13",
                                      "Normo.Pair16", "Hypo.Pair16",
                                      "Normo.Pair17", "Hypo.Pair17"))
Vec.Pos = c(12,18,32,38,52,58)
DF.morphology_3 = cbind(DF.Morphology_2,
                          Position = case_when(DF.Morphology_2$Int=="Normo.Pair13" ~ Vec.Pos[1],
                                               DF.Morphology_2$Int=="Hypo.Pair13" ~ Vec.Pos[2],
                                               DF.Morphology_2$Int=="Normo.Pair16" ~ Vec.Pos[3],
                                               DF.Morphology_2$Int=="Hypo.Pair16" ~ Vec.Pos[4],
                                               DF.Morphology_2$Int=="Normo.Pair17" ~ Vec.Pos[5],
                                               DF.Morphology_2$Int=="Hypo.Pair17" ~ Vec.Pos[6]))
Plot.CAFmorphology = 
  ggplot(DF.morphology_3, aes(x=Position, y=AR, group=Int)) +
  geom_boxplot(aes(fill=Oxygen2, color=Oxygen2),
               width=5.0, staplewidth=0.9, linewidth=1.1, 
               outliers=F, outlier.shape=NA,
               median.linewidth=0.5)+
  geom_dotplot(binaxis ="y", color="black",fill="black", stroke=2,
               stackdir="center", binwidth=0.4, alpha=0.8,
               dotsize = 1.0) +
  #stat_summary(fun=mean, geom="point") +
  labs(subtitle=NULL,
       caption="30 cells per CAF culture. Statistics:Wilcoxon rank sum test",
       y="Aspect ratio \n(Long / short axis)", x=NULL, fill="Oxygen") +
  scale_x_continuous(
    breaks=c(Vec.Pos[1],Vec.Pos[2],
             Vec.Pos[3],Vec.Pos[4],
             Vec.Pos[5],Vec.Pos[6]),
    labels=c(" Normo-CAF.Pair #13", " Hypo-CAF.Pair #13",
             " Normo-CAF.Pair #16", " Hypo-CAF.Pair #16",
             " Normo-CAF.Pair #17", " Hypo-CAF.Pair #17")) +
  scale_y_continuous(limits=c(0, 25), 
                     breaks=c(0, 5, 10, 15, 20),
                     expand=expansion(mult=c(0, 0))) +
  scale_fill_manual(values=c(Hypo="#91BFFA", Normo="#FEA0A0"), #"#B2182B" "#D6604D" "#F4A582" "#FDDBC7" "#D1E5F0" "#92C5DE" "#4393C3" "#2166AC"
                    labels=c(Hypo="H-CAF", Normo="N-CAF")) +
  scale_color_manual(values=c(Hypo="#2166AC", Normo="#B2182B")) +
  guides(fill = guide_legend(override.aes=list(size=18))) +
  geom_segment(x=Vec.Pos[1]-0.1,xend=Vec.Pos[2]+0.1,y=20.25,yend=20.25,linewidth=1.2) +
  geom_segment(x=Vec.Pos[1],xend=Vec.Pos[1],y=19.75,yend=20.25,linewidth=1.2) +
  geom_segment(x=Vec.Pos[2],xend=Vec.Pos[2],y=19.75,yend=20.25,linewidth=1.2) +
  annotate(geom="text", x=(Vec.Pos[1]+Vec.Pos[2])/2, y=20.75,fontface="bold", size=10, label="***") + 
  geom_segment(x=Vec.Pos[3]-0.1,xend=Vec.Pos[4]+0.1,y=20.83,yend=20.83,linewidth=1.2) +
  geom_segment(x=Vec.Pos[3],xend=Vec.Pos[3],y=20.33,yend=20.83,linewidth=1.2) +
  geom_segment(x=Vec.Pos[4],xend=Vec.Pos[4],y=20.33,yend=20.83,linewidth=1.2) +
  annotate(geom="text", x=(Vec.Pos[3]+Vec.Pos[4])/2, y=21.33,fontface="bold", size=10, label="***") +
  geom_segment(x=Vec.Pos[5]-0.1,xend=Vec.Pos[6]+0.1,y=18.5,yend=18.5,linewidth=1.25) +
  geom_segment(x=Vec.Pos[5],xend=Vec.Pos[5],y=18.0,yend=18.5,linewidth=1.2) +
  geom_segment(x=Vec.Pos[6],xend=Vec.Pos[6],y=18.0,yend=18.5,linewidth=1.2) +
  annotate(geom="text", x=(Vec.Pos[5]+Vec.Pos[6])/2, y=19.0,fontface="bold", size=10, label="***") +
  theme(text = element_text(size = 30),
        legend.position = "none",
        axis.text = element_text(face = "bold"),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, color="black"),
        axis.text.y = element_text(color="black"),
        axis.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        plot.background = element_rect(fill="transparent", color=NA),
        panel.background = element_rect(fill="white", color=NA),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"))
Plot.CAFmorphology
ggsave(plot=Plot.CAFmorphology,
       file="[3]Figures/2[Figures]CAFmorphology_Plot.png",
       width=9, height=8.5, dpi=400, bg="transparent")
df.No13 = subset(DF.morphology_3, subset=Pair=="Pair13") %>% group_by(Oxygen2)
df.No16 = subset(DF.morphology_3, subset=Pair=="Pair16") %>% group_by(Oxygen2)
df.No17 = subset(DF.morphology_3, subset=Pair=="Pair17") %>% group_by(Oxygen2)
summarize(df.No13, Mean = mean(AR), Median = median(AR), SD = sd(AR))
summarize(df.No16, Mean = mean(AR), Median = median(AR), SD = sd(AR))
summarize(df.No17, Mean = mean(AR), Median = median(AR), SD = sd(AR))

saveRDS(Plot.CAFmorphology, "[3]Figures/[Figures]CAFmorphology_Plot.rds")




ggplot(DF.morphology_3, aes(x=Position, y=AR, group=Int)) +
  geom_boxplot(aes(fill=Oxygen2), width=4, staplewidth=1.1, linewidth=1.5) +
  ggbeeswarm::geom_quasirandom(
    color="black",
    groupOnX = TRUE,
    shape = 4,
    size = 3,
    stroke = 1,
    alpha = 0.7,
    width = 1.8
  ) +
  stat_summary(fun=mean, geom="point", size=4) +
  scale_fill_manual(values=c(Hypo="#91BFFA", Normo="#FEA0A0"),
                    labels=c(Hypo="H-CAF", Normo="N-CAF")) +
  scale_color_manual(values=c(Hypo="#91BFFA", Normo="#FEA0A0"),
                     guide="none")