library(ggplot2)
library(dplyr)
library(magrittr)
library(cowplot)
library(tidyverse)
library(ggsignif)
#library(psych)
library(exactRankTests)

Data.Outgrowth = read.csv(file = "/Volumes/PortableSSD/[4]myR/[1]dataset/CAFoutgrowth.csv", header = TRUE)
Data_total = read.csv("[1]dataset/CAFlist.csv", header=TRUE, row.names = 1, stringsAsFactors=FALSE)
Data.Outgrowth = Data.Outgrowth %>% 
  dplyr::mutate(Site = as.character(Data_total["Site", ]),
                NAT = as.character(Data_total["NAT", ]))
str(Data.Outgrowth)
#psych::describeBy(Data.Outgrowth, group = Data.Outgrowth$Oxygen)

# Shapiro test
# if p<0.05  --> use non-parametric test (i.e. wilcoxon)
Shapiro_Outgrowth = shapiro.test(Data.Outgrowth$Outgrowth)
Shapiro_Outgrowth

# non-parametiric test
p.value_Oxygen_0 = 
             wilcox.exact(x=subset(Data.Outgrowth, subset=Oxygen=="21%")$Outgrowth,
                          y=subset(Data.Outgrowth, subset=Oxygen=="1%")$Outgrowth,
                          paired=T)$p.value #p=0.001953
p.value_Oxygen_1 = round(p.value_Oxygen_0, 3)
DF.summary_Oxygen = summarize((group_by(Data.Outgrowth, Oxygen)), 
                                Mean = mean(Outgrowth),
                                SD = sd(Outgrowth)) %>% as.data.frame() %>% 
  dplyr::mutate(text = paste0(Oxygen,":",formatC(Mean,digits=3),"(±",formatC(SD,digits=3),")"))
######### Location #########
p.value_Location_0 = 
             kruskal.test(Outgrowth ~ Site,
                          data=Data.Outgrowth)$p.value 
p.value_Location_1 = round(p.value_Location_0, 3) #p=0.722
DF.summary_Location = summarize((group_by(Data.Outgrowth, Site)), 
                           Mean = mean(Outgrowth),
                           SD = sd(Outgrowth)) %>% as.data.frame() %>% 
                      dplyr::mutate(text = paste0(Site,":",formatC(Mean,digits=3),"(±",formatC(SD,digits=3),")"))
######### NAT #########
p.value_NAT_0 = 
             kruskal.test(Outgrowth ~ NAT,
                          data=Data.Outgrowth)$p.value #p=0.26
p.value_NAT_1 = round(p.value_NAT_0, 3) #p=0.269
DF.summary_NAT = summarize((group_by(Data.Outgrowth,NAT)), 
                            Mean = mean(Outgrowth),
                            SD = sd(Outgrowth)) %>% as.data.frame() %>% 
                 dplyr::mutate(text = paste0(NAT,":",formatC(Mean,digits=3),"(±",formatC(SD,digits=3),")"))
# box plot
Data.Outgrowth_2 = 
  Data.Outgrowth %>% 
  dplyr::mutate(Oxygen = str_replace(Oxygen, pattern="21%", replacement="N_CAF"),
                Oxygen = str_replace(Oxygen, pattern="1%", replacement="H_CAF"), 
                Oxygen = factor(Oxygen, levels = c("N_CAF", "H_CAF")),
                Site = factor(Site, levels = c("Ph","Pb","Pt")),
                NAT = factor(NAT, levels = c("None", "GS", "GnP")))
ScaleYcontinu = 
  scale_y_continuous(
    limits=c(0, max(Data.Outgrowth$Outgrowth)+4),
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
Plot.CAFoutgrowth=
  ggplot(Data.Outgrowth_2, 
         aes(x=Oxygen, 
             y=Outgrowth, 
             fill=Oxygen)) +
  geom_boxplot(width = 0.45, fill = "white", staplewidth = 0.50, linewidth = 0.8, outliers = F) +
  #stat_summary(fun = "median", geom="crossbar", width=0.35, linewidth=0.7, color = "black") +
  geom_line(aes(group = No.), color = "gray60", linewidth=0.4) +
  geom_point(shape=21, 
             size=TextSize*0.4, 
             position=position_jitter(width=0.10, height=0)) +
  #geom_dotplot(binaxis="x", stackdir="center",
  #             binwidth=0.6,
  #             alpha=0.8) +
  geom_signif(comparisons = list(c("H_CAF","N_CAF")),
              map_signif_level = FALSE,
              y_position = 16, 
              annotation=paste0("**"),
              textsize=TextSize*1.5, fontface="bold", vjust=0.3, size=1.5) +
  labs(subtitle=NULL, caption=NULL) + CommonLabs + #paste0("p=",p.value_Oxygen_1," (Wilcoxon signed-rank test)"),
  scale_x_discrete(labels=c(H_CAF="Hypo-\nCAFs", N_CAF="Normo-\nCAFs")) +
  coord_fixed(ratio = 1 / 7.0) +
  scale_fill_manual(
    values=c(H_CAF="#91BFFA", N_CAF="#FEA0A0"),
    labels=c(H_CAF=expression(bold("H-CAF")), N_CAF=expression(bold("N-CAF"))))　+
  ScaleYcontinu + CommonTheme + theme(axis.text.x = element_text(color="transparent"))
#Plot.CAFoutgrowth
ggsave(plot=Plot.CAFoutgrowth,
       file="[3]Figures/[Figures]CAFoutgrowth_Plot4.png",
       width=7, height=7, dpi=400, bg="transparent")
saveRDS(Plot.CAFoutgrowth,
        "[3]Figures/[Figures]CAFoutgrowth_Plot.rds")

Plot.CAFoutgrowth_Site=
  Data.Outgrowth_2 %>% 
  ggplot(aes(y=Outgrowth, x=Site)) +
  geom_boxplot(width = .6, fill = "white", staplewidth = .6, linewidth = 0.8) +
  geom_dotplot(aes(fill=Site),
               binaxis="y", stackdir="center",
               binwidth=0.6,
               alpha=0.6) +
  geom_text(label=paste0("Kruskal–Wallis, p = ", p.value_Location_1), 
            x=Inf, y=18, hjust=1.05, size=TextSize*0.8, fontface="bold") +
  labs(subtitle="Cancer Location",
       caption=paste0("p=",p.value_Location_1," (Kruskal-Wallis test)\n",
                      paste(DF.summary_Location$text, collapse=", "))) + CommonLabs + 
  scale_fill_manual(values=c(Ph = "#1019BF",
                             Pb = "#38BDFF",
                             Pt = "#8FFFEB")) + 
  coord_fixed(ratio = 3 / 15) +
  ScaleYcontinu + CommonTheme
Plot.CAFoutgrowth_NAT =
  Data.Outgrowth_2 %>% 
  ggplot(aes(y=Outgrowth, x=NAT)) +
  geom_text(label=paste0("Kruskal–Wallis, p = ", p.value_NAT_1), 
            x=Inf, y=18, hjust=1.05, size=TextSize*0.8, fontface="bold") +
  geom_boxplot(width = .6, fill = "white", staplewidth = .6, linewidth = 0.8) +
  geom_dotplot(aes(fill=NAT), binaxis="y", stackdir="center",
               binwidth=0.6,
               alpha=0.6) +
  labs(subtitle="NAT",
       caption=paste0("p=",p.value_NAT_1," (Kruskal-Wallis test)\n",
                      paste(DF.summary_NAT$text, collapse=", "))) + CommonLabs +
  scale_fill_manual(values=c(None = "#BDBDBD",
                             GS = "pink",
                             GnP = "orange")) +
  coord_fixed(ratio = 3 / 15) +
  ScaleYcontinu + CommonTheme
Plot.Supp = 
  (Plot.CAFoutgrowth_NAT | Plot.CAFoutgrowth_Site) &
  theme(plot.background = element_rect(fill="transparent", color=NA),
        panel.background = element_rect(fill="white", color=NA))
ggsave(plot=Plot.Supp,
       file="[3]Figures/[Figures]CAFoutgrowth_Plot5_NATandLocation.png",
       width=16, height=7, dpi=400)
