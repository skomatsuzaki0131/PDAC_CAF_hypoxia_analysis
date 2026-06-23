# ggplot2,Hmisc,cowplot,RColorBrewerパッケージの読み込み
library(ggplot2)
library(dplyr)
library(ggpubr)

# Parameter
Dir.OutProlife = "/Volumes/PortableSSD/[4]myR/[3]figures/"
List.title = list(
  Pair9 = "Pair #9",
  Pair12 = "Pair #12",
  Pair13 = "Pair #13")
Param = list(
  Ylab = expression(bold(""*log[2]*"(NL)") ) ,
  Errorbar_linewidth = 1.0,
  Errorbar_width = 8,
  Line_width = 1.2,
  ColCode_H = "#607DF0",
  ColCode_HN = "#7FBFFF",
  ColCode_NH = "#FFA3A3", 
               #"#FFDDDC"もいい
  ColCode_N = "#FF7F7F",
  PtSize = 3.5,
  Alpha = 0.9,
  pd = position_dodge( width = 2.0 ),
  XLim = 218,
  YLim = c(0, 4.5),
  YBreak = seq(0, 4, by = 1),    # breaks=seq([軸の始点]:[軸の終点],by[目盛幅]
  GridCol = "gray90",
  SigMkRight = 178,
  SigMkXwidth = 3,
  AsteriskDist1 = 3,
  AsteriskDist2 = 9,
  AsteriskSize = 10,
  NSSize = 5)
List.Col = list(Simple = c(h="#607DF0", n="#FF7F7F"),
                Under1 = c(h="#607DF0", nh="#FFA3A3"),
                Under21 = c(hn="#7FBFFF", n="#FF7F7F"),
                SwitchH = c(h="#607DF0", hn="#7FBFFF"),
                SwitchN = c(nh="#FFA3A3", n="#FF7F7F"))
List.LgdLab = list(Simple = c(h="H-CAF", n="N-CAF"),
                   Under1 = c(h="H-CAF", nh="N-CAF\n  under 1% O2"),
                   Under21 = c(hn="H-CAF\n under  21% O2", n="N-CAF"),
                   SwitchH = c(h="H-CAF", hn="H-CAF\n  under 21% O2"),
                   SwitchN = c(nh="N-CAF\n  under 1% O2", n="N-CAF"))
List.BGCol = list(Simple = adjustcolor("white", alpha.f=0.8),
                  Under1 = adjustcolor("#EAF7FF", alpha.f=0.8),
                  Under21 = adjustcolor("#FFEFEF", alpha.f=0.8),
                  SwitchH = adjustcolor("white", alpha.f=0.8),
                  SwitchN = adjustcolor("white", alpha.f=0.8))
Scale.Settings.X = scale_x_continuous(limits=c(0, Param[["XLim"]]),
                                      breaks=seq(0, 168, by=24),
                                      minor_breaks=NULL, 
                                      expand=expansion(mult = c(0,0)))
Scale.Settings.Y = scale_y_continuous(limits=Param[["YLim"]],
                                      breaks=Param[["YBreak"]],
                                      oob=scales::squish,
                                      expand=expansion(mult=c(0,0)))



##########################################################
##########################################################
##                                                      ##
##            1.  Simple comparison                     ##
##            2.  Under 1% O2 comparison                ##
##            3.  Under 21% O2 comparison               ##
##                                                      ##
##########################################################
##########################################################
# dataの読み込み
data_all=read.csv(file="/Volumes/PortableSSD/[4]myR/[1]dataset/proliferation_all.csv", header=T, sep=",", stringsAsFactors=F) %>% 
  rowwise %>% 
  mutate(log2nl = log(nl, 2)) %>% ungroup() %>%
  mutate(Day = rep(c("Day1","Day1","Day1","Day2","Day2","Day2",
                     "Day3","Day3","Day3","Day4","Day4","Day4",
                     "Day5","Day5","Day5","Day6","Day6","Day6",
                     "Day7","Day7","Day7"), times=12))

Func.ProPlot = function(Pair, Comp){
  O2Condition = case_when(Comp=="Simple" ~ c("h","n"),
                          Comp=="Under1" ~ c("h","nh"),
                          Comp=="Under21" ~ c("hn","n"),
                          Comp=="SwitchH" ~ c("h", "hn"),
                          Comp=="SwitchN" ~ c("nh", "n"))
  List.Data = list(
    Pair9 = subset(data_all,subset= patient=="23-1122" & oxygen%in%O2Condition),
    Pair12 = subset(data_all,subset= patient=="24-0110" & oxygen%in%O2Condition),
    Pair13 = subset(data_all,subset= patient=="24-0207" & oxygen%in%O2Condition))
  DataFrame = List.Data[[Pair]]
  Signif = t.test(x = DataFrame$log2nl[19:21],
                  y = DataFrame$log2nl[40:42])$p.value
  ggplot(DataFrame, aes(x=time, y=log2nl, color=oxygen)) +
    stat_summary(geom="line", fun="mean",
                 alpha=Param[["Alpha"]], 
                 linewidth=Param[["Line_width"]], 
                 position=Param[["pd"]] ) +
    stat_summary(geom="errorbar", 
                 fun.data=mean_cl_normal, 
                 alpha=Param[["Alpha"]],
                 linewidth=Param[["Errorbar_linewidth"]],
                 width=Param[["Errorbar_width"]], 
                 position=Param[["pd"]] ) + 
    stat_summary(geom="point",
                 fun="mean",
                 size=Param[["PtSize"]],
                 alpha=Param[["Alpha"]],
                 position=Param[["pd"]] ) +
    labs(title = paste0(List.title[[Pair]],"  ", case_when(Signif<0.001 ~ paste0("***p<0.001"), 
                                                           Signif<0.01 ~ paste0("**p=",round(Signif, 4)), 
                                                           Signif<0.05 ~ paste0("*p=",round(Signif, 4)), 
                                                           TRUE ~ paste0("p=",round(Signif, 4),"(n.s.)"))),
         x="Elapsed time(hour)", y=Param[["Ylab"]], color=NULL)+
    Scale.Settings.X + Scale.Settings.Y +
    scale_colour_manual(values=List.Col[[Comp]],
                        labels=List.LgdLab[[Comp]]) + 
    theme(text = element_text(face="bold"),
          legend.key = element_blank(),
          panel.background = element_rect(fill = List.BGCol[[Comp]]),
          panel.grid = element_blank(), #element_line(color = Param[["GridCol"]]),
          #panel.grid.minor.y = element_blank(),
          axis.line = element_line())+
    # vertical line
    annotate("segment", 
             x = Param[["SigMkRight"]], 
             y = mean(DataFrame$log2nl[19:21]), 
             yend = mean(DataFrame$log2nl[40:42])) +
    # upper horizontal line
    annotate("segment",
             x = Param[["SigMkRight"]],
             xend = Param[["SigMkRight"]]-Param[["SigMkXwidth"]], 
             y = mean(DataFrame$log2nl[19:21])) +
    # lower horizontal line
    annotate("segment", 
             x = Param[["SigMkRight"]], 
             xend = Param[["SigMkRight"]]-Param[["SigMkXwidth"]],
             y = mean(DataFrame$log2nl[40:42])) +
    # significance
    geom_text(x = Param[["SigMkRight"]]+Param[["AsteriskDist1"]], 
              label = case_when(Signif<0.001 ~ "***", 
                                Signif<0.01 ~ "**", 
                                Signif<0.05 ~ "*", 
                                TRUE ~ "ns"),
              hjust=0,
              color="black",
              size=ifelse(Signif<0.05, Param[["AsteriskSize"]], Param[["NSSize"]]),
              y=(mean(DataFrame$log2nl[19:21]) + mean(DataFrame$log2nl[40:42]))/2-0.08 ) 
}

Plot1.Pair9 = Func.ProPlot(Pair="Pair9", Comp="Simple")
Plot1.Pair12 = Func.ProPlot(Pair="Pair12", Comp="Simple")
Plot1.Pair13 = Func.ProPlot(Pair="Pair13", Comp="Simple")
Plot2.Pair9 = Func.ProPlot(Pair="Pair9", Comp="Under1")
Plot2.Pair12 = Func.ProPlot(Pair="Pair12", Comp="Under1")
Plot2.Pair13 = Func.ProPlot(Pair="Pair13", Comp="Under1")
Plot3.Pair9 = Func.ProPlot(Pair="Pair9", Comp="Under21")
Plot3.Pair12 = Func.ProPlot(Pair="Pair12", Comp="Under21")
Plot3.Pair13 = Func.ProPlot(Pair="Pair13", Comp="Under21")
Plot4.Pair9 = Func.ProPlot(Pair="Pair9", Comp="SwitchH")
Plot4.Pair12 = Func.ProPlot(Pair="Pair12", Comp="SwitchH")
Plot4.Pair13 = Func.ProPlot(Pair="Pair13", Comp="SwitchH")
Plot5.Pair9 = Func.ProPlot(Pair="Pair9", Comp="SwitchN")
Plot5.Pair12 = Func.ProPlot(Pair="Pair12", Comp="SwitchN")
Plot5.Pair13 = Func.ProPlot(Pair="Pair13", Comp="SwitchN")
NL = theme(legend.position="none")
Lgd1 = ggpubr::get_legend(Plot1.Pair9+theme(legend.position = c(0,0.5), legend.justification = c(0, 0.5)))
Lgd2 = ggpubr::get_legend(Plot2.Pair9+theme(legend.position = c(0,0.5), legend.justification = c(0, 0.5)))
Lgd3 = ggpubr::get_legend(Plot3.Pair9+theme(legend.position = c(0,0.5), legend.justification = c(0, 0.5)))
Lgd4 = ggpubr::get_legend(Plot4.Pair9+theme(legend.position = c(0,0.5), legend.justification = c(0, 0.5)))
Lgd5 = ggpubr::get_legend(Plot5.Pair9+theme(legend.position = c(0,0.5), legend.justification = c(0, 0.5)))

Plot = ggarrange(Plot1.Pair9+NL, Plot1.Pair12+NL, Plot1.Pair13+NL, Lgd1,
                    Plot3.Pair9+NL, Plot3.Pair12+NL, Plot3.Pair13+NL, Lgd3,
                    Plot2.Pair9+NL, Plot2.Pair12+NL, Plot2.Pair13+NL, Lgd2,
                    Plot5.Pair9+NL, Plot5.Pair12+NL, Plot5.Pair13+NL, Lgd5,
                    Plot4.Pair9+NL, Plot4.Pair12+NL, Plot4.Pair13+NL, Lgd4,
                    ncol=4, nrow=5, widths=c(3, 3, 3, 1.5)) %>% 
       annotate_figure(bottom=text_grob("NL: Luminescence normalized to that at 24 hours."))
ggsave(plot=Plot, file=paste0(Dir.OutProlife, "[Figure][Proliferation].png"),
       width=12, height=15, dpi=300, bg="white")

data_day7_0 = subset(data_all, subset=Day=="Day7")
data_day7 = mutate(data_day7_0,
                   Sample = paste0(data_day7_0$patient, data_day7_0$oxygen))
DF.summary = summarize(data_day7 %>% group_by(Sample),
                       Mean = mean(log2nl),
                       SD = sd(log2nl))
#   Sample     Mean  SD
#1  23-1122h   2.59  0.109 
#2  23-1122hn  3.01  0.0649
#3  23-1122n   1.65  0.172 
#4  23-1122nh  1.20  0.0611
#5  24-0110h   2.93  0.0198
#6  24-0110hn  3.46  0.0555
#7  24-0110n   1.66  0.153 
#8  24-0110nh  1.52  0.252 
#9  24-0207h   3.62  0.140 
#10 24-0207hn  4.05  0.0448
#11 24-0207n   3.29  0.117 
#12 24-0207nh  2.88  0.207 
Func.SignifText = function(SAMPLE){
    p_value_simple = t.test(x = subset(data_all, subset=patient==SAMPLE & oxygen=="h" & Day=="Day7")$log2nl,
                            y = subset(data_all, subset=patient==SAMPLE & oxygen=="n" & Day=="Day7")$log2nl)$p.value
    p_vaue_under21 = t.test(x = subset(data_all, subset=patient==SAMPLE & oxygen=="hn" & Day=="Day7")$log2nl,
                            y = subset(data_all, subset=patient==SAMPLE & oxygen=="n" & Day=="Day7")$log2nl)$p.value
    p_value_under1 = t.test(x = subset(data_all, subset=patient==SAMPLE & oxygen=="h" & Day=="Day7")$log2nl,
                            y = subset(data_all, subset=patient==SAMPLE & oxygen=="nh" & Day=="Day7")$log2nl)$p.value
    p_value_SwitchN = t.test(x = subset(data_all, subset=patient==SAMPLE & oxygen=="n" & Day=="Day7")$log2nl,
                             y = subset(data_all, subset=patient==SAMPLE & oxygen=="nh" & Day=="Day7")$log2nl)$p.value
    p_value_SwitchH = t.test(x = subset(data_all, subset=patient==SAMPLE & oxygen=="h" & Day=="Day7")$log2nl,
                             y = subset(data_all, subset=patient==SAMPLE & oxygen=="hn" & Day=="Day7")$log2nl)$p.value
    signif_simple = case_when(p_value_simple<0.001 ~ "***", p_value_simple<0.01 ~ "**", p_value_simple<0.05 ~ "*", TRUE ~ "(n.s.)")
    signif_under21 = case_when(p_vaue_under21<0.001 ~ "***", p_vaue_under21<0.01 ~ "**", p_vaue_under21<0.05 ~ "*", TRUE ~ "(n.s.)")
    signif_under1 = case_when(p_value_under1<0.001 ~ "***", p_value_under1<0.01 ~ "**", p_value_under1<0.05 ~ "*", TRUE ~ "(n.s.)")
    signif_SwitchN = case_when(p_value_SwitchN<0.001 ~ "***", p_value_SwitchN<0.01 ~ "**", p_value_SwitchN<0.05 ~ "*", TRUE ~ "(n.s.)")
    signif_SwitchH = case_when(p_value_SwitchH<0.001 ~ "***", p_value_SwitchH<0.01 ~ "**", p_value_SwitchH<0.05 ~ "*", TRUE ~ "(n.s.)")
    paste0(
      "p(Simple)=", formatC(p_value_simple,digits=2),signif_simple,",\n",
      "p(Under21%)=", formatC(p_vaue_under21,digits=2),signif_under21,", ",
      "p(Under1%)=", formatC(p_value_under1,digits=2),signif_under1,",\n",
      "p(nvsnh)=", formatC(p_value_SwitchN,digits=2),signif_SwitchN,", ",
      "p(hvshn)=", formatC(p_value_SwitchH,digits=2),signif_SwitchH )
    }
List.Summary = list("23-1122" = paste0("H:Mean",formatC(DF.summary$Mean[1],digits=3),"(SD", formatC(DF.summary$SD[1],digits=3),"), ",
                                      "HN:Mean",formatC(DF.summary$Mean[2],digits=3),"(SD", formatC(DF.summary$SD[2],digits=3),")\n",
                                      "N:Mean",formatC(DF.summary$Mean[3],digits=3),"(SD", formatC(DF.summary$SD[3],digits=3),"), ",
                                      "NH:Mean",formatC(DF.summary$Mean[4],digits=3),"(SD", formatC(DF.summary$SD[4],digits=3),")\n",
                                      Func.SignifText(SAMPLE="23-1122") ), 
                    "24-0110" = paste0("H:Mean",formatC(DF.summary$Mean[5],digits=3),"(SD", formatC(DF.summary$SD[5],digits=3),"), ",
                                      "HN:Mean",formatC(DF.summary$Mean[6],digits=3),"(SD", formatC(DF.summary$SD[6],digits=3),")\n",
                                      "N:Mean",formatC(DF.summary$Mean[7],digits=3),"(SD", formatC(DF.summary$SD[7],digits=3),"), ",
                                      "NH:Mean",formatC(DF.summary$Mean[8],digits=3),"(SD", formatC(DF.summary$SD[8],digits=3),")\n",
                                      Func.SignifText(SAMPLE="24-0110") ),
                    "24-0207" = paste0("H:Mean",formatC(DF.summary$Mean[9],digits=3),"(SD", formatC(DF.summary$SD[9],digits=3),"), ",
                                      "HN:Mean",formatC(DF.summary$Mean[10],digits=3),"(SD", formatC(DF.summary$SD[10],digits=3),")\n",
                                      "N:Mean",formatC(DF.summary$Mean[11],digits=3),"(SD", formatC(DF.summary$SD[11],digits=3),"), ",
                                      "NH:Mean",formatC(DF.summary$Mean[12],digits=3),"(SD", formatC(DF.summary$SD[12],digits=3),")\n",
                                      Func.SignifText(SAMPLE="24-0207") ) )
###############################
data_all_2 = data_all
data_all_2$oxygen = factor(data_all_2$oxygen, levels=c("hn","h","n","nh"))
Func.ProPlot_2 = function(SAMPLE){
  ggplot(subset(data_all_2, subset=patient==SAMPLE), 
         aes(x=time, y=log2nl, color=oxygen)) +
  stat_summary(geom="line", fun="mean",
               alpha=Param[["Alpha"]], 
               linewidth=Param[["Line_width"]], 
               position=Param[["pd"]] ) +
  stat_summary(geom="errorbar", 
               fun.data=mean_cl_normal, 
               alpha=Param[["Alpha"]],
               linewidth=Param[["Errorbar_linewidth"]],
               width=Param[["Errorbar_width"]], 
               position=Param[["pd"]] ) + 
  stat_summary(geom="point",
               fun="mean",
               size=Param[["PtSize"]],
               alpha=Param[["Alpha"]],
               position=Param[["pd"]] ) +
  labs(x="Elapsed time(hour)", y=Param[["Ylab"]], color=NULL,
       caption=List.Summary[[SAMPLE]]) +
  guides(color = guide_legend(override.aes = list(linetype = NA))) +
  facet_wrap( ~ patient) +
  Scale.Settings.X + Scale.Settings.Y +
  scale_color_manual(values=c(h = "#607DF0",
                              hn = "#7FBFFF",
                              nh = "#FFA3A3", 
                              n = "#FF7F7F"),
                     labels = c(h="H-CAF", n="N-CAF",
                                nh="N-CAF\n  under 1% O2",
                                hn="H-CAF\n under  21% O2")) +
  coord_fixed(ratio=48) +
  theme(text = element_text(face="bold", size=15),
        legend.key = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(), #element_line(color = Param[["GridCol"]]),
        strip.text = element_text(face = "bold"),
        strip.background = element_rect(color=NULL, fill="white"),
        axis.line = element_line()) }
NL = theme(legend.position="none")
ggsave(plot=ggarrange(Func.ProPlot_2(SAMPLE="23-1122")+NL,
                      Func.ProPlot_2(SAMPLE="24-0110")+NL,
                      Func.ProPlot_2(SAMPLE="24-0207")+NL,
                      ggpubr::get_legend(Func.ProPlot_2(SAMPLE="23-1122")+
                                         theme(legend.position = c(0,0.5), legend.justification = c(0, 0.5))),
                      ncol=4, widths=c(3,3,3,2), nrow=1),
       file=paste0(Dir.OutProlife, "[Figure][Proliferation]_4group_3.png"),
       width=15, height=5, dpi=400, bg="white")




################

ggplot(subset(data_all, subset=patient=="23-1122"), 
       aes(x=time, y=log2nl, color=oxygen)) +
  stat_summary(geom="line", fun="mean",
               alpha=Param[["Alpha"]], 
               linewidth=Param[["Line_width"]], 
               position=Param[["pd"]] ) 
  stat_summary(geom="errorbar", 
               fun.data=mean_cl_normal, 
               alpha=Param[["Alpha"]],
               linewidth=Param[["Errorbar_linewidth"]],
               width=Param[["Errorbar_width"]], 
               position=Param[["pd"]] ) + 
  stat_summary(geom="point",
               fun="mean",
               size=Param[["PtSize"]],
               alpha=Param[["Alpha"]],
               position=Param[["pd"]] ) +
  labs(x="Elapsed time(hour)", y=Param[["Ylab"]], color=NULL,
       caption=List.Summary[["23-1122"]]) +
  guides(color = guide_legend(override.aes = list(linetype = NA))) +
  coord_fixed(ratio=40) +
  Scale.Settings.X + Scale.Settings.Y +
  scale_color_manual(values=c(h = "#607DF0",
                              hn = "#7FBFFF",
                              nh = "#FFA3A3", 
                              n = "#FF7F7F"),
                     labels = c(h="H-CAF", n="N-CAF",
                                nh="N-CAF\n  under 1% O2",
                                hn="H-CAF\n under  21% O2")) +
  theme(text = element_text(face="bold"),
        legend.key = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(), #element_line(color = Param[["GridCol"]]),
        axis.line = element_line())