# ggplot2,Hmisc,cowplot,RColorBrewerパッケージの読み込み
library(ggplot2)
library(Hmisc)

# dataの読み込み
data_all=read.table(file="[1]dataset/proliferation_all.csv", header=T, sep=",", stringsAsFactors=F)

data_1HvsN_0207 = subset(data_all,subset= patient=="24-0207" & oxygen%in%c("h","n"))
data_1HvsN_0110 = subset(data_all,subset= patient=="24-0110" & oxygen%in%c("h","n"))
data_1HvsN_1122 = subset(data_all,subset= patient=="24-1122" & oxygen%in%c("h","n"))


# Parameter
Errorbar_linewidth = 1.0
Errorbar_width = 8
Line_width = 1.2
BackGroundColor_1pct = "#EAF7FF"
BackGroundColor_21pct = "#FFEFEF"
ColorCode_H = "#607DF0"
ColorCode_HN = "#7FBFFF"
ColorCode_NH = "#FFA3A3" #FFDDDCもいい
ColorCode_N = "#FF7F7F"
PointSize = 3.5
Alpha = 0.9
pd = position_dodge( width = 3.0 )
ScaleYLimits = c(0,20)
ScaleYBreaks = seq(0,20,by=5)     # breaks=seq([軸の始点]:[軸の終点],by[目盛幅]
ScaleYMinorBreaks = seq(0,20,by=2.5)
GridColor = "gray90"

# graphの作成

graph =
  ggplot(data_1HvsN_0207, aes(x=time, y=nl, color=oxygen, shape=oxygen, group=oxygen)) +
  stat_summary(geom = "line", 
               fun = "mean",
               alpha = Alpha,
               linewidth = Line_width, 
               position = pd ) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_normal, 
               alpha = Alpha,
               linewidth = Errorbar_linewidth, 
               width = Errorbar_width, 
               position = pd ) + 
  stat_summary(fun="mean", geom="point", size = PointSize, alpha = Alpha, position = pd ) +
  scale_x_continuous(limits = c(0, 192),
                     breaks = seq(0, 168, by=24),
                     minor_breaks = NULL, expand = c(0,0))+
  scale_y_continuous(limits = ScaleYLimits, breaks = ScaleYBreaks, 
                     minor_breaks = ScaleYMinorBreaks, expand = c(0,0)) +
  scale_shape_manual(values = c("circle", "square")) +
  scale_colour_manual(values = c(ColorCode_H, ColorCode_N)) +
  labs(title = "CAF24-0207",
       x = "Elapsed time(hr)", y = "Normalized luminescence")+
  theme(axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = GridColor),
        axis.line = element_line(),
        legend.title=element_blank()) #凡例titleの非表示

graph
