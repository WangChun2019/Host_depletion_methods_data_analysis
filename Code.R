########## library packages ##########
library(readxl)# necessary to import the data from Excel file
library(readr)
library(dplyr)
library(tidyverse)
library(reshape2)

library(pheatmap)
library(ggplot2)
library(grid) # for combining plots
library(cowplot)
library(ggstatsplot)
library(ggpubr)
library(RColorBrewer)
library(paletteer) # color

library(scales) # to access break formatting functions

library(phyloseq)
library(vegan)

#install.packages("picante")
library(picante)
library(ggplot2)
library(readxl)# necessary to import the data from Excel file
library(data.table)
library(readr)
library(dplyr)
library(tidyverse)
library(stringr)
library(reshape2)
library(pheatmap)
library(ggplot2)
library(ggpmisc)#拟合曲线加方程和p值
library(ggbiplot)
library(ggrepel)
library(grid) # for combining plots
library(cowplot)
theme_set(theme_cowplot())
library(ggstatsplot)
library(ggpubr)
library(ggalluvial)
library(RColorBrewer)
library(paletteer) # color
library(ggsci)
library(scales) # to access break formatting functions
library(phyloseq)
library(vegan)
library(philentropy)#JSD距离
getDistMethods()#展示这46种距离
library(ade4)
#加载α多样性需要的包
library(agricolae)
library(corrplot)#speasman rank correlation作图 https://mp.weixin.qq.com/s/qzlMMaN8JwzHlNC5EPhC0g
library(decontam)#去污染软件decontam
library(picante)
#icc重复性相关R包
#install.packages("irr")
#install.packages("lpSolve")
library(lpSolve)
library(irr)
library(car)#正态分布函数
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(decontam); packageVersion("decontam")
library(patchwork)

library("plyr")
library("lattice")
library("ggplot2")
library("dplyr")
library("readr")
library("rmarkdown")
# install.packages("Rmisc")
library("Rmisc")
library("devtools")
# install.packages("rlang")
library(rlang)
# install.packages("gghalves")
library("gghalves")
library(ggplot2)
library(gghalves)
# Fig. 1---------------------
setwd("F:/Code")
plotdf <- read_excel("data.xlsx",sheet = "mass") %>% as.data.frame()
sig <- compare_means(bac_load~sample_type1, plotdf, 
                             method = "wilcox.test",
                             p.adjust.method = "fdr",
                             paired = FALSE)
set.seed(1)
plotdf$sample_type.jitter <- jitter(as.numeric(plotdf$sample_type1 %>% as.factor()),amount =0.05)
# 设置扰动点在图形中的位置，以绘制连接线
plotdf$sample_type.jitter.line <- 0
plotdf$sample_type.jitter.line[which(plotdf$sample_type1 =="BALF")] <- plotdf$sample_type.jitter[which(plotdf$sample_type1 == "BALF")] +0.15
plotdf$sample_type.jitter.line[which(plotdf$sample_type1 =="OP")] <- plotdf$sample_type.jitter[which(plotdf$sample_type1 == "OP")] - 0.15
dat_1 <- plotdf[which(plotdf$sample_type1 =="BALF"), ]
dat_2 <- plotdf[which(plotdf$sample_type1 =="OP"), ]
# Fig 1A
p1 <- ggplot(plotdf, aes(x = sample_type1, y = bac_load/1000, 
                   fill = sample_type1))+
  geom_half_violin(data = dat_1,
                   aes(x = sample_type1, y = bac_load/1000),
                   side ="l", trim =F, color =NA, alpha =0.75,
                   position = position_nudge(x = -0.15, y =0))+
  geom_half_violin(data = dat_2,
                   aes(x = sample_type1, y = bac_load/1000),
                   side ="r", trim =F, color =NA, alpha =0.75,
                   position = position_nudge(x =0.15, y =0))+
  geom_boxplot(width =0.15,size=0.3, alpha =0.75,outlier.colour = NA)+
  geom_line(aes(x = sample_type.jitter.line, group = sampleNo1), 
            color ="lightgrey",linewidth=0.2)+
    geom_point(data = dat_1,
               aes(x = sample_type.jitter +0.15, y = bac_load/1000, 
                   color = sample_type1),size =0.5, shape =19)+
    geom_point(data = dat_2,
               aes(x = sample_type.jitter -0.15, y = bac_load/1000, 
                   color = sample_type1),size =0.5, shape =19)+
    scale_fill_manual(limits = base::unique(plotdf$sample_type1),
                      values = c(BALF="#2B3A55",OP="#CE7777"))+
    scale_color_manual(limits = base::unique(plotdf$sample_type1),
                       values = c(BALF="#2B3A55",OP="#CE7777"))+
  theme_bw()+
  scale_y_log10(limits=c(0.00001,1000000),
                breaks=c(0.00001,0.001,0.1,10,1000,100000,1000000),
                labels=trans_format("log10",math_format(10^.x)))+
    labs(x='',y="Bacterial DNA (ng/ml BALF or ng/swab)");p1 
# Fig 1B
sig <- compare_means(host_load~sample_type1, plotdf, 
                             method = "wilcox.test",
                             p.adjust.method = "fdr",
                             paired = FALSE) 
p2 <- ggplot(plotdf, aes(x = sample_type1, y = host_load/1000, 
                         fill = sample_type1))+
  geom_half_violin(data = dat_1,
                   aes(x = sample_type1, y = host_load/1000),
                   side ="l", trim =F, color =NA, alpha =0.75,
                   position = position_nudge(x = -0.15, y =0))+
  geom_half_violin(data = dat_2,
                   aes(x = sample_type1, y = host_load/1000),
                   side ="r", trim =F, color =NA, alpha =0.75,
                   position = position_nudge(x =0.15, y =0))+
  geom_boxplot(width =0.15,size=0.3, alpha =0.75,outlier.colour = NA)+
  geom_line(aes(x = sample_type.jitter.line, group = sampleNo1), 
            color ="lightgrey",linewidth=0.2)+
    geom_point(data = dat_1,
               aes(x = sample_type.jitter +0.15, y = host_load/1000, 
                   color = sample_type1),size =0.5, shape =19)+
    geom_point(data = dat_2,
               aes(x = sample_type.jitter -0.15, y = host_load/1000, 
                   color = sample_type1),size =0.5, shape =19)+
  scale_fill_manual(limits = base::unique(plotdf$sample_type1),
                      values = c(BALF="#2B3A55",OP="#CE7777"))+
  scale_color_manual(limits = base::unique(plotdf$sample_type1),
                       values = c(BALF="#2B3A55",OP="#CE7777"))+
  theme_bw()+
  scale_y_log10(limits=c(0.01,10000000),
                breaks=c(0.01,1,100,10000,1000000,10000000),
                labels=trans_format("log10",math_format(10^.x)))+
  labs(x='',y="Human DNA (ng/ml BALF or ng/swab)");p2
# Fig 1C
plotdf <- read_excel("data.xlsx",sheet = "microbe_host_ratio") %>% as.data.frame()
sig <- compare_means(ratio~sample_type1, plotdf, 
                     method = "wilcox.test",
                     p.adjust.method = "fdr",
                     paired = FALSE)
set.seed(1)
plotdf$sample_type.jitter <- jitter(as.numeric(plotdf$sample_type1 %>% as.factor()),amount =0.05)
# 设置扰动点在图形中的位置，以绘制连接线
plotdf$sample_type.jitter.line <- 0
plotdf$sample_type.jitter.line[which(plotdf$sample_type1 =="BALF")] <- 
  plotdf$sample_type.jitter[which(plotdf$sample_type1 == "BALF")] +0.15
plotdf$sample_type.jitter.line[which(plotdf$sample_type1 =="OP")] <- 
  plotdf$sample_type.jitter[which(plotdf$sample_type1 == "OP")] - 0.15
dat_1 <- plotdf[which(plotdf$sample_type1 =="BALF"), ]
dat_2 <- plotdf[which(plotdf$sample_type1 =="OP"), ]

p3 <- ggplot(plotdf, aes(x = sample_type1, y = ratio, 
                         fill = sample_type1))+
  geom_half_violin(data = dat_1,
                   aes(x = sample_type1, y = ratio),
                   side ="l", trim =F, color =NA, alpha =0.75,
                   position = position_nudge(x = -0.15, y =0))+
  geom_half_violin(data = dat_2,
                   aes(x = sample_type1, y = ratio),
                   side ="r", trim =F, color =NA, alpha =0.75,
                   position = position_nudge(x =0.15, y =0))+
  geom_boxplot(width =0.15, alpha =0.75,outlier.colour = NA,
               size=0.3)+
  geom_line(aes(x = sample_type.jitter.line, group = sampleNo1), 
            color ="lightgrey",linewidth=0.2)+
  geom_point(data = dat_1,
             aes(x = sample_type.jitter +0.15, y = ratio, 
                 color = sample_type1),size =0.5, shape =19)+
  geom_point(data = dat_2,
             aes(x = sample_type.jitter -0.15, y = ratio, 
                 color = sample_type1),size =0.5, shape =19)+
  scale_fill_manual(values = c(BALF="#2B3A55",OP="#CE7777"))+
  scale_color_manual(values = c(BALF="#2B3A55",OP="#CE7777"))+
  theme_bw()+
  scale_y_log10()+
  labs(x='',y="Microbe-to-host read ratio")+
  scale_y_log10(limits=c(0.000001,100),
                breaks=c(0.000001,0.00001,0.0001,0.001,0.01,0.1,1,10,100,1000));p3
# Fig S1-----------------------------------------------------------
setwd("F:/Code")
#Fig S1A
df <- read_excel("data.xlsx",sheet = "optimise1") %>% as.data.frame() %>% 
  filter(exp=="PMA_con"&metrics=="human") %>% .[order(.$sampleID)&grepl("BA",.$sampleID),]
wilcox <- compare_means(value~group, df, method = "wilcox.test",
                        p.adjust.method = "fdr",paired = TRUE) %>% 
  mutate(y.position = seq(log10(max(df$value))*1.05,log10(max(df$value)), 
                          length.out = nrow(.)))
df$group <- factor(df$group,levels = c("10uM","50uM"))
plot <- ggplot(df,aes(group,value,group=sampleID))+
  geom_point(size=1)+
  geom_line(linewidth=0.5)+theme_classic()+
  labs(x='',y='pg/ml BALF',title='Human DNA')+
  theme(plot.title = element_text(hjust = 0.5))+
  annotation_logticks(sides = "l",outside = TRUE)+
  coord_cartesian(clip = "off")+
  scale_y_log10(limits=c(100,10000000),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels=trans_format("log10",math_format(10^.x)))+
  theme(axis.text = element_text(size = 7.5,color = "black"),
        axis.line = element_line(linewidth=0.5),
        axis.ticks=element_line(color="black",linewidth=0.5),
        axis.title = element_text(siz=9,color = "black"),
        plot.title = element_text(size = 9,color = "black"),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA));plot
#Fig S1B
df <- read_excel("data.xlsx",sheet = "optimise1") %>% as.data.frame() %>% 
  filter(exp=="PMA_con"&metrics=="bac") %>% 
  .[order(.$sampleID)&grepl("BA",.$sampleID),]
wilcox <- compare_means(value~group, df, method = "wilcox.test",
                        p.adjust.method = "fdr",paired = TRUE) %>% 
  filter(p.adj < 0.05) %>% 
  mutate(y.position = seq(log10(max(df$value))*1.05,log10(max(df$value)), 
                          length.out = nrow(.)))
df$group <- factor(df$group,levels = c("10uM","50uM"))
plot <- ggplot(df,aes(group,value,group=sampleID))+
  geom_point(size=0.5)+
  geom_line(linewidth=0.3)+theme_classic()+
  labs(x='',y='pg/ml BALF',title='Bacterial DNA')+
  theme(plot.title = element_text(hjust = 0.5))+
  annotation_logticks(sides = "l",outside = TRUE)+
  coord_cartesian(clip = "off")+
  scale_y_log10(limits=c(1,10000),breaks = c(1,10,100,1000,10000),
                labels=trans_format("log10",math_format(10^.x)))+
  stat_pvalue_manual(wilcox, label = "p.signif", tip.length = 0,
                     vjust=0, size=5)+
  theme(axis.text = element_text(size = 7,color = "black"),
        axis.line = element_line(linewidth=0.5),
        axis.ticks=element_line(color="black",linewidth=0.5),
        axis.title = element_text(siz=9,color = "black"),
        plot.title = element_text(size = 9,color = "black"),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA));plot
#Fig S1C
df <- read_excel("data.xlsx",sheet = "optimise1") %>% as.data.frame() %>% 
  filter(exp=="PMA_con"&metrics=="bac_prop") %>% .[order(.$sampleID)&grepl("BA",.$sampleID),]
wilcox <- compare_means(value~group, df, method = "wilcox.test",
                        p.adjust.method = "fdr",paired = TRUE) %>% 
  mutate(y.position = seq(log10(max(df$value))*1.05,log10(max(df$value)), 
                          length.out = nrow(.)))
df$group <- factor(df$group,levels = c("10uM","50uM"))
plot <- ggplot(df,aes(group,value,group=sampleID))+
  geom_point(size=1)+
  geom_line(linewidth=0.5)+theme_classic()+
  labs(x='',y='',title='Bacterial DNA percentage (%)')+
  theme(plot.title = element_text(hjust = 0.5))+
  annotation_logticks(sides = "l",outside = TRUE)+
  coord_cartesian(clip = "off")+
  scale_y_log10(limits=c(0.01,100),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels=trans_format("log10",math_format(10^.x)))+
  theme(axis.text = element_text(size = 7.5,color = "black"),
        axis.line = element_line(linewidth=0.5),
        axis.ticks=element_line(color="black",linewidth=0.5),
        axis.title = element_text(siz=9,color = "black"),
        plot.title = element_text(size = 9,color = "black"),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA));plot
#Fig S1D 
df <- read_excel("optimise.xlsx",sheet = "optimise1") %>% as.data.frame() %>% 
  filter(exp=="PMA_con"&metrics=="human") %>% .[order(.$sampleID)&grepl("OP",.$sampleID),]
wilcox <- compare_means(value~group, df, method = "wilcox.test",
                        p.adjust.method = "fdr",paired = TRUE) %>% 
  mutate(y.position = seq(log10(max(df$value))*1.05,log10(max(df$value)), 
                          length.out = nrow(.)))
df$group <- factor(df$group,levels = c("10uM","50uM"))
plot <- ggplot(df,aes(group,value,group=sampleID))+
  geom_point(size=1)+
  geom_line(size=0.5)+theme_classic()+
  labs(x='',y='pg/OP swab',title='Human DNA')+
  theme(plot.title = element_text(hjust = 0.5))+
  annotation_logticks(sides = "l",outside = TRUE)+
  coord_cartesian(clip = "off")+
  scale_y_log10(limits=c(100,100000),
                breaks = c(100,1000,10000,100000),
                labels=trans_format("log10",math_format(10^.x)))+
  theme(axis.text = element_text(size = 7.5,color = "black"),
        axis.line = element_line(linewidth=0.5),
        axis.ticks=element_line(color="black",linewidth=0.5),
        axis.title = element_text(siz=9,color = "black"),
        plot.title = element_text(size = 9,color = "black"),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA));plot
#Fig S1E 
df <- read_excel("data.xlsx",sheet = "optimise1") %>% as.data.frame() %>% 
  filter(exp=="PMA_con"&metrics=="bac") %>% .[order(.$sampleID)&grepl("OP",.$sampleID),]
wilcox <- compare_means(value~group, df, method = "wilcox.test",
                        p.adjust.method = "fdr",paired = TRUE) %>% 
  mutate(y.position = seq(log10(max(df$value))*1.05,log10(max(df$value)), 
                          length.out = nrow(.)))
df$group <- factor(df$group,levels = c("10uM","50uM"))
plot <- ggplot(df,aes(group,value,group=sampleID))+
  geom_point(size=1)+
  geom_line(size=0.5)+theme_classic()+
  labs(x='',y='pg/OP swab',title='Bacterial DNA')+
  theme(plot.title = element_text(hjust = 0.5))+
  annotation_logticks(sides = "l",outside = TRUE)+
  coord_cartesian(clip = "off")+
  scale_y_log10(limits=c(1,10000),breaks = c(1,10,100,1000,10000),
                labels=trans_format("log10",math_format(10^.x)))+
  theme(axis.text = element_text(size = 7,color = "black"),
        axis.line = element_line(linewidth=0.5),
        axis.ticks=element_line(color="black",linewidth=0.5),
        axis.title = element_text(siz=9,color = "black"),
        plot.title = element_text(size = 9,color = "black"),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA));plot
#Fig S1F
df <- read_excel("data.xlsx",sheet = "optimise1") %>% as.data.frame() %>% 
  filter(exp=="PMA_con"&metrics=="bac_prop") %>% .[order(.$sampleID)&grepl("OP",.$sampleID),]
wilcox <- compare_means(value~group, df, method = "wilcox.test",
                        p.adjust.method = "fdr",paired = TRUE) %>% 
  mutate(y.position = seq(log10(max(df$value))*1.05,log10(max(df$value)), 
                          length.out = nrow(.)))
df$group <- factor(df$group,levels = c("10uM","50uM"))
plot <- ggplot(df,aes(group,value,group=sampleID))+
  geom_point(size=1)+
  geom_line(size=0.5)+theme_classic()+
  labs(x='',y='',title='Bacterial DNA percentage (%)')+
  theme(plot.title = element_text(hjust = 0.5))+
  annotation_logticks(sides = "l",outside = TRUE)+
  coord_cartesian(clip = "off")+
  scale_y_log10(limits=c(1,100),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels=trans_format("log10",math_format(10^.x)))+
  theme(axis.text = element_text(size = 7.5,color = "black"),
        axis.line = element_line(linewidth=0.5),
        axis.ticks=element_line(color="black",linewidth=0.5),
        axis.title = element_text(siz=9,color = "black"),
        plot.title = element_text(size = 9,color = "black"),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA));plot
#Fig S1G
df <- read_excel("data.xlsx",sheet = "optimise1") %>% as.data.frame() %>% 
  filter(exp=="sap_con"&metrics=="human") %>% .[order(.$sampleID)&grepl("BA",.$sampleID),] %>% 
  filter(group!=0)
wilcox <- compare_means(value~group, df, method = "wilcox.test",
                        p.adjust.method = "fdr",paired = TRUE)
df$group <- factor(df$group,levels = c("0.00025","0.001","0.005"))
plot <- ggplot(df,aes(group,value,group=sampleID))+
  geom_point(size=1)+
  geom_line(size=0.5)+theme_classic()+
  labs(x='',y='pg/ml BALF',title='Human DNA')+
  theme(plot.title = element_text(hjust = 0.5))+
  annotation_logticks(sides = "l",outside = TRUE)+
  coord_cartesian(clip = "off")+
  scale_y_log10(limits=c(100,100000),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels=trans_format("log10",math_format(10^.x)))+
  theme(axis.text = element_text(size = 7.5,color = "black"),
        axis.line = element_line(linewidth=0.5),
        axis.ticks=element_line(color="black",linewidth=0.5),
        axis.title = element_text(siz=9,color = "black"),
        plot.title = element_text(size = 9,color = "black"),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA));plot
#Fig S1H
df <- read_excel("data.xlsx",sheet = "optimise1") %>% as.data.frame() %>% 
  filter(exp=="sap_con"&metrics=="bac") %>% .[order(.$sampleID)&grepl("BA",.$sampleID),] %>% 
  filter(group!=0)
wilcox <- compare_means(value~group, df, method = "wilcox.test",
                        p.adjust.method = "fdr",paired = TRUE)
df$group <- factor(df$group,levels = c("0.00025","0.001","0.005"))
plot <- ggplot(df,aes(group,value,group=sampleID))+
  geom_point(size=1)+
  geom_line(size=0.5)+theme_classic()+
  labs(x='',y='pg/ml BALF',title='Bacterial DNA')+
  theme(plot.title = element_text(hjust = 0.5))+
  annotation_logticks(sides = "l",outside = TRUE)+
  coord_cartesian(clip = "off")+
  scale_y_log10(limits=c(1,10000),breaks = c(1,10,100,1000,10000),
                labels=trans_format("log10",math_format(10^.x)))+
  theme(axis.text = element_text(size = 7,color = "black"),
        axis.line = element_line(linewidth=0.5),
        axis.ticks=element_line(color="black",linewidth=0.5),
        axis.title = element_text(siz=9,color = "black"),
        plot.title = element_text(size = 9,color = "black"),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA));plot
#Fig S1I
df <- read_excel("data.xlsx",sheet = "optimise1") %>% as.data.frame() %>% 
  filter(exp=="sap_con"&metrics=="bac_prop") %>% .[order(.$sampleID)&grepl("BA",.$sampleID),] %>% 
  filter(group!=0)
wilcox <- compare_means(value~group, df, method = "wilcox.test",
                        p.adjust.method = "fdr",paired = TRUE) 
df$group <- factor(df$group,levels = c("0.00025","0.001","0.005"))
plot <- ggplot(df,aes(group,value,group=sampleID))+
  geom_point(size=1)+
  geom_line(size=0.5)+theme_classic()+
  labs(x='',y='',title='Bacterial DNA percentage (%)')+
  theme(plot.title = element_text(hjust = 0.5))+
  annotation_logticks(sides = "l",outside = TRUE)+
  coord_cartesian(clip = "off")+
  scale_y_log10(limits=c(0.1,100),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels=trans_format("log10",math_format(10^.x)))+
  theme(axis.text = element_text(size = 7.5,color = "black"),
        axis.line = element_line(linewidth=0.5),
        axis.ticks=element_line(color="black",linewidth=0.5),
        axis.title = element_text(siz=9,color = "black"),
        plot.title = element_text(size = 9,color = "black"),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA));plot
#Fig S1J
df <- read_excel("data.xlsx",sheet = "optimise1") %>% as.data.frame() %>% 
  filter(exp=="sap_con"&metrics=="human") %>% .[order(.$sampleID)&grepl("OP",.$sampleID),] %>% 
  filter(group!=0)
wilcox <- compare_means(value~group, df, method = "wilcox.test",
                        p.adjust.method = "fdr",paired = TRUE)
df$group <- factor(df$group,levels = c("0.00025","0.001","0.005"))
plot <- ggplot(df,aes(group,value,group=sampleID))+
  geom_point(size=1)+
  geom_line(size=0.5)+theme_classic()+
  labs(x='',y='pg/OP swab',title='Human DNA')+
  theme(plot.title = element_text(hjust = 0.5))+
  annotation_logticks(sides = "l",outside = TRUE)+
  coord_cartesian(clip = "off")+
  scale_y_log10(limits=c(100,300000),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels=trans_format("log10",math_format(10^.x)))+
  theme(axis.text = element_text(size = 7.5,color = "black"),
        axis.line = element_line(linewidth=0.5),
        axis.ticks=element_line(color="black",linewidth=0.5),
        axis.title = element_text(siz=9,color = "black"),
        plot.title = element_text(size = 9,color = "black"),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA));plot
#Fig S1K
df <- read_excel("data.xlsx",sheet = "optimise1") %>% as.data.frame() %>% 
  filter(exp=="sap_con"&metrics=="bac") %>% .[order(.$sampleID)&grepl("OP",.$sampleID),] %>% 
  filter(group!=0)
wilcox <- compare_means(value~group, df, method = "wilcox.test",
                        p.adjust.method = "fdr",paired = TRUE)
df$group <- factor(df$group,levels = c("0.00025","0.001","0.005"))
plot <- ggplot(df,aes(group,value,group=sampleID))+
  geom_point(size=1)+
  geom_line(size=0.5)+theme_classic()+
  labs(x='',y='pg/OP swab',title='Bacterial DNA')+
  theme(plot.title = element_text(hjust = 0.5))+
  annotation_logticks(sides = "l",outside = TRUE)+
  coord_cartesian(clip = "off")+
  scale_y_log10(limits=c(100,13000),breaks = c(100,1000,10000,13000),
                labels=trans_format("log10",math_format(10^.x)))+
  theme(axis.text = element_text(size = 7,color = "black"),
        axis.line = element_line(linewidth=0.5),
        axis.ticks=element_line(color="black",linewidth=0.5),
        axis.title = element_text(siz=9,color = "black"),
        plot.title = element_text(size = 9,color = "black"),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA));plot
##Fig S1L
df <- read_excel("data.xlsx",sheet = "optimise1") %>% as.data.frame() %>% 
  filter(exp=="sap_con"&metrics=="bac_prop") %>% .[order(.$sampleID)&grepl("OP",.$sampleID),] %>% 
  filter(group!=0)
wilcox <- compare_means(value~group, df, method = "wilcox.test",
                        p.adjust.method = "fdr",paired = TRUE) 
df$group <- factor(df$group,levels = c("0.00025","0.001","0.005"))
plot <- ggplot(df,aes(group,value,group=sampleID))+
  geom_point(size=1)+
  geom_line(size=0.5)+theme_classic()+
  labs(x='',y='',title='Bacterial DNA percentage (%)')+
  theme(plot.title = element_text(hjust = 0.5))+
  annotation_logticks(sides = "l",outside = TRUE)+
  coord_cartesian(clip = "off")+
  scale_y_log10(limits=c(1,100),
                breaks = c(1,10,100),
                labels=trans_format("log10",math_format(10^.x)))+
  theme(axis.text = element_text(size = 7.5,color = "black"),
        axis.line = element_line(linewidth=0.5),
        axis.ticks=element_line(color="black",linewidth=0.5),
        axis.title = element_text(siz=9,color = "black"),
        plot.title = element_text(size = 9,color = "black"),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA));plot
##Fig S1N
df <- read_excel("data.xlsx",sheet = "optimise1") %>% as.data.frame() %>% 
  filter(exp=="gly"&metrics=="bac") %>% .[order(.$sampleID),]
wilcox <- compare_means(value~group, df, method = "wilcox.test",
                        p.adjust.method = "fdr",paired = TRUE) %>% 
  filter(p.adj < 0.05) %>% 
  mutate(y.position = seq(log10(max(df$value))*1.05,log10(max(df$value)), 
                          length.out = nrow(.)))
df$group <- factor(df$group,levels = c("negative","positive"))
plot <- ggplot(df,aes(group,value,group=sampleID))+
  geom_point(size=1)+
  geom_line(size=0.5)+theme_classic()+
  stat_pvalue_manual(wilcox, label = "p.signif", tip.length = 0.01,
                     vjust=0.8, size=5)+
  labs(x='',y='pg/ml BALF',title='Bacterial DNA')+
  theme(plot.title = element_text(hjust = 0.5))+
  annotation_logticks(sides = "l",outside = TRUE)+
  coord_cartesian(clip = "off")+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels=trans_format("log10",math_format(10^.x)))+
  theme(axis.text = element_text(size = 7,color = "black"),
        axis.line = element_line(linewidth=0.5),
        axis.ticks=element_line(color="black",linewidth=0.5),
        axis.title = element_text(siz=9,color = "black"),
        plot.title = element_text(size = 9,color = "black"),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA));plot
##Fig S1O
df <- read_excel("data.xlsx",sheet = "optimise1") %>% as.data.frame() %>% 
  filter(exp=="gly"&metrics=="bac_prop") %>% .[order(.$sampleID),]
wilcox <- compare_means(value~group, df, method = "wilcox.test",
                        p.adjust.method = "fdr",paired = TRUE) %>% 
  mutate(y.position = seq(log10(max(df$value))*1.05,log10(max(df$value)), 
                          length.out = nrow(.)))
df$group <- factor(df$group,levels = c("negative","positive"))
plot <- ggplot(df,aes(group,value,group=sampleID))+
  geom_point(size=1)+
  geom_line(size=0.5)+theme_classic()+
  labs(x='',y='',title='Bacterial DNA percentage (%)')+
  theme(plot.title = element_text(hjust = 0.5))+
  annotation_logticks(sides = "l",outside = TRUE)+
  coord_cartesian(clip = "off")+
  scale_y_log10(limits=c(0.001,100),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels=trans_format("log10",math_format(10^.x)))+
  theme(axis.text = element_text(size = 7.5,color = "black"),
        axis.line = element_line(linewidth=0.5),
        axis.ticks=element_line(color="black",linewidth=0.5),
        axis.title = element_text(siz=9,color = "black"),
        plot.title = element_text(size = 9,color = "black"),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA));plot
# Fig. 2 ------------------------------------------------------------------
#Fig 2A
setwd("F:/Code")
plotdf <- read_excel("data.xlsx",sheet = "mass.enrich") %>% as.data.frame() %>% filter(sample_type1=="BALF") %>% 
  .[order(.$sampleNo1),]
plotdf$Methods <- factor(plotdf$Methods,levels = c("Raw","R_ase","O_pma","O_ase","F_ase","S_ase","K_qia","K_zym"))
stat <- compare_means(host.amount~Methods, plotdf, 
                      method = "wilcox.test",
                      p.adjust.method = "fdr",
                      paired = TRUE)
plotdf$host.amount[plotdf$host.amount==0] <- 0.1*min(plotdf$host.amount[plotdf$host.amount!=0])
anova <- aov(host.amount~Methods, data = plotdf)
tukey <- TukeyHSD(anova)
tukey[["Methods"]][,4] <- stat$p.adj
cld <- multcompLetters4(anova, tukey)
aa=cld$Methods$Letters %>% as.data.frame() %>% rownames_to_column("Methods")
colnames(aa)[2] <- "sign"
plot <- ggplot(plotdf,
               aes(Methods,host.amount,fill=Methods))+
  geom_boxplot(outlier.color = NA,size=0.35,color="black")+#theme_bw()+
  geom_jitter(height = 0,width = 0.1,size=0.3,color="black")+
  scale_fill_manual(values = c(Raw="black",R_ase="#696969",O_pma="#7876B1",O_ase="#FFDC91",F_ase="#E18727",
                               S_ase="#BC3C29",K_qia="#20854E",K_zym="#0072B5"))+
  labs(x='',y="Human DNA load (pg/ml BALF)")+
  coord_cartesian(clip = "off")+
  annotation_logticks(sides = "l",outside = TRUE)+
  scale_y_log10(limits=c(min(plotdf$host.amount),max(plotdf$host.amount)*2),
                labels=trans_format("log10",math_format(10^.x)),
                breaks=c(1,10,100,1000,10000,100000,1000000,
                         10000000,100000000))+
  geom_text(aes(label = sign, y = max(plotdf$host.amount)), vjust = -0.4,
            size=3.5,data = aa,show.legend = F)+
  theme(axis.ticks=element_line(color="black",linewidth=0.1),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1));plot
#Fig 2B
plotdf <- read_excel("data.xlsx",sheet = "mass.enrich") %>% as.data.frame() %>% filter(sample_type1=="OP") %>% 
  .[order(.$sampleNo1),]
plotdf$Methods <- factor(plotdf$Methods,levels = c("Raw","R_ase","O_pma","O_ase","F_ase","S_ase","K_qia","K_zym"))
stat <- compare_means(host.amount~Methods, plotdf, 
                      method = "wilcox.test",
                      p.adjust.method = "fdr",
                      paired = TRUE)
plotdf$host.amount[plotdf$host.amount==0] <- 0.1*min(plotdf$host.amount[plotdf$host.amount!=0])
anova <- aov(host.amount~Methods, data = plotdf)
tukey <- TukeyHSD(anova)
tukey[["Methods"]][,4] <- stat$p.adj
cld <- multcompLetters4(anova, tukey)
aa=cld$Methods$Letters %>% as.data.frame() %>% rownames_to_column("Methods")
colnames(aa)[2] <- "sign"
plot <- ggplot(plotdf,
               aes(Methods,host.amount,fill=Methods))+
  geom_boxplot(outlier.color = NA,size=0.35,color="black")+#theme_bw()+
  geom_jitter(height = 0,width = 0.1,size=0.3,color="black")+
  scale_fill_manual(values = c(Raw="black",R_ase="#696969",O_pma="#7876B1",O_ase="#FFDC91",F_ase="#E18727",
                               S_ase="#BC3C29",K_qia="#20854E",K_zym="#0072B5"))+
  labs(x='',y="Host DNA load (pg/swab)")+
  coord_cartesian(clip = "off")+
  annotation_logticks(sides = "l",outside = TRUE)+
  scale_y_log10(limits=c(min(plotdf$host.amount),max(plotdf$host.amount)*2),
                labels=trans_format("log10",math_format(10^.x)),
                breaks=c(1,10,100,1000,10000,100000,1000000,
                         10000000,100000000))+
  geom_text(aes(label = sign, y = max(plotdf$host.amount)), vjust = -0.4,
            size=3.5,data = aa,show.legend = F)+
  theme(axis.ticks=element_line(color="black",linewidth=0.1),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1));plot
##Fig 2C,D--Bacterial retention rate
setwd("F:/code")
df1 <- read_excel("data.xlsx",sheet = "qpcr") %>% as.data.frame()
df2 <- df1 %>% filter(methods3=="R") %>% 
  .[,c("sampleNo1","sample_type1","bac.amount")]
colnames(df2)[3] <- str_c(colnames(df2)[3],".R",sep = "")
df1 <- left_join(df1,df2) %>% mutate(bac.remain=bac.amount/bac.amount.R)
df1$bac.remain[df1$bac.remain>1] <- 1
summary <- df1 %>% dplyr::group_by(sample_type1,methods3) %>% 
  dplyr::summarise(med.retention=median(bac.remain*100),
                   Q1.bac.retention=quantile(bac.remain*100,0.25),
                   Q3.bac.retention=quantile(bac.remain*100,0.75))
df1$Methods <- "Raw"
df1$Methods[df1$methods3=="RM"] <- "R_ase"
df1$Methods[df1$methods3=="P"] <- "O_pma"
df1$Methods[df1$methods3=="W"] <- "O_ase"
df1$Methods[df1$methods3=="Fil"] <- "F_ase"
df1$Methods[df1$methods3=="S"] <- "S_ase"
df1$Methods[df1$methods3=="M"] <- "K_qia"
df1$Methods[df1$methods3=="H"] <- "K_zym"

library(multcompView)
for (i in c("BALF","OP")) {
  plotdf <- df1 %>% filter(sample_type1==i) %>% .[order(.$sampleNo1),]
  plotdf$Methods <- factor(plotdf$Methods,
                           levels = c("Raw","R_ase","O_pma","O_ase","F_ase","S_ase","K_qia","K_zym"))
  #bac.remain
  stat <- compare_means(bac.remain~Methods, plotdf, 
                        method = "wilcox.test",
                        p.adjust.method = "fdr",
                        paired = TRUE)
  anova <- aov(bac.remain~Methods, data = plotdf)
  tukey <- TukeyHSD(anova)
  tukey[["Methods"]][,4] <- stat$p.adj
  cld <- multcompLetters4(anova, tukey)
  aa=cld$Methods$Letters %>% as.data.frame() %>% rownames_to_column("Methods")
  colnames(aa)[2] <- "sign"
  plot <- ggplot(plotdf,
                 aes(Methods,bac.remain*100,fill=Methods))+
    geom_boxplot(outlier.color = NA,size=0.35,color="black")+#theme_bw()+
    geom_jitter(height = 0,width = 0.1,size=0.3,color="black")+
    scale_fill_manual(values = c(Raw="black",R_ase="#696969",
                                 O_pma="#7876B1",O_ase="#FFDC91",F_ase="#E18727",
                                 S_ase="#BC3C29",K_qia="#20854E",K_zym="#0072B5"))+
    labs(x='',y="Bacterial retention rate (%)")+
    scale_y_continuous(limits = c(0,110),
                       breaks = seq(0,100,20))+
    geom_text(aes(label = sign, y = max(plotdf$bac.remain)*100), vjust = -0.4,
              size=3.5,data = aa,show.legend = F)+
    theme(axis.ticks=element_line(color="black"),
          axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1));plot
  ggsave(paste0("result/Bac.retention.rate.",i,".pdf"),
         plot,height = 3.5,width = 4)
}  

#Fig 2E,F
df1 <- read_excel("data.xlsx",sheet = "seq.prop") %>% as.data.frame() %>% 
  .[,c("sampleID","low_quality_prop","abfv.contam","abfv.taxa","host_prop","nonabfv_prop")] %>% 
  reshape2::melt(variable.name="class") %>% 
  left_join(meta[,c("sampleID","sampleNo1","sample_type1","methods3")])
df1$Methods <- "Raw"
df1$Methods[df1$methods3=="RM"] <- "R_ase"
df1$Methods[df1$methods3=="P"] <- "O_pma"
df1$Methods[df1$methods3=="W"] <- "O_ase"
df1$Methods[df1$methods3=="Fil"] <- "F_ase"
df1$Methods[df1$methods3=="S"] <- "S_ase"
df1$Methods[df1$methods3=="M"] <- "K_qia"
df1$Methods[df1$methods3=="H"] <- "K_zym"
#Fig 2E -balf
plotdf <- df1 %>% filter(class=="abfv.taxa"&sample_type1=="BALF") %>% 
  .[order(.$sampleNo1),]
plotdf$Methods <- factor(plotdf$Methods,
                         levels = c("Raw","R_ase","O_pma","O_ase","F_ase","S_ase","K_qia","K_zym"))
stat <- compare_means(value~Methods, plotdf, 
                      method = "wilcox.test",
                      p.adjust.method = "fdr",
                      paired = TRUE)
library(multcompView)
anova <- aov(value~Methods, data = plotdf)
tukey <- TukeyHSD(anova)
tukey[["Methods"]][,4] <- stat$p.adj
cld <- multcompLetters4(anova, tukey)
aa=cld$Methods$Letters %>% as.data.frame() %>% rownames_to_column("Methods")
colnames(aa)[2] <- "sign"
plot <- ggplot(plotdf,
               aes(Methods,value*100,fill=Methods))+
  geom_boxplot(outlier.color = NA,size=0.5)+
  geom_jitter(height = 0,width = 0.1,size=0.3,color="black")+
  scale_fill_manual(values = c(Raw="black",R_ase="#696969",O_pma="#7876B1",O_ase="#FFDC91",F_ase="#E18727",
                               S_ase="#BC3C29",K_qia="#20854E",K_zym="#0072B5"))+
  labs(x='',y="Proportion of microbial reads (%)")+
  coord_cartesian(clip = "off")+
  annotation_logticks(sides = "l",outside = TRUE)+
  scale_y_log10(limits=c(min(plotdf$value)*100,max(plotdf$value)*100*2),
                labels=trans_format("log10",math_format(10^.x)),
                breaks=c(0.001,0.0001,0.001,0.01,0.1,1,10,100))+
  geom_text(aes(label = sign, y = max(plotdf$value)*100), vjust = -0.4,size=3.5,
            data = aa,show.legend = F)+
  theme(axis.ticks=element_line(color="black",linewidth=0.1),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1));plot
#Fig 2F - OP
plotdf <- df1 %>% filter(class=="abfv.taxa"&sample_type1=="OP") %>% 
  .[order(.$sampleNo1),]
plotdf$Methods <- factor(plotdf$Methods,
                         levels = c("Raw","R_ase","O_pma","O_ase","F_ase","S_ase","K_qia","K_zym"))
stat <- compare_means(value~Methods, plotdf, 
                      method = "wilcox.test",
                      p.adjust.method = "fdr",
                      paired = TRUE)
library(multcompView)
anova <- aov(value~Methods, data = plotdf)
tukey <- TukeyHSD(anova)
tukey[["Methods"]][,4] <- stat$p.adj
cld <- multcompLetters4(anova, tukey)
aa=cld$Methods$Letters %>% as.data.frame() %>% rownames_to_column("Methods")
colnames(aa)[2] <- "sign"
plot <- ggplot(plotdf,
               aes(Methods,value*100,fill=Methods))+
  geom_boxplot(outlier.color = NA,size=0.5)+
  geom_jitter(height = 0,width = 0.1,size=0.3,color="black")+
  scale_fill_manual(values = c(Raw="black",R_ase="#696969",O_pma="#7876B1",O_ase="#FFDC91",F_ase="#E18727",
                               S_ase="#BC3C29",K_qia="#20854E",K_zym="#0072B5"))+
  labs(x='',y="Proportion of microbial reads (%)")+
  coord_cartesian(clip = "off")+
  annotation_logticks(sides = "l",outside = TRUE)+
  scale_y_log10(limits=c(min(plotdf$value)*100,max(plotdf$value)*100*2),
                labels=trans_format("log10",math_format(10^.x)),
                breaks=c(0.001,0.0001,0.001,0.01,0.1,1,10,100))+
  geom_text(aes(label = sign, y = max(plotdf$value)*100), vjust = -0.4,size=3.5,
            data = aa,show.legend = F)+
  theme(axis.ticks=element_line(color="black",linewidth=0.1),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1));plot


# Fig. 3 ------------------------------------------------------------------
#Fig3A,B
setwd("F:/Code")
JSD2 <- read_excel("data.xlsx",sheet = "JSD.contam") %>% as.data.frame() %>% mutate(contam=contam*100)
#Fig. 3A
plotdata <- JSD2 %>% filter(sample_type1=="BALF")
summary <- plotdata %>% dplyr::group_by(methods3,sample_type1) %>% 
  dplyr::summarise(
    mean.rm=mean(value),
    sd.rm=sd(value),
    mean.nc=mean(contam),
    sd.nc=sd(contam)
  )
plotdata <- left_join(plotdata,summary)
plot <- ggplot(plotdata,aes(value,contam,color=methods3))+
  geom_point(alpha=0.5,size=0.2)+theme_bw()+
  geom_rect(aes(xmin=0.09,xmax=0.33,ymin=-Inf,ymax=Inf),
            fill="#e2e2e2",color=NA)+
  geom_rect(aes(xmin=0.52,xmax=0.83,ymin=-Inf,ymax=Inf),
            fill="#bababa",alpha=0.02,color=NA)+
  scale_color_manual(values = c(P="#7876B1",H="#0072B5",M="#20854E",R="black",
                                W="#FFDC91",Fil="#E18727",S="#BC3C29",RM="#696969"))+
  geom_errorbar(aes(x=mean.rm,
                    ymin=mean.nc-sd.nc,
                    ymax=mean.nc+sd.nc),
                width=0.02,size=0.1)+
  geom_errorbarh(aes(y=mean.nc,
                     xmin=mean.rm-sd.rm,
                     xmax=mean.rm+sd.rm),height=0.02,size=0.1)+
  geom_point(aes(mean.rm,mean.nc),size=1)+
  labs(x='JSD to R_ase',y='Contamination proportion (%)');plot 
ggsave("JSD2RM vs contam.b.pdf",plot,height = 3.5,width = 4.5)
#Fig. 3B
plotdata <- JSD2 %>% filter(sample_type1=="OP")
summary <- plotdata %>% dplyr::group_by(methods3,sample_type1) %>% 
  dplyr::summarise(
    mean.rm=mean(value),
    sd.rm=sd(value),
    mean.nc=mean(contam),
    sd.nc=sd(contam)
  )
plotdata <- left_join(plotdata,summary)
plot <- ggplot(plotdata,aes(value,contam,color=methods3))+
  geom_point(alpha=0.5,size=0.2)+theme_bw()+
  geom_rect(aes(xmin=0.09,xmax=0.33,ymin=-Inf,ymax=Inf),
            fill="#e2e2e2",color=NA)+
  geom_rect(aes(xmin=0.44,xmax=0.83,ymin=-Inf,ymax=Inf),
            fill="#bababa",alpha=0.02,color=NA)+
  scale_color_manual(values = c(P="#7876B1",H="#0072B5",M="#20854E",R="black",
                                W="#FFDC91",Fil="#E18727",S="#BC3C29",RM="#696969"))+
  geom_errorbar(aes(x=mean.rm,
                    ymin=mean.nc-sd.nc,
                    ymax=mean.nc+sd.nc),
                width=0.02,size=0.1)+
  geom_errorbarh(aes(y=mean.nc,
                     xmin=mean.rm-sd.rm,
                     xmax=mean.rm+sd.rm),height=0.02,size=0.1)+
  geom_point(aes(mean.rm,mean.nc),size=1)+
  scale_y_continuous(breaks = seq(0,100,25),labels = seq(0,100,25),limits = c(-8.3,100))+
  labs(x='JSD to R_ase',y='Contamination proportion');plot #,color="black"
ggsave("JSD2RM vs contam.o.pdf",plot,height = 3.5,width = 4.5)


#####Figure 6A######

pacman::p_load(philentropy,dplyr,ggplot2,phyloseq,vegan,tidyverse,ggsci,ggpubr,ggstatsplot,reshape2)
options(stringsAsFactors = F)

df_raw <- read.csv("reads.decontam.tbl2zjx.csv") %>% column_to_rownames("tax_name")
meta <- read.table("meta.tsv",sep="\t")
input_data <- df_raw[,c(meta_pair$sampleID_BA,meta_pair$sampleID_OP)]
###rarefy
input_data_rarefy=apply(input_data,2,function(x){x=rep(NA,nrow(input_data))})
rownames(input_data_rarefy)=rownames(input_data)
input_data_rarefy=as.matrix(t(rrarefy(t(input_data),10000)))
input_data_rarefy=apply(input_data_rarefy,2,function(x){x/sum(x)})
###abundance filter
input_data_rarefy <- input_data_rarefy[,c(meta_pair$sampleID_BA,meta_pair$sampleID_OP)] %>% .[rowSums(.>0.01)>=1,]

###JSD compare
ll=c(meta_pair$sampleID_BA,meta_pair$sampleID_OP)
jsd <- JSD(input_data_rarefy[,ll] %>% as.matrix() %>% prop.table(.,2) %>% t())
colnames(jsd) <- ll
rownames(jsd) <- ll
jsd[upper.tri(jsd)]=0
jsd <- jsd %>% as.data.frame() %>% rownames_to_column("id") %>% melt %>% filter(value!=0) %>% mutate(id=as.character(id),variable=as.character(variable))

jsd$ii1 <-  meta[jsd$id,"subject"]
jsd$ii2 <-  meta[jsd$variable,"subject"]
jsd$group=ifelse(jsd$ii1==jsd$ii2,"intra","inter")
jsd$group2 <-  ifelse(str_sub(jsd$id,start = 1,end = 2)==str_sub(jsd$variable,start = 1,end = 2),"Same","Diff")
jsd$group3=str_sub(jsd$id,start = 1,end = 2)
jsd$group4 <-  paste(jsd$group,jsd$group2,jsd$group3,sep = "-")
jsd$group4 <- gsub("intra-Diff-OP","intra-paired",jsd$group4)
jsd$group4 <- gsub("inter-Diff-OP","inter-unpaired",jsd$group4)
jsd$group4 <- gsub("inter-Diff-BA","inter-unpaired",jsd$group4)
jsd$group4 <- factor(jsd$group4,
                     levels = c("intra-paired","inter-Same-OP","inter-Same-BA","inter-unpaired"),
                     labels = c("paired-OP-BALF","inter-OP","inter-BALF","unpaired-OP-BALF"))

fig6a_1 <- ggboxplot(jsd,"group4","value",fill="grey",width = 0.5,legend = "none")+
  ylab("JSD distance")+xlab("")+
  stat_compare_means(comparisons = list(c("paired-OP-BALF","inter-OP"),
                                        c("paired-OP-BALF","inter-BALF"),
                                        c("paired-OP-BALF","unpaired-OP-BALF"),
                                        c("inter-OP","inter-BALF"),
                                        c("inter-OP","unpaired-OP-BALF"),
                                        c("inter-BALF","unpaired-OP-BALF")),
                     label = "p.signif")+
  theme(axis.text.x = element_text(angle = 90))
###pcoa
dist = phyloseq::distance(phyloseq(phyloseq::otu_table(input_data_rarefy,taxa_are_rows = T)),method = "jsd")


adonis2(dist~sample_type,data=meta_unpair,permutations = 1000,strata = meta_unpair$subject)

pcoa <- cmdscale(dist, k = 67 , eig = TRUE)
pcoa_eig <- pcoa$eig
pcoa_exp <- pcoa$eig/sum(pcoa$eig)
pcoa_site <- data.frame(pcoa$point)[1:2]
pcoa_site$group <- substr(rownames(pcoa_site),1,2)
pcoa_site$SampleID <-  rownames(pcoa_site)
pcoa_site$subject <- meta_unpair$subject
colnames(pcoa_site)
pcoa1 <- paste('PCoA axis1 :', round(100*pcoa_exp[1], 2), '%')
pcoa2 <- paste('PCoA axis2 :', round(100*pcoa_exp[2], 2), '%')
fig6a_2 <-   ggplot(pcoa_site, aes(X1,X2,color=group))+
  geom_point(size=4,alpha=0.9)+
  geom_line(aes(group=subject),color="grey")+
  theme_bw()+
  theme(axis.title = element_text(size = rel(1.5), 
                                  color = "black",))+
  scale_color_manual(values = c('#2B3A55','#CE7777')) +
  scale_x_continuous(limits = c(-0.4, 0.4))+
  scale_y_continuous(limits = c(-0.4, 0.4))+
  
  xlab(paste("PCo1"," (",round(pcoa_exp[1]*100,1),"%",")",sep = ""))+
  ylab(paste("PCo2"," (",round(pcoa_exp[2]*100,1),"%",")",sep = ""))


#####Figure 6B#####
input_data_rarefy_pair=apply(input_data,2,function(x){x=rep(NA,nrow(input_data))})
rownames(input_data_rarefy_pair)=rownames(input_data)
for(i in seq(1,68,2)){
  tmp=input_data[,c(i,i+1)]
  tmp=as.matrix(t(rrarefy(t(tmp),min(colSums(tmp)))))
  input_data_rarefy_pair[,colnames(tmp)]=tmp
}

input_data_rarefy_pair=apply(input_data_rarefy_pair,2,function(x){x/sum(x)})
input_data_rarefy_pair <- input_data_rarefy_pair[,meta_unpair$sampleID] %>% .[rowSums(.>0.01)>=1,]
input_data_rarefy_pair=apply(input_data_rarefy_pair,2,function(x){x/sum(x)})
input_BA <- input_data_rarefy_pair[,1:34]
input_OP<- input_data_rarefy_pair[,35:68]
result_diff_species <- NULL
for(i in unique(rownames(input_data_rarefy_pair))){
  input_tmp_BA_i <- input_BA[i,]
  input_tmp_OP_i <- input_OP[i,]
  input_tmp_i <- input_data_rarefy_pair[i,]
  fc <- mean(input_tmp_OP_i )/mean(input_tmp_BA_i )
  r <- wilcox.test(input_tmp_BA_i,input_tmp_OP_i,paired = T) 
  abundance <- mean(input_tmp_i)
  result_diff_species <- rbind(result_diff_species,c(i,r$p.value,fc,abundance))}

result_diff_species <- as.data.frame(result_diff_species)
colnames(result_diff_species) <- c("species","p","fc","abundance")
result_diff_species[,2:4] <- as.numeric(unlist(result_diff_species[,2:4]))   
result_diff_species$P_fdr <- p.adjust(result_diff_species$p,method = "fdr")
result_diff_species <- result_diff_species %>%
  mutate(color = case_when(
    fc > 1 & P_fdr < 0.05 ~ "OP",
    fc < 1 & P_fdr < 0.05 ~ "BA",
    TRUE ~ "Other"
  ))
library(ggrepel)
fig6b <- ggplot(result_diff_species,aes(log2(fc), -log10(P_fdr),color = color))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "#999999")+
  geom_point(aes(size=abundance),  alpha = 0.5) +
  scale_color_manual(values = c(OP = "#CE7777", BA = "#2B3A55", Other = "grey")) +
  theme_bw(base_size = 12)+
  geom_text_repel(data = filter(result_diff_species,  -log10(P_fdr) > 1.30103),
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                  aes(label = species),color="black",
                  size = 4) +
  xlab("Log2 (fold change)")+
  ylab("-Log10(p.adj)")
rm(tmp)
rm(r)

#####Figure 6C######
cor_all_subject=data.frame(BALF=NA,OP=NA,species=NA,subject=NA)
for(i in 1:nrow(input_data_rarefy_pair)){
  tmp=data.frame(BALF=input_data_rarefy_pair[i,meta_pair$sampleID_BA],
                 OP=input_data_rarefy_pair[i,meta_pair$sampleID_OP],
                 species=rep(rownames(input_data_rarefy_pair)[i],dim(meta_pair)[1]),
                 subject=meta_pair$subject)
  cor_all_subject=rbind(cor_all_subject,tmp)
}
cor_all_subject=cor_all_subject[-1,]

cor_all_subject$BALF=log10(cor_all_subject$BALF+0.0001)
cor_all_subject$OP=log10(cor_all_subject$OP+0.0001)
length(which(cor_all_subject$BALF==-4 & cor_all_subject$OP==-4))
cor.test(cor_all_subject$BALF,cor_all_subject$OP)
fig6c <- 
  ggplot(data=cor_all_subject[!(cor_all_subject$BALF==-4 & cor_all_subject$OP==-4),],
         aes(x=BALF,y=OP))+
  xlab("Log10 relative abundance of species in BALF")+
  ylab("Log10 relative abundance of species in OP")+
  geom_bin2d(bins = 40) +
  geom_abline(intercept = 0, slope = 1, color = "darkred")+
  scico::scale_fill_scico(palette = "vik",limits=c(0, 80), breaks=seq(0,80,by=10))+
  theme_bw()
rm(tmp)
#####Figure 6D######
meta_pair$share_species_RAinBA=apply(meta_pair,1,function(x){
  ba=x[2];op=x[3]
  ba_name <- names(which(input_data_rarefy_pair[,ba]>0.01))
  ba_abun <- input_data_rarefy_pair[ba_name,ba]/sum(input_data_rarefy_pair[ba_name,ba])
  op_name <- names(which(input_data_rarefy_pair[,op]>0.01))
  op_abun <- input_data_rarefy_pair[op_name,op]/sum(input_data_rarefy_pair[op_name,op])
  share_species=intersect( ba_name,op_name)
  return(sum(ba_abun[share_species]))
})
meta_pair$share_species_RAinOP=apply(meta_pair,1,function(x){
  ba=x[2];op=x[3]
  ba_name <- names(which(input_data_rarefy_pair[,ba]>0.01))
  ba_abun <- input_data_rarefy_pair[ba_name,ba]/sum(input_data_rarefy_pair[ba_name,ba])
  
  op_name <- names(which(input_data_rarefy_pair[,op]>0.01))
  op_abun <- input_data_rarefy_pair[op_name,op]/sum(input_data_rarefy_pair[op_name,op])
  share_species=intersect( ba_name,op_name)
  return(sum(op_abun[share_species]))
})
meta_pair[is.na(meta_pair$share_species_RAinBA),"share_species_RAinBA"] <- 1
meta_pair$op_top1_taxa=apply(meta_pair,1,function(x){rownames(input_data_rarefy_pair)[which.max(input_data_rarefy_pair[,x[3]] %>% as.numeric())]})
meta_pair$balf_top1_taxa=apply(meta_pair,1,function(x){rownames(input_data_rarefy_pair)[which.max(input_data_rarefy_pair[,x[2]] %>% as.numeric())]})
meta_pair$sameDominantTaxa=ifelse(meta_pair$op_top1_taxa==meta_pair$balf_top1_taxa,"same","not_same")

proportion_all_subject <- data.frame(share_species_RAinBA=c(meta_pair$share_species_RAinBA,meta_pair$share_species_RAinBA),
                                     share_species_RAinOP=c(meta_pair$share_species_RAinOP,meta_pair$share_species_RAinOP),
                                     pneumonia=c(meta_pair$pneumonia,meta_pair$pneumonia),
                                     sample_type=c(rep("BA",34),rep("OP",34)),
                                     bac_load=c(meta_pair$bac_load_BA,meta_pair$bac_load_OP),
                                     shape=c(rep(16,34),rep(1,34)))
fig6d <- 
  ggplot(proportion_all_subject,aes(x=share_species_RAinBA,y=share_species_RAinOP,color=pneumonia,shape=sample_type1,size=log10(bac_load)))+
  geom_point(alpha=0.5,shape=proportion_all_subject$shape)+
  scale_color_manual(values = c("#8B658B","#CD661D"),name="Pneumonia")+ 
  theme_bw()+
  xlab("Relative abundance of shared species in BALF")+
  ylab("Relative abundance of shared species in OP")

#####Figure 6E###############
fi6e_1 <- 
  ggplot(meta_pair, aes(x=pneumonia,y=100*share_species_RAinBA,color = pneumonia)) +
  geom_boxplot(outlier.size = 0.5,width=0.5) +
  geom_point(size = 0.5) +
  stat_compare_means(comparisons = list(c("no", "yes")), method = "wilcox.test",
                     size = 4,label = "p.format")+
  scale_color_manual(values = c("#2B3A55","#2B3A55")) +
  
  xlab("")+ylab("Proportion of shared species (abundance%)") +theme_test()
fig6e_2 <-
  ggplot(meta_pair, aes(x=pneumonia,y=100*share_species_RAinOP,color = pneumonia)) +
  geom_boxplot(outlier.size = 0.5,width=0.5) +
  geom_point(size = 0.5) +
  stat_compare_means(comparisons = list(c("no", "yes")), 
                     size = 4,label = "p.format")+
  scale_color_manual(values = c("#CE7777","#CE7777")) +
  xlab("")+ylab("Proportion of shared species (abundance%)") +theme_test()

#####Figure 6F##########
meta_pair$share_species_num_BA <- NULL
meta_pair$share_species_num_OP<- NULL
for(i in unique(cor_all_subject$subject)){
  p_i <- cor_all_subject[cor_all_subject$subject==i,]
  ba_name <-  p_i[p_i$BALF>log10(0.01+0.0001),"species"]
  op_name <-  p_i[p_i$OP>log10(0.01+0.0001),"species"]
  share_species=intersect( ba_name,op_name)
  ba_pro=0
  op_pro=0
  if(length(share_species)!=0){
    ba_pro  <- length(share_species)/length(ba_name)
    op_pro  <- length(share_species)/length(op_name)
  }
  meta_pair[meta_pair$subject==i,"share_species_num_BA"] <- ba_pro
  meta_pair[meta_pair$subject==i,"share_species_num_OP"] <- op_pro
}
fig6f_1 <-
  ggplot(meta_pair, aes(x=pneumonia,y=100*share_species_num_BA,color = pneumonia)) +
  geom_boxplot(outlier.size = 0.5,width=0.5) +
  geom_point(size = 0.5) +
  stat_compare_means(comparisons = list(c("no", "yes")), 
                     size = 4,label = "p.format")+
  scale_color_manual(values = c("#2B3A55","#2B3A55")) +
  xlab("")+ylab("Proportion of shared species (number%)") +theme_test()
fig6f_2 <-
  ggplot(meta_pair, aes(x=pneumonia,y=100*share_species_num_OP,color = pneumonia)) +
  geom_boxplot(outlier.size = 0.5,width=0.5) +
  geom_point(size = 0.5) +
  stat_compare_means(comparisons = list(c("no", "yes")), 
                     size = 4,label = "p.format")+
  scale_color_manual(values = c("#CE7777","#CE7777")) +
  xlab("")+ylab("Proportion of shared species (number%)") +theme_test()


#####Figure 6G#####
strain_res <- read.table("strain_res.tsv",sep="\t")
strain_res_freq <- read.table("strain_res_freq.tsv",sep="\t")
strain_res_compare <- read.table("strain_res_compare.tsv",sep="\t")
colnames(strain_res_compare)[6:7] <- c("OP>BA","OP<BA")
strain_res_note=strain_res_compare[strain_res_compare$`OP<BA`==1,]

tt1=as.data.frame(table(strain_res_compare[strain_res_compare$strain_num_BA==0,"species"]))
names(tt1)=c("species","only_OP")

tt2=as.data.frame(table(strain_res_compare[strain_res_compare$strain_num_BA==strain_res_compare$strain_num_OP & strain_res_compare$strain_num_OP==strain_res_compare$strain_overlap_num,"species"]))
names(tt2)=c("species","OP==BALF")

tt3=as.data.frame(table(strain_res_compare[strain_res_compare$strain_num_BA<strain_res_compare$strain_num_OP & strain_res_compare$strain_num_BA!=0 &strain_res_compare$strain_num_BA<=strain_res_compare$strain_overlap_num,"species"]))
names(tt3)=c("species","OP_all_include_BALF")

tt4=as.data.frame(table(strain_res_compare[strain_res_compare$strain_num_BA<=strain_res_compare$strain_num_OP & strain_res_compare$strain_num_BA!=0 & strain_res_compare$strain_overlap_num!=0 &strain_res_compare$strain_num_BA>strain_res_compare$strain_overlap_num,"species"]))
names(tt4)=c("species","OP_part_include_BALF")


tt5=as.data.frame(table(strain_res_compare[strain_res_compare$strain_num_BA>0 & strain_res_compare$strain_num_OP>0 & strain_res_compare$strain_overlap_num==0,"species"]))
names(tt5)=c("species","OP_not_include_BALF")

tt6=as.data.frame(table(strain_res_note[strain_res_note$strain_num_OP==0,"species"]))
names(tt6)=c("species","only_BALF")

tt7=as.data.frame(table(strain_res_note[strain_res_note$strain_num_OP!=0 & strain_res_note$strain_num_OP==strain_res_note$strain_overlap_num,"species"]))
names(tt7)=c("species","BALF_all_include_OP")
tt_sum=Reduce(function(x, y) merge(x, y, by = "species", all = TRUE),list(tt1,tt2,tt3,tt4,tt5,tt6,tt7))
tt_sum=tt_sum[,c(grep("species",names(tt_sum))[1],grep("species",names(tt_sum),invert = T))]
tt_sum[is.na(tt_sum)] <- 0
apply(tt_sum,1,function(x){sum(x[-1] %>% as.numeric())}) %>% sum

tmp=tt_sum$species
tt_sum=tt_sum[,-1]
rownames(tt_sum)=tmp
tt_sum$`OP>BALF`=tt_sum$only_OP+tt_sum$OP_all_include_BALF
tt_sum$`OP<BALF`=tt_sum$BALF_all_include_OP+tt_sum$only_BALF
tt_sum=select(tt_sum,-only_OP,-OP_all_include_BALF,-BALF_all_include_OP,-only_BALF)

tt_sum=tt_sum[order(tt_sum$`OP>BALF`,decreasing = T),]
tt_sum=tt_sum[,c("OP>BALF","OP==BALF","OP_part_include_BALF","OP_not_include_BALF","OP<BALF")]
fig6g_1 <- 
  corrplot::corrplot(tt_sum %>% as.matrix(),is.corr = F,method = c("circle"),
                     bg = "white",
                     shade.col = "black",
                     na.label.col="black",
                     tl.pos="lt", tl.cex=0.7, tl.col="black",
                     cl.pos = "r",cl.length = 8,
                     cl.ratio = 0.3,cl.offset=1,cl.cex = 1
  )

colnames(strain_res_compare)
strain_reads <- data.frame(species=c(strain_res_compare$species,strain_res_compare$species),
                           reads=c(strain_res_compare$reads_BALF,strain_res_compare$reads_OP),
                           sample_type=c(rep("BALF",207),rep("OP",207)))
fig6g_2 <- 
  ggplot(strain_reads, aes(x = species, y = log10(reads+1), color = sample_type)) +
  geom_boxplot(outlier.size = 0.5,width=0.5) +
  
  stat_compare_means(comparisons = list(c("BALF", "OP")), 
                     paired = T, size = 4,label = "p.format")+
  scale_color_manual(values = c(BALF="#2B3A55",OP="#CE7777"),name="") +
  xlab("")+ylab("Reads (log10)") +theme_test()+
  ylim(0,8)


#####Figure S6A#####
meta_unpair$shannon_rarefy=diversity(t(input_data_rarefy),index = "shannon")
ggplot(meta_unpair, aes(x = sample_type, y = shannon_rarefy, color = sample_type)) +
  geom_boxplot(outlier.size = 0.5,width=0.5) +
  geom_line(aes(group = subject), color = "grey80", size = 0.5) +
  geom_point(size = 0.5) +
  stat_compare_means(comparisons = list(c("BALF", "OP")), 
                     paired = T, size = 4,label = "p.format")+
  scale_color_manual(values = c(BALF="#2B3A55",OP="#CE7777"),name="") +
  xlab("")+ylab("Shannon index") +theme_test()

#####Figure S6B#####
cluster_abundance <- input_data_rarefy_pair %>% as.data.frame() %>% rownames_to_column("taxname")
cluster_abundance$taxname[!cluster_abundance$taxname %in% unique(c(meta_pair$op_top1_taxa,meta_pair$balf_top1_taxa))] <- "Others"#16:nrow(cluster_abundance)
cluster_abundance <- cluster_abundance %>% 
  group_by(taxname) %>%
  summarise_each(funs = sum) %>%
  as.data.frame()

cluster_abundance <- cluster_abundance[!cluster_abundance$taxname=="Others",]
cluster_abundance[dim(cluster_abundance)[1]+1,] <- c("Others",as.numeric(1-colSums(cluster_abundance[,2:69])))
cluster_abundance[dim(cluster_abundance)[1],2:69] <- as.numeric(cluster_abundance[dim(cluster_abundance)[1],2:69])

cluster_abundance$taxname <- factor(cluster_abundance$taxname, levels = c("Others",unique(c(meta_pair$op_top1_taxa,meta_pair$balf_top1_taxa))))
cluster_abundance <- cluster_abundance %>% melt(id="taxname",variable.name="sampleID") %>% right_join(meta_unpair,.)
cluster_abundance$subject <- factor(cluster_abundance$subject,levels = c(meta_pair[order(meta_pair$share_species_RAinBA,decreasing = T),"subject"]))
library(pals)
ggplot(cluster_abundance) +
  geom_bar(aes(x = subject, y = as.numeric(value), fill = taxname),width=1, stat = "identity") +
  scale_fill_manual(values = c("grey90",pals::alphabet2()[c(1:4,6,8,10:12,14:17)]%>% as.character(),
                               pals::kelly() %>% as.character() %>% .[c(4:8,11:13,16:22)]),
                    name = "Taxonomy") +
  theme_bw() +theme(axis.text.x = element_text(size = 12), axis.ticks.x = element_blank(), 
                    plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  #theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  scale_y_continuous(expand = c(0, 0.01), name = "Relative abundance") +
  xlab("")+facet_grid(sample_type~pneumonia,scales = "free_x")

#####Figure S6C#####
input_data_rarefy_pair_pneumonia <- input_data_rarefy_pair[,meta_unpair[meta_unpair$pneumonia=="yes","sampleID"]]
input_data_rarefy_pair_non_pneumonia <- input_data_rarefy_pair[,meta_unpair[meta_unpair$pneumonia=="no","sampleID"]]
meta_pair_pneumonia <- meta_pair[meta_pair$pneumonia=="yes",]
meta_pair_non_pneumonia <- meta_pair[meta_pair$pneumonia=="no",]
cor_all_subject_pneumonia=data.frame(BALF=NA,OP=NA,species=NA,subject=NA)
for(i in 1:nrow(input_data_rarefy_pair_pneumonia)){
  tmp=data.frame(BALF=input_data_rarefy_pair_pneumonia[i,meta_pair_pneumonia $sampleID_BA],
                 OP=input_data_rarefy_pair_pneumonia[i,meta_pair_pneumonia $sampleID_OP],
                 species=rep(rownames(input_data_rarefy_pair_pneumonia)[i],dim(meta_pair_pneumonia )[1]),
                 subject=meta_pair_pneumonia $subject)
  cor_all_subject_pneumonia=rbind(cor_all_subject_pneumonia,tmp)
}
cor_all_subject_pneumonia=cor_all_subject_pneumonia[-1,]

cor_all_subject_pneumonia$BALF=log10(cor_all_subject_pneumonia$BALF+0.0001)
cor_all_subject_pneumonia$OP=log10(cor_all_subject_pneumonia$OP+0.0001)
length(which(cor_all_subject_pneumonia$BALF==-4 & cor_all_subject_pneumonia$OP==-4))
cor.test(cor_all_subject_pneumonia$BALF,cor_all_subject_pneumonia$OP)

ggplot(data=cor_all_subject_pneumonia[!(cor_all_subject_pneumonia$BALF==-4 & cor_all_subject_pneumonia$OP==-4),],
       aes(x=BALF,y=OP))+
  xlab("Log10 relative abundance of species in BALF")+
  ylab("Log10 relative abundance of species in OP")+
  geom_bin2d(bins = 40) +
  geom_abline(intercept = 0, slope = 1, color = "darkred")+
  scico::scale_fill_scico(palette = "vik",limits=c(0, 80), breaks=seq(0,80,by=10))+
  theme_bw()

cor_all_subject_non_pneumonia=data.frame(BALF=NA,OP=NA,species=NA,subject=NA)
for(i in 1:nrow(input_data_rarefy_pair_non_pneumonia)){
  tmp=data.frame(BALF=input_data_rarefy_pair_non_pneumonia[i,meta_pair_non_pneumonia $sampleID_BA],
                 OP=input_data_rarefy_pair_non_pneumonia[i,meta_pair_non_pneumonia $sampleID_OP],
                 species=rep(rownames(input_data_rarefy_pair_non_pneumonia)[i],dim(meta_pair_non_pneumonia )[1]),
                 subject=meta_pair_non_pneumonia $subject)
  cor_all_subject_non_pneumonia=rbind(cor_all_subject_non_pneumonia,tmp)
}
cor_all_subject_non_pneumonia=cor_all_subject_non_pneumonia[-1,]

cor_all_subject_non_pneumonia$BALF=log10(cor_all_subject_non_pneumonia$BALF+0.0001)
cor_all_subject_non_pneumonia$OP=log10(cor_all_subject_non_pneumonia$OP+0.0001)
length(which(cor_all_subject_non_pneumonia$BALF==-4 & cor_all_subject_non_pneumonia$OP==-4))
cor.test(cor_all_subject_non_pneumonia$BALF,cor_all_subject_non_pneumonia$OP)

ggplot(data=cor_all_subject_non_pneumonia[!(cor_all_subject_non_pneumonia$BALF==-4 & cor_all_subject_non_pneumonia$OP==-4),],
       aes(x=BALF,y=OP))+
  xlab("Log10 relative abundance of species in BALF")+
  ylab("Log10 relative abundance of species in OP")+
  geom_bin2d(bins = 40) +
  geom_abline(intercept = 0, slope = 1, color = "darkred")+
  scico::scale_fill_scico(palette = "vik",limits=c(0, 80), breaks=seq(0,80,by=10))+
  theme_bw()
#####Figure S6D#####
ggplot(meta_unpair[meta_unpair$pneumonia=="yes",], aes(x=sample_type,y=shannon_rarefy,color = pneumonia)) +
  geom_boxplot(outlier.size = 0.5,width=0.5) +
  geom_point(size = 0.5) +
  geom_line(aes(group=subject),color='grey',lwd=0.35)+
  stat_compare_means(comparisons = list(c("BALF", "OP")), method = "wilcox.test",
                     size = 4,label = "p.format")+
  scale_color_manual(values = c("#CD661D","#CD661D")) +
  
  xlab("")+ylab("Shanno index") +theme_test()

ggplot(meta_unpair[meta_unpair$pneumonia=="no",], aes(x=sample_type,y=shannon_rarefy,color = pneumonia)) +
  geom_boxplot(outlier.size = 0.5,width=0.5) +
  geom_point(size = 0.5) +
  geom_line(aes(group=subject),color='grey',lwd=0.35)+
  stat_compare_means(comparisons = list(c("BALF", "OP")), method = "wilcox.test",
                     size = 4,label = "p.format")+
  scale_color_manual(values = c("#8B658B","#8B658B")) +
  
  xlab("")+ylab("Shanno index") +theme_test()
#####Figure S6E#######
jsd_pair <- jsd[jsd$group4=="paired-OP-BALF",]
colnames(jsd_pair)[4] <- "subject"
jsd_pair <- merge(jsd_pair,meta_unpair,by="subject")
ggplot(jsd_pair, aes(x=pneumonia,y=value,color = pneumonia)) +
  geom_boxplot(outlier.size = 0.5,width=0.5) +
  geom_point(size = 0.5) +
  stat_compare_means(comparisons = list(c("yes", "no")), method = "wilcox.test",
                     size = 4,label = "p.format")+
  scale_color_manual(values = c("#8B658B","#CD661D")) +
  
  xlab("")+ylab("JSD") +theme_test()

#####Figure S7A######

input_data_before <- df_raw[,c(meta_pair$sampleID_BA_before,meta_pair$sampleID_OP_before)]
###rarefy
input_data_before_rarefy=apply(input_data_before,2,function(x){x=rep(NA,nrow(input_data_before))})
rownames(input_data_before_rarefy)=rownames(input_data_before)
input_data_before_rarefy=as.matrix(t(rrarefy(t(input_data_before),10000)))
input_data_before_rarefy=apply(input_data_before_rarefy,2,function(x){x/sum(x)})
###abundance filter
input_data_before_rarefy <- input_data_before_rarefy[,c(meta_pair$sampleID_BA_before,meta_pair$sampleID_OP_before)] %>% .[rowSums(.>0.01)>=1,]

###JSD compare
ll=c(meta_pair$sampleID_BA_before,meta_pair$sampleID_OP_before)
jsd <- JSD(input_data_before_rarefy[,ll] %>% as.matrix() %>% prop.table(.,2) %>% t())
colnames(jsd) <- ll
rownames(jsd) <- ll
jsd[upper.tri(jsd)]=0
jsd <- jsd %>% as.data.frame() %>% rownames_to_column("id") %>% melt %>% filter(value!=0) %>% mutate(id=as.character(id),variable=as.character(variable))

jsd$ii1 <-  meta[jsd$id,"subject"]
jsd$ii2 <-  meta[jsd$variable,"subject"]
jsd$group=ifelse(jsd$ii1==jsd$ii2,"intra","inter")
jsd$group2 <-  ifelse(str_sub(jsd$id,start = 1,end = 2)==str_sub(jsd$variable,start = 1,end = 2),"Same","Diff")
jsd$group3=str_sub(jsd$id,start = 1,end = 2)
jsd$group4 <-  paste(jsd$group,jsd$group2,jsd$group3,sep = "-")
jsd$group4 <- gsub("intra-Diff-OP","intra-paired",jsd$group4)
jsd$group4 <- gsub("inter-Diff-OP","inter-unpaired",jsd$group4)
jsd$group4 <- gsub("inter-Diff-BA","inter-unpaired",jsd$group4)
jsd$group4 <- factor(jsd$group4,
                     levels = c("intra-paired","inter-Same-OP","inter-Same-BA","inter-unpaired"),
                     labels = c("paired-OP-BALF","inter-OP","inter-BALF","unpaired-OP-BALF"))

ggboxplot(jsd,"group4","value",fill="grey",width = 0.5,legend = "none")+
  ylab("JSD distance")+xlab("")+
  stat_compare_means(comparisons = list(c("paired-OP-BALF","inter-OP"),
                                        c("paired-OP-BALF","inter-BALF"),
                                        c("paired-OP-BALF","unpaired-OP-BALF"),
                                        c("inter-OP","inter-BALF"),
                                        c("inter-OP","unpaired-OP-BALF"),
                                        c("inter-BALF","unpaired-OP-BALF")),
                     label = "p.signif")+
  theme(axis.text.x = element_text(angle = 90))
###pcoa
dist = phyloseq::distance(phyloseq(phyloseq::otu_table(input_data_before_rarefy,taxa_are_rows = T)),method = "jsd")

adonis2(dist~sample_type,data=meta_unpair,permutations = 1000,strata = meta_unpair$subject)

pcoa <- cmdscale(dist, k = 67 , eig = TRUE)
pcoa_eig <- pcoa$eig
pcoa_exp <- pcoa$eig/sum(pcoa$eig)
pcoa_site <- data.frame(pcoa$point)[1:2]
pcoa_site$group <- substr(rownames(pcoa_site),1,2)
pcoa_site$SampleID <-  rownames(pcoa_site)
pcoa_site$subject <- meta_unpair$subject
colnames(pcoa_site)
pcoa1 <- paste('PCoA axis1 :', round(100*pcoa_exp[1], 2), '%')
pcoa2 <- paste('PCoA axis2 :', round(100*pcoa_exp[2], 2), '%')
ggplot(pcoa_site, aes(X1,X2,color=group))+
  geom_point(size=4,alpha=0.9)+
  geom_line(aes(group=subject),color="grey")+
  theme_bw()+
  theme(axis.title = element_text(size = rel(1.5), 
                                  color = "black",))+
  scale_color_manual(values = c('#2B3A55','#CE7777')) +
  scale_x_continuous(limits = c(-0.4, 0.4))+
  scale_y_continuous(limits = c(-0.4, 0.4))+
  
  xlab(paste("PCo1"," (",round(pcoa_exp[1]*100,1),"%",")",sep = ""))+
  ylab(paste("PCo2"," (",round(pcoa_exp[2]*100,1),"%",")",sep = ""))


#####Figure S7B#####
input_data_before_rarefy_pair=apply(input_data_before,2,function(x){x=rep(NA,nrow(input_data_before))})
rownames(input_data_before_rarefy_pair)=rownames(input_data_before)
for(i in seq(1,68,2)){
  tmp=input_data_before[,c(i,i+1)]
  tmp=as.matrix(t(rrarefy(t(tmp),min(colSums(tmp)))))
  input_data_before_rarefy_pair[,colnames(tmp)]=tmp
}

input_data_before_rarefy_pair=apply(input_data_before_rarefy_pair,2,function(x){x/sum(x)})
input_data_before_rarefy_pair <- input_data_before_rarefy_pair[,meta_unpair$sampleID_before] %>% .[rowSums(.>0.01)>=1,]
input_data_before_rarefy_pair=apply(input_data_before_rarefy_pair,2,function(x){x/sum(x)})
input_BA <- input_data_before_rarefy_pair[,1:34]
input_OP<- input_data_before_rarefy_pair[,35:68]
result_diff_species <- NULL
for(i in unique(rownames(input_data_before_rarefy_pair))){
  input_tmp_BA_i <- input_BA[i,]
  input_tmp_OP_i <- input_OP[i,]
  input_tmp_i <- input_data_before_rarefy_pair[i,]
  fc <- mean(input_tmp_OP_i )/mean(input_tmp_BA_i )
  r <- wilcox.test(input_tmp_BA_i,input_tmp_OP_i,paired = T) 
  abundance <- mean(input_tmp_i)
  result_diff_species <- rbind(result_diff_species,c(i,r$p.value,fc,abundance))}

result_diff_species <- as.data.frame(result_diff_species)
colnames(result_diff_species) <- c("species","p","fc","abundance")
result_diff_species[,2:4] <- as.numeric(unlist(result_diff_species[,2:4]))   
result_diff_species$P_fdr <- p.adjust(result_diff_species$p,method = "fdr")
result_diff_species <- result_diff_species %>%
  mutate(color = case_when(
    fc > 1 & P_fdr < 0.05 ~ "OP",
    fc < 1 & P_fdr < 0.05 ~ "BA",
    TRUE ~ "Other"
  ))
library(ggrepel)
ggplot(result_diff_species,aes(log2(fc), -log10(P_fdr),color = color))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "#999999")+
  geom_point(aes(size=abundance),  alpha = 0.5) +
  scale_color_manual(values = c(OP = "#CE7777", BA = "#2B3A55", Other = "grey")) +
  theme_bw(base_size = 12)+
  geom_text_repel(data = filter(result_diff_species,  -log10(P_fdr) > 1.30103),
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                  aes(label = species),color="black",
                  size = 4) +
  xlab("Log2 (fold change)")+
  ylab("-Log10(p.adj)")
rm(tmp)
rm(r)

#####Figure S7C######
cor_all_subject=data.frame(BALF=NA,OP=NA,species=NA,subject=NA)
for(i in 1:nrow(input_data_before_rarefy_pair)){
  tmp=data.frame(BALF=input_data_before_rarefy_pair[i,meta_pair$sampleID_BA_before],
                 OP=input_data_before_rarefy_pair[i,meta_pair$sampleID_OP_before],
                 species=rep(rownames(input_data_before_rarefy_pair)[i],dim(meta_pair)[1]),
                 subject=meta_pair$subject)
  cor_all_subject=rbind(cor_all_subject,tmp)
}
cor_all_subject=cor_all_subject[-1,]

cor_all_subject$BALF=log10(cor_all_subject$BALF+0.0001)
cor_all_subject$OP=log10(cor_all_subject$OP+0.0001)
length(which(cor_all_subject$BALF==-4 & cor_all_subject$OP==-4))
cor.test(cor_all_subject$BALF,cor_all_subject$OP)

ggplot(data=cor_all_subject[!(cor_all_subject$BALF==-4 & cor_all_subject$OP==-4),],
       aes(x=BALF,y=OP))+
  xlab("Log10 relative abundance of species in BALF")+
  ylab("Log10 relative abundance of species in OP")+
  geom_bin2d(bins = 40) +
  geom_abline(intercept = 0, slope = 1, color = "darkred")+
  scico::scale_fill_scico(palette = "vik",limits=c(0, 80), breaks=seq(0,80,by=10))+
  theme_bw()
rm(tmp)
#####Figure S7D######
meta_pair$share_species_RAinBA=apply(meta_pair,1,function(x){
  ba=x[2];op=x[3]
  ba_name <- names(which(input_data_before_rarefy_pair[,ba]>0.01))
  ba_abun <- input_data_before_rarefy_pair[ba_name,ba]/sum(input_data_before_rarefy_pair[ba_name,ba])
  op_name <- names(which(input_data_before_rarefy_pair[,op]>0.01))
  op_abun <- input_data_before_rarefy_pair[op_name,op]/sum(input_data_before_rarefy_pair[op_name,op])
  share_species=intersect( ba_name,op_name)
  return(sum(ba_abun[share_species]))
})
meta_pair$share_species_RAinOP=apply(meta_pair,1,function(x){
  ba=x[2];op=x[3]
  ba_name <- names(which(input_data_before_rarefy_pair[,ba]>0.01))
  ba_abun <- input_data_before_rarefy_pair[ba_name,ba]/sum(input_data_before_rarefy_pair[ba_name,ba])
  
  op_name <- names(which(input_data_before_rarefy_pair[,op]>0.01))
  op_abun <- input_data_before_rarefy_pair[op_name,op]/sum(input_data_before_rarefy_pair[op_name,op])
  share_species=intersect( ba_name,op_name)
  return(sum(op_abun[share_species]))
})
meta_pair[is.na(meta_pair$share_species_RAinBA),"share_species_RAinBA"] <- 1
meta_pair$op_top1_taxa=apply(meta_pair,1,function(x){rownames(input_data_before_rarefy_pair)[which.max(input_data_before_rarefy_pair[,x[3]] %>% as.numeric())]})
meta_pair$balf_top1_taxa=apply(meta_pair,1,function(x){rownames(input_data_before_rarefy_pair)[which.max(input_data_before_rarefy_pair[,x[2]] %>% as.numeric())]})
meta_pair$sameDominantTaxa=ifelse(meta_pair$op_top1_taxa==meta_pair$balf_top1_taxa,"same","not_same")

proportion_all_subject <- data.frame(share_species_RAinBA=c(meta_pair$share_species_RAinBA,meta_pair$share_species_RAinBA),
                                     share_species_RAinOP=c(meta_pair$share_species_RAinOP,meta_pair$share_species_RAinOP),
                                     pneumonia=c(meta_pair$pneumonia,meta_pair$pneumonia),
                                     sample_type=c(rep("BA",34),rep("OP",34)),
                                     bac_load=c(meta_pair$bac_load_BA,meta_pair$bac_load_OP),
                                     shape=c(rep(16,34),rep(1,34)))

ggplot(proportion_all_subject,aes(x=share_species_RAinBA,y=share_species_RAinOP,color=pneumonia,shape=sample_type1,size=log10(bac_load)))+
  geom_point(alpha=0.5,shape=proportion_all_subject$shape)+
  scale_color_manual(values = c("#8B658B","#CD661D"),name="Pneumonia")+ 
  theme_bw()+
  xlab("Relative abundance of shared species in BALF")+
  ylab("Relative abundance of shared species in OP")

#####Figure S7E###############
fi6e_1 <- 
  ggplot(meta_pair, aes(x=pneumonia,y=100*share_species_RAinBA,color = pneumonia)) +
  geom_boxplot(outlier.size = 0.5,width=0.5) +
  geom_point(size = 0.5) +
  stat_compare_means(comparisons = list(c("no", "yes")), method = "wilcox.test",
                     size = 4,label = "p.format")+
  scale_color_manual(values = c("#2B3A55","#2B3A55")) +
  
  xlab("")+ylab("Proportion of shared species (abundance%)") +theme_test()

ggplot(meta_pair, aes(x=pneumonia,y=100*share_species_RAinOP,color = pneumonia)) +
  geom_boxplot(outlier.size = 0.5,width=0.5) +
  geom_point(size = 0.5) +
  stat_compare_means(comparisons = list(c("no", "yes")), 
                     size = 4,label = "p.format")+
  scale_color_manual(values = c("#CE7777","#CE7777")) +
  xlab("")+ylab("Proportion of shared species (abundance%)") +theme_test()

#####Figure S7F##########
meta_pair$share_species_num_BA <- NULL
meta_pair$share_species_num_OP<- NULL
for(i in unique(cor_all_subject$subject)){
  p_i <- cor_all_subject[cor_all_subject$subject==i,]
  ba_name <-  p_i[p_i$BALF>log10(0.01+0.0001),"species"]
  op_name <-  p_i[p_i$OP>log10(0.01+0.0001),"species"]
  share_species=intersect( ba_name,op_name)
  ba_pro=0
  op_pro=0
  if(length(share_species)!=0){
    ba_pro  <- length(share_species)/length(ba_name)
    op_pro  <- length(share_species)/length(op_name)
  }
  meta_pair[meta_pair$subject==i,"share_species_num_BA"] <- ba_pro
  meta_pair[meta_pair$subject==i,"share_species_num_OP"] <- op_pro
}

ggplot(meta_pair, aes(x=pneumonia,y=100*share_species_num_BA,color = pneumonia)) +
  geom_boxplot(outlier.size = 0.5,width=0.5) +
  geom_point(size = 0.5) +
  stat_compare_means(comparisons = list(c("no", "yes")), 
                     size = 4,label = "p.format")+
  scale_color_manual(values = c("#2B3A55","#2B3A55")) +
  xlab("")+ylab("Proportion of shared species (number%)") +theme_test()

ggplot(meta_pair, aes(x=pneumonia,y=100*share_species_num_OP,color = pneumonia)) +
  geom_boxplot(outlier.size = 0.5,width=0.5) +
  geom_point(size = 0.5) +
  stat_compare_means(comparisons = list(c("no", "yes")), 
                     size = 4,label = "p.format")+
  scale_color_manual(values = c("#CE7777","#CE7777")) +
  xlab("")+ylab("Proportion of shared species (number%)") +theme_test()


#####Figure S7G#####
colnames(strain_res_compare)[6:7] <- c("OP>BA","OP<BA")
strain_res_note=strain_res_compare[strain_res_compare$`OP<BA`==1,]

tt1=as.data.frame(table(strain_res_compare[strain_res_compare$strain_num_BA==0,"species"]))
names(tt1)=c("species","only_OP")

tt2=as.data.frame(table(strain_res_compare[strain_res_compare$strain_num_BA==strain_res_compare$strain_num_OP & strain_res_compare$strain_num_OP==strain_res_compare$strain_overlap_num,"species"]))
names(tt2)=c("species","OP==BALF")

tt3=as.data.frame(table(strain_res_compare[strain_res_compare$strain_num_BA<strain_res_compare$strain_num_OP & strain_res_compare$strain_num_BA!=0 &strain_res_compare$strain_num_BA<=strain_res_compare$strain_overlap_num,"species"]))
names(tt3)=c("species","OP_all_include_BALF")

tt4=as.data.frame(table(strain_res_compare[strain_res_compare$strain_num_BA<=strain_res_compare$strain_num_OP & strain_res_compare$strain_num_BA!=0 & strain_res_compare$strain_overlap_num!=0 &strain_res_compare$strain_num_BA>strain_res_compare$strain_overlap_num,"species"]))
names(tt4)=c("species","OP_part_include_BALF")


tt5=as.data.frame(table(strain_res_compare[strain_res_compare$strain_num_BA>0 & strain_res_compare$strain_num_OP>0 & strain_res_compare$strain_overlap_num==0,"species"]))
names(tt5)=c("species","OP_not_include_BALF")

tt6=as.data.frame(table(strain_res_note[strain_res_note$strain_num_OP==0,"species"]))
names(tt6)=c("species","only_BALF")

tt7=as.data.frame(table(strain_res_note[strain_res_note$strain_num_OP!=0 & strain_res_note$strain_num_OP==strain_res_note$strain_overlap_num,"species"]))
names(tt7)=c("species","BALF_all_include_OP")
tt_sum=Reduce(function(x, y) merge(x, y, by = "species", all = TRUE),list(tt1,tt2,tt3,tt4,tt5,tt6,tt7))
tt_sum=tt_sum[,c(grep("species",names(tt_sum))[1],grep("species",names(tt_sum),invert = T))]
tt_sum[is.na(tt_sum)] <- 0
apply(tt_sum,1,function(x){sum(x[-1] %>% as.numeric())}) %>% sum

tmp=tt_sum$species
tt_sum=tt_sum[,-1]
rownames(tt_sum)=tmp
tt_sum$`OP>BALF`=tt_sum$only_OP+tt_sum$OP_all_include_BALF
tt_sum$`OP<BALF`=tt_sum$BALF_all_include_OP+tt_sum$only_BALF
tt_sum=select(tt_sum,-only_OP,-OP_all_include_BALF,-BALF_all_include_OP,-only_BALF)

tt_sum=tt_sum[order(tt_sum$`OP>BALF`,decreasing = T),]
tt_sum=tt_sum[,c("OP>BALF","OP==BALF","OP_part_include_BALF","OP_not_include_BALF","OP<BALF")]

corrplot::corrplot(tt_sum %>% as.matrix(),is.corr = F,method = c("circle"),
                   bg = "white",
                   shade.col = "black",
                   na.label.col="black",
                   tl.pos="lt", tl.cex=0.7, tl.col="black",
                   cl.pos = "r",cl.length = 8,
                   cl.ratio = 0.3,cl.offset=1,cl.cex = 1
)

colnames(strain_res_compare)
strain_reads <- data.frame(species=c(strain_res_compare$species,strain_res_compare$species),
                           reads=c(strain_res_compare$reads_BALF,strain_res_compare$reads_OP),
                           sample_type=c(rep("BALF",207),rep("OP",207)))

ggplot(strain_reads, aes(x = species, y = log10(reads+1), color = sample_type)) +
  geom_boxplot(outlier.size = 0.5,width=0.5) +
  
  stat_compare_means(comparisons = list(c("BALF", "OP")), 
                     paired = T, size = 4,label = "p.format")+
  scale_color_manual(values = c(BALF="#2B3A55",OP="#CE7777"),name="") +
  xlab("")+ylab("Reads (log10)") +theme_test()+
  ylim(0,8)
