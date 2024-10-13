<<<<<<< HEAD
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
# Fig. 2A-B ------------------------------------------------------------------
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

# Fig 2C-D--Bacterial retention rate --------------------------------------
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

# Fig 2E-F ----------------------------------------------------------------
setwd("F:/code")
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

# Fig. 2G-H, K ---------------------------------------------------------------
setwd("F:/code")
de.tbl1 <- read_excel("data.xlsx",sheet = "de.species") %>% 
  as.data.frame() %>% column_to_rownames("tax_name") #taxonomy table after removing contaminants
qc <- read_excel("data.xlsx",sheet = "qc") %>% 
  as.data.frame() %>% .[,c("sampleID2","total_reads","sequins","ABFV")] %>% 
  mutate(bac.prop=ABFV/(total_reads-sequins),
         total=total_reads-sequins)
meta <- read_excel("data.xlsx",sheet = "meta") %>% 
  as.data.frame() %>% 
  .[,c("sampleID2","sample_type1","methods3","sampleNo1")] %>% 
  .[!grepl("_2",.$sampleID),] %>% unique() %>% filter(methods3!="NC")
qc2 <- left_join(qc,meta) %>% filter(!is.na(methods3))
sum <- qc2 %>% dplyr::group_by(sample_type1,sampleNo1) %>% 
  dplyr::summarise(rare=min(total))
qc2 <- left_join(qc2,sum)
rare <- data.frame(tax_name=rownames(de.tbl1))
rich <- data.frame()
bac.raw <- qc2 %>% filter(methods3=="R")
for (i in qc2$sampleID2) {
  qc1.i <- filter(qc2,sampleID2==i)
  df.i <- data.frame(tax_name=rownames(de.tbl1),
                     sample=round(de.tbl1[,i])) %>% 
    column_to_rownames("tax_name")
  if (nrow(df.i)>0) {
    set.seed(100)
    rare.i <- rrarefy(df.i,qc1.i$rare) %>% t() %>% 
      as.data.frame() %>% 
      rownames_to_column("tax_name")
    
    rich.i <- data.frame(sampleID=i,
                         richniss=length(rare.i$sample[rare.i$sample>0]),
                         rare.theory=qc1.i$rare,
                         rare.real=sum(rare.i$sample))
    rich <- rbind(rich,rich.i)
    
    colnames(rare.i)[2] <- i
    rare <- left_join(rare,rare.i) 
  }
}
rich1 <- rich %>% mutate(diff=rare.theory-rare.real)
rare[is.na(rare)] <- 0
colnames(rich1)[1] <- "sampleID2"
df <- left_join(rich1,qc2)
df <- left_join(rich1,qc) %>% left_join(meta)
df.R <- df %>% filter(methods3=="R") %>% 
  .[,c("sample_type1","sampleNo1","richniss")]
colnames(df.R)[3] <- "rich.R"
df1 <- left_join(df,df.R) %>% mutate(fold=richniss/rich.R)
#BALF
plotdf <- df1 %>% filter(sample_type1=="BALF") %>% 
  .[order(.$sampleNo1),]
plotdf$Methods <- "Raw"
plotdf$Methods[plotdf$methods3=="RM"] <- "R_ase"
plotdf$Methods[plotdf$methods3=="P"] <- "O_pma"
plotdf$Methods[plotdf$methods3=="W"] <- "O_ase"
plotdf$Methods[plotdf$methods3=="Fil"] <- "F_ase"
plotdf$Methods[plotdf$methods3=="S"] <- "S_ase"
plotdf$Methods[plotdf$methods3=="M"] <- "K_qia"
plotdf$Methods[plotdf$methods3=="H"] <- "K_zym"
plotdf$Methods <- factor(plotdf$Methods,
                         levels = c("Raw","R_ase","O_pma","O_ase","F_ase",
                                    "S_ase","K_qia","K_zym"))
stat <- compare_means(richniss~Methods, plotdf, 
                      method = "wilcox.test",
                      p.adjust.method = "fdr",
                      paired = TRUE)
library(multcompView)
anova <- aov(richniss~Methods, data = plotdf)
tukey <- TukeyHSD(anova)
tukey[["Methods"]][,4] <- stat$p.adj
cld <- multcompLetters4(anova, tukey)
aa=cld$Methods$Letters %>% as.data.frame() %>% rownames_to_column("Methods")
colnames(aa)[2] <- "sign"
# fold change
plotdf <- df1 %>% filter(sample_type1=="BALF") %>% 
  .[order(.$sampleNo1),]
plotdf$Methods <- "Raw"
plotdf$Methods[plotdf$methods3=="RM"] <- "R_ase"
plotdf$Methods[plotdf$methods3=="P"] <- "O_pma"
plotdf$Methods[plotdf$methods3=="W"] <- "O_ase"
plotdf$Methods[plotdf$methods3=="Fil"] <- "F_ase"
plotdf$Methods[plotdf$methods3=="S"] <- "S_ase"
plotdf$Methods[plotdf$methods3=="M"] <- "K_qia"
plotdf$Methods[plotdf$methods3=="H"] <- "K_zym"
plotdf$Methods <- factor(plotdf$Methods,
                         levels = c("Raw","R_ase","O_pma","O_ase","F_ase",
                                    "S_ase","K_qia","K_zym"))
stat <- compare_means(fold~Methods, plotdf, 
                      method = "wilcox.test",
                      p.adjust.method = "fdr",
                      paired = TRUE)
library(multcompView)
anova <- aov(richniss~Methods, data = plotdf)
tukey <- TukeyHSD(anova)
tukey[["Methods"]][,4] <- stat$p.adj
cld <- multcompLetters4(anova, tukey)
aa=cld$Methods$Letters %>% as.data.frame() %>% rownames_to_column("Methods")
colnames(aa)[2] <- "sign"
plotdf$fold[plotdf$fold>10] <- 10
plot <- ggplot(plotdf,aes(Methods,fold,fill=Methods))+
  geom_boxplot(outlier.color = NA,size=0.5)+
  geom_jitter(width = 0.1,size=0.3,color="black")+
  scale_fill_manual(values = c(Raw="black",R_ase="#696969",
                               O_pma="#7876B1",O_ase="#FFDC91",F_ase="#E18727",
                               S_ase="#BC3C29",K_qia="#20854E",K_zym="#0072B5"))+
  labs(x='',y="Fold change of species detected")+
  coord_cartesian(clip = "off")+
  annotation_logticks(sides = "l",outside = TRUE)+
  scale_y_log10(labels=trans_format("log10",math_format(10^.x)),
                breaks=c(1,10,100,1000))+
  geom_text(aes(label = sign, y = 9.9), vjust = -0.4,size=3.5,
            data = aa %>% filter(Methods!="raw"),show.legend = F)+
  theme(axis.ticks=element_line(color="black"),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1));plot
ggsave("result/Figure2.G.rich_fold.species.ba.pdf",plot,height = 3.5,width = 4)

#OP
plotdf <- df1 %>% filter(sample_type1=="OP") %>% .[order(.$sampleNo1),]
plotdf$Methods <- "Raw"
plotdf$Methods[plotdf$methods3=="RM"] <- "R_ase"
plotdf$Methods[plotdf$methods3=="P"] <- "O_pma"
plotdf$Methods[plotdf$methods3=="W"] <- "O_ase"
plotdf$Methods[plotdf$methods3=="Fil"] <- "F_ase"
plotdf$Methods[plotdf$methods3=="S"] <- "S_ase"
plotdf$Methods[plotdf$methods3=="M"] <- "K_qia"
plotdf$Methods[plotdf$methods3=="H"] <- "K_zym"

plotdf$Methods <- factor(plotdf$Methods,
                         levels = c("Raw","R_ase","O_pma","O_ase","F_ase",
                                    "S_ase","K_qia","K_zym"))
stat <- compare_means(fold~Methods, plotdf, 
                      method = "wilcox.test",
                      p.adjust.method = "fdr",
                      paired = TRUE)
library(multcompView)
anova <- aov(richniss~Methods, data = plotdf)
tukey <- TukeyHSD(anova)
tukey[["Methods"]][,4] <- stat$p.adj
cld <- multcompLetters4(anova, tukey)
aa=cld$Methods$Letters %>% as.data.frame() %>% rownames_to_column("Methods")
colnames(aa)[2] <- "sign"
plotdf$fold[plotdf$fold>2] <- 2
plot <- ggplot(plotdf,aes(Methods,fold,fill=Methods))+
  geom_boxplot(outlier.color = NA,size=0.5)+
  geom_jitter(height = 0,width = 0.1,size=0.3,color="black")+
  scale_fill_manual(values = c(raw="black",R_ase="#696969",
                               O_pma="#7876B1",O_ase="#FFDC91",F_ase="#E18727",
                               S_ase="#BC3C29",K_qia="#20854E",K_zym="#0072B5"))+
  labs(x='',y="Fold change of species detected")+
  coord_cartesian(clip = "off")+
  geom_text(aes(label = sign, y = 2), vjust = -0.4,size=3.5,
            data = aa %>% filter(Methods!="raw"),show.legend = F)+
  theme(axis.ticks=element_line(color="black"),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1));plot
ggsave("result/Figure2.H.rich_fold.species.op.pdf",plot,height = 3.5,width = 4)
# Fig. 2K
meta <- read_excel("data.xlsx",sheet = "meta") %>% as.data.frame() %>% filter(methods3!="NC") %>% 
  .[!grepl("_2",.$sampleID),] %>% .[,c("sampleID2","sampleNo1","sample_type1","methods3")]
df1 <- read_excel("data.xlsx",sheet = "JSD.contam") %>% as.data.frame() %>% .[,c(1:3,5)] %>% 
  mutate(taxa.prop=1-contam) %>% left_join(meta)
de.tbl1 <- read_excel("data.xlsx",sheet = "de.species") %>% 
  as.data.frame() %>% column_to_rownames("tax_name")
qc <- read_excel("data.xlsx",sheet = "qc") %>% as.data.frame() %>% .[,c("sampleID2","total_reads","sequins","ABFV")] %>%
  left_join(df1) %>% mutate(bac.prop=ABFV*taxa.prop/(total_reads-sequins),total=total_reads-sequins)
meta <- read_excel("data.xlsx",sheet = "meta") %>% 
  as.data.frame() %>% 
  .[,c("sampleID2","sample_type1","methods3","sampleNo1")] %>% 
  .[!grepl("_2",.$sampleID),] %>% unique() %>% filter(methods3!="NC")
qc2 <- left_join(qc,meta) %>% filter(!is.na(methods3))
sum <- qc2 %>% dplyr::group_by(sample_type1,sampleNo1) %>% 
  dplyr::summarise(rare=min(total))
qc2 <- left_join(qc2,sum)
rare <- data.frame(tax_name=rownames(de.tbl1))
rich <- data.frame()
bac.raw <- qc2 %>% filter(methods3=="R")
for (i in qc2$sampleID2) {
  qc1.i <- filter(qc2,sampleID2==i)
  df.i <- data.frame(tax_name=rownames(de.tbl1),
                     sample=round(de.tbl1[,i])) %>% 
    column_to_rownames("tax_name")
  if (nrow(df.i)>0) {
    set.seed(100)
    rare.i <- rrarefy(df.i,qc1.i$rare) %>% t() %>% 
      as.data.frame() %>% 
      rownames_to_column("tax_name")
    
    rich.i <- data.frame(sampleID=i,
                         richniss=length(rare.i$sample[rare.i$sample>0]),
                         rare.theory=qc1.i$rare,
                         rare.real=sum(rare.i$sample))
    rich <- rbind(rich,rich.i)
    
    colnames(rare.i)[2] <- i
    rare <- left_join(rare,rare.i) 
  }
}
rich1 <- rich %>% mutate(diff=rare.theory-rare.real)
rare[is.na(rare)] <- 0
colnames(rich1)[1] <- "sampleID2"
df <- left_join(rich1,qc2)
df <- left_join(rich1,qc) %>% left_join(meta)
df.R <- df %>% filter(methods3=="R") %>% 
  .[,c("sample_type1","sampleNo1","richniss")]
colnames(df.R)[3] <- "rich.R"
df1 <- left_join(df,df.R) %>% mutate(change=richniss-rich.R)

df2 <- left_join(df1[,c("sampleID2","richniss","sample_type1","methods3","sampleNo1","rich.R")],
                 bac.raw[,c("sample_type1","sampleNo1","bac.prop")])
d2 <- df2 %>% mutate(change=richniss-rich.R)
d2$Methods <- "Raw"
d2$Methods[d2$methods3=="RM"] <- "R_ase"
d2$Methods[d2$methods3=="P"] <- "O_pma"
d2$Methods[d2$methods3=="W"] <- "O_ase"
d2$Methods[d2$methods3=="Fil"] <- "F_ase"
d2$Methods[d2$methods3=="S"] <- "S_ase"
d2$Methods[d2$methods3=="M"] <- "K_qia"
d2$Methods[d2$methods3=="H"] <- "K_zym"
plot <- ggplot(d2)+
  theme_bw()+
  geom_point(aes(bac.prop*100,change,
                 fill=Methods,color=Methods),size=0.3,alpha=0.3)+
  geom_smooth(aes(bac.prop*100,change,
                  fill=Methods,color=Methods),se=FALSE,linewidth=0.3,linetype=1,alpha=0.5)+
  scale_color_manual(values = c(R="black",R_ase="#696969",
                                O_pma="#7876B1",K_zym="#0072B5",K_qia="#20854E",O_ase="#FFDC91",
                                F_ase="#E18727",S_ase="#BC3C29"))+
  scale_fill_manual(values = c(raw="black",RM="#696969",
                               P="#7876B1",H="#0072B5",M="#20854E",W="#FFDC91",
                               Fil="#E18727",S="#BC3C29"))+
  geom_hline(yintercept = 1,linetype="dashed")+
  labs(x="Microbial proportion (%)",y='Increased number of species')+scale_x_log10()+
  scale_y_continuous(limits = c(-100,250),breaks = c(-100,-50,0,50,100,150,200,250))+
  theme(axis.text.x = element_text(hjust = 0.5,vjust = 0.5));plot
ggsave(filename = "result/Figure2.K.species.limit.pdf",plot,height = 3,width = 4)

# Fig. 2I-J, L ---------------------------------------------------------------
#Fig. 2I BALF-gene family
meta <- read_excel("data.xlsx",sheet = "meta") %>% as.data.frame() %>% 
  .[,c("sampleID2","sample_type1","methods3","sampleNo1")] %>% 
  .[!grepl("_2",.$sampleID),] %>% filter(sample_type1=="BALF"&methods3!="NC")
d <- fread("balf.gene.csv") %>% as.data.frame() %>% 
  column_to_rownames("V1") %>% .[,meta$sampleID2]
reads <- data.frame(sampleID2=colnames(d),species=colSums(d))
qc <- read_excel("data.xlsx",sheet = "qc") %>% 
  as.data.frame() %>% 
  .[,c("sampleID2","total_reads","sequins","ABFV","human")] %>%
  left_join(reads,.) %>% 
  mutate(total=total_reads-sequins,
         s.specific.prop=species/total,
         bac.prop=ABFV/total) %>% 
  .[!grepl("_2",.$sampleID2),]
all(qc$sampleID2==colnames(d))
qc <- left_join(qc,meta)
sum <- qc %>% dplyr::group_by(sample_type1,sampleNo1) %>% 
  dplyr::summarise(rare=min(total))
qc <- left_join(qc,sum) %>% 
  mutate(rare.s=round(rare*s.specific.prop))
qc1 <- qc
rare <- data.frame(tax_name=rownames(d))
rich <- data.frame()
for (i in qc1$sampleID[qc1$sampleID!="BA044"]) {
  qc1.i <- filter(qc1,sampleID2==i)
  df.i <- data.frame(tax_name=rownames(d),
                     sample=round(d[,i])) %>% 
    filter(.[,2]>0) %>% 
    column_to_rownames("tax_name")
  set.seed(100)
  rare.i <- rrarefy(df.i,qc1.i$rare.s) %>% t() %>% 
    as.data.frame() %>% 
    rownames_to_column("tax_name")
  
  rich.i <- data.frame(sampleID=i,
                       richniss=length(rare.i$sample[rare.i$sample>0]))
  rich <- rbind(rich,rich.i)
  
  colnames(rare.i)[2] <- i
  rare <- left_join(rare,rare.i)
}
rare[is.na(rare)] <- 0
colnames(rich)[1] <- "sampleID2"
df <- left_join(qc,rich)
df[is.na(df)] <- 0
df.R <- df %>% filter(methods3=="R") %>% 
  .[,c("sample_type1","sampleNo1","richniss","bac.prop")]
colnames(df.R)[3:4] <- c("rich.R","bac.prop.r")
df1 <- left_join(df,df.R) %>% mutate(fold=richniss/rich.R)
sum <- df1 %>% filter(sampleNo1!="LG02") %>% dplyr::group_by(sample_type1,methods3) %>% 
  dplyr::summarise(med.fold=median(fold),
                   Q1=quantile(fold,0.25),
                   Q3=quantile(fold,0.75))
plotdf <- df1 %>% filter(sample_type1=="BALF") %>% .[order(.$sampleNo1),]
plotdf$Methods <- "Raw"
plotdf$Methods[plotdf$methods3=="RM"] <- "R_ase"
plotdf$Methods[plotdf$methods3=="P"] <- "O_pma"
plotdf$Methods[plotdf$methods3=="W"] <- "O_ase"
plotdf$Methods[plotdf$methods3=="Fil"] <- "F_ase"
plotdf$Methods[plotdf$methods3=="S"] <- "S_ase"
plotdf$Methods[plotdf$methods3=="M"] <- "K_qia"
plotdf$Methods[plotdf$methods3=="H"] <- "K_zym"
plotdf$Methods <- factor(plotdf$Methods,
                         levels = c("Raw","R_ase","O_pma","O_ase","F_ase",
                                    "S_ase","K_qia","K_zym"))
plotdf$fold[is.nan(plotdf$fold)] <- 1
plotdf$fold[plotdf$fold>1000] <- 1000
plotdf <- plotdf %>% filter(Methods!="raw")
stat <- compare_means(fold~Methods, plotdf, 
                      method = "wilcox.test",
                      p.adjust.method = "fdr",
                      paired = TRUE)
library(multcompView)
anova <- aov(fold~Methods, data = plotdf)
tukey <- TukeyHSD(anova)
tukey[["Methods"]][,4] <- stat$p.adj
cld <- multcompLetters4(anova, tukey)
aa=cld$Methods$Letters %>% as.data.frame() %>% rownames_to_column("Methods")
colnames(aa)[2] <- "sign"

plot <- ggplot(plotdf,
               aes(Methods,fold,fill=Methods))+
  geom_boxplot(outlier.color = NA,size=0.5)+#theme_bw()+
  geom_jitter(width = 0.1,size=0.3,color="black")+
  scale_fill_manual(values = c(Raw="black",R_ase="#696969",
                               O_pma="#7876B1",O_ase="#FFDC91",F_ase="#E18727",
                               S_ase="#BC3C29",K_qia="#20854E",K_zym="#0072B5"))+
  labs(x='',y="Fold change of genefamily detected")+
  coord_cartesian(clip = "off")+
  annotation_logticks(sides = "l",outside = TRUE)+
  scale_y_log10(labels=trans_format("log10",math_format(10^.x)))+
  geom_text(aes(label = sign, y = max(plotdf$fold)), vjust = -0.4,size=3.5,
            data = aa,show.legend = F)+
  theme(axis.ticks=element_line(color="black"),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1));plot
ggsave("result/Figure2.I.rich.gene.ba.fold.change.pdf",plot,height = 3.5,width = 4)
#Fig. 2J OP-gene family
meta <- read_excel("data.xlsx",sheet = "meta") %>% as.data.frame() %>% 
  .[,c("sampleID2","sample_type1","methods3","sampleNo1")] %>% 
  .[!grepl("_2",.$sampleID),] %>% filter(sample_type1=="OP"&methods3!="NC")
d <- fread("OP.gene.csv") %>% as.data.frame() %>% column_to_rownames("V1") %>% .[,meta$sampleID2]
reads <- data.frame(sampleID2=colnames(d),species=colSums(d))
qc <- read_excel("data.xlsx",sheet = "qc") %>% as.data.frame() %>% 
  .[,c("sampleID2","total_reads","sequins","ABFV","human")] %>% left_join(reads,.) %>% 
  mutate(total=total_reads-sequins,
         s.specific.prop=species/total,
         bac.prop=ABFV/total) %>% .[!grepl("_2",.$sampleID2),]
all(qc$sampleID2==colnames(d))
qc <- left_join(qc,meta)
sum <- qc %>% dplyr::group_by(sample_type1,sampleNo1) %>% 
  dplyr::summarise(rare=min(total))
qc <- left_join(qc,sum) %>% 
  mutate(rare.s=round(rare*s.specific.prop))
qc1 <- qc
rare <- data.frame(tax_name=rownames(d))
rich <- data.frame()
for (i in qc1$sampleID) {
  qc1.i <- filter(qc1,sampleID2==i)
  df.i <- data.frame(tax_name=rownames(d),
                     sample=round(d[,i])) %>% 
    filter(.[,2]>0) %>% 
    column_to_rownames("tax_name")
  set.seed(100)
  rare.i <- rrarefy(df.i,qc1.i$rare.s) %>% t() %>% 
    as.data.frame() %>% 
    rownames_to_column("tax_name")
  
  rich.i <- data.frame(sampleID=i,
                       richniss=length(rare.i$sample[rare.i$sample>0]))
  rich <- rbind(rich,rich.i)
  
  colnames(rare.i)[2] <- i
  rare <- left_join(rare,rare.i)
}
rare[is.na(rare)] <- 0
colnames(rich)[1] <- "sampleID2"
df <- left_join(qc,rich)
df.R <- df %>% filter(methods3=="R") %>% 
  .[,c("sample_type1","sampleNo1","richniss","bac.prop")]
colnames(df.R)[3:4] <- c("rich.R","bac.prop.r")
df1 <- left_join(df,df.R) %>% mutate(fold=richniss/rich.R)
sum <- df1 %>% filter(sampleNo1!="LG02") %>% dplyr::group_by(sample_type1,methods3) %>% 
  dplyr::summarise(med.fold=median(fold),
                   Q1=quantile(fold,0.25),
                   Q3=quantile(fold,0.75))
plotdf <- df1 %>% filter(sample_type1=="OP") %>% .[order(.$sampleNo1),]
plotdf$Methods <- "Raw"
plotdf$Methods[plotdf$methods3=="RM"] <- "R_ase"
plotdf$Methods[plotdf$methods3=="P"] <- "O_pma"
plotdf$Methods[plotdf$methods3=="W"] <- "O_ase"
plotdf$Methods[plotdf$methods3=="Fil"] <- "F_ase"
plotdf$Methods[plotdf$methods3=="S"] <- "S_ase"
plotdf$Methods[plotdf$methods3=="M"] <- "K_qia"
plotdf$Methods[plotdf$methods3=="H"] <- "K_zym"
plotdf$Methods <- factor(plotdf$Methods,
                         levels = c("Raw","R_ase","O_pma","O_ase","F_ase",
                                    "S_ase","K_qia","K_zym"))
stat <- compare_means(fold~Methods, plotdf, 
                      method = "wilcox.test",
                      p.adjust.method = "fdr",
                      paired = TRUE)
library(multcompView)
anova <- aov(fold~Methods, data = plotdf)
tukey <- TukeyHSD(anova)
tukey[["Methods"]][,4] <- stat$p.adj
cld <- multcompLetters4(anova, tukey)
aa=cld$Methods$Letters %>% as.data.frame() %>% rownames_to_column("Methods")
colnames(aa)[2] <- "sign"
plot <- ggplot(plotdf,
               aes(Methods,fold,fill=Methods))+
  geom_boxplot(outlier.color = NA,size=0.5)+
  geom_jitter(height = 0,width = 0.1,size=0.3,color="black")+
  scale_fill_manual(values = c(Raw="black",R_ase="#696969",
                               O_pma="#7876B1",O_ase="#FFDC91",F_ase="#E18727",
                               S_ase="#BC3C29",K_qia="#20854E",K_zym="#0072B5"))+
  labs(x='',y="Fold change of gene family detected")+
  coord_cartesian(clip = "off")+
  annotation_logticks(sides = "l",outside = TRUE)+
  scale_y_log10(breaks=c(0.1,1,10),labels=trans_format("log10",math_format(10^.x)))+#
  geom_text(aes(label = sign, y = max(plotdf$fold)), vjust = -0.4,size=3.5,
            data = aa,show.legend = F)+
  theme(axis.ticks=element_line(color="black"),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1));plot
ggsave("result/Figure2.J.rich.gene.fold.OP.pdf",plot,height = 3.5,width = 4)

#fIG. 2L
meta <- read_excel("data.xlsx",sheet = "meta") %>% as.data.frame() %>% filter(methods3!="NC") %>% 
  .[!grepl("_2",.$sampleID),] %>% .[,c("sampleID2","sampleNo1","sample_type1","methods3")]
df1 <- read_excel("data.xlsx",sheet = "JSD.contam") %>% as.data.frame() %>% .[,c(1:3,5)] %>% 
  mutate(taxa.prop=1-contam) %>% left_join(meta)
d <- fread("balf.gene.csv") %>% 
  as.data.frame() %>% 
  column_to_rownames("V1") %>% .[,meta$sampleID2[grepl("BA",meta$sampleID2)]]
reads <- data.frame(sampleID2=colnames(d),species=colSums(d))
qc <- read_excel("data.xlsx",sheet = "qc") %>% 
  as.data.frame() %>% 
  .[,c("sampleID2","total_reads","sequins","ABFV","human")] %>%
  left_join(reads,.) %>% left_join(df1) %>% 
  mutate(total=total_reads-sequins,
         s.specific.prop=species/total,
         bac.prop=ABFV*taxa.prop/total) %>% 
  .[!grepl("_2",.$sampleID2),]
all(qc$sampleID2==colnames(d))
qc <- left_join(qc,meta)
sum <- qc %>% dplyr::group_by(sample_type1,sampleNo1) %>% 
  dplyr::summarise(rare=min(total))
qc <- left_join(qc,sum) %>% 
  mutate(rare.s=round(rare*s.specific.prop))
qc1 <- qc
rare <- data.frame(tax_name=rownames(d))
rich <- data.frame()
for (i in qc1$sampleID[qc1$sampleID!="BA044"]) {
  qc1.i <- filter(qc1,sampleID2==i)
  df.i <- data.frame(tax_name=rownames(d),
                     sample=round(d[,i])) %>% 
    filter(.[,2]>0) %>% 
    column_to_rownames("tax_name")
  set.seed(100)
  rare.i <- rrarefy(df.i,qc1.i$rare.s) %>% t() %>% 
    as.data.frame() %>% 
    rownames_to_column("tax_name")
  
  rich.i <- data.frame(sampleID=i,
                       richniss=length(rare.i$sample[rare.i$sample>0]))
  rich <- rbind(rich,rich.i)
  
  colnames(rare.i)[2] <- i
  rare <- left_join(rare,rare.i)
}
rare[is.na(rare)] <- 0
colnames(rich)[1] <- "sampleID2"
df <- left_join(qc,rich)
df[is.na(df)] <- 0
df.R <- df %>% filter(methods3=="R") %>% 
  .[,c("sample_type1","sampleNo1","richniss","bac.prop")]
colnames(df.R)[3:4] <- c("rich.R","bac.prop.r")
df1 <- left_join(df,df.R) %>% mutate(change=richniss-rich.R)
data.ba <- df1[,c("sampleID2","sample_type1","methods3","sampleNo1","change","bac.prop.r")]
data.ba <- data.ba %>% as.data.frame() %>% .[,2:6] %>% left_join(qc %>% filter(methods3=="R") %>% 
            .[,c("sample_type1","sampleNo1","bac.prop")])
colnames(data.ba)[6] <- "bac.prop.r"
#OP
meta <- read_excel("data.xlsx",sheet = "meta") %>% as.data.frame() %>% filter(methods3!="NC") %>% 
  .[!grepl("_2",.$sampleID),] %>% .[,c("sampleID2","sampleNo1","sample_type1","methods3")]
df1 <- read_excel("data.xlsx",sheet = "JSD.contam") %>% as.data.frame() %>% .[,c(1:3,5)] %>% 
  mutate(taxa.prop=1-contam) %>% left_join(meta)
meta <- read_excel("data.xlsx",sheet = "meta") %>% as.data.frame() %>% 
  .[,c("sampleID2","sample_type1","methods3","sampleNo1")] %>% 
  .[!grepl("_2",.$sampleID),] %>% filter(sample_type1=="OP"&methods3!="NC")
d <- fread("op.gene.csv") %>% 
  as.data.frame() %>% 
  column_to_rownames("V1") %>% .[,meta$sampleID2]
reads <- data.frame(sampleID2=colnames(d),species=colSums(d))
qc <- read_excel("data.xlsx",sheet = "qc") %>% as.data.frame() %>% 
  .[,c("sampleID2","total_reads","sequins","ABFV","human")] %>% left_join(df1) %>% left_join(reads,.) %>% 
  mutate(total=total_reads-sequins,
         s.specific.prop=species/total,
         bac.prop=ABFV*taxa.prop/total) %>% .[!grepl("_2",.$sampleID2),]
all(qc$sampleID2==colnames(d))

qc <- left_join(qc,meta)
sum <- qc %>% dplyr::group_by(sample_type1,sampleNo1) %>% 
  dplyr::summarise(rare=min(total))
qc <- left_join(qc,sum) %>% 
  mutate(rare.s=round(rare*s.specific.prop))
qc1 <- qc
rare <- data.frame(tax_name=rownames(d))
rich <- data.frame()
for (i in qc1$sampleID) {
  qc1.i <- filter(qc1,sampleID2==i)
  df.i <- data.frame(tax_name=rownames(d),
                     sample=round(d[,i])) %>% 
    filter(.[,2]>0) %>% 
    column_to_rownames("tax_name")
  set.seed(100)
  rare.i <- rrarefy(df.i,qc1.i$rare.s) %>% t() %>% 
    as.data.frame() %>% 
    rownames_to_column("tax_name")
  
  rich.i <- data.frame(sampleID=i,
                       richniss=length(rare.i$sample[rare.i$sample>0]))
  rich <- rbind(rich,rich.i)
  
  colnames(rare.i)[2] <- i
  rare <- left_join(rare,rare.i)
}
rare[is.na(rare)] <- 0
colnames(rich)[1] <- "sampleID2"
df <- left_join(qc,rich)
df.R <- df %>% filter(methods3=="R") %>% 
  .[,c("sample_type1","sampleNo1","richniss","bac.prop")]
colnames(df.R)[3:4] <- c("rich.R","bac.prop.r")
df1 <- left_join(df,df.R) %>% mutate(change=richniss-rich.R)
data.op <- df1[,c("sampleID2","sample_type1","methods3","sampleNo1","change","bac.prop.r")]
plotdf <- rbind(data.ba,data.op)
plotdf$Methods <- "Raw"
plotdf$Methods[plotdf$methods3=="RM"] <- "R_ase"
plotdf$Methods[plotdf$methods3=="P"] <- "O_pma"
plotdf$Methods[plotdf$methods3=="W"] <- "O_ase"
plotdf$Methods[plotdf$methods3=="Fil"] <- "F_ase"
plotdf$Methods[plotdf$methods3=="S"] <- "S_ase"
plotdf$Methods[plotdf$methods3=="M"] <- "K_qia"
plotdf$Methods[plotdf$methods3=="H"] <- "K_zym"
plotdf$Methods <- factor(plotdf$Methods,
                         levels = c("Raw","R_ase","O_pma","O_ase","F_ase",
                                    "S_ase","K_qia","K_zym"))
plotdf1 <- plotdf
plotdf1$bac.prop.r2 <- plotdf1$bac.prop.r
plot <- ggplot(plotdf1)+theme_bw()+
  geom_point(aes(bac.prop.r*100,change,
                 fill=Methods,color=Methods),size=0.3,alpha=0.3)+
  geom_smooth(aes(bac.prop.r*100,change,
                  fill=Methods,color=Methods),se=FALSE,size=0.3,linetype=1,alpha=0.5)+
  scale_color_manual(values = c(Raw="black",R_ase="#696969",
                                O_pma="#7876B1",K_zym="#0072B5",K_qia="#20854E",O_ase="#FFDC91",
                                F_ase="#E18727",S_ase="#BC3C29"))+
  scale_fill_manual(values = c(Raw="black",R_ase="#696969",
                               O_pma="#7876B1",K_zym="#0072B5",K_qia="#20854E",O_ase="#FFDC91",
                               F_ase="#E18727",S_ase="#BC3C29"))+
  geom_hline(yintercept = 1,linetype="dashed")+
  scale_y_continuous(limits = c(-100000,100000),breaks = c(-100000,-50000,0,50000,100000))+
  labs(x="Microbial proportion (%)",y='Increased number of gene family')+
  scale_x_log10(breaks=c(0.001,0.01,0.1,1,10,30,100))+
  theme(axis.text.x = element_text(hjust = 0.5,vjust = 0.5));plot
ggsave(filename = "result/Figure2.L.gene.limit.pdf",plot,height = 3,width = 4.1)
# Fig. S2 -----------------------------------------------------------------
meta <- read_excel("data.xlsx",sheet = "meta") %>% as.data.frame() %>% filter(methods3!="NC") %>% 
  .[!grepl("_2",.$sampleID),]
df <- read_excel("data.xlsx",sheet = "qc") %>%
  as.data.frame() %>% 
  .[,c("sampleID2","total_reads","sequins","low_quality","human","ABFV","nonABFV")] %>% filter(sampleID2%in%meta$sampleID)
colnames(df)[1] <- "sampleID"
df1 <- read_excel("data.xlsx",sheet = "JSD.contam") %>% as.data.frame() %>% .[,c(1:3,5)] %>% 
  left_join(meta[,c("sampleNo1","sample_type1","methods3","sampleID")])
df2 <- left_join(df,df1) %>% 
  mutate(abfv_prop=ABFV/(total_reads-sequins),
         host_prop=human/(total_reads-sequins),
         nonabfv_prop=nonABFV/(total_reads-sequins),
         low_quality_prop=low_quality/(total_reads-sequins),
         abfv.taxa=abfv_prop*(1-contam),
         abfv.contam=abfv_prop*contam) %>%
  .[,c("sampleID","low_quality_prop","abfv.contam","abfv.taxa",
       "host_prop","nonabfv_prop")] %>% reshape2::melt(variable.name="class") %>% 
  left_join(meta[,c("sampleID","sampleNo1","sample_type1","methods3")])

summary <- df2 %>% 
  dplyr::group_by(sample_type1,methods3,class) %>% 
  dplyr::summarise(Q1=quantile(value, 0.25),
                   med=median(value),
                   Q3=quantile(value, 0.75),
                   mean=mean(value),
                   min=min(value),
                   max=max(value))
df2$Methods <- "Raw"
df2$Methods[df2$methods3=="RM"] <- "R_ase"
df2$Methods[df2$methods3=="P"] <- "O_pma"
df2$Methods[df2$methods3=="W"] <- "O_ase"
df2$Methods[df2$methods3=="Fil"] <- "F_ase"
df2$Methods[df2$methods3=="S"] <- "S_ase"
df2$Methods[df2$methods3=="M"] <- "K_qia"
df2$Methods[df2$methods3=="H"] <- "K_zym"
#bacterial proportion
plotdf <- df2 %>% 
  filter(sample_type1=="BALF") %>% 
  .[order(.$sampleNo1),]

plotdf$Methods <- factor(plotdf$Methods,
                         levels = c("Raw","R_ase","O_pma","O_ase","F_ase",
                                    "S_ase","K_qia","K_zym"))
##overview
sum <- plotdf %>% dplyr::group_by(sample_type1,Methods,class) %>% 
  dplyr::summarise(mean=mean(value))
sum$class <- factor(sum$class,
                    levels = rev(c("abfv.taxa","abfv.contam","nonabfv_prop",
                                   "low_quality_prop","host_prop")))
sum$Methods <- factor(sum$Methods,
                      levels = c("Raw","R_ase","O_pma","O_ase","F_ase",
                                 "S_ase","K_qia","K_zym"))
p1 <- ggplot(sum,aes(Methods,mean,color=class,
                     fill=class))+
  geom_col(position = "stack")+
  scale_color_manual(values = c("grey90","grey75","grey60","grey45","grey30"))+
  scale_fill_manual(values = c("grey90","grey75","grey60","grey45","grey30"))+
  labs(x='',y='Proportion')+
  theme(axis.ticks=element_line(color="black"),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1));p1
ggsave("result/Fig.S2A.pdf",p1,height = 3.5,width = 5.5)

#OP
plotdf <- df2 %>% filter(sample_type1=="OP") %>% .[order(.$sampleNo1),]
plotdf$Methods <- factor(plotdf$Methods,
                         levels = c("Raw","R_ase","O_pma","O_ase","F_ase",
                                    "S_ase","K_qia","K_zym"))
##overview
sum <- plotdf %>% dplyr::group_by(sample_type1,Methods,class) %>% 
  dplyr::summarise(mean=mean(value))
sum$class <- factor(sum$class,
                    levels = rev(c("abfv.taxa","abfv.contam","nonabfv_prop",
                                   "low_quality_prop","host_prop")))
sum$Methods <- factor(sum$Methods,
                      levels = c("Raw","R_ase","O_pma","O_ase","F_ase",
                                 "S_ase","K_qia","K_zym"))
p1 <- ggplot(sum,aes(Methods,mean,color=class,
                     fill=class))+
  geom_col(position = "stack")+
  scale_color_manual(values = c("grey90","grey75","grey60","grey45","grey30"))+
  scale_fill_manual(values = c("grey90","grey75","grey60","grey45","grey30"))+
  labs(x='',y='Proportion')+
  theme(axis.ticks=element_line(color="black"),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1));p1
ggsave("result/Fig.S2B.pdf",p1,height = 3.5,width = 5.5)

# Fig. S3 -----------------------------------------------------------------
setwd("G:/Lilab/7.enrichment prject/2022/data/coverage/depth0214")
df <- read.csv("summary.comp2.all.csv") %>% .[,2:7]
setwd("G:/Lilab/7.enrichment prject/2022/data/coverage")
taxa <- read_excel("coverage.xlsx",sheet = "decontam.tbl") %>% as.data.frame() %>% 
  column_to_rownames("tax_name") %>% apply(2,function(x){x/sum(x)}) %>% 
  as.data.frame() %>% rownames_to_column("tax_name") %>% 
  reshape2::melt(variable.name="sampleID") %>% 
  left_join(read_excel("coverage.xlsx",sheet = "meta") %>% as.data.frame() %>% 
              .[,c("sampleID","oldName","sample_type1","methods3","sampleNo1")]) %>% 
  .[!grepl("_2",.$sampleID),]
taxa.rm <- taxa %>% filter(methods3=="RM"&value>0)
bias.tax <- read_excel("coverage.xlsx",sheet = "bias.taxa") %>% left_join(taxa.rm)

coverage <- left_join(df,read_excel("coverage.xlsx",sheet = "genome") %>% 
                        as.data.frame() %>% .[,c("taxid","tax_name")]) %>% 
  left_join(read_excel("coverage.xlsx",sheet = "meta") %>% as.data.frame() %>% 
              .[,c("sampleID","sample_type1","sampleNo1","methods3")]) %>% 
  left_join(bias.tax[,c("tax_name","sample_type1","sampleNo1","value")]) %>% 
  filter(!is.na(.$value)) %>% .[!grepl("_2",.$sampleID),]

merge1 <- coverage %>% filter(sample_type1=="BALF")
b <- data.frame()
for (i in unique(merge1$tax_name)) {
  for (j in c("RM","P","W","Fil","S","M","H")) {
    d1.1 <- merge1 %>% filter(tax_name==i&methods3%in%c(j,"R")) %>% 
      .[,c("tax_name","sampleNo1","sample_type1","methods3",
           "cover.prop")] %>% spread(methods3,cover.prop)
    d1.1[is.na(d1.1)] <- 0
    d1.2 <- d1.1[,1:5] %>% reshape2::melt(variable.name="methods3")
    wilcox1 <- compare_means(value~methods3, d1.2, 
                             method = "wilcox.test",
                             p.adjust.method = "fdr",paired = TRUE)
    
    d1.1 <- d1.1 %>% mutate(fold=.[,j]/R) %>% filter(!is.nan(fold))
    d1.1$fold[is.infinite(d1.1$fold)] <- max(d1.1$fold[!is.infinite(d1.1$fold)])
    
    b.ij <- data.frame(tax_name=i,methods3=j,sample_type1="BALF",
                       med.cover.prop=median(d1.1$fold),
                       p.cover.prop=wilcox1$p)
    b <- rbind(b,b.ij)
  }
}

b1 <- spread(b[,c("tax_name","methods3","med.cover.prop")],methods3,med.cover.prop) %>% 
  .[,c("tax_name","RM","P","W","Fil","S","M","H")]
# write.csv(b1,"result/b.fold.cover.prop.csv")
b2 <- spread(b[,c("tax_name","methods3","p.cover.prop")],methods3,p.cover.prop) %>% 
  .[,c("tax_name","RM","P","W","Fil","S","M","H")]
b3 <- data.frame(tax_name=b2$tax_name,
                 p.adj.RM=p.adjust(b2$RM,method = "fdr"),
                 p.adj.P=p.adjust(b2$P,method = "fdr"),
                 p.adj.W=p.adjust(b2$W,method = "fdr"),
                 p.adj.Fil=p.adjust(b2$Fil,method = "fdr"),
                 p.adj.S=p.adjust(b2$S,method = "fdr"),
                 p.adj.M=p.adjust(b2$M,method = "fdr"),
                 p.adj.H=p.adjust(b2$H,method = "fdr"))
# write.csv(b1,"result/b.fold.cover.prop.pvalue.csv")
result <- left_join(b1,b3)
result2 <- read_excel("coverage.xlsx",sheet = "bias.taxa") %>% filter(sample_type1=="BALF") %>% 
  left_join(result) %>% filter(!is.na(RM))
write.csv(result2,"result/b.coverage.pvalue.csv")
#OP
setwd("F:/1.working file/7.enrichment prject/2022/data/coverage/depth0214")
df <- read.csv("summary.comp2.all.csv") %>% .[,2:7]
setwd("F:/1.working file/7.enrichment prject/2022/data/coverage")
taxa <- read_excel("coverage.xlsx",sheet = "decontam.tbl") %>% as.data.frame() %>% 
  column_to_rownames("tax_name") %>% apply(2,function(x){x/sum(x)}) %>% 
  as.data.frame() %>% rownames_to_column("tax_name") %>% 
  reshape2::melt(variable.name="sampleID") %>% 
  left_join(read_excel("coverage.xlsx",sheet = "meta") %>% as.data.frame() %>% 
              .[,c("sampleID","oldName","sample_type1","methods3","sampleNo1")]) %>% 
  .[!grepl("_2",.$sampleID),]
taxa.rm <- taxa %>% filter(methods3=="RM"&value>0)
bias.tax <- read_excel("coverage.xlsx",sheet = "bias.taxa") %>% left_join(taxa.rm)

coverage <- left_join(df,read_excel("coverage.xlsx",sheet = "genome") %>% 
                        as.data.frame() %>% .[,c("taxid","tax_name")]) %>% 
  left_join(read_excel("coverage.xlsx",sheet = "meta") %>% as.data.frame() %>% 
              .[,c("sampleID","sample_type1","sampleNo1","methods3")]) %>% 
  left_join(bias.tax[,c("tax_name","sample_type1","sampleNo1","value")]) %>% 
  filter(!is.na(.$value)) %>% .[!grepl("_2",.$sampleID),]

merge1 <- coverage %>% filter(sample_type1=="OP")

# summary <- merge1 %>% dplyr::group_by(tax_name,methods3) %>% 
#     dplyr::summarise(med=median(cover.prop)) %>% 
#     spread(methods3,med) %>% .[,c("tax_name","R","RM","P","W","Fil","S","M","H")]
# write.csv(summary,"result/coverage.med.csv")
b <- data.frame()
for (i in unique(merge1$tax_name)) {
  for (j in c("RM","P","W","Fil","S","M","H")) {
    d1.1 <- merge1 %>% filter(tax_name==i&methods3%in%c(j,"R")) %>% 
      .[,c("tax_name","sampleNo1","sample_type1","methods3",
           "cover.prop")] %>% spread(methods3,cover.prop)
    d1.1[is.na(d1.1)] <- 0
    d1.2 <- d1.1[,1:5] %>% reshape2::melt(variable.name="methods3")
    wilcox1 <- compare_means(value~methods3, d1.2, 
                             method = "wilcox.test",
                             p.adjust.method = "fdr",paired = TRUE)
    
    d1.1 <- d1.1 %>% mutate(fold=.[,j]/R) %>% filter(!is.nan(fold))
    d1.1$fold[is.infinite(d1.1$fold)] <- max(d1.1$fold[!is.infinite(d1.1$fold)])
    
    b.ij <- data.frame(tax_name=i,methods3=j,sample_type1="OP",
                       med.cover.prop=median(d1.1$fold),
                       p.cover.prop=wilcox1$p)
    b <- rbind(b,b.ij)
  }
}

b1 <- spread(b[,c("tax_name","methods3","med.cover.prop")],methods3,med.cover.prop) %>% 
  .[,c("tax_name","RM","P","W","Fil","S","M","H")]
# write.csv(b1,"result/b.fold.cover.prop.csv")
b2 <- spread(b[,c("tax_name","methods3","p.cover.prop")],methods3,p.cover.prop) %>% 
  .[,c("tax_name","RM","P","W","Fil","S","M","H")]
b3 <- data.frame(tax_name=b2$tax_name,
                 p.adj.RM=p.adjust(b2$RM,method = "fdr"),
                 p.adj.P=p.adjust(b2$P,method = "fdr"),
                 p.adj.W=p.adjust(b2$W,method = "fdr"),
                 p.adj.Fil=p.adjust(b2$Fil,method = "fdr"),
                 p.adj.S=p.adjust(b2$S,method = "fdr"),
                 p.adj.M=p.adjust(b2$M,method = "fdr"),
                 p.adj.H=p.adjust(b2$H,method = "fdr"))
# write.csv(b1,"result/b.fold.cover.prop.pvalue.csv")
result <- left_join(b1,b3)
result2 <- read_excel("coverage.xlsx",sheet = "bias.taxa") %>% filter(sample_type1=="OP") %>% 
  left_join(result) #%>% filter(!is.na(RM))
write.csv(result2,"result/o.coverage.pvalue.csv")
# heatmap -----------------------------------------------------------------
#balf
plotdf <- read_excel("coverage.xlsx",sheet = "cover.1011",
                     col_types = c(rep("text",3),rep("numeric",14))) %>% 
  filter(sample_type1=="BALF"&class=="cover.prop") %>% 
  .[order(.$tax_name),]
rownames(plotdf) <- plotdf$tax_name
plotdf.f <- plotdf[,11:17]
rownames(plotdf.f) <- plotdf$tax_name
# plotdf.f[plotdf.f>150] <- 150
plotdf.f1 <- log(plotdf.f,2)
rownames(plotdf.f1) <- plotdf$tax_name
range(plotdf.f1)
plotdf.p <- plotdf[,4:10] %>% apply(2,as.numeric)
plotdf.p[plotdf.p<0.001] <- "***"
plotdf.p[plotdf.p>0.001&plotdf.p<0.01] <- "**"
plotdf.p[plotdf.p>0.01&plotdf.p<0.05] <- "*"
plotdf.p[plotdf.p>0.05] <- ""
colnames(plotdf.p) <- colnames(plotdf.f)
col_fun = circlize::colorRamp2(c(-3, 0, 10), c("#006837", "white", "#AB0626"))
# col_fun = circlize::colorRamp2(c(-2.5, 0, 8.4), c("#006837", "white", "#AB0626"))
# col_fun = circlize::colorRamp2(c(0, 160, 324), c("#006837", "white", "#AB0626"))
# col_fun = circlize::colorRamp2(c(0,50,100,150), 
#                                     c("white","#ffe4da","#ff8f70","#ff0000"))
library("ComplexHeatmap")
Heatmap(plotdf.f1, "log2(fold.change)",col = col_fun,
        cluster_rows = FALSE,cluster_columns = FALSE,
        row_names_side = "left",column_names_side = "top",
        column_names_rot = 90,
        rect_gp = gpar(col = "white", lwd = 1),
        cell_fun = function(j, i, x, y, width, height, fill) 
        {if( !is.na(plotdf.p[i,j]))
          grid.text(sprintf("%s", plotdf.p[i, j]), x, y, 
                    gp = gpar(fontsize = 10))
        },row_title_rot = 0,border = T)
#OP 1011
plotdf <- read_excel("coverage.xlsx",sheet = "cover.1011",
                     col_types = c(rep("text",3),rep("numeric",14))) %>% 
  filter(sample_type1=="OP"&class=="cover.prop") %>% .[order(.$tax_name),]
rownames(plotdf) <- plotdf$tax_name
plotdf.f <- plotdf[,11:17]
rownames(plotdf.f) <- plotdf$tax_name
# plotdf.f[plotdf.f>150] <- 150
plotdf.f1 <- log(plotdf.f,2)
range(plotdf.f1)
rownames(plotdf.f1) <- plotdf$tax_name
plotdf.p <- plotdf[,4:10] %>% apply(2,as.numeric)
plotdf.p[plotdf.p<0.001] <- "***"
plotdf.p[plotdf.p>0.001&plotdf.p<0.01] <- "**"
plotdf.p[plotdf.p>0.01&plotdf.p<0.05] <- "*"
plotdf.p[plotdf.p>0.05] <- ""
colnames(plotdf.p) <- colnames(plotdf.f)
col_fun = circlize::colorRamp2(c(-3, 0, 9.5), c("#006837", "white", "#AB0626"))
# col_fun = circlize::colorRamp2(c(-2.5, 0, 8.4), c("#006837", "white", "#AB0626"))
# col_fun = circlize::colorRamp2(c(0, 160, 324), c("#006837", "white", "#AB0626"))
# col_fun = circlize::colorRamp2(c(0,50,100,150), 
#                                     c("white","#ffe4da","#ff8f70","#ff0000"))
library("ComplexHeatmap")
Heatmap(plotdf.f1, "log2(fold.change)",col = col_fun,
        cluster_rows = FALSE,cluster_columns = FALSE,
        row_names_side = "left",column_names_side = "top",
        column_names_rot = 90,
        rect_gp = gpar(col = "white", lwd = 1),
        cell_fun = function(j, i, x, y, width, height, fill) 
        {if( !is.na(plotdf.p[i,j]))
          grid.text(sprintf("%s", plotdf.p[i, j]), x, y, 
                    gp = gpar(fontsize = 10))
        },row_title_rot = 0,border = T)
##op样本与balf使用相同的量程，使OP和BALF之间的趋势可比
#OP
plotdf <- read_excel("coverage.xlsx",sheet = "cover.depth",
                     col_types = c(rep("text",10),rep("numeric",7))) %>% 
  filter(sample_type1=="OP"&class=="cover.prop")
rownames(plotdf) <- plotdf$tax_name
plotdf.f <- plotdf[,11:17]
rownames(plotdf.f) <- plotdf$tax_name
# plotdf.f[plotdf.f>150] <- 150
plotdf.f1 <- log(plotdf.f,2)
rownames(plotdf.f1) <- plotdf$tax_name
plotdf.p <- plotdf[,4:10]
colnames(plotdf.p) <- colnames(plotdf.f)
col_fun = circlize::colorRamp2(c(-3, 0, 8.5), c("#006837", "white", "#AB0626"))
# col_fun = circlize::colorRamp2(c(-2.5, 0, 8.4), c("#006837", "white", "#AB0626"))
# col_fun = circlize::colorRamp2(c(0, 160, 324), c("#006837", "white", "#AB0626"))
# col_fun = circlize::colorRamp2(c(0,50,100,150), 
#                                     c("white","#ffe4da","#ff8f70","#ff0000"))
library("ComplexHeatmap")
Heatmap(plotdf.f1, "log2(fold.change)",col = col_fun,
        cluster_rows = FALSE,cluster_columns = FALSE,
        row_names_side = "left",column_names_side = "top",
        column_names_rot = 90,
        rect_gp = gpar(col = "white", lwd = 1),
        cell_fun = function(j, i, x, y, width, height, fill) 
        {if( !is.na(plotdf.p[i,j]))
          grid.text(sprintf("%s", plotdf.p[i, j]), x, y, 
                    gp = gpar(fontsize = 10))
        },row_title_rot = 0,border = T)
##depth-OP
plotdf <- read_excel("coverage.xlsx",sheet = "cover.depth",
                     col_types = c(rep("text",10),rep("numeric",7))) %>% 
  filter(sample_type1=="OP"&class=="mean.depth")
rownames(plotdf) <- plotdf$tax_name
plotdf.f <- plotdf[,11:17]
rownames(plotdf.f) <- plotdf$tax_name
# plotdf.f[plotdf.f>150] <- 150
plotdf.f1 <- log(plotdf.f,2)
rownames(plotdf.f1) <- plotdf$tax_name
plotdf.p <- plotdf[,4:10]
colnames(plotdf.p) <- colnames(plotdf.f)
col_fun = circlize::colorRamp2(c(-1, 0, 5.5), c("#006837", "white", "#AB0626"))
# col_fun = circlize::colorRamp2(c(-2.5, 0, 8.4), c("#006837", "white", "#AB0626"))
# col_fun = circlize::colorRamp2(c(0, 160, 324), c("#006837", "white", "#AB0626"))
# col_fun = circlize::colorRamp2(c(0,50,100,150), 
#                                     c("white","#ffe4da","#ff8f70","#ff0000"))
library("ComplexHeatmap")
Heatmap(plotdf.f1, "log2(fold.change)",col = col_fun,
        cluster_rows = FALSE,cluster_columns = FALSE,
        row_names_side = "left",column_names_side = "top",
        column_names_rot = 90,
        rect_gp = gpar(col = "white", lwd = 1),
        cell_fun = function(j, i, x, y, width, height, fill) 
        {if( !is.na(plotdf.p[i,j]))
          grid.text(sprintf("%s", plotdf.p[i, j]), x, y, 
                    gp = gpar(fontsize = 10))
        },row_title_rot = 0,border = T)

# Fig3A,B ------------------------------------------------------------------
setwd("F:/Code")
JSD2 <- read_excel("data.xlsx",sheet = "JSD.contam") %>% as.data.frame() %>% mutate(contam=contam*100)
#Fig. 3A
plotdata <- JSD2 %>% filter(sample_type1=="BALF")
summary <- plotdata %>% dplyr::group_by(methods3,sample_type1) %>% 
  dplyr::summarise(
    mean.rm=median(value),
    sd.rm.Q1=quantile(value,0.25),
    sd.rm.Q3=quantile(value,0.75),
    mean.nc=median(contam),
    sd.nc.Q1=quantile(contam,0.25),
    sd.nc.Q3=quantile(contam,0.75)
  )
write.csv(summary,"result/med.ba.accuracy.csv")
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
                    ymin=sd.nc.Q1,
                    ymax=sd.nc.Q3),
                width=0.02,size=0.1)+
  geom_errorbarh(aes(y=mean.nc,
                     xmin=sd.rm.Q1,
                     xmax=sd.rm.Q3),height=1,size=0.1)+
  geom_point(aes(mean.rm,mean.nc),size=1)+
  scale_y_continuous(breaks = seq(0,50,10),labels = seq(0,50,10),limits = c(0,50))+
  labs(x='JSD to R_ase',y='Contamination proportion (%)');plot  
ggsave("result/JSD2RM vs contam.b.pdf",plot,height = 3.5,width = 4.5)
#Fig. 3B
plotdata <- JSD2 %>% filter(sample_type1=="OP")
summary2 <- plotdata %>% dplyr::group_by(methods3,sample_type1) %>% 
  dplyr::summarise(
    mean.rm=median(value),
    sd.rm.Q1=quantile(value,0.25),
    sd.rm.Q3=quantile(value,0.75),
    mean.nc=median(contam),
    sd.nc.Q1=quantile(contam,0.25),
    sd.nc.Q3=quantile(contam,0.75)
  )
write.csv(summary,"result/med.op.accuracy.csv")
plotdata <- left_join(plotdata,summary2)
plot <- ggplot(plotdata,aes(value,contam,color=methods3))+
  geom_point(alpha=0.5,size=0.2)+theme_bw()+
  geom_rect(aes(xmin=0.09,xmax=0.33,ymin=-Inf,ymax=Inf),
            fill="#e2e2e2",color=NA)+
  geom_rect(aes(xmin=0.44,xmax=0.83,ymin=-Inf,ymax=Inf),
            fill="#bababa",alpha=0.02,color=NA)+
  scale_color_manual(values = c(P="#7876B1",H="#0072B5",M="#20854E",R="black",
                                W="#FFDC91",Fil="#E18727",S="#BC3C29",RM="#696969"))+
  geom_errorbar(aes(x=mean.rm,
                    ymin=sd.nc.Q1,
                    ymax=sd.nc.Q3),
                width=0.02,size=0.1)+
  geom_errorbarh(aes(y=mean.nc,
                     xmin=sd.rm.Q1,
                     xmax=sd.rm.Q3),height=1,size=0.1)+
  geom_point(aes(mean.rm,mean.nc),size=1)+
  scale_y_continuous(breaks = seq(0,50,10),labels = seq(0,50,10),limits = c(0,50))+
  labs(x='JSD to R_ase',y='Contamination proportion');plot #,color="black"
ggsave("result/JSD2RM vs contam.o.pdf",plot,height = 3.5,width = 4.5)

# Fig3A,B(2) ------------------------------------------------------------------
setwd("F:/Code")
JSD2 <- read_excel("data.xlsx",sheet = "JSD.contam") %>% as.data.frame() %>% mutate(contam=contam*100)
#Fig. 3A
plotdata <- JSD2 %>% filter(sample_type1=="BALF")
summary <- plotdata %>% dplyr::group_by(methods3,sample_type1) %>% 
  dplyr::summarise(
    mean.rm=median(value),
    sd.rm.Q1=quantile(value,0.25),
    sd.rm.Q3=quantile(value,0.75),
    sd.rm=sd(value),
    mean.nc=median(contam),
    sd.nc.Q1=quantile(contam,0.25),
    sd.nc.Q3=quantile(contam,0.75),
    sd.nc=sd(contam),
  )
write.csv(summary,"result/med.ba.accuracy.csv")
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
                     xmax=mean.rm+sd.rm),height=1,size=0.1)+
  geom_point(aes(mean.rm,mean.nc),size=1)+
  scale_y_continuous(breaks = seq(-20,60,10),labels = seq(-20,60,10),limits = c(-20,60))+
  labs(x='JSD to R_ase',y='Contamination proportion (%)');plot  
ggsave("result/JSD2RM vs contam.b.pdf",plot,height = 3.5,width = 4.5)
#Fig. 3B
plotdata <- JSD2 %>% filter(sample_type1=="OP")
summary <- plotdata %>% dplyr::group_by(methods3,sample_type1) %>% 
  dplyr::summarise(
    mean.rm=median(value),
    sd.rm.Q1=quantile(value,0.25),
    sd.rm.Q3=quantile(value,0.75),
    sd.rm=sd(value),
    mean.nc=median(contam),
    sd.nc.Q1=quantile(contam,0.25),
    sd.nc.Q3=quantile(contam,0.75),
    sd.nc=sd(contam),
  )
write.csv(summary,"result/med.ba.accuracy.csv")
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
                     xmax=mean.rm+sd.rm),height=1,size=0.1)+
  geom_point(aes(mean.rm,mean.nc),size=1)+
  scale_y_continuous(breaks = seq(-20,60,10),labels = seq(-20,60,10),limits = c(-20,60))+
  labs(x='JSD to R_ase',y='Contamination proportion (%)');plot  
ggsave("result/JSD2RM vs contam.o.pdf",plot,height = 3.5,width = 4.5)

# Fig3C,D ---------------------------------------------------------
#C
setwd("F:/Code")
plotdf <- read_excel("data.xlsx",sheet = "bias",
                     col_types = c("text",rep("numeric",31),rep("text",4))) %>% 
  as.data.frame() %>% filter(sample_type1=="BALF")
rownames(plotdf) <- plotdf$tax_name
plotdf.f <- plotdf[,c(2:8,16:22)]
plotdf.p <- plotdf[,c(9:15,23:29)]
colnames(plotdf.p) <- colnames(plotdf.f)
plotdf.p[plotdf.p<0.001] <- "***"
plotdf.p[plotdf.p<0.01&plotdf.p>0.001] <- "**"
plotdf.p[plotdf.p<0.05&plotdf.p>0.01] <- "*"
plotdf.p[plotdf.p>0.05] <- NA
plotdf.f <- plotdf.f %>% log(.,2) %>% as.matrix()
range(plotdf.f)
plotdf.f[plotdf.f>5] <- 5
plotdf.f[plotdf.f< -5] <- -5
col_fun = circlize::colorRamp2(c(-5, 0, 5), c("#006837", "white", "#AB0626"))
range(plotdf$abun.RM)
col_fun_prop = circlize::colorRamp2(c(-2.2,-1.8, -1.4,-1), 
                                    c("white","#ffe4da","#ff8f70","#ff0000"))
# c(-6,,0), c("white","#ffbdae","#fa7963","#E41A1C")
library("ComplexHeatmap")
lgd = Legend(col_fun = col_fun_prop, title = "Prop",
             break_dist = 1)

# "#E64B35FF","#4DBBD5FF","#8491B4FF","#00A087FF","#3C5488FF","#91D1C2FF","#F39B7FFF"


Heatmap(plotdf.f, "log2(fold.chage)",col = col_fun,cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_side = "left",column_names_side = "top",
        column_names_rot = 90,
        rect_gp = gpar(col = "white", lwd = 1),
        cell_fun = function(j, i, x, y, width, height, fill) 
        {if( !is.na(plotdf.p[i,j]))
          grid.text(sprintf("%s", plotdf.p[i, j]), x, y, gp = gpar(fontsize = 10))
        },column_split = rep(c("abs","rel"),each=7),
        row_split = c(rep("Actinobacteria",5),rep("Bacteroidetes",8),"Euryarchaeota",
                      rep("Firmicutes",14),
                      rep("Proteobacteria",5),"Tenericutes","Uroviricota"),
        row_title_rot = 0,border = T)+
  Heatmap(plotdf$phylum,name = "Phylum",column_names_side = "top",
          col=c("#E64B35FF","#4DBBD5FF","#8491B4FF","#00A087FF","#3C5488FF","#91D1C2FF","#F39B7FFF"),
          width = unit(5.5,"mm"))+
  Heatmap(plotdf$Gram,name = "Gram",column_names_side = "top",
          col=c("#6666cc","#ffff99"),width = unit(5.5,"mm"),
          rect_gp = gpar(col = "white", lwd = 0.7))+
  Heatmap(plotdf$O2_consuming,name = "O2",column_names_side = "top",
          col=c("#003f5c","#58508d","#bc5090","#ff6361"),width = unit(5.5,"mm"),
          rect_gp = gpar(col = "white", lwd = 0.7))+
  Heatmap(log10(plotdf$abun.RM), name = "abun.RM",col = col_fun_prop,
          width = unit(5.5, "mm"),
          #lgd = Legend(col_fun = col_fun_prop, title = "Prop", break_dist = 1),
          column_names_side = "top",
          rect_gp = gpar(col = "white", lwd = 0.7)) +
  Heatmap(log10(plotdf$abun.NC), name = "abun.NC", width = unit(5.5, "mm"),
          col = col_fun_prop,
          column_names_side = "top",
          rect_gp = gpar(col = "white", lwd = 0.7)) # + ha2
draw(lgd, x = unit(26, "cm"), y = unit(0.7, "cm"), just = c("left", "bottom"))

#D
plotdf <- read_excel("data.xlsx",sheet = "bias",
                     col_types = c("text",rep("numeric",31),rep("text",4))) %>% 
  as.data.frame() %>% filter(sample_type1=="OP")
rownames(plotdf) <- plotdf$tax_name
plotdf.f <- plotdf[,c(2:8,16:22)]
plotdf.p <- plotdf[,c(9:15,23:29)]
colnames(plotdf.p) <- colnames(plotdf.f)
plotdf.p[plotdf.p<0.001] <- "***"
plotdf.p[plotdf.p<0.01&plotdf.p>0.001] <- "**"
plotdf.p[plotdf.p<0.05&plotdf.p>0.01] <- "*"
plotdf.p[plotdf.p>0.05] <- NA
plotdf.f <- plotdf.f %>% log(.,2) %>% as.matrix()
range(plotdf.f)
plotdf.f[plotdf.f>5] <- 5
plotdf.f[plotdf.f< -5] <- -5
col_fun = circlize::colorRamp2(c(-5, 0, 5), c("#006837", "white", "#AB0626"))
range(plotdf$abun.RM)
col_fun_prop = circlize::colorRamp2(c(-2.2,-1.8, -1.4,-1), 
                                    c("white","#ffe4da","#ff8f70","#ff0000"))
# c(-6,,0), c("white","#ffbdae","#fa7963","#E41A1C")
lgd = Legend(col_fun = col_fun_prop, title = "Prop",
             break_dist = 1)
library("ComplexHeatmap")

# "#E64B35FF","#4DBBD5FF","#00A087FF","#7E6148FF","#3C5488FF","#F39B7FFF"


Heatmap(plotdf.f, "log2(fold.chage)",col = col_fun,cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_side = "left",column_names_side = "top",
        column_names_rot = 90,
        rect_gp = gpar(col = "white", lwd = 0.7),
        cell_fun = function(j, i, x, y, width, height, fill) 
        {if( !is.na(plotdf.p[i,j]))
          grid.text(sprintf("%s", plotdf.p[i, j]), x, y, gp = gpar(fontsize = 10))
        },column_split = rep(c("abs","rel"),each=7),
        row_split = c(rep("Actinobacteria",4),rep("Bacteroidetes",8),
                      rep("Firmicutes",16),"Fusobacteria",
                      rep("Proteobacteria",8),rep("Uroviricota",2)),
        row_title_rot = 0,border = T)+
  Heatmap(plotdf$phylum,name = "Phylum",column_names_side = "top",
          col=c("#E64B35FF","#4DBBD5FF","#00A087FF","#7E6148FF","#3C5488FF","#F39B7FFF"),
          width = unit(5.5,"mm"))+
  Heatmap(plotdf$Gram,name = "Gram",column_names_side = "top",
          col=c("#6666cc","#ffff99"),width = unit(5.5,"mm"),
          rect_gp = gpar(col = "white", lwd = 0.7))+
  Heatmap(plotdf$O2_consuming,name = "O2",column_names_side = "top",
          col=c("#003f5c","#58508d","#bc5090","#ff6361"),width = unit(5.5,"mm"),
          rect_gp = gpar(col = "white", lwd = 0.7))+
  Heatmap(log10(plotdf$abun.RM), name = "abun.RM",col = col_fun_prop,
          width = unit(5.5, "mm"),
          #lgd = Legend(col_fun = col_fun_prop, title = "Prop", break_dist = 1),
          column_names_side = "top",
          rect_gp = gpar(col = "white", lwd = 0.7)) +
  Heatmap(log10(plotdf$abun.NC), name = "abun.NC", width = unit(5.5, "mm"),
          col = col_fun_prop,
          column_names_side = "top",
          rect_gp = gpar(col = "white", lwd = 0.7)) # + ha2
draw(lgd, x = unit(26, "cm"), y = unit(0.7, "cm"), just = c("left", "bottom"))



# Figure 4 ---------------------------------------------------------------
setwd("F:/Code")
df <- read_excel("data.xlsx",sheet = "radar") %>% as.data.frame()
#balf Figure 4A
df1 <- df %>% filter(sample_type1=="BALF") %>% 
  spread(methods3,value) %>% 
  dplyr::select(-"sample_type1") %>% 
  column_to_rownames("index") %>% as.data.frame() %>% t() %>% 
  .[,c("microbial.prop","richness","bacterial.load","JSD2RM","contamination")] %>% 
  as.data.frame()

maxmin <- data.frame(
  microbial.prop=c(max(df1$microbial.prop),min(df1$microbial.prop)),
  richness=c(max(df1$richness),min(df1$richness)),
  bacterial.load=c(max(df1$bacterial.load),min(df1$bacterial.load)),
  JSD2RM=c(max(df1$JSD2RM),min(df1$JSD2RM)),
  contamination=c(max(df1$contamination),min(df1$contamination)))
dat <- rbind(maxmin,df1)
library(fmsb)
radarchart(dat, 
           # axistype=1, #设定axes的类型,1 means center axis label only
           seg=5, #设定网格的数目
           plty=1, #设定point连线的线型
           # plwd = 2,
           pcol=c("#E18727","#0072B5","#20854E","#7876B1","#696969","#BC3C29","#FFDC91"),
           vlcex=1, #设置标签的字体粗细大小
           cglcol = "#696969",
           cglwd = 0.6
)

#OP
df1 <- df %>% filter(sample_type1=="OP") %>% 
  spread(methods3,value) %>% 
  dplyr::select(-"sample_type1") %>% 
  column_to_rownames("index") %>% as.data.frame() %>% t() %>% 
  .[,c("microbial.prop","richness","bacterial.load","JSD2RM","contamination")] %>% 
  as.data.frame()

maxmin <- data.frame(
  microbial.prop=c(max(df1$microbial.prop),min(df1$microbial.prop)),
  richness=c(max(df1$richness),min(df1$richness)),
  bacterial.load=c(max(df1$bacterial.load),min(df1$bacterial.load)),
  JSD2RM=c(max(df1$JSD2RM),min(df1$JSD2RM)),
  contamination=c(max(df1$contamination),min(df1$contamination)))
dat <- rbind(maxmin,df1)
library(fmsb)
radarchart(dat, 
           # axistype=1, #设定axes的类型,1 means center axis label only
           seg=5, #设定网格的数目
           plty=1, #设定point连线的线型
           # plwd = 2,
           pcol=c("#E18727","#0072B5","#20854E","#7876B1","#696969","#BC3C29","#FFDC91"),
           vlcex=1, #设置标签的字体粗细大小
           cglcol = "#696969",
           cglwd = 0.6
)

# Figure 5  ---------------------------------------------------------------
#Figure 5A
setwd("F:/Code")
meta <- read_excel("data.xlsx",sheet = "meta.mock") %>% as.data.frame() %>% 
  .[order(.$order),]
species <- read_excel("data.xlsx",sheet = "species.mock") %>% as.data.frame() %>% 
  column_to_rownames("tax_name") %>% apply(2,function(x){x/sum(x)}) %>% 
  as.data.frame() %>% .[,meta$oldID] %>% filter(rowSums(.)>0)
tax_name <- read_excel("data.xlsx",sheet = "tax_name.mock") %>% as.data.frame()
sample <- species[tax_name$tax_name,] %>% cbind(tax_name) %>% 
  dplyr::select(-"tax_name")
tax.order <- rownames(sample)[order(sample$ID18,decreasing = F)]
sample["others",] <- 1-colSums(sample)
sample1 <- sample %>% as.data.frame() %>% rownames_to_column("tax_name") %>% 
  reshape2::melt(id.vars="tax_name",variable.name="oldID") %>% left_join(meta)

JSD_df <- JSD(t(species),unit = "log") %>% sqrt() %>% as.data.frame()  
row.names(JSD_df) <- colnames(JSD_df) <- meta$sampleID #添加行名
JSD_df1 <- JSD_df
JSD_df <- as.dist(JSD_df)
hc <- hclust(JSD_df,method = 'average')
plot(hc)

sample1$sampleID <- factor(sample1$sampleID,levels = meta$sampleID[hc$order])
sample1$tax_name <- factor(sample1$tax_name,
                           levels = c("others",tax.order))
plot <- ggplot(sample1,aes(sampleID,fill=tax_name,y=value*100))+
  geom_col(position='fill')+
  theme_bw() +
  scale_y_continuous(expand=c(0, 0),breaks = seq(0,1,0.25),
                     labels = seq(0,100,25))+
  labs(x='', y='Percentage (%)')+
  scale_fill_igv()+
  theme(axis.text.x = element_text(angle = 90,hjust = 0.5,vjust = 0.5,
                                   color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.title = element_text(size = rel(1.2), color = "black"),
        legend.background = element_blank(),
        legend.text = element_text(color = "black"),
        legend.box.spacing=unit(0,"cm"),
        legend.title = element_blank(),
        legend.box = "horizontal",
        strip.background = element_blank(),
        strip.text = element_text(color = "black"));plot
#Figure5B
meta <- read_excel("data.xlsx",sheet = "meta.mock") %>% as.data.frame()
species <- read_excel("data.xlsx",sheet = "species.mock") %>% as.data.frame() %>% 
  column_to_rownames("tax_name") %>% apply(2,function(x){x/sum(x)}) %>% 
  as.data.frame() %>% .[,meta$oldID] %>% filter(rowSums(.)>0) %>% .[,meta$oldID]
tax_name <- read_excel("data.xlsx",sheet = "tax_name.mock") %>% as.data.frame()
sample <- species %>% filter(rownames(.)%in%tax_name$tax_name) %>% 
  as.data.frame() %>% 
  rownames_to_column("tax_name") %>% 
  reshape2::melt(id.vars="tax_name",variable.name="oldID") %>% left_join(meta) %>% 
  mutate(abs=con_ng_ul*value*50)
fold.change <- data.frame()
for (i in unique(sample$tax_name)) {
  sample.i <- sample %>% filter(tax_name==i)
  for (j in c("R","P","W","Fil","S","M","H")) {
    sample.ij <- sample.i %>% filter(methods3==j)
    sample.ij.rm <- sample.i %>% filter(methods3=="RM")
    fold.ij <- data.frame(methods3=rep(j,9),
                          tax_name=rep(i,9),
                          fold=c(sample.ij$abs[1]/sample.ij.rm$abs,
                                 sample.ij$abs[2]/sample.ij.rm$abs,
                                 sample.ij$abs[3]/sample.ij.rm$abs))
    fold.change <- rbind(fold.change,fold.ij)
  }
}
sum <- fold.change %>% dplyr::group_by(tax_name,methods3) %>% 
  dplyr::summarise(med=median(fold),
                   Q1=quantile(fold,0.25),
                   Q3=quantile(fold,0.75))
sum$methods3 <- factor(sum$methods3,levels = c("R","P","W","Fil","S","M","H"))
order <- read_excel("data.xlsx",sheet = "bias.mock") 
sum$tax_name <- factor(sum$tax_name,levels = rev(order$tax_name))

plot <- ggplot(sum,aes(log2(med),tax_name,color=methods3))+
  geom_point()+theme_bw()+
  #cale_x_log10()+#geom_path(aes(group=methods3))+
  geom_vline(xintercept = 0)+
  theme(axis.text.x = element_text(angle = 90,hjust = 0.5,vjust = 0.5))+
  facet_wrap(~methods3,ncol = 7)+
  geom_segment(aes(y=tax_name,x=log2(Q1),yend=tax_name,xend=log2(Q3)),color="black")+
  scale_color_manual(values = c(R="black",RM="#696969",P="#7876B1",W="#FFDC91",
                                Fil="#E18727",S="#BC3C29",M="#20854E",H="#0072B5"));plot

annotate <- read_excel("data.xlsx",sheet = "annotate") %>% as.data.frame() %>% .[,1:2]
annotate$tax_name <- factor(annotate$tax_name,levels = rev(order$tax_name))
p2 <- ggplot(annotate, aes(x = 1, y = tax_name, fill = phylum)) +
  geom_tile() +
  # geom_text(data = label_pos,
  #           aes(x = midpoint, y = 1, label = phylum),
  #           fontface = "bold",
  #           size = 4) +
  theme_bw(base_size = 12) +
  scale_fill_npg() +
  scale_x_continuous(expand = c(0, 0)) +
  xlab("Feature");p2

p <- plot_grid(p2,plot, ncol = 2, rel_widths = c(0.3, 0.8), align = "h", axis = "lr")

annotate <- read_excel("data.xlsx",sheet = "annotate") %>% as.data.frame() %>% .[,c(1,3)]
annotate$tax_name <- factor(annotate$tax_name,levels = rev(order$tax_name))
annotate$Gram <- factor(annotate$Gram,levels = c("negative","positive"))
p2 <- ggplot(annotate, aes(x = 1, y = tax_name, fill = Gram)) +
  geom_tile() +
  # geom_text(data = label_pos,
  #           aes(x = midpoint, y = 1, label = phylum),
  #           fontface = "bold",
  #           size = 4) +
  theme_bw(base_size = 12) +
  scale_fill_manual(values = c("#6666cc","#ffff99")) +
  scale_x_continuous(expand = c(0, 0)) +
  xlab("Feature");p2

p <- plot_grid(p2,plot, ncol = 2, rel_widths = c(0.3, 0.8), align = "h", axis = "lr")
ggsave("result/gram.pdf",p,height = 4,width = 18)

annotate <- read_excel("data.xlsx",sheet = "annotate") %>% as.data.frame() %>% .[,c(1,4)]
annotate$tax_name <- factor(annotate$tax_name,levels = rev(order$tax_name))
annotate$O2_consuming <- factor(annotate$O2_consuming,levels = c("Aerobic","Anaerobic","Facultatively_Anaerobic"))
p2 <- ggplot(annotate, aes(x = 1, y = tax_name, fill = O2_consuming)) +
  geom_tile() +
  # geom_text(data = label_pos,
  #           aes(x = midpoint, y = 1, label = phylum),
  #           fontface = "bold",
  #           size = 4) +
  theme_bw(base_size = 12) +
  scale_fill_manual(values = c("#003f5c","#58508d","#bc5090","#ff6361")) +
  scale_x_continuous(expand = c(0, 0)) +
  xlab("Feature");p2

p <- plot_grid(p2,plot, ncol = 2, rel_widths = c(0.3, 0.8), align = "h", axis = "lr")

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

# Figure S9A Bacterial retention rate ----------------------------------------------------------
setwd("F:/Code")
data <- read_excel("data.xlsx",sheet = "d1") %>% as.data.frame()
data$Bacterial_retention_rate[data$Bacterial_retention_rate>1] <- 1
data$Methods <- factor(data$Methods, levels = c("Raw","R_ase","F_ase","S_ase","MEM"))

stat.test <- compare_means(Bacterial_retention_rate~Methods, data, 
                           method = "t.test",
                           p.adjust.method = "fdr",
                           paired = TRUE) 
library(multcompView)
anova <- aov(Bacterial_retention_rate~Methods, data = data)
tukey <- TukeyHSD(anova)
tukey[["Methods"]][,4] <- stat.test$p.adj
cld <- multcompLetters4(anova, tukey)
aa=cld$Methods$Letters %>% as.data.frame() %>% rownames_to_column("Methods")
colnames(aa)[2] <- "sign"
#作图
p1 <- ggplot(data)+ 
  geom_boxplot(aes(x = Methods, y =Bacterial_retention_rate*100 , fill = Methods),width = 0.5)+  
  geom_point(aes(x = Methods, y = Bacterial_retention_rate*100, fill = Methods))+
  labs(x = "", y = "Bacterial_retention_rate")+  
  theme_bw()+ 
  scale_fill_manual(values = c(Raw = "black", R_ase = "#696969",F_ase="#E18727",
                               S_ase="#BC3C29",MEM = "#4f3"))+
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  annotation_logticks(sides = "l",outside = TRUE)+ylim(0,110)+
  geom_text(aes(label = sign, y = 105,x=Methods), vjust = -0.4,size=3.5,
            data = aa,show.legend = F);p1
ggsave("result/bac.retention.rate.pdf",p1,height = 3,width = 3.5)

# Figure S9B Microbial proportion ----------------------------------------------------------
setwd("F:/Code")
data <- read_excel("data.xlsx",sheet = "d2") %>% as.data.frame() %>% 
  mutate(ABFV_taxa.prop=abfv/(total_reads-low_quality)*taxa_prop)
data$Methods <- factor(data$Methods, levels = c("Raw","R_ase","F_ase","S_ase","MEM"))
stat.test <- compare_means(ABFV_taxa.prop~Methods, data, 
                           method = "wilcox.test",
                           p.adjust.method = "fdr",
                           paired = TRUE) 
library(multcompView)
anova <- aov(ABFV_taxa.prop~Methods, data = data)
tukey <- TukeyHSD(anova)
tukey[["Methods"]][,4] <- stat.test$p.adj
cld <- multcompLetters4(anova, tukey)
aa=cld$Methods$Letters %>% as.data.frame() %>% rownames_to_column("Methods")
colnames(aa)[2] <- "sign"
#作图
p1 <- ggplot(data)+ 
  geom_boxplot(aes(x = Methods, y =ABFV_taxa.prop*100 , fill = Methods),width = 0.5)+  
  geom_point(aes(x = Methods, y = ABFV_taxa.prop*100, fill = Methods))+
  labs(x = "", y = "Proportion of microbial reads(%)")+  
  theme_bw()+ 
  scale_fill_manual(values = c(Raw = "black", R_ase = "#696969",F_ase="#E18727",
                               S_ase="#BC3C29",MEM = "#4f3"))+
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  coord_cartesian(clip = "off")+
  annotation_logticks(sides = "l",outside = TRUE)+
  scale_y_log10(limits=c(min(data$ABFV_taxa.prop)*100,max(data$ABFV_taxa.prop)*100),
                labels=trans_format("log10",math_format(10^.x)),
                breaks=c(0.00001,0.0001,0.001,0.01,0.1,1,10,100))+
  geom_text(aes(label = sign, y = max(data$ABFV_taxa.prop)*100,x=Methods), vjust = -0.4,size=3.5,
            data = aa,show.legend = F);p1
ggsave("result/microbial.prop.pdf",p1,height = 3,width = 3.5)

# Figure S9C species richness ----------------------------------------------------------
setwd("F:/code")
de.tbl1 <- read_excel("data.xlsx",sheet = "reads") %>% 
  as.data.frame() %>% column_to_rownames("tax_name") #taxonomy table after removing contaminants
qc <- read_excel("data.xlsx",sheet = "d2") %>% 
  as.data.frame() %>% .[,c("ID","total_reads","low_quality","reads.decontam")] %>% 
  mutate(bac.prop=reads.decontam/(total_reads-low_quality),
         total=total_reads-low_quality)
meta <- read_excel("data.xlsx",sheet = "d2") %>% 
  as.data.frame() %>% 
  .[,c("ID","Methods","sampleNo1")]
qc2 <- left_join(qc,meta)
sum <- qc2 %>% dplyr::group_by(sampleNo1) %>% 
  dplyr::summarise(rare=min(total))
qc2 <- left_join(qc2,sum)
rare <- data.frame(tax_name=rownames(de.tbl1))
rich <- data.frame()
# bac.raw <- qc2 %>% filter(methods3=="R")
for (i in qc2$ID) {
  qc1.i <- filter(qc2,ID==i)
  df.i <- data.frame(tax_name=rownames(de.tbl1),
                     sample=round(de.tbl1[,i])) %>% 
    column_to_rownames("tax_name")
  if (nrow(df.i)>0) {
    set.seed(100)
    rare.i <- rrarefy(df.i,qc1.i$rare) %>% t() %>% 
      as.data.frame() %>% 
      rownames_to_column("tax_name")
    
    rich.i <- data.frame(ID=i,
                         richniss=length(rare.i$sample[rare.i$sample>0]),
                         rare.theory=qc1.i$rare,
                         rare.real=sum(rare.i$sample))
    rich <- rbind(rich,rich.i)
    
    colnames(rare.i)[2] <- i
    rare <- left_join(rare,rare.i) 
  }
}
rich1 <- rich
write.csv(rich1,"rich1.csv")

data <- read_excel("data.xlsx",sheet = "d3") %>% as.data.frame()
data$Methods <- factor(data$Methods, levels = c("Raw","R_ase","F_ase","S_ase","MEM"))
stat.test <- compare_means(num_s~Methods, data, 
                           method = "wilcox.test",
                           p.adjust.method = "fdr",
                           paired = TRUE) 
library(multcompView)
anova <- aov(num_s~Methods, data = data)
tukey <- TukeyHSD(anova)
tukey[["Methods"]][,4] <- stat.test$p.adj
cld <- multcompLetters4(anova, tukey)
aa=cld$Methods$Letters %>% as.data.frame() %>% rownames_to_column("Methods")
colnames(aa)[2] <- "sign"
#作图
p1 <- ggplot(data)+ 
  geom_boxplot(aes(x = Methods, y =num_s , fill = Methods),width = 0.5)+  
  geom_point(aes(x = Methods, y = num_s, fill = Methods))+
  labs(x = "", y = "Species richness")+  
  theme_bw()+ 
  scale_fill_manual(values = c(Raw = "black", R_ase = "#696969",F_ase="#E18727",
                               S_ase="#BC3C29",MEM = "#4f3"))+
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  annotation_logticks(sides = "l",outside = TRUE)+
  scale_y_continuous(limits=c(min(data$num_s),max(data$num_s)))+
  geom_text(aes(label = sign, y = max(data$num_s),x=Methods), vjust = -0.4,size=3.5,
            data = aa,show.legend = F);p1
ggsave("result/species.richness.pdf",p1,height = 3,width = 3.5)

# Figure S9D accuracy ----------------------------------------------------------------
setwd("F:/Code")
data <- read_excel("data.xlsx",sheet = "d4") %>% as.data.frame()
data$Methods <- factor(data$Methods, levels = c("Raw","R_ase","F_ase","S_ase","MEM"))
sum <- data %>% dplyr::group_by(Methods) %>% dplyr::summarise(mean.nc=median(contam_level)*100,
                                                              sd.nc.Q1=quantile(contam_level*100,0.25),
                                                              sd.nc.Q3=quantile(contam_level*100,0.75),
                                                              mean.rm=median(JSD),
                                                              sd.rm.Q1=quantile(JSD,0.25),
                                                              sd.rm.Q3=quantile(JSD,0.75),
                                                              sd.rm=sd(JSD),
                                                              sd.nc=sd(contam_level*100))
plotdata <- left_join(data,sum,by="Methods") 
plot <- ggplot(plotdata,aes(JSD,contam_level*100,color=Methods))+
  geom_point(size=0.2)+theme_bw()+
  scale_color_manual(values = c(Raw = "black", R_ase = "#696969",F_ase="#E18727",
                                S_ase="#BC3C29",MEM = "#4f3"))+
  geom_errorbar(aes(x=mean.rm,
                    ymin=mean.nc-sd.nc,
                    ymax=mean.nc+sd.nc),
                width=0.02,size=0.1)+
  geom_errorbarh(aes(y=mean.nc,
                     xmin=mean.rm-sd.rm,
                     xmax=mean.rm+sd.rm),height=1,size=0.1)+
  geom_point(aes(mean.rm,mean.nc),size=1)+
  scale_y_continuous(limits = c(-25,50))+
  labs(x='JSD to R_ase',y='Contamination proportion (%)');plot 
ggsave("result/accuracy.pdf",plot,height = 3,width = 3.9)

# Figure S9E hetamap bias ------------------------------------------------------------
setwd("F:/Code")
plotdf <- read_excel("data.xlsx",sheet = "d5",col_types = c("text",rep("numeric",16)))
plotdf.f <- plotdf[,grepl(".f.|tax_name",colnames(plotdf))] %>% as.data.frame() %>% column_to_rownames("tax_name")
plotdf.p <- plotdf[,grepl(".p.|tax_name",colnames(plotdf))] %>% as.data.frame() %>% column_to_rownames("tax_name")
colnames(plotdf.p) <- colnames(plotdf.f)
plotdf.p[plotdf.p<0.001] <- "***"
plotdf.p[plotdf.p<0.01&plotdf.p>0.001] <- "**"
plotdf.p[plotdf.p<0.05&plotdf.p>0.01] <- "*"
plotdf.p[plotdf.p>0.05] <- NA
plotdf.f <- plotdf.f %>% log(.,2) %>% as.matrix()
plotdf.f[plotdf.f>5] <- 5
plotdf.f[plotdf.f< -5] <- -5
col_fun = circlize::colorRamp2(c(-5, 0, 5), c("#006837", "white", "#AB0626"))
col_fun_prop = circlize::colorRamp2(c(-3.5,-2.7, -1.9,-1), 
                                    c("white","#ffe4da","#ff8f70","#ff0000"))
lgd = Legend(col_fun = col_fun_prop, title = "Prop",
             break_dist = 1)
colnames(plotdf.f) <- c("Raw","F_ase","S_ase","MEM","Raw","F_ase","S_ase","MEM")
library(ComplexHeatmap)
Heatmap(plotdf.f, "log2(fold.chage)",col = col_fun,cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_side = "left",column_names_side = "top",
        column_names_rot = 90,
        rect_gp = gpar(col = "white", lwd = 0.7),
        cell_fun = function(j, i, x, y, width, height, fill) 
        {if( !is.na(plotdf.p[i,j]))
          grid.text(sprintf("%s", plotdf.p[i, j]), x, y, gp = gpar(fontsize = 10))
        },
        column_split = rep(c("Relative\nabundance","Absolute\nabundance"),each=4),
        row_title_rot = 0,border = T,row_names_gp = gpar(fontsize = 10, fontface = 3))


# Figure S9F radar -------------------------------------------------------------------
setwd("F:/Code")
data <- read_excel("data.xlsx",sheet = "d7") %>% as.data.frame() %>% filter(Methods!="Raw")
df1 <- data %>% mutate(decontam.prop=1-contam_level,distance=1-JSD) %>% 
  .[,c("ID","Bacterial_retention_rate", "ABFV_taxa.prop","num_s","decontam.prop","distance","Methods")]
median_df <- df1 %>% dplyr::group_by(Methods) %>%  
  dplyr::summarise(ABFV_taxa.prop=median(ABFV_taxa.prop),num_s=median(num_s),
                   Bacterial_retention_rate=median(Bacterial_retention_rate),
                   distance=median(distance),
                   decontam.prop=median(decontam.prop)) %>% 
  column_to_rownames(var="Methods")
maxmin <- data.frame(
  `Microbial proportion`=c(max(median_df$ABFV_taxa.prop),min(median_df$ABFV_taxa.prop)),
  `Species richness`=c(max(median_df$num_s),min(median_df$num_s)),
  `Bacterial retention rate`=c(max(median_df$Bacterial_retention_rate),min(median_df$Bacterial_retention_rate)),
  Accuracy=c(max(median_df$distance),min((median_df$distance))),
  `Decontamination level`=c(max(median_df$decontam.prop),min(median_df$decontam.prop))
)
names(median_df) <- names(maxmin)
dat <- rbind(maxmin, median_df )
library(fmsb)
radarchart(dat, 
           seg=5, #设定网格的数目
           plty=1, #设定point连线的线型
           plwd = 2,
           pcol=c("#E18727","#4f3","#696969","#BC3C29"),#"#FFDC91"),
           vlcex=0.8, #设置标签的字体粗细大小
           cglcol = "#696969",
           cglwd = 0.6, 
)

# decontam ----------------------------------------------------------------
setwd("E:/富集/MEM/MEM/data")
meta <- read_excel("abfv_s_rel.xlsx",sheet = "meta1") %>% as.data.frame() 
d1 <- read_excel("abfv_s_rel.xlsx",sheet = "species") %>% as.data.frame() %>% 
  column_to_rownames("tax_name") %>% .[,meta$seq]
d1[is.na(d1)] <- 0
d1 <- d1 %>% apply(2,function(x){x/sum(x)}) %>% as.data.frame() %>% filter(rowSums(.)>0)
colnames(d1) <- meta$ID
NC <- d1[,meta$ID[meta$Methods=="NC"]] #%>% filter(rowSums(.)>0)

meta <- read_excel("abfv_s_rel.xlsx",sheet = "meta1") %>% as.data.frame() %>% filter(Methods=="NC")
rel.abun <- d1[,meta$ID] #%>% filter(rowSums(.)>0)

df <- rel.abun %>% as.matrix()
JSD_df <- JSD(t(df),unit = "log") %>% sqrt() %>% as.data.frame()  
colnames(JSD_df) <- rownames(JSD_df) <- colnames(df)
JSD_df <- as.dist(JSD_df)

set.seed(123)  # 设置随机种子，以确保结果可重复
k_values <- 1:10  # 尝试的类别数量范围
within_cluster_variability <- sapply(k_values, function(k) {
  kmeans_result <- kmeans(JSD_df, centers = k)
  sum(kmeans_result$withinss)  # 计算类内变异性的总和
})
plot(k_values, within_cluster_variability, type = "b", xlab = "Number of clusters", 
     ylab = "Within-cluster variability")
hc <- hclust(JSD_df,method = "ward.D")
plot(hc,hang=-1, ylab = "Height") 
rect.hclust(hc, k = 2)
cluster <- cutree(hc,k=2) %>% as.data.frame() %>% rownames_to_column("seqID")
colnames(cluster)[2] <- "cluster"
tax.ord <- colnames(df)[hc$order]
NC1 <- rel.abun[,grepl("-N|-F|S",colnames(rel.abun))] %>% mutate(NC1=rowMeans(.)) %>% rownames_to_column("tax_name")
NC2 <- rel.abun[,!grepl("-N|-F|S",colnames(rel.abun))] %>% mutate(NC2=rowMeans(.)) %>% rownames_to_column("tax_name")
meta1 <- read_excel("abfv_s_rel.xlsx",sheet = "meta1") %>% as.data.frame() %>% filter(Methods2!="MEM")
sample1 <- d1[,meta1$ID] %>% rownames_to_column("tax_name") %>% reshape2::melt() %>% left_join(NC1[,c("tax_name","NC1")]) %>% 
  mutate(fold=value/NC1) %>% filter(fold>5)
sum1 <- sample1 %>% dplyr::group_by(variable) %>% dplyr::summarise(taxa=sum(value))
meta1 <- read_excel("abfv_s_rel.xlsx",sheet = "meta1") %>% as.data.frame() %>% filter(Methods2=="MEM")
sample2 <- d1[,meta1$ID] %>% rownames_to_column("tax_name") %>% reshape2::melt() %>% left_join(NC2[,c("tax_name","NC2")]) %>% 
  mutate(fold=value/NC2) %>% filter(fold>5)
sum2 <- sample2 %>% dplyr::group_by(variable) %>% dplyr::summarise(taxa=sum(value))
meta <- read_excel("abfv_s_rel.xlsx",sheet = "meta1") %>% as.data.frame() 
colnames(meta)[2] <- "variable"
sum <- rbind(sum1,sum2) %>% left_join(meta) %>% filter(Methods!="NC")
write.csv(sum,"taxa.csv")
sum <- rbind(sum1,sum2) %>% left_join(meta) %>% filter(Methods!="NC") %>% dplyr::group_by(Methods) %>% 
  dplyr::summarise(med.taxa=median(taxa))
sample <- rbind(sample1[,c("tax_name","variable","value")],sample2[,c("tax_name","variable","value")]) %>% unique()
reads1 <- read_excel("abfv_s_rel.xlsx",sheet = "species") %>% as.data.frame() %>% 
  .[,c("tax_name",meta$seq)]
colnames(reads1) <- c("tax_name",meta$variable)
reads1 <- reads1 %>% reshape2::melt() %>% filter(!is.na(value))
colnames(reads1)[3] <- "reads"
reads2 <- left_join(sample,reads1) %>% .[,c("tax_name","variable","reads")] %>% spread(variable,reads)
write.csv(reads2,"reads.decontam.csv")
sum.reads <- left_join(sample,reads1) %>% dplyr::group_by(variable) %>% dplyr::summarise(reads=sum(reads))
write.csv(sum.reads,"sum.reads.csv")
# meta1 <- read_excel("abfv_s_rel.xlsx",sheet = "meta1") %>% as.data.frame() %>% filter(Methods%in%c("R_ase","F_ase","S_ase"))
# sample1 <- d1[,meta1$ID] %>% rownames_to_column("tax_name") %>% reshape2::melt() %>% left_join(NC1[,c("tax_name","NC1")]) %>% 
#   mutate(fold=value/NC1) %>% filter(fold>5)
# sum1 <- sample1 %>% dplyr::group_by(variable) %>% dplyr::summarise(taxa=sum(value))
# 
# meta1 <- read_excel("abfv_s_rel.xlsx",sheet = "meta1") %>% as.data.frame() %>% filter(!Methods%in%c("R_ase","F_ase","S_ase"))
# sample1 <- d1[,meta1$ID] %>% rownames_to_column("tax_name") %>% reshape2::melt() %>% left_join(NC2[,c("tax_name","NC2")]) %>% 
#   mutate(fold=value/NC2) %>% filter(fold>5)
# sum2 <- sample1 %>% dplyr::group_by(variable) %>% dplyr::summarise(taxa=sum(value))
# 
# meta <- read_excel("abfv_s_rel.xlsx",sheet = "meta1") %>% as.data.frame() 
# colnames(meta)[2] <- "variable"
# sum <- rbind(sum1,sum2) %>% left_join(meta) %>% filter(Methods!="NC") %>% dplyr::group_by(Methods) %>% 
#   dplyr::summarise(med.taxa=median(taxa))
meta <- read_excel("abfv_s_rel.xlsx",sheet = "meta1") %>% as.data.frame() %>% .[,c("ID","Methods","sampleNo1")]
sample <- rbind(sample1[,c("tax_name","variable","value")],sample2[,c("tax_name","variable","value")]) %>% 
  spread(variable,value) %>% column_to_rownames("tax_name")
sample[is.na(sample)] <- 0
sample <- apply(sample, 2, function(x){x/sum(x)})
df <- sample %>% as.matrix()
JSD_df <- JSD(t(df),unit = "log") %>% sqrt() %>% as.data.frame()  
colnames(JSD_df) <- rownames(JSD_df) <- colnames(df)
JSD_df <- JSD_df %>% rownames_to_column("ID")
JSD_df[upper.tri(JSD_df)] <- NA
JSD_df1 <- JSD_df %>% reshape2::melt(variable.name="ID.2") %>% filter(!is.na(value)) %>% left_join(meta)
colnames(JSD_df1) <- c("ID.1","ID.2","value","Methods.1","sampleNo1.1")
colnames(meta)[1] <- "ID.2"
JSD_df2 <- left_join(JSD_df1,meta) %>% filter(sampleNo1.1==sampleNo1) %>% filter(Methods.1=="R_ase"|Methods=="R_ase")
JSD_df2$Method <- JSD_df2$Methods
JSD_df2$Method[JSD_df2$Method=="R_ase"] <- JSD_df2$Methods.1[JSD_df2$Method=="R_ase"]
sum <- JSD_df2 %>% dplyr::group_by(Method) %>% dplyr::summarise(med=median(value))
write.csv(JSD_df2,"JSD_df2.csv")


=======
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
>>>>>>> 02ea453d2363a40f79f06b71cd4ace6c9af28af3
