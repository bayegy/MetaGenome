#!/usr/bin/env Rscript

# 画热图
library(pheatmap)
setwd("Bin_all/Bin_cazy/")

data=read.table("cazy_group.txt", header=T, row.names=1, sep="\t")  # 读取

h=length(data[, 1])%/%10
w=length(data[1,])%/%10
Height=ifelse(length(data[, 1])>10, 7+h, 7)
Width=ifelse(length(data[1,])>10, 7+w, 7)

pdf("cazy_group_pheatmap.pdf", heigh=Height, width=Width)
pheatmap(data, scale="column")
dev.off()

# 画箱图
library(ggplot2)
library(reshape)

pdf("cazy_group_boxplot.pdf")
ggplot(melt(t(data)), aes(x=X2, y=value, fill=X2))+
# 添加数据、xy值、 颜色
geom_boxplot()+
# 盒图
labs(x="", y="CAZy number")+
# 标题：xy轴标签
theme(plot.title=element_text(hjust=0.5), legend.title=element_blank())+
# 图例：
theme(panel.grid=element_blank(), panel.background=element_rect(color='black', fill='transparent'))
# 背景
dev.off()

