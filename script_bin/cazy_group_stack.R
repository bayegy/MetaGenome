#!/usr/bin/env Rscript

# 画堆叠图
library(ggplot2)
library(reshape)
setwd("Bin_all/Bin_cazy/")

data=read.table("cazy_group.txt", header=T, sep="\t", row.names=1)
data=data.frame(t(data))
data$sum=rowSums(data[, 1:length(data[1,])])  # 求和
data2=data[order(data$sum, decreasing=T),]  # 排序
data2=data2[,-length(data2[1,])]
data2=data.frame(t(data2))
data2$CAZyme=rownames(data2)

data_stack=melt(data2, id='CAZyme')

wide = round(length(data[,1])/40)*7

result = ggplot(data_stack, aes(variable, fill=CAZyme, value))+
geom_col(position='stack')+
# stack：堆叠图
labs(x="", y='CAZyme number') + labs(fill="CAZyme")+
# 给xy轴取名；给图例取名
scale_y_continuous(expand=c(0, 0))+
# 调整y轴属性
theme(axis.text.x=element_text(angle=45, hjust=1))+
# angle：调整横轴标签倾斜角度
# hjust：上下移动横轴标签
theme(panel.grid=element_blank(), panel.background=element_rect(color='black', fill='transparent'))
# 调整背景

ggsave(result, filename="cazy_group_stack.pdf", width=wide)
