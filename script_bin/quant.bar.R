#!/usr/bin/env Rscript

# 画柱形图
library(ggplot2)

data=read.table("Bin_all/Bin_quant/bin_abundance_table.txt", header=T, sep="\t")
data$sum=apply(data[,2:length(data[1,])], 1, FUN=sum)  # 求和，加到列尾
data=data[, c("Genomic.bins", "sum")]  # 取sum信息
data=data[order(data$sum, decreasing=T),]  # 倒序

write.table(data, file="Bin_all/Bin_quant/bin_abundance_sum.txt", quote=F, sep="\t", row.names=F)

data$color = paste("color", 1:length(data[,1]), sep="_")

result = ggplot(data, aes(x=data[,1], y=data[,2], fill=color)) +
geom_bar(stat="identity", width=0.5) +
# 柱图
labs(x="", y="Bin abundance") +
scale_y_continuous(expand=c(0, 0))+
# 调整y轴属性
scale_x_discrete(limits=factor(data[,1])) +
# 设第一列为因子，不排序
theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1)) +
# angle：调整横轴标签倾斜角度
theme(panel.grid=element_blank(), panel.background=element_rect(color='black', fill='transparent'))
# 调整背景
wide = round(length(data[,1])/40)*7

ggsave(result, filename="Bin_all/Bin_quant/bin_abundance_bar.pdf", width=wide)
