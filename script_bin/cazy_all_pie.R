#!/usr/bin/env Rscript

# 画饼图（一个cazy一个饼，top 20 all）
library(RColorBrewer)
setwd("Bin_all/Bin_cazy/")

data=read.table("cazy_all.txt", header=T, sep="\t", row.names=1)
data=data.frame(t(data[1:length(data[,1]), 1:(length(data[1,])-1)]))
data$bin=rownames(data)

id = colnames(data)

dir.create("cazy_all_pie")

for(i in 1:(length(data[1,])-1))
{
    name = paste("cazy_all_pie", id[i], sep = "/")
    input = data.frame(data[, i])
    input$id = data$bin
    rownames(input) = data$bin
    colnames(input) = colnames(data)[i]
    input = input[order(input[, 1], decreasing = T), ]
    
    pdf(paste(name, "pdf", sep = "."))
    #percent <- round(data[,i]/sum(data[,i])*100, 1)
    #label <- paste(data$bin, "(", percent, "% )")
    label <- paste(input[,2][1:10], "(", input[,1][1:10], ")")
    pie(input[,1][1:10], main=colnames(data)[i], label=label, border="white", col=brewer.pal(12, "Set3"))
    dev.off()
}
