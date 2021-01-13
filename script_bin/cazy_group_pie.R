#!/usr/bin/env Rscript

# 画饼图（一个bin一个饼，group）

library(RColorBrewer)
setwd("Bin_all/Bin_cazy/")

data=read.table("cazy_group.txt", header=T, sep="\t")

dir.create("cazy_group_pie")

id = colnames(data)

for(i in 2:(length(data[1,])))
{
    name = paste("cazy_group_pie", id[i], sep = "/")
    pdf(paste(name, "pdf", sep="."))
    percent <- round(data[,i]/sum(data[,i])*100, 1)
    label <- paste(data$CAZy_ID, "(", percent, "% )")
    pie(data[,i], main=colnames(data)[i], label=label, border="white", col=brewer.pal(6, "Set3"))
    dev.off()
}

