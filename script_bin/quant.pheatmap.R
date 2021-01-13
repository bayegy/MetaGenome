#!/usr/bin/env Rscript

library("pheatmap")
library(ggplot2)

data=read.table("Bin_all/Bin_quant/bin_abundance_table.txt", header=T, row.names=1, sep="\t")

pheatmap(data, filename="Bin_all/Bin_quant/bin_abundance_pheatmap.pdf", scale="column", cellwidth=20, cellheight=20, fontsize_number=18)
