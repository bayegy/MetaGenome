#!/usr/bin/env Rscript

library(ggplot2) # 加载ggplot2
library(ggtree) # 加载ggtree

tree=read.tree(list.files()[grepl("denovo.tree.nwk", list.files())]) # 读取nwk文件
map=read.table("map.txt", header=T, sep="\t", comment.char="")
data=fortify(tree)

tregraph=ggtree(tree, layout="circular", ladderize=FALSE, size=0.8, branch.length="none", aes(col=Phylum)) %<+% map +
# 树体：树文件、树形、粗细、颜色
geom_tiplab2(hjust=-0.3) +
# 枝名：大小、颜色、高度
geom_tippoint(size=3)+
# 端点颜色、大小
theme(legend.title=element_text(face="bold"), legend.position="bottom", legend.box="horizontal", legend.text=element_text(size=rel(0.8)))
# 图例位置、文字大小

ggsave(tregraph, file="tree_bins.pdf", width=9, height=9)
