data=read.table("Bin_all/Bin_blobology/final.contigs.binned.blobplot", header=T, sep="\t")
colnames(data)[4]="abundance"  # 重命名一列
colnames

data$Bin=as.character(strsplit(as.character(data[, 'bin']), split=".fa"))  # 字符处理旧列，建新列
data=data[, c("gc", "abundance", "binned_phylum", "Bin")]
colnames(data)[3]="Phylum"

library(ggplot2)
pdf("Bin_all/Bin_blobology/bin_contig_bin.pdf")

ggplot(data, aes(x=gc, y=abundance, colour=Bin)) + geom_point() +
labs(x="GC percent", y="Contig abundance") +
theme(panel.grid=element_blank(), panel.background=element_rect(color='black', fill='transparent'))
# 去背景，加黑框

dev.off()

pdf("Bin_all/Bin_blobology/bin_contig_phylum.pdf")

ggplot(data, aes(x=gc, y=abundance, colour=Phylum)) + geom_point() +
labs(x="GC percent", y="Contig abundance") +
theme(panel.grid=element_blank(), panel.background=element_rect(color='black', fill='transparent'))
# 去背景，加黑框

dev.off()

