#!/usr/bin/env Rscript

data=read.table("bin_taxonomy.txt", header=T, sep="\t")
id=data[,"id"]
Phylum=rep("Bin", length(data[,1]))
lab=data[,"id"]
data2=data.frame(id, tax, lab)
write.table(data2, file="map_part_bin.txt", quote=F, row.names=F, sep="\t")

