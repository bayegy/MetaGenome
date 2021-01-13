#!/usr/bin/env Rscript

setwd("Bin_all/Bin_summary")

checkm=read.table("../Bin_quality/checkm_pick.txt", header=T, sep="\t")
gc=read.table("bin.gc_out.txt", header=T, sep="\t")
n50n90=read.table("N50.N90.out.txt", header=T, sep="\t")

merge=merge(merge(checkm, n50n90, by="BinID"), gc, by="BinID")

write.table(merge, file="bin_summary.txt", row.names=F, quote=F, sep="\t")
