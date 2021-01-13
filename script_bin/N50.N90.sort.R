#!/usr/bin/env Rscript

setwd("Bin_all/Bin_summary/")

data=read.table("N50.N90.txt", header=F)

BinID=vector()
N50=vector()
N90=vector()

index=seq(from=1, to=length(data[, 1]), by=3)
a=1
for(i in index)
{
    BinID[a]=as.character(data[i, "V1"])
    N50[a]=as.character(data[i+1, "V1"])
    N90[a]=as.character(data[i+2, "V1"])
    a=a+1
}

N50_N90_out=data.frame(BinID, N50, N90)
write.table(N50_N90_out, file="N50.N90.out.txt", sep="\t", quote=F, row.names=F)
