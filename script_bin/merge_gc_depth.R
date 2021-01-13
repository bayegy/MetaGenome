#!/usr/bin/env Rscript
args = commandArgs(T)

route_file = unlist(strsplit(args[1], "/"))
route = paste(route_file[1:(length(route_file)-1)], collapse="/")
setwd(route)
file_name = route_file[length(route_file)]

base = as.character(strsplit(file_name, split=".gc.txt"))

data1 = read.table(file_name, sep="\t", header=T)
data2 = read.table(args[2], sep="\t", header=T)

result = merge(data1, data2, by.x="id", by.y="contigName", all.x=T)

result$bin_id = rep(base, length(data1[,1]))

write.table(result, file=paste(base, "gc.depth.txt", sep="."), sep="\t", quote=F, row.names=F)

