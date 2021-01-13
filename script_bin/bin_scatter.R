#!/usr/bin/env Rscript
args = commandArgs(T)

route_file = unlist(strsplit(args[1], "/"))
route = paste(route_file[1:(length(route_file)-1)], collapse="/")
setwd(route)
file_name = route_file[length(route_file)]

data=read.table(file_name, header=T, sep="\t")

library(ggplot2)

result = ggplot(data, aes(x=GC_percent, y=Depth, color=Bin)) +
geom_point(size = 0.5) +
labs(x="GC Percent", y="Average Depth") +
theme(legend.title=element_text(face="bold"), panel.grid=element_blank(), panel.background=element_rect(color='black', fill='transparent'))  # 去背景，加黑框

ggsave(result, filename=args[2])
