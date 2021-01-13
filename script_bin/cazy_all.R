#!/usr/bin/env Rscript

setwd("Bin_all/Bin_cazy/")
files=list.files(pattern="bin.*.sum")  # 读取所有文件名

Bin_ID=vector()
for(i in 1:length(files))
{
    Bin_ID[i]=as.character(strsplit(files[i], split=".sum"))
    # 提取所有文件名
}

ml=list()
for(i in 1:length(files))
{
    ml[[i]]=read.table(files[i], sep='\t', na.string="", stringsAsFactors=F, header=T, quote="", comment.char="")
    # 读取所有数据框
}

CAZy_ID=vector()
for(i in 1:length(files))
{
    CAZy_ID=c(CAZy_ID, ml[[i]][, "entries"])
    # 合并所有cazy id
}
CAZy_ID=unique(CAZy_ID)  # 去重

origin=data.frame(CAZy_ID)  # 新建数据框

part=list()
for(i in 1:length(files))
{
    part[[i]]=data.frame(ml[[i]][,"entries"], ml[[i]][,"count"])
    colnames(part[[i]])=c("CAZy_ID", Bin_ID[i])
    origin=merge(origin, part[[i]], by="CAZy_ID", all=T)
    # 合并count列，并重命名表头
}
origin[is.na(origin)]=0  # 0替换NA

origin$sum=rowSums(origin[, 2:length(origin[1,])])  # 求和

result=origin[order(origin$sum, decreasing=T),]  # 排序

write.table(result, file="cazy_all.txt", sep="\t", row.names=F, quote=F) 

