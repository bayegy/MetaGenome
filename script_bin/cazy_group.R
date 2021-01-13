#!/usr/bin/env Rscript

setwd("Bin_all/Bin_cazy/")
files=list.files(pattern="bin.*.clean")  # 读取所有文件名
print(files)
Bin_ID=vector()
for(i in 1:length(files))
{
    Bin_ID[i]=as.character(strsplit(files[i], split=".clean"))
    # 提取所有文件名
}

ml=list()
for(i in 1:length(files))
{
    print(files[i])
    ml[[i]]=read.table(files[i], sep='\t', na.string="", stringsAsFactors=F, header=T, quote="", comment.char="")
    # 读取所有数据框
}

origin=data.frame(CAZy_ID=c("AA", "CBM", "CE", "GH", "GT", "PL"))  # 新建数据框

part=list()
for(i in 1:length(files))
{
    part[[i]]=data.frame(table(ml[[i]]))
    colnames(part[[i]])=c("CAZy_ID", Bin_ID[i])
    origin=merge(origin, part[[i]], by="CAZy_ID", all=T)
    # 合并Freq列，并重命名表头，空为NA
}

origin[is.na(origin)]=0  # 0替换NA

write.table(origin, file="cazy_group.txt", sep="\t", row.names=F, quote=F)  # 保存
