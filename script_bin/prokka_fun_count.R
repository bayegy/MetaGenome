setwd("Bin_all/Bin_prokka/prokka_out_table")
files=list.files(pattern="bin.*.tsv")  # 读取所有文件名

Bin_ID=vector()
for(i in 1:length(files))
{
    Bin_ID[i]=as.character(strsplit(files[i], split=".tsv"))
    # 提取所有文件名
}

ml=list()
for(i in 1:length(files))
{
    ml[[i]]=read.table(files[i], sep='\t', na.string="", stringsAsFactors=F, header=T, quote="", comment.char="")
    # 读取所有数据框
}

EC_num=vector()
for(i in 1:length(files))
{
    EC_num[i]=lengths(ml[[i]]["EC_number"])-sum(is.na(ml[[i]]["EC_number"]))
    # EC总数/BIN
}

COG_num=vector()
for(i in 1:length(files))
{
    COG_num[i]=lengths(ml[[i]]["COG"])-sum(is.na(ml[[i]]["COG"]))
    # COG总数/BIN
}

GENE_num=vector()
for(i in 1:length(files))
{
    GENE_num[i]=lengths(ml[[i]]["gene"])-sum(is.na(ml[[i]]["gene"]))
    # 基因总数/BIN
}

Total_num=vector()
for(i in 1:length(files))
{
    Total_num[i]=lengths(ml[[i]]["ftype"])
    # 功能预测总数/BIN
}

CDS_num=vector()
rRNA_num=vector()
tRNA_num=vector()
tmRNA_num=vector()

for(i in 1:length(files))
{
    CDS_num[i]=0
    rRNA_num[i]=0
    tRNA_num[i]=0
    tmRNA_num[i]=0
    
    for(j in 1:length(ml[[i]][,"ftype"]))
    {
        if(ml[[i]][j, "ftype"]=="CDS")
        {
            CDS_num[i]=CDS_num[i]+1
        }
         else if(ml[[i]][j, "ftype"]=="rRNA")
        {
            rRNA_num[i]=rRNA_num[i]+1
        }
        else if(ml[[i]][j, "ftype"]=="tRNA")
        {
            tRNA_num[i]=tRNA_num[i]+1
        }
        else if(ml[[i]][j, "ftype"]=="tmRNA")
        {
            tmRNA_num[i]=tmRNA_num[i]+1
        }
    }
}

prokka_result=data.frame(Bin_ID, Total_num, EC_num, COG_num, GENE_num, CDS_num, tRNA_num, rRNA_num, tmRNA_num)

write.table(prokka_result, file="../prokka_out_table.txt", sep="\t", quote=F, row.names=F)
