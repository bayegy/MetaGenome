#!/usr/bin/env bash
SCRIPT_BIN="$( cd "$(dirname "$0")" ; pwd -P )"
python3_path=/home/bayegy/pipelines/metagenome/miniconda2/bin/python3.8
R_path=/home/bayegy/pipelines/metagenome/miniconda2/bin/Rscript
perl_path=/home/bayegy/pipelines/metagenome/miniconda2/bin/perl
metagenome_home=/home/bayegy/pipelines/metagenome/MetaGenome
bayegy_home=/home/bayegy/pipelines/metagenome/Bayegy
result_root=$1
metadata=$2
all_group=$3

cd $result_root/../
<<skip
mkdir -p $result_root/Bin_summary

n90_file=$result_root/Bin_summary/N50.N90.txt
echo -e "\033[32m统计N50和N90: \033[0m"
echo -e "BinID\tN50\tN90" > $n90_file
for i in $result_root/Bin_pick/*.fa; do
    base=${i##*/}
    name=${base%.fa}
    echo -ne "$name\t" >> $n90_file
    $perl_path $SCRIPT_BIN/N50.N90.pl $i | \
    sed 's/N50: //g' | \
    sed 's/N90: //g' | \
    tr -s '\n' '\t' | \
    awk '{printf"%s\t%s\n",$1,$2}' >> $n90_file
    echo -e "\033[32m$i Done...\033[0m"
done


cat $result_root/bin_summary/bin_drep/data/checkM/checkM_outdir/results.tsv | \
 sed 's/\.fa\t/\t/g' | awk -F '\t' 'BEGIN{OFS="\t"}{print $1,$12,$13}' \
 > $result_root/Bin_summary/completeness.txt

$python3_path $metagenome_home/stat_bin_gc.py $result_root/Bin_pick/ \
 $result_root/Bin_summary/bin.gc.txt

$python3_path $bayegy_home/merge_tables.py  \
 $result_root/Bin_summary/N50.N90.txt \
 $result_root/Bin_summary/bin.gc.txt \
 $result_root/Bin_summary/completeness.txt \
 - | awk -F "\t" 'BEGIN{OFS="\t"}{print $1,$6,$7,$2,$3,$4,$5}' \
 > $result_root/Bin_summary/bin_summary.txt


mkdir $result_root/Bin_plot
echo -e "\033[32m计算每个Bin中每条contig的GC含量：\033[0m"
for i in $result_root/Bin_pick/bin.*.fa; do
    file=${i##*/}
    fold=${file%.fa}
    mkdir $result_root/Bin_plot/$fold
    $python3_path $SCRIPT_BIN/fasta_gc_percent.py $i $result_root/Bin_plot/$fold/${fold}.gc.txt
    echo -e "\033[32m\t $i gc percent Done...\033[0m"
done

echo -e "\033[32m获取contig深度数据：\033[0m"
cat $result_root/bin_out/*/depth.txt | sed -e '/totalAvgDepth/d' | \
 sed '1i contigName\tcontigLen\ttotalAvgDepth\talign.sorted.bam\talign.sorted.bam-var' | \
 awk -F"\t" 'BEGIN{OFS="\t"}{print $1, $3}' > $result_root/Bin_plot/all_contig_avgdepth.txt

echo -e "\033[32mmerge contig gc depth：\033[0m"
for i in $result_root/Bin_plot/bin.*; do
    fold=${i##*/}
    $R_path $SCRIPT_BIN/merge_gc_depth.R $i/${fold}.gc.txt ../all_contig_avgdepth.txt
    echo -e "\033[32mmerge $i gc depth Done...\033[0m"
done

echo -e "\033[32mmerge all bin gc depth data: \033[0m"
touch $result_root/Bin_plot/all_bin_gc_depth.txt
for i in $result_root/Bin_plot/bin*; do
    fold=${i##*/}
    cat $i/${fold}.gc.depth.txt | sed '1d' >> $result_root/Bin_plot/all_bin_gc_depth.txt
    echo -e "\033[32m\tadd $i gc depth data: \033[0m"
done
sed -i '1 icontig\tGC_percent\tDepth\tBin' $result_root/Bin_plot/all_bin_gc_depth.txt

echo -e "\033[32mbin gc depth 可视化: \033[0m"
$R_path $SCRIPT_BIN/bin_scatter.R $result_root/Bin_plot/all_bin_gc_depth.txt all_bin_gc_depth.pdf

convert -density 400 -quality 200 $result_root/Bin_plot/all_bin_gc_depth.pdf $result_root/Bin_plot/all_bin_gc_depth.png


#######################
##### 计算bin丰度 #####
#######################

mv $result_root/Bin_quant/bin_abundance_table.tab $result_root/Bin_quant/bin_abundance_table.txt



####################################
##### prokka：type gene EC COG #####
####################################


if [ ! -d "$result_root/Bin_prokka" ]; then
    echo -e "\033[31merror: 文件夹不存在\033[0m"
else
    mkdir -p $result_root/Bin_prokka/prokka_out_table \
    $result_root/Bin_prokka/prokka_map \
    $result_root/Bin_prokka/prokka_map_table
fi

echo -e "\033[32m集合prokka注释tsv文件：\033[0m"
for i in $result_root/Bin_prokka/bin*; do
    fold=${i##*/};
    mv $i/$fold.tsv $result_root/Bin_prokka/prokka_out_table;  # prokka注释结果
    mv $i/$fold.gff $result_root/Bin_prokka/prokka_map; # prokka比对文件
done


echo -e "\033[32m集合提取prokka gff比对信息：\033[0m"
for i in $result_root/Bin_prokka/prokka_map/bin.*.gff; do
    base=${i##*/}
    grep '^contig' $i > $result_root/Bin_prokka/prokka_map_table/${base}.txt
    sed -i '1 iseqid\tsource\ttype\tstart\tend\tscore\tstrand\tphase\tattributes' $result_root/Bin_prokka/prokka_map_table/${base}.txt
done

echo -e "\033[32m统计prokka注释tsv文件：\033[0m"
$R_path $SCRIPT_BIN/prokka_fun_count.R

mv $result_root/*.Rout tmp
mv *.Rout tmp
skip


############################
##### cazyme数据库注释 #####
############################
if [ ! -f "$result_root/Bin_prokka/prokka_out_table.txt" ]; then
    echo -e "\033[31merror: 文件不存在\033[0m"
fi

echo -e "\033[32m统计cazy level 1 (group): \033[0m"
for i in $result_root/Bin_cazy/bin.*.raw; do
    base=${i##*/}
    name=${base%.raw}
    awk '{print $1,$2}' $i | uniq | awk '{print $2}' | sed 's/|/ /g' | awk '{print $2}' | sed 's/[|_]/ /g' | awk '{print $1}' | sed 's/[0123456789]//g' > $result_root/Bin_cazy/$name.clean
    echo -e "\033[32m$i Done...\033[0m"
done

echo -e "\033[32mcazyme group 统计: \033[0m"
$R_path $SCRIPT_BIN/cazy_group.R
echo -e "\033[32mcazyme group 热图、箱图: \033[0m"
$R_path $SCRIPT_BIN/cazy_group_pheatmap_boxplot.R

echo -e "\033[32mcazyme group 饼图: \033[0m"
$R_path $SCRIPT_BIN/cazy_group_pie.R
echo -e "\033[32mcazyme group 堆叠图: \033[0m"
$R_path $SCRIPT_BIN/cazy_group_stack.R

echo -e "\033[32m统计cazy level 2 (all) 统计: \033[0m"
for i in $result_root/Bin_cazy/bin.*.raw; do
    base=${i##*/}
    name=${base%.raw}
    count_cazy.py $i > $result_root/Bin_cazy/$name.sum
    echo -e "\033[32m$i Done...\033[0m"
done
echo -e "\033[32m统计cazy level 2 (all) 统计: \033[0m"
$R_path $SCRIPT_BIN/cazy_all.R
echo -e "\033[32m统计cazy level 2 (all) 饼图: \033[0m"
$R_path $SCRIPT_BIN/cazy_all_pie.R

echo -e "\033[32mcazyme pdf to png: \033[0m"
convert -density 400 -quality 200 $result_root/Bin_cazy/cazy_group_pheatmap.pdf $result_root/Bin_cazy/cazy_group_pheatmap.png
convert -density 400 -quality 200 $result_root/Bin_cazy/cazy_group_boxplot.pdf $result_root/Bin_cazy/cazy_group_boxplot.png
convert -density 400 -quality 200 $result_root/Bin_cazy/cazy_group_stack.pdf $result_root/Bin_cazy/cazy_group_stack.png

for i in $result_root/Bin_cazy/cazy_all_pie/*.pdf; do
    base=${i%.pdf}
    convert -density 400 -quality 200 $i ${base}.png
    echo -e "\033[32m\t $i all pie Done...\033[0m"
done

for i in $result_root/Bin_cazy/cazy_group_pie/*.pdf; do
    base=${i%.pdf}
    convert -density 400 -quality 200 $i ${base}.png
    echo -e "\033[32m\t $i group pie Done...\033[0m"
done

mv *.Rout tmp

<<skip
#######################
##### bin进化分析 #####
#######################
if [ ! -f "$result_root/Bin_cazy/cazy_group_stack.png" ]; then
    echo -e "\033[31merror: 文件不存在\033[0m"
else
    mkdir $result_root/Bin_phylo
fi

echo -e "\033[32mphylophlan进化分析: \033[0m"
mkdir /home/cheng/huty/softwares/phylophlan/input/$project.denovo
mkdir /home/cheng/huty/softwares/phylophlan/input/$project.insert

for i in $result_root/Bin_prokka/bin.*; do
    fold=${i##*/}
    cp $i/$fold.faa /home/cheng/huty/softwares/phylophlan/input/$project.denovo
    cp $i/$fold.faa /home/cheng/huty/softwares/phylophlan/input/$project.insert
done

temp=`pwd`
cd /home/cheng/huty/softwares/phylophlan
echo -e "\033[32mde novo phylophlan进化分析: \033[0m"
python phylophlan.py -u $project.denovo --nproc $thread
echo -e "\033[32minsert phylophlan进化分析: \033[0m"
python phylophlan.py -i -t $project.insert --nproc $thread
cd $temp

temp=`pwd`
cd /home/cheng/huty/softwares/phylophlan/output/$project.insert
echo -e "\033[32m合并insert分析物种注释结果: \033[0m"
cat imputed_conf_* | sort -t $'\t' -k 2r | sed '1 iid\ttax' > bin_taxonomy.txt
cd ../
mv /home/cheng/huty/softwares/phylophlan/output/$project.insert $temp/$result_root/Bin_phylo 
mv /home/cheng/huty/softwares/phylophlan/output/$project.denovo $temp/$result_root/Bin_phylo
cd $temp

rm -r /home/cheng/huty/softwares/phylophlan/input/$project.denovo
rm -r /home/cheng/huty/softwares/phylophlan/input/$project.insert
rm -r /home/cheng/huty/softwares/phylophlan/data/$project.denovo
rm -r /home/cheng/huty/softwares/phylophlan/data/$project.insert

echo -e "\033[32m绘制de novo phylophlan进化树: \033[0m"
temp=`pwd`
cd $result_root/Bin_phylo/$project.denovo
sed '1d' ../$project.insert/bin_taxonomy.txt | sed 's/.p__/\t/g' | sed 's/.c__/\t/g' | awk '{print $1,$3}' | sed 's/ /\t/g' | sed '1 ibin\tPhylum' > map.txt
$R_path $SCRIPT_BIN/tree_bins_circ.R  # 画图
convert -density 400 -quality 200 tree_bins.pdf tree_bins.png
cd $temp

mv *.Rout tmp

##########################
##### 绘制circos圈图 #####
##########################
if [ ! -f "$result_root/Bin_phylo/$project.insert/bin_taxonomy.txt" ]; then
    echo -e "\033[31merror: 文件不存在\033[0m"
else
    mkdir $result_root/Bin_circos
    mkdir $result_root/Bin_circos/Figure
fi

bash /home/cheng/huty/softwares/script_circos/circos.sh

echo -e "\033[32mpdf to png: \033[0m"
for i in $result_root/Bin_circos/Figure/*.circos.pdf; do
    base=${i%.pdf}
    convert -density 400 -quality 200 $i ${base}.png
done 

mv *.Rout tmp
mv .RData tmp

##########################
##### kegg注释、绘图 #####
##########################
if [ ! -f "$result_root/Bin_prokka/prokka_out_table.txt" ]; then
    echo -e "\033[31merror: 文件不存在\033[0m"
else
    mkdir $result_root/kegg
fi

echo -e "\033[32m\nkofamscan注释KEGG：\033[0m"
for i in $result_root/Bin_prokka/bin.*; do
    fold=${i##*/}
    mkdir $result_root/kegg/$fold
    exec_annotation -f mapper -o $result_root/kegg/$fold/${fold}_kegg_raw.txt $i/${fold}.faa
    echo -e "\033[32m\t$i KEGG Done...\033[0m"
done

for i in $result_root/kegg/bin*; do
    fold=${i##*/}
    cat $i/${fold}_kegg_raw.txt | awk -F"\t" '{if($2 != "") print $2,$1}' | sed 's/ /\t/' | sed '1 ik_id\tgene_id' > $i/${fold}_kegg.txt
    echo -e "\033[32m\t$i KEGG注释结果整理 Done...\033[0m"
done

for i in $result_root/kegg/bin*; do
    fold=${i##*/}
    $R_path /home/cheng/huty/softwares/script_genome/kofamscan_annotation.R $i ${fold}_kegg.txt ${fold}_kegg_pathway.txt /home/cheng/huty/databases/kofamkoala/KEGG_orthology2pathway.txt
    echo -e "\033[32m\t$i KEGG to pathway Done...\033[0m"
done

for i in $result_root/kegg/bin*; do
    fold=${i##*/}
    cat $i/${fold}_kegg_pathway.txt | sed '1d' | awk -F"\t" '{print $3}' | sort | uniq -c | awk '{printf("%s\t%s\n", $2, $1)}' | sort -t $'\t' -k 2nr | sed '1 ipathway_id\tpathway_num' > $i/${fold}_kegg_pathway_num.txt
    echo -e "\033[32m\t$i KEGG pathway 统计 Done...\033[0m"
    $R_path $SCRIPT_BIN/kegg_level_annotation.R $i ${fold}_kegg_pathway_num.txt ${fold}_kegg_pathway_num_annotation.txt /home/cheng/huty/databases/KEGG/pathway_annotation.txt    
    echo -e "\033[32m\t$i KEGG pathway 分类注释 Done...\033[0m"
    $R_path $SCRIPT_BIN/kegg_pathway_level_bar.R ${i}/${fold}_kegg_pathway_num_annotation.txt ${fold}_kegg_pathway_num_annotation
    echo -e "\033[32m\t$i KEGG pathway barplot绘图 Done...\033[0m"
done

for i in $result_root/kegg/bin*; do
    fold=${i##*/}
    $R_path /home/cheng/huty/softwares/script_genome/kofamscan_annotation.R $i ${fold}_kegg.txt ${fold}_kegg_annotation.txt /home/cheng/huty/databases/kofamkoala/ko_mapping.txt
    echo -e "\033[32m\t$i KEGG pathway 详细注释 Done...\033[0m"
done

for i in $result_root/kegg/bin*; do
    fold=${i##*/}
    mkdir $i/Figure
    map=`cat $i/${fold}_kegg_pathway.txt | sed '1d' | awk -F"\t" '{print $3}' | sort | uniq`
    for j in $map; do
        cp /home/cheng/Databases/map/${j}.png $i/Figure/
    done
    echo -e "\033[32m\t$i Figure Done...\033[0m"
done

########################
##### GO注释、绘图 #####
########################
if [ ! -f "$result_root/Bin_prokka/prokka_out_table.txt" ]; then
    echo -e "\033[31merror: 文件不存在\033[0m"
else
    mkdir $result_root/emapper
fi

echo -e "\033[32m\n emapper注释GO：\033[0m"

for i in $result_root/Bin_prokka/bin.*; do
    fold=${i##*/}
    mkdir $result_root/emapper/$fold
    emapper.py -i $i/${fold}.ffn \
--output_dir $result_root/emapper/$fold \
-o $fold \
--seed_ortholog_evalue 0.00001 \
--data_dir /media/cheng/disk2/Databases/eggnog \
-m diamond --cpu $thread
    echo -e "\033[32m\t$i emapper Done...\033[0m"
done

for i in $result_root/emapper/bin.*; do
    fold=${i##*/}
    #cat $i/${fold}.emapper.annotations | awk -F"\t" '{print $7}' | sed '/^$/d' | sed 's/,/\n/g' | sed '1d' | sort | uniq | sed '1 iKEGG_KOs' > $i/${fold}_kegg.txt
    cat $i/${fold}.emapper.annotations | sed '/^#/d' | awk -F"\t" '{print $1,$6}' | awk '{if($2 != "") print $0}' | sed '1d'| sed '1 igene_id\tGO_terms' > $i/${fold}_go_map.txt
    cat $i/${fold}_go_map.txt | sed '1d' | awk '{print $2}' | sed 's/,/\n/g' | sort | uniq -c | awk '{print $2,$1}' | sed 's/ /\t/' | sed '1 igo_id\tgo_gene_num' > $i/${fold}_go.txt
    $R_path /home/cheng/huty/softwares/script_genome/emapper_go_annotation.R $i ${fold}_go.txt ${fold}_go_annotation.txt /home/cheng/huty/databases/GO/ncbi/go.annotation.txt
    $R_path /home/cheng/huty/softwares/script_genome/emapper_go_bar.R $i/${fold}_go_annotation.txt ${fold}_go_annotation.pdf
    convert -density 400 -quality 200 $i/${fold}_go_annotation.pdf $i/${fold}_go_annotation.png
    echo -e "\033[32m\t$i GO Done...\033[0m"
done

mv *.Rout tmp
mv .RData tmp

########################
##### cog注释、绘图 #####
########################
if [ ! -f "$result_root/Bin_prokka/prokka_out_table.txt" ]; then
    echo -e "\033[31merror: 文件不存在\033[0m"
else
    mkdir $result_root/cog
fi

echo -e "\033[32mextract cog: \033[0m"
for i in $result_root/Bin_prokka/prokka_out_table/bin.*.tsv; do
    file=${i##*/}
    fold=${file%.tsv}
    mkdir $result_root/cog/$fold
    cat $i | awk -F"\t" '{if($6 != "") printf("%s\t%s\n", $1, $6)}' > $result_root/cog/$fold/${fold}_cog.txt
    echo -e "\033[32m\t $i cog extract Done...\033[0m"
done

echo -e "\033[32mannotate cog: \033[0m"
for i in $result_root/cog/bin.*; do
    fold=${i##*/}
    $R_path /home/cheng/huty/softwares/script_genome/prokka_cog_annotation.R $i/${fold}_cog.txt /home/cheng/huty/databases/cog/cog_anno.txt ${fold}_cog_annotation.txt
    echo -e "\033[32m\t $i annotate cog Done...\033[0m"
done

echo -e "\033[32msummary annotated cog: \033[0m"
for i in $result_root/cog/bin.*; do
    fold=${i##*/}
    cat $i/${fold}_cog_annotation.txt | awk -F"\t" 'BEGIN{OFS="\t"}{if($6 != "NA") print $3}' | sed '1d' | sort | uniq -c | awk 'BEGIN{OFS="\t"}{print $2,$1}' | sed '1 ifunc\tCount' > $i/${fold}_cog_sum.txt
    echo -e "\033[32m\t $i summary annotated cog Done...\033[0m"
done

echo -e "\033[32m分类 annotated cog: \033[0m"
for i in $result_root/cog/bin.*; do
    fold=${i##*/}
    $R_path /home/cheng/huty/softwares/script_genome/prokka_cog_count_annotation.R $i/${fold}_cog_sum.txt /home/cheng/huty/databases/cog/cog_category.txt ${fold}_cog_sum_annotation.txt
    echo -e "\033[32m\t $i 分类 annotated cog Done...\033[0m"
done

echo -e "\033[32m分类可视化 annotated cog: \033[0m"
for i in $result_root/cog/bin.*; do
    fold=${i##*/}
    $R_path /home/cheng/huty/softwares/script_genome/prokka_cog_count_annotation_bar.R $i/${fold}_cog_sum_annotation.txt ${fold}_cog_sum_annotation.pdf
    convert -density 400 -quality 200 $i/${fold}_cog_sum_annotation.pdf $i/${fold}_cog_sum_annotation.png
    echo -e "\033[32m\t $i 分类可视化 annotated cog Done...\033[0m"
done

########################
##### 整理结果文件 #####
########################
if [ ! -f "$result_root/Bin_prokka/prokka_out_table.txt" ]; then
    echo -e "\033[31merror: 文件不存在\033[0m"
else
    mkdir $result_root/Results
fi

echo -e "\033[32mcopy Report.html image: \n\033[0m"
cp $SCRIPT_BIN/Report.html $result_root/Results/
cp -r $SCRIPT_BIN/image $result_root/Results/

echo -e "\033[32m结果文件整理: \n\033[0m"
# 0 1 质控
echo -e "\033[32m\t1 QC\n\033[0m"
mkdir $result_root/Results/00-RawDataQC
mkdir $result_root/Results/01-CleanDataQC
cp $result_root/kneaddata/fastqc/reformatted_* $result_root/Results/00-RawDataQC
cp $result_root/kneaddata/fastqc/*_R1_kneaddata.trimmed.[12]_fastqc* $result_root/Results/01-CleanDataQC

cp $result_root/summary_rawdata.txt $result_root/Results/00-RawDataQC
cp $result_root/summary_cleandata.txt $result_root/Results/01-CleanDataQC

$python3_path $SCRIPT_BIN/table_insert_html.py $result_root/Results/Report.html table4 $result_root/summary_rawdata_html.txt $result_root/Results/temp.html

$python3_path $SCRIPT_BIN/table_insert_html.py $result_root/Results/Report.html table5 $result_root/summary_cleandata_html.txt $result_root/Results/temp.html

# 2 分箱
echo -e "\033[32m\t2 binning\n\033[0m"
mkdir $result_root/Results/02-Bin
cp $result_root/Bin_quality/checkm_pick.txt $result_root/Results/02-Bin
cp $result_root/Bin_quality/checkm_raw.txt $result_root/Results/02-Bin
cp -r $result_root/Bin_pick $result_root/Results/02-Bin
cp -r $result_root/Bin $result_root/Results/02-Bin
cp $result_root/Bin_summary/bin_summary.txt $result_root/Results/02-Bin
echo -e "\033[32m\ttable to heml\n\033[0m"
$python3_path $SCRIPT_BIN/table_to_html_bin_summary.py $result_root/Bin_summary/bin_summary.txt $result_root/Bin_summary/bin_summary_html.txt

$python3_path $SCRIPT_BIN/table_insert_html.py $result_root/Results/Report.html table1 $result_root/Bin_summary/bin_summary_html.txt $result_root/Results/temp.html

# 3 可视化
echo -e "\033[32m\t3 visualize\n\033[0m"
mkdir $result_root/Results/03-Bin_Plot
cp $result_root/Bin_plot/all_bin_gc_depth* $result_root/Results/03-Bin_Plot

# 4 丰度
echo -e "\033[32m\t4 abundance\n\033[0m"
mkdir $result_root/Results/04-Bin_Abundance
cp $result_root/Bin_quant/bin* $result_root/Results/04-Bin_Abundance
cp -r $result_root/Bin_quant/lefse $result_root/Results/04-Bin_Abundance/

cp $result_root/Results/04-Bin_Abundance/lefse/Category1/lefse.png $result_root/Results/image/example_Category1_lefse.png

# 5 功能
echo -e "\033[32m\t5 Function\n\033[0m"
mkdir $result_root/Results/05-Bin_Function
cp $result_root/Bin_prokka/prokka_out_table.txt $result_root/Results/05-Bin_Function/
cp -r $result_root/Bin_prokka/prokka_out_table $result_root/Results/05-Bin_Function/
cp -r $result_root/cog $result_root/Results/05-Bin_Function/  # 举例图

example_cog=`ls $result_root/Results/05-Bin_Function/cog/ | head -n 1`
cp $result_root/Results/05-Bin_Function/cog/$example_cog/${example_cog}_cog_sum_annotation.png $result_root/Results/image/example_cog.png

echo -e "\033[32m\ttable to heml\n\033[0m"
$python3_path $SCRIPT_BIN/table_to_html_Function.py $result_root/Bin_prokka/prokka_out_table.txt $result_root/Bin_prokka/prokka_out_table_html.txt

$python3_path $SCRIPT_BIN/table_insert_html.py $result_root/Results/Report.html table2 $result_root/Bin_prokka/prokka_out_table_html.txt $result_root/Results/temp.html

# 6 kegg
echo -e "\033[32m\t6 KEGG\n\033[0m"
mkdir $result_root/Results/06-Bin_KEGG
cp -r $result_root/kegg/* $result_root/Results/06-Bin_KEGG/  # 举例图

example_bin=`ls $result_root/Results/06-Bin_KEGG/ | head -n 1`
example_kegg=`ls $result_root/Results/06-Bin_KEGG/$example_bin/Figure/ | head -n 1`
cp $result_root/Results/06-Bin_KEGG/$example_bin/Figure/$example_kegg $result_root/Results/image/example_kegg.png
cp $result_root/Results/06-Bin_KEGG/$example_bin/${example_bin}_kegg_pathway_num_annotation.png $result_root/Results/image/example_kegg_bar.png

# 7 GO
echo -e "\033[32m\t7 GO\n\033[0m"
mkdir $result_root/Results/07-Bin_GO

for i in $result_root/emapper/bin.*; do
    fold=${i##*/}
    mkdir $result_root/Results/07-Bin_GO/$fold
    # cp $i/${fold}_kegg.txt $result_root/Results/07-Bin_GO/$fold/
    cp $i/${fold}_go* $result_root/Results/07-Bin_GO/$fold/ 
done

example_go=`ls $result_root/Results/07-Bin_GO/ | head -n 1`
cp $result_root/Results/07-Bin_GO/$example_go/${example_go}_go_annotation.png $result_root/Results/image/example_go.png

# 8 cazyme
echo -e "\033[32m\t8 CAZyme\n\033[0m"
mkdir $result_root/Results/08-Bin_CAZyme
mkdir $result_root/Results/08-Bin_CAZyme/Data
cp -r $result_root/Bin_cazy/cazy_* $result_root/Results/08-Bin_CAZyme/
cp $result_root/Bin_cazy/bin* $result_root/Results/08-Bin_CAZyme/Data/

example_group_pie=`ls $result_root/Results/08-Bin_CAZyme/cazy_group_pie | head -n 2 | grep '.png'`
cp $result_root/Results/08-Bin_CAZyme/cazy_group_pie/$example_group_pie $result_root/Results/image/example_group_pie.png

example_all_pie=`ls $result_root/Results/08-Bin_CAZyme/cazy_all_pie | head -n 2 | grep '.png'`
cp $result_root/Results/08-Bin_CAZyme/cazy_all_pie/$example_all_pie $result_root/Results/image/example_all_pie.png

# 9 进化
echo -e "\033[32m\t5 phylophlan\n\033[0m"
mkdir $result_root/Results/09-Bin_Taxonomy
cp $result_root/Bin_phylo/$project.denovo/tree_bins.* $result_root/Results/09-Bin_Taxonomy
cp $result_root/Bin_phylo/$project.denovo/*tree.nwk $result_root/Results/09-Bin_Taxonomy
cp $result_root/Bin_phylo/$project.denovo/map.txt $result_root/Results/09-Bin_Taxonomy
cp $result_root/Bin_phylo/$project.insert/bin_taxonomy.txt $result_root/Results/09-Bin_Taxonomy

echo -e "\033[32m\ttable to heml\n\033[0m"
$python3_path $SCRIPT_BIN/table_to_html_taxonomy.py $result_root/Bin_phylo/$project.insert/bin_taxonomy.txt $result_root/Bin_phylo/$project.insert/bin_taxonomy_html.txt

$python3_path $SCRIPT_BIN/table_insert_html.py $result_root/Results/Report.html table3 $result_root/Bin_phylo/$project.insert/bin_taxonomy_html.txt $result_root/Results/temp.html

# 10 circos
echo -e "\033[32m\t10 Circos\n\033[0m"
mkdir $result_root/Results/10-Bin_Circos
mkdir $result_root/Results/10-Bin_Circos/Data
cp -r $result_root/Bin_circos/Figure $result_root/Results/10-Bin_Circos/
cp -r $result_root/Bin_circos/bin.* $result_root/Results/10-Bin_Circos/Data/  # 举例图

example_circos=`ls $result_root/Results/10-Bin_Circos/Figure | head -n 3 | grep 'circos.png'`
cp $result_root/Results/10-Bin_Circos/Figure/$example_circos $result_root/Results/image/example_circos.png
skip