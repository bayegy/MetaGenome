
mapping_file=$1
category_1=$2
#prefix=$3

if [ -z "$2" ]; then
    echo "##########

          Please prepare the directory structure before starting the program like below:
          raw/fastq_files ...
          mapping_file
          manifest_file
          \n\n"

    echo "Please provide following input parameters
        1) Full path of the mapping file. (Accept both .csv or txt format.)
        2) The name of the first category in the mapping file. 

        Sample Usage:
        sh $0 M231_Mapping_2.tsv ${category_1} readme.pdf
        "
    exit 0
else
    echo "################
    Running: sh $0 $1 $2"
fi

check_file() {
    echo "Checking file for $1 ..."
    file_name=$1
    if [ -f $file_name ]; then
        echo "File $file_name exists"
    else
        echo "File $file_name does not exist"
        exit
    fi
}



echo "##############################################################\n#Organize the Result folder"

SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"
if [ -d "./Result_Metagenomics" ];then
    rm -r ./Result_Metagenomics;
fi;


mkdir -p Result_Metagenomics/FiguresTablesForReport \
    Result_Metagenomics/2-TaxaAundanceAnalysis/1-AbundanceSummary/1-RelativeAbundance \
    Result_Metagenomics/2-TaxaAundanceAnalysis/1-AbundanceSummary/2-Barplots \
    Result_Metagenomics/2-TaxaAundanceAnalysis/1-AbundanceSummary/4-Abundance_Metaphlan2 \
    Result_Metagenomics/2-TaxaAundanceAnalysis/1-AbundanceSummary/3-Heatmaps \
    Result_Metagenomics/2-TaxaAundanceAnalysis/2-AbundanceComparison/VennAndFlower \
    Result_Metagenomics/2-TaxaAundanceAnalysis/2-AbundanceComparison/ANCOM \
    Result_Metagenomics/2-TaxaAundanceAnalysis/2-AbundanceComparison \
    Result_Metagenomics/2-TaxaAundanceAnalysis/2-AbundanceComparison/LEfSe \
    Result_Metagenomics/2-TaxaAundanceAnalysis/3-DiversityAnalysis \
    Result_Metagenomics/1-QCStats/2-QC_report_Filtered \
    Result_Metagenomics/1-QCStats/1-QC_report_Rawfastq \
    Result_Metagenomics/3-FuctionAnalysis/2-Metacyc_Humann2 \
    Result_Metagenomics/3-FuctionAnalysis/1-KEGG/4-SignificanceAnalysis \
    Result_Metagenomics/4-AMRAnalysis/4-SignificanceAnalysis \
    Result_Metagenomics/FiguresTablesForReport

cp $mapping_file Result_Metagenomics/
cp QC_report/*html Result_Metagenomics/1-QCStats/1-QC_report_Rawfastq/
rm Result_Metagenomics/1-QCStats/1-QC_report_Rawfastq/*unmapped*
mv Result_Metagenomics/1-QCStats/1-QC_report_Rawfastq/*good* Result_Metagenomics/1-QCStats/2-QC_report_Filtered/
cp Report/reads_summary.txt Result_Metagenomics/1-QCStats/


cp Kraken2/*/Relative/Classified_stat_relative.png Kraken2/*/All.Taxa.OTU.taxonomy.biom Kraken2/*/All.Taxa.OTU.txt Result_Metagenomics/2-TaxaAundanceAnalysis/1-AbundanceSummary/
cp Kraken2/*/Relative/*relative.txt Result_Metagenomics/2-TaxaAundanceAnalysis/1-AbundanceSummary/1-RelativeAbundance/
cp -rp Kraken2/*/All.Taxa.OTU.taxa-bar-plots Kraken2/*/All.Taxa.OTU.taxa-bar-plots.1000  Kraken2/*/Barplot-of-Group-Mean Kraken2/*/Taxa-bar-plots-top20  Result_Metagenomics/2-TaxaAundanceAnalysis/1-AbundanceSummary/2-Barplots/
cp -r Kraken2/*/Heatmap/* Result_Metagenomics/2-TaxaAundanceAnalysis/1-AbundanceSummary/3-Heatmaps/
cp -r Metagenome/Metaphlan/All*txt Result_Metagenomics/2-TaxaAundanceAnalysis/1-AbundanceSummary/4-Abundance_Metaphlan2/


cp -r Kraken2/*/ANCOM/*ANCOM* Result_Metagenomics/2-TaxaAundanceAnalysis/2-AbundanceComparison/ANCOM/
cp -r Kraken2/*/Lefse/* Result_Metagenomics/2-TaxaAundanceAnalysis/2-AbundanceComparison/LEfSe/

cp -r Kraken2/*/VennAndFlower/*  Result_Metagenomics/2-TaxaAundanceAnalysis/2-AbundanceComparison/VennAndFlower/
# Kraken2/*/core-metrics/alpha
cp -rp  Kraken2/*/core-metrics/bray_curtis_emperor Kraken2/*/PCoA-NMDS/* Result_Metagenomics/2-TaxaAundanceAnalysis/3-DiversityAnalysis/
cp -rp Kraken2/*/CorrelationAnalysis Result_Metagenomics/2-TaxaAundanceAnalysis/4-CorrelationAnalysis
#rm -r Result_Metagenomics/2-TaxaAundanceAnalysis/3-DiversityAnalysis/alpha/chao1 \
#    Result_Metagenomics/2-TaxaAundanceAnalysis/3-DiversityAnalysis/alpha/observed_otus \
#    Result_Metagenomics/2-TaxaAundanceAnalysis/3-DiversityAnalysis/alpha/shannon





cp AMR/All.AMR.abundance.txt Result_Metagenomics/4-AMRAnalysis/
  
cp -rp AMR/Barplots    Result_Metagenomics/4-AMRAnalysis/1-Barplots 
cp -rp AMR/Heatmaps    Result_Metagenomics/4-AMRAnalysis/2-Heatmaps
cp -rp AMR/Circos   Result_Metagenomics/4-AMRAnalysis/3-Circos
rm -r Result_Metagenomics/4-AMRAnalysis/3-Circos/circos_conf
cp AMR/LEfSe/* Result_Metagenomics/4-AMRAnalysis/4-SignificanceAnalysis/
cp -rp AMR/CorrelationAnalysis Result_Metagenomics/4-AMRAnalysis/5-CorrelationAnalysis







cp FMAP/All*.txt Result_Metagenomics/3-FuctionAnalysis/1-KEGG/
cp -rp FMAP/Barplots    Result_Metagenomics/3-FuctionAnalysis/1-KEGG/1-Barplots
cp -rp FMAP/Heatmaps    Result_Metagenomics/3-FuctionAnalysis/1-KEGG/2-Heatmaps
cp -rp FMAP/Circos    Result_Metagenomics/3-FuctionAnalysis/1-KEGG/3-Circos
rm -r Result_Metagenomics/3-FuctionAnalysis/1-KEGG/3-Circos/circos_conf
cp  FMAP/LEfSe/*    Result_Metagenomics/3-FuctionAnalysis/1-KEGG/4-SignificanceAnalysis/
cp -rp FMAP/ColoredMaps    Result_Metagenomics/3-FuctionAnalysis/1-KEGG/5-ColoredMaps
cp -rp FMAP/CorrelationAnalysis    Result_Metagenomics/3-FuctionAnalysis/1-KEGG/6-CorrelationAnalysis


#mv Result_Metagenomics/3-FuctionAnalysis/1-KEGG/*lefse* Result_Metagenomics/3-FuctionAnalysis/1-KEGG/LEfSe/
cp Metagenome/Humann/All.*tsv Result_Metagenomics/3-FuctionAnalysis/2-Metacyc_Humann2/




################################################make FiguresTablesForReport
cp -rp ${SCRIPTPATH}/Report/src Result_Metagenomics/FiguresTablesForReport/
cp ${SCRIPTPATH}/Report/结题报告.html Result_Metagenomics/

cd Result_Metagenomics/FiguresTablesForReport
cp ../2-TaxaAundanceAnalysis/1-AbundanceSummary/Classified_stat_relative.png Figure4-1.png
cp ../2-TaxaAundanceAnalysis/1-AbundanceSummary/2-Barplots/Taxa-bar-plots-top20/Phylum_${category_1}_barplot.pdf Figure4-2.pdf
cp ../2-TaxaAundanceAnalysis/1-AbundanceSummary/3-Heatmaps/Heatmap_top20_clustered/Phylum/${category_1}_heatmap.pdf  Figure4-3.pdf
cp ../2-TaxaAundanceAnalysis/2-AbundanceComparison/LEfSe/Genus/${category_1}_Genus_lefse_LDA2.pdf  Figure4-4.pdf
cp ../2-TaxaAundanceAnalysis/2-AbundanceComparison/VennAndFlower/${category_1}_Venn_plot.png Figure4-5.png

cp ../3-FuctionAnalysis/1-KEGG/1-Barplots/Pathway.Level1_${category_1}_barplot.pdf Figure5-1.pdf
cp ../3-FuctionAnalysis/1-KEGG/4-SignificanceAnalysis/Pathway_${category_1}_lefse_LDA2.pdf   Figure5-2.pdf
cp ../3-FuctionAnalysis/1-KEGG/5-ColoredMaps/${category_1}/${category_1}_map00010.png   Figure5-3.png

cp ../4-AMRAnalysis/1-Barplots/AMR_${category_1}_barplot.pdf Figure6-1.pdf
cp ../4-AMRAnalysis/2-Heatmaps/AMR_${category_1}_nocluster_heatmap.pdf Figure6-2.pdf
# cp ../4-AMRAnalysis/2-Heatmaps/AMR_${category_1}_clustered_heatmap.pdf Figure6-2.pdf
cp ../4-AMRAnalysis/3-Circos/AMR_${category_1}_circos.png Figure6-3.png


cp ../2-TaxaAundanceAnalysis/4-CorrelationAnalysis/RDA/Genus/${category_1}_RDA_features_location_plot.pdf Figure7-1.pdf
cp ../4-AMRAnalysis/5-CorrelationAnalysis/CorrelationHeatmap/AMR_Correlation_heatmap.pdf Figure7-2.pdf


#cp -r ../2-TaxaAundanceAnalysis/1-AbundanceSummary/2-Barplots/All.Taxa.OTU.taxa-bar-plots page4-2
#cp -r ../2-TaxaAundanceAnalysis/3-DiversityAnalysis/bray_curtis_emperor page4-5


# python3 ${SCRIPTPATH}/convert_to_html_table.py -i ../1-QCStats/reads_summary.txt -o src/pages/main_cleaned.html -t txt -k '{{table1}}'


if [ -f Figure4-2.pdf ];then echo "Converting pdf to png"; for pdfs in *.pdf; do echo $pdfs; base=$(basename $pdfs .pdf); convert  -density 300 -quality 80 $pdfs ${base}.png; rm $pdfs;done;fi;
