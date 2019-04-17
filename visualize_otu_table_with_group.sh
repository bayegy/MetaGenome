#/bin/sh -S

#########
#Please address any bugs to Cheng. 
#Date 2017.11.31
#########
SCRIPTPATH=$1
otu_table=$2
mapping_file=$3
category_set=${4//,/ }
output_prefix=$5
taxa_filtered=$6

frequency=1000
depth=$(Rscript $SCRIPTPATH/min.R $otu_table)

echo "Selected depth was $depth"

echo "Check wheather your categories are the following:"
for i in $category_set;do echo $i;done


declare -A tax_aa;
tax_aa=([k]=Kingdom [p]=Phylum [c]=Class [o]=Order [f]=Family [g]=Genus [s]=Species);
tax_levels["1"]="Kingdom"
tax_levels["2"]="Phylum"
tax_levels["3"]="Class"
tax_levels["4"]="Order"
tax_levels["5"]="Family"
tax_levels["6"]="Genus"
tax_levels["7"]="Species"




if [ -z "$3" ]; then
	echo "\n\n"

	echo "Please provide following input parameters
		1) path of the OTU table file. 
		2) path of the meta data file.
		3) Group of interest from the meta file to be shown in the heatmaps.
		4) Output directory prefix name. (Required when filtered taxonomy will be specified)
		5) Taxonomy to be filtered, in the format of "tax1,tax2,tax3". (Optional)

		Sample Usage:
		sh $0 /share/data/IlluminaData/LYM/OTU_tables/Report/LYM_RTS.megablast.OTU_taxonomySummaryCounts.minHitCount_5.Bacteria.txt /share/data/IlluminaData/LYM/OTU_tables/Report/LYM_RTS_metadata.txt Group FilteredwithP_F_A Proteobacteria,Firmicutes,Actinobacteria
		"
	exit 0
else
	echo "################
	Running: sh $0 $1 $2 $3 $4 $5"
fi



if [ -z ${4+x} ]; then 
        echo "output directory prefix is unset yet but will be auto-valued"
        output_prefix="NoPrefix"
        else
        echo "output directory prefix is set to '$4'"
fi

if [ -z ${5+x} ]; then 
	echo "taxa_filtered is unset yet but will be auto-valued"
	taxa_filtered="UNREALTAX"
	else 
	echo "taxa_filtered is set to '$5'"
fi

check_file() {
	echo "Checking file for $1 ..."
	file_name=$1
	if [ -f $file_name ]; then
		echo "File $file_name exists."
	else
		echo "File $file_name does not exist."
		exit
	fi
}

check_dir() {
	echo "Cheking the directory for $1 ..."
	dir_name=$1
	if [ ! -d "$dir_name" ]; then
		echo "Creating directory $dir_name"
		mkdir -pv $dir_name
	else
		echo "Directory $dir_name exists."
	fi
}

MAIN() {

	echo -e "\n#Enter Qiime2 enviroment"
	echo -e "\n#Checking files"
	check_file $otu_table
	check_file $mapping_file

	echo -e "\n#Setting up the directory structure"
	echo -e "\n#The output directory prefix is $output_prefix"
	sample_name=`echo $(basename $otu_table) | sed -e 's/\.txt$//'`
	output_dir=$(dirname $otu_table)/${sample_name}_Qiime2_output_${output_prefix}/
	check_dir $output_dir




	echo "##############################################################\n#Generate the figure for the percentage of annotated level"
	cp $otu_table ${output_dir}/
	perl ${SCRIPTPATH}/stat_otu_tab.pl -unif min $otu_table -prefix ${output_dir}/Relative/otu_table --even ${output_dir}/Relative/otu_table.even.txt -spestat ${output_dir}/Relative/classified_stat_relative.xls
	perl ${SCRIPTPATH}/stat_otu_tab.pl -unif min $otu_table -prefix ${output_dir}/Absolute/otu_table -nomat -abs -spestat exported/Absolute/classified_stat.xls
	perl ${SCRIPTPATH}/bar_diagram.pl -table ${output_dir}/Relative/classified_stat_relative.xls -style 1 -x_title "Sample Name" -y_title "Sequence Number Percent" -right -textup -rotate='-45' --y_mun 1,7 > ${output_dir}/Relative/Classified_stat_relative.svg



	echo -e "\n#Convert OTU table to biom format"
	biom convert -i $otu_table -o ${output_dir}/${sample_name}.taxonomy.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy
	#biom convert -i $otu_table -o ${output_dir}/${sample_name}.taxonomy.biom --table-type="OTU table" --process-obs-metadata taxonomy

	source activate qiime2-2019.1
	echo -e "\n#Generate Qiime2 artifacts"
	qiime tools import   --input-path ${output_dir}/${sample_name}.taxonomy.biom --type 'FeatureTable[Frequency]'   --input-format BIOMV210Format   --output-path ${output_dir}/${sample_name}.count.qza
	qiime tools import   --input-path ${output_dir}/${sample_name}.taxonomy.biom --type "FeatureData[Taxonomy]"   --input-format BIOMV210Format   --output-path ${output_dir}/${sample_name}.taxonomy.qza

	echo -e "\n#Filter OTU table by taxonomy"
	qiime taxa filter-table   --i-table ${output_dir}/${sample_name}.count.qza --i-taxonomy ${output_dir}/${sample_name}.taxonomy.qza --p-exclude $taxa_filtered --o-filtered-table ${output_dir}/${sample_name}.count.filtered.temp.qza
	qiime feature-table filter-features --i-table ${output_dir}/${sample_name}.count.filtered.temp.qza --p-min-frequency 10 --o-filtered-table ${output_dir}/${sample_name}.count.filtered.qza


	echo -e "\n#Generate barplot"
	qiime taxa barplot --i-table ${output_dir}/${sample_name}.count.filtered.qza --i-taxonomy ${output_dir}/${sample_name}.taxonomy.qza  --m-metadata-file $mapping_file --o-visualization ${output_dir}/${sample_name}.taxa-bar-plots.qzv
	qiime feature-table filter-features --i-table ${output_dir}/${sample_name}.count.filtered.qza   --p-min-frequency $frequency  --o-filtered-table ${output_dir}/${sample_name}.count.filtered.${frequency}.qza
	qiime taxa barplot --i-table ${output_dir}/${sample_name}.count.filtered.${frequency}.qza --i-taxonomy ${output_dir}/${sample_name}.taxonomy.qza  --m-metadata-file $mapping_file --o-visualization ${output_dir}/${sample_name}.taxa-bar-plots.${frequency}.qzv

	echo -e "Conduct non-phylogenetic diversity analysis"
	qiime diversity core-metrics --i-table ${output_dir}/${sample_name}.count.qza  --p-sampling-depth $depth --m-metadata-file $mapping_file  --output-dir ${output_dir}/core-metrics
	

	echo -e "\n#Generate the heatmaps with the OTU (>= $frequency read) at different levels after collapsing."

	
	mkdir ${output_dir}/collapsed
	#mkdir ${output_dir}/Heatmap
	for n in 2 3 4 5 6 7; 
		do echo $n; qiime taxa collapse   --i-table ${output_dir}/${sample_name}.count.filtered.qza  --i-taxonomy ${output_dir}/${sample_name}.taxonomy.qza  --p-level $n  --o-collapsed-table ${output_dir}/collapsed/${sample_name}-l${n}.qza;
		#qiime feature-table filter-features   --i-table ${output_dir}/collapsed/${sample_name}-l${n}.qza   --p-min-frequency $frequency  --o-filtered-table ${output_dir}/collapsed/${sample_name}-l${n}.${frequency}.qza;
		#for category_1 in $category_set;
		#	do echo $category_1;
		#	Rscript ${SCRIPTPATH}/clean_na_of_inputs.R -m $mapping_file --group $category_1 -t ${output_dir}/collapsed/${sample_name}-l${n}.${frequency}.qza -o media_files;
		#   This heatmap was deprecated as various problems reflect by clients
		#	qiime feature-table heatmap --i-table media_files/filtered_feature_table.qza --m-metadata-file $mapping_file --m-metadata-column $category_1 --o-visualization ${output_dir}/Heatmap/${category_1}-${tax_levels[${n}]}-${frequency}-heatmap.qzv;
		#done;
	done;

	for n in 2 3 4 5 6 7;
		do echo $n; 
		for category_1 in $category_set;
			do echo $category_1;
			Rscript ${SCRIPTPATH}/abundance_heatmap.R  -m $mapping_file -c $category_1 -n 20 -i ${output_dir}/Absolute/otu_table.${n}.absolute.txt -o ${output_dir}/Heatmap/Heatmap_top20/${n}/ -l T -t F;
			Rscript ${SCRIPTPATH}/abundance_heatmap.R -m $mapping_file -c $category_1  -n 20 -i ${output_dir}/Absolute/otu_table.${n}.absolute.txt -o ${output_dir}/Heatmap/Heatmap_top20_clustered/${n}/ -l T -t F -u T;
			Rscript ${SCRIPTPATH}/abundance_heatmap.R  -m $mapping_file -c $category_1 -n 20 -i ${output_dir}/Absolute/otu_table.${n}.absolute.txt -o ${output_dir}/Heatmap/Heatmap_top20/${n}/ -b T -l T -p 'group_mean_' -t T;
		done;
	done;



	echo "ANCOM analaysis for differential OTU"
	mkdir ${output_dir}/ANCOM
	for n2 in 2 3 4 5 6 7;
		do echo $n2;
		for category_1 in $category_set;
			do echo $category_1;
				Rscript ${SCRIPTPATH}/clean_na_of_inputs.R -m $mapping_file --group $category_1 -t ${output_dir}/collapsed/${sample_name}-l${n}.qza -o media_files
				qiime composition add-pseudocount   --i-table media_files/filtered_feature_table.qza  --o-composition-table ${output_dir}/ANCOM/composition.${tax_levels[${n2}]}.qza;
				qiime composition ancom  --i-table ${output_dir}/ANCOM/composition.${tax_levels[${n2}]}.qza --m-metadata-file media_files/cleaned_map.txt --m-metadata-column $category_1 --o-visualization ${output_dir}/ANCOM/${category_1}-ANCOM-${tax_levels[${n2}]}.qzv;
			done;
			#qiime composition ancom  --i-table ${output_dir}/ANCOM/composition.${tax_levels[${n2}]}.qza --m-metadata-file $mapping_file --m-metadata-column $category_2 --o-visualization ${output_dir}/ANCOM/SecondaryGroup/ANCOM.${tax_levels[${n2}]}.qzv;
	done;


	echo -e "\n############################################Converting qzv files to html"
	for f in $(find ${output_dir}/ -type f -name "*.qzv"); do echo $f; base=$(basename $f .qzv); dir=$(dirname $f); new=${dir}/${base}; qiime tools export --input-path $f --output-path ${new}.qzv.exported; done
	for f in $(find ${output_dir}/ -type d -name "*qzv.exported"); do echo $f; base=$(basename $f .qzv.exported); dir=$(dirname $f); mv $f ${dir}/${base}; done
	for f in $(find ${output_dir}/ -type f -name "index.html") ; do echo $f; base=$(basename $f .html); dir=$(dirname $f); new=${dir}/Summary_请点此文件查看.html; mv $f $new; done	


	echo -e "\n############################################Generate the distance matrix and visulization artifact (rarefraction depth = $depth)."
	#qiime diversity core-metrics --i-table ${output_dir}/${sample_name}.count.filtered.qza --p-sampling-depth $depth --m-metadata-file $mapping_file --output-dir ${output_dir}/diversity
	echo -e "\n############################################Exit Qiime2 enviroment"
	conda deactivate

	mv ${output_dir}/Relative/otu_table.p.relative.mat ${output_dir}/Relative/otu_table.Phylum.relative.txt
	mv ${output_dir}/Relative/otu_table.c.relative.mat ${output_dir}/Relative/otu_table.Class.relative.txt
	mv ${output_dir}/Relative/otu_table.o.relative.mat ${output_dir}/Relative/otu_table.Order.relative.txt
	mv ${output_dir}/Relative/otu_table.f.relative.mat ${output_dir}/Relative/otu_table.Family.relative.txt
	mv ${output_dir}/Relative/otu_table.g.relative.mat ${output_dir}/Relative/otu_table.Genus.relative.txt
	mv ${output_dir}/Relative/otu_table.s.relative.mat ${output_dir}/Relative/otu_table.Species.relative.txt

	mv ${output_dir}/Absolute/otu_table.p.absolute.mat ${output_dir}/Absolute/otu_table.Phylum.absolute.txt
	mv ${output_dir}/Absolute/otu_table.c.absolute.mat ${output_dir}/Absolute/otu_table.Class.absolute.txt
	mv ${output_dir}/Absolute/otu_table.o.absolute.mat ${output_dir}/Absolute/otu_table.Order.absolute.txt
	mv ${output_dir}/Absolute/otu_table.f.absolute.mat ${output_dir}/Absolute/otu_table.Family.absolute.txt
	mv ${output_dir}/Absolute/otu_table.g.absolute.mat ${output_dir}/Absolute/otu_table.Genus.absolute.txt
	mv ${output_dir}/Absolute/otu_table.s.absolute.mat ${output_dir}/Absolute/otu_table.Species.absolute.txt


	echo "##############################################################\n#Barplot according to group mean"
	for n7 in "Phylum" "Class" "Order" "Family" "Genus" "Species"; 
		do echo $n7; 
		perl -lane '$,="\t";pop(@F);print(@F)' ${output_dir}/Relative/otu_table.${n7}.relative.txt > ${output_dir}/Relative/otu_table.${n7}.relative.lastcolumn.txt; 
		perl ${SCRIPTPATH}/get_table_head2.pl ${output_dir}/Relative/otu_table.${n7}.relative.lastcolumn.txt 20 -trantab > ${output_dir}/Relative/otu_table.${n7}.relative.lastcolumn.trans; 
		perl ${SCRIPTPATH}/bar_diagram.pl -table ${output_dir}/Relative/otu_table.${n7}.relative.lastcolumn.trans -style 1 -x_title "Sample Name" -y_title "Sequence Number Percent (%)" -right -textup -rotate='-45' --y_mun 0.2,5 --micro_scale --percentage > ${output_dir}/Relative/otu_table.${n7}.relative.svg
	done;

	for svg_file in ${output_dir}/Relative/*svg; do echo $svg_file; n=$(basename "$svg_file" .svg); echo $n; rsvg-convert -h 3200 -b white $svg_file > ${output_dir}/Relative/${n}.png; done;


	for category_1 in $category_set;
	do echo $category_1;
		for n7 in "Phylum" "Class" "Order" "Family" "Genus" "Species"; 
			do echo $n7;
			/usr/bin/Rscript ${SCRIPTPATH}/abundance_barplot.R -n 20 -m $mapping_file -c $category_1 -i ${output_dir}/Relative/otu_table.${n7}.relative.txt -o ${output_dir}/taxa-bar-plots-top20-group-ordered/ -p ${n7}_${category_1}_ordered_ -b F;
			/usr/bin/Rscript ${SCRIPTPATH}/abundance_barplot.R -n 20 -m $mapping_file -c $category_1 -i ${output_dir}/Relative/otu_table.${n7}.relative.txt -o ${output_dir}/Barplot-of-Group-Mean/ -p ${category_1}_${n7}_mean_ -b T;
		done;
	done;




	echo "#####################Run Lefse for taxa abundance";
	source activate lefse
	cd ${output_dir}/
	mkdir Lefse/
	for n7 in "Phylum" "Class" "Order" "Family" "Genus" "Species";
		do echo $n7;
			mkdir Lefse/${n7}
			cp ${output_dir}/Relative/otu_table.${n7}.relative.txt Lefse/${n7}/
			cd Lefse/${n7}
			for category_1 in $category_set;
				do echo $category_1;
					Rscript ${SCRIPTPATH}/write_data_for_lefse.R -i  otu_table.${n7}.relative.txt -m  $mapping_file -c  $category_1 -o  ${category_1}_${n7}_lefse.txt -u l;
					base="${category_1}_${n7}_lefse_LDA2"; lefse-format_input.py ${category_1}_${n7}_lefse.txt ${base}.lefseinput.txt -c 2 -u 1 -o 1000000; run_lefse.py ${base}.lefseinput.txt ${base}.LDA.txt -l 2;  
					lefse-plot_res.py  --max_feature_len 200 --orientation h --format pdf --left_space 0.3 --dpi 300 ${base}.LDA.txt ${base}.pdf; lefse-plot_cladogram.py ${base}.LDA.txt --dpi 300 ${base}.cladogram.pdf --clade_sep 1.8 --format pdf --right_space_prop 0.45 --label_font_size 10;
					base="${category_1}_${n7}_lefse_LDA4"; lefse-format_input.py ${category_1}_${n7}_lefse.txt ${base}.lefseinput.txt -c 2 -u 1 -o 1000000; run_lefse.py ${base}.lefseinput.txt ${base}.LDA.txt -l 4;  
					lefse-plot_res.py   --max_feature_len 200 --orientation h --format pdf --left_space 0.3 --dpi 300 ${base}.LDA.txt ${base}.pdf; lefse-plot_cladogram.py ${base}.LDA.txt --dpi 300 ${base}.cladogram.pdf --clade_sep 1.8 --format pdf --right_space_prop 0.45 --label_font_size 10;
				done;
			cd ../../
		done;
	conda deactivate


	echo "############################################Additional R related plot "
	sed 's/taxonomy/Consensus Lineage/' < $otu_table | sed 's/ConsensusLineage/Consensus Lineage/' > ${output_dir}/feature-table.ConsensusLineage.txt
	for category_1 in $category_set;
		do echo $category_1;
		Rscript ${SCRIPTPATH}/clean_na_of_inputs.R -m $mapping_file --group $category_1 -o media_files
		map="./media_files/cleaned_map.txt"
		Rscript ${SCRIPTPATH}/venn_and_flower_plot.R  -i $otu_table -m $mapping_file -c $category_1 -o ${output_dir}/VennAndFlower;
		Rscript ${SCRIPTPATH}/pcoa_and_nmds.R  -i ${output_dir}/feature-table.ConsensusLineage.txt -m $map -c $category_1 -o ${output_dir}/PCoA-NMDS;
		done;

}

MAIN;
