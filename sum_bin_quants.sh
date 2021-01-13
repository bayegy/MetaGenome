SOFT=$1
out=$2
bin_folder=$3
assembly=$4

comm () { ${SOFT}/print_comment.py "$1" "-"; }
error () { ${SOFT}/print_comment.py "$1" "*"; exit 1; }
warning () { ${SOFT}/print_comment.py "$1" "*"; }
announcement () { ${SOFT}/print_comment.py "$1" "#"; }


comm "summarize salmon files..."
home=$(pwd)
cd ${out}/alignment_files/
${SOFT}/summarize_salmon_files.py
if [[ $? -ne 0 ]]; then error "something went wrong with summarizing salmon output with python script! Exiting..."; fi
cd $home
mkdir ${out}/quant_files
for f in $(ls ${out}/alignment_files/ | grep .quant.counts); do mv ${out}/alignment_files/$f ${out}/quant_files/; done



########################################################################################################
########################        EXTRACTING AVERAGE ABUNDANCE OF EACH BIN        ########################
########################################################################################################
announcement "EXTRACTING AVERAGE ABUNDANCE OF EACH BIN"

n=$(ls ${out}/quant_files/ | grep counts | wc -l)
if [[ $n -lt 1 ]]; then error "There were no files found in ${out}/quant_files/"; fi
comm "There were $n samples detected. Making abundance table!"

${SOFT}/split_salmon_out_into_bins.py ${out}/quant_files/ $bin_folder $assembly > ${out}/bin_abundance_table.tab
if [[ $? -ne 0 ]]; then error "something went wrong with making summary abundance table. Exiting..."; fi
comm "Average bin abundance table stored in ${out}/abundance_table.tab"




########################################################################################################
########################            MAKING GENOME ABUNDANCE HEATMAP             ########################
########################################################################################################
### ugly heatmap
# if [[ $n -gt 1 ]]; then
#     announcement "MAKING GENOME ABUNDANCE HEATMAP WITH SEABORN"

#     comm "making heatmap with Seaborn"
#     ${SOFT}/make_heatmap.py ${out}/bin_abundance_table.tab ${out}/bin_abundance_heatmap.png
#     if [[ $? -ne 0 ]]; then error "something went wrong with making the heatmap. Exiting..."; fi

#     comm "cleaning up..."
#     rm -r ${out}/alignment_files/ 
# else
#     warning "Cannot make clustered heatmap with just one sample... Skipping heatmap"
# fi

########################################################################################################
########################     QUANT_BINS PIPELINE SUCCESSFULLY FINISHED!!!       ########################
########################################################################################################
announcement "QUANT_BINS PIPELINE SUCCESSFULLY FINISHED!!!"

