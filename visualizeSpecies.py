import os
from visualize import Visualize
from osEnv import OSEnv


class VisualizeSpecies(Visualize):
    """docstring for VisualizeSpecies"""

    def __init__(self, abundance_table, mapping_file, categories=False, prefix=False, out_dir=False, filter_species="exclude:Environmentalsamples", tmp_dir='./'):
        super(VisualizeSpecies, self).__init__(
            abundance_table, mapping_file, categories, prefix, out_dir, filter_abc_description=filter_species)

        self.set_attr(
            sample_name=os.path.basename(self.abundance_table).rstrip('.txt'),
            # exclude_species=exclude_species
        )

    def __visualize_with_group__(self, exclude='none'):
        with OSEnv(path=self.qiime2_home + 'bin/', pythonpath=self.qiime2_home + 'lib/python3.6/site-packages/'):
            self.system("""
    echo "Check wheather your categories are the following:"
    for i in {category_set};do echo $i;done

    not_rda='{exclude}'
    echo "Check CorrelationAnalysis config:"
    for i in $not_rda;do echo $i;done

    tax_levels["1"]="Kingdom"
    tax_levels["2"]="Phylum"
    tax_levels["3"]="Class"
    tax_levels["4"]="Order"
    tax_levels["5"]="Family"
    tax_levels["6"]="Genus"
    tax_levels["7"]="Species"

    echo "##############################################################\n#Generate the figure for the percentage of annotated level"
    {perl_path} {bayegy_home}/stat_otu_tab.pl -unif min {abundance_table} -prefix {tmp_dir}/Relative/otu_table --even {tmp_dir}/Relative/otu_table.even.txt -spestat {tmp_dir}/Relative/classified_stat_relative.xls
    {perl_path} {bayegy_home}/stat_otu_tab.pl -unif min {abundance_table} -prefix {tmp_dir}/Absolute/otu_table -nomat -abs -spestat exported/Absolute/classified_stat.xls
    {perl_path} {bayegy_home}/bar_diagram.pl -table {tmp_dir}/Relative/classified_stat_relative.xls -style 1 -x_title "Sample Name" -y_title "Sequence Number Percent" -right -textup -rotate='-45' --y_mun 1,7 > {tmp_dir}/Relative/Classified_stat_relative.svg

    # source {base_dir}/path/activate_qiime2.sh

    echo -e "\n#Convert OTU table to biom format"
    biom convert -i {abundance_table} -o {tmp_dir}/{sample_name}.taxonomy.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy


    echo -e "\n#Generate Qiime2 artifacts"
    qiime tools import   --input-path {tmp_dir}/{sample_name}.taxonomy.biom --type 'FeatureTable[Frequency]'   --input-format BIOMV210Format   --output-path {tmp_dir}/{sample_name}.count.qza
    qiime tools import   --input-path {tmp_dir}/{sample_name}.taxonomy.biom --type "FeatureData[Taxonomy]"   --input-format BIOMV210Format   --output-path {tmp_dir}/{sample_name}.taxonomy.qza

    echo -e "\n#Filter OTU table by taxonomy"
    qiime feature-table filter-features --i-table {tmp_dir}/{sample_name}.count.qza --p-min-frequency 10 --o-filtered-table {tmp_dir}/{sample_name}.count.filtered.qza

    echo -e "\n#Generate barplot"
    qiime taxa barplot --i-table {tmp_dir}/{sample_name}.count.filtered.qza --i-taxonomy {tmp_dir}/{sample_name}.taxonomy.qza  --m-metadata-file {mapping_file} --o-visualization {tmp_dir}/{sample_name}.taxa-bar-plots.qzv
    qiime feature-table filter-features --i-table {tmp_dir}/{sample_name}.count.filtered.qza   --p-min-frequency {frequency}  --o-filtered-table {tmp_dir}/{sample_name}.count.filtered.{frequency}.qza
    qiime taxa barplot --i-table {tmp_dir}/{sample_name}.count.filtered.{frequency}.qza --i-taxonomy {tmp_dir}/{sample_name}.taxonomy.qza  --m-metadata-file {mapping_file} --o-visualization {tmp_dir}/{sample_name}.taxa-bar-plots.{frequency}.qzv

    echo -e "Conduct non-phylogenetic diversity analysis"
    depth=$({R_path} {bayegy_home}/min.R {abundance_table})
    echo "depth is $depth"
    qiime diversity core-metrics --i-table {tmp_dir}/{sample_name}.count.qza  --p-sampling-depth $depth --m-metadata-file {mapping_file}  --output-dir {tmp_dir}/core-metrics


    mkdir {tmp_dir}/collapsed
    for n2 in 2 3 4 5 6 7;
        do echo $n2; qiime taxa collapse   --i-table {tmp_dir}/{sample_name}.count.filtered.qza  --i-taxonomy {tmp_dir}/{sample_name}.taxonomy.qza  --p-level $n2  --o-collapsed-table {tmp_dir}/collapsed/{sample_name}-${{tax_levels[${{n2}}]}}.qza;
    done;


    echo "ANCOM analaysis for differential OTU"
    mkdir {tmp_dir}/ANCOM
    for n2 in "Phylum" "Class" "Order" "Family" "Genus" "Species";
        do echo $n2;
        for category_1 in {category_set};
            do echo $category_1;
                {R_path} {bayegy_home}/clean_na_of_inputs.R -m {mapping_file} --group $category_1 -t {tmp_dir}/collapsed/{sample_name}-${{n2}}.qza -o {tmp_dir}
                qiime composition add-pseudocount   --i-table {tmp_dir}/filtered_feature_table.qza  --o-composition-table {tmp_dir}/ANCOM/composition.${{n2}}.qza;
                qiime composition ancom  --i-table {tmp_dir}/ANCOM/composition.${{n2}}.qza --m-metadata-file {tmp_dir}/cleaned_map.txt --m-metadata-column $category_1 --o-visualization {tmp_dir}/ANCOM/${{category_1}}-ANCOM-${{n2}}.qzv;
            done;
            #qiime composition ancom  --i-table {tmp_dir}/ANCOM/composition.${{n2}}.qza --m-metadata-file {mapping_file} --m-metadata-column $category_2 --o-visualization {tmp_dir}/ANCOM/SecondaryGroup/ANCOM.${{n2}}.qzv;
    done;


    echo -e "\n############################################Converting qzv files to html"
    for f in $(find {tmp_dir}/ -type f -name "*.qzv"); do echo $f; base=$(basename $f .qzv); dir=$(dirname $f); new=${{dir}}/${{base}}; qiime tools export --input-path $f --output-path ${{new}}.qzv.exported; done
    for f in $(find {tmp_dir}/ -type d -name "*qzv.exported"); do echo $f; base=$(basename $f .qzv.exported); dir=$(dirname $f); mv $f ${{dir}}/${{base}}; done
    for f in $(find {tmp_dir}/ -type f -name "index.html") ; do echo $f; base=$(basename $f .html); dir=$(dirname $f); new=${{dir}}/Summary_请点此文件查看.html; mv $f $new; done



    mv {tmp_dir}/Relative/otu_table.p.relative.mat {tmp_dir}/Relative/otu_table.Phylum.relative.txt
    mv {tmp_dir}/Relative/otu_table.c.relative.mat {tmp_dir}/Relative/otu_table.Class.relative.txt
    mv {tmp_dir}/Relative/otu_table.o.relative.mat {tmp_dir}/Relative/otu_table.Order.relative.txt
    mv {tmp_dir}/Relative/otu_table.f.relative.mat {tmp_dir}/Relative/otu_table.Family.relative.txt
    mv {tmp_dir}/Relative/otu_table.g.relative.mat {tmp_dir}/Relative/otu_table.Genus.relative.txt
    mv {tmp_dir}/Relative/otu_table.s.relative.mat {tmp_dir}/Relative/otu_table.Species.relative.txt

    mv {tmp_dir}/Absolute/otu_table.p.absolute.mat {tmp_dir}/Absolute/otu_table.Phylum.absolute.txt
    mv {tmp_dir}/Absolute/otu_table.c.absolute.mat {tmp_dir}/Absolute/otu_table.Class.absolute.txt
    mv {tmp_dir}/Absolute/otu_table.o.absolute.mat {tmp_dir}/Absolute/otu_table.Order.absolute.txt
    mv {tmp_dir}/Absolute/otu_table.f.absolute.mat {tmp_dir}/Absolute/otu_table.Family.absolute.txt
    mv {tmp_dir}/Absolute/otu_table.g.absolute.mat {tmp_dir}/Absolute/otu_table.Genus.absolute.txt
    mv {tmp_dir}/Absolute/otu_table.s.absolute.mat {tmp_dir}/Absolute/otu_table.Species.absolute.txt

    # source {base_dir}/path/deactivate_qiime2.sh


    for n in "Phylum" "Class" "Order" "Family" "Genus" "Species";
        do echo $n;
        for category_1 in {category_set};
            do echo $category_1;
            {R_path} {bayegy_home}/abundance_heatmap.R  -m {mapping_file} -c $category_1 -n 20 -i {tmp_dir}/Absolute/otu_table.${{n}}.absolute.txt -o {tmp_dir}/Heatmap/Heatmap_top20/${{n}}/ -l T -t F -p "${{category_1}}_";
            {R_path} {bayegy_home}/abundance_heatmap.R -m {mapping_file} -c $category_1  -n 20 -i {tmp_dir}/Absolute/otu_table.${{n}}.absolute.txt -o {tmp_dir}/Heatmap/Heatmap_top20_clustered/${{n}}/ -l T -t F -u T -p "${{category_1}}_";
            {R_path} {bayegy_home}/abundance_heatmap.R  -m {mapping_file} -c $category_1 -n 20 -i {tmp_dir}/Absolute/otu_table.${{n}}.absolute.txt -o {tmp_dir}/Heatmap/Heatmap_top20/${{n}}/ -b T -l T -p "${{category_1}}_group_mean_" -t T;
        done;
    done;



    test=${{not_rda// */}}
    if [ ! $test == "all" ];then
        echo "##############################################################\nCorrelation heatmap analysis"
        for nrda in $not_rda;
            do echo $nrda;
            prefix=${{nrda//,/_}}_excluded_;
            prefix=${{prefix//none_excluded_/}};
            prefix=${{prefix//\//-}};
            prefix=${{prefix//\\/-}};
            prefix=${{prefix//\(/}};
            prefix=${{prefix//\)/}};
            prefix=${{prefix//\%/}};
            for n7 in "Phylum" "Class" "Order" "Family" "Genus" "Species";
                do echo $n7;
                {R_path} {bayegy_home}/cor_heatmap.R -i {tmp_dir}/Relative/otu_table.${{n7}}.relative.txt -o {tmp_dir}/CorrelationAnalysis/CorrelationHeatmap/${{n7}}/ -n 25 -m {mapping_file} -e $nrda -p "$prefix";
                for category_1 in {category_set};do echo $category_1;{R_path} {bayegy_home}/RDA.R -i {tmp_dir}/Absolute/otu_table.${{n7}}.absolute.txt -m {mapping_file} -c $category_1 -o {tmp_dir}/CorrelationAnalysis/RDA/${{n7}} -n 25 -e $nrda -p "$prefix";done;
            done;
        done;
    fi;


    echo "##############################################################\n#Barplot according to group mean"
    for n7 in "Phylum" "Class" "Order" "Family" "Genus" "Species";
        do echo $n7;
        {perl_path} -lane '$,="\t";pop(@F);print(@F)' {tmp_dir}/Relative/otu_table.${{n7}}.relative.txt > {tmp_dir}/Relative/otu_table.${{n7}}.relative.lastcolumn.txt;
        {perl_path} {bayegy_home}/get_table_head2.pl {tmp_dir}/Relative/otu_table.${{n7}}.relative.lastcolumn.txt 20 -trantab > {tmp_dir}/Relative/otu_table.${{n7}}.relative.lastcolumn.trans;
        {perl_path} {bayegy_home}/bar_diagram.pl -table {tmp_dir}/Relative/otu_table.${{n7}}.relative.lastcolumn.trans -style 1 -x_title "Sample Name" -y_title "Sequence Number Percent (%)" -right -textup -rotate='-45' --y_mun 0.2,5 --micro_scale --percentage > {tmp_dir}/Relative/otu_table.${{n7}}.relative.svg
    done;

    for svg_file in {tmp_dir}/Relative/*svg; do echo $svg_file; n=$(basename "$svg_file" .svg); echo $n; convert $svg_file {tmp_dir}/Relative/${{n}}.png; done;


    for category_1 in {category_set};
    do echo $category_1;
        for n7 in "Phylum" "Class" "Order" "Family" "Genus" "Species";
            do echo $n7;
            {R_path} {bayegy_home}/abundance_barplot.R -n 20 -m {mapping_file} -c $category_1 -i {tmp_dir}/Relative/otu_table.${{n7}}.relative.txt -o {tmp_dir}/Taxa-bar-plots-top20/ -p ${{n7}}_${{category_1}}_ -b F;
            {R_path} {bayegy_home}/abundance_barplot.R -n 20 -m {mapping_file} -c $category_1 -i {tmp_dir}/Relative/otu_table.${{n7}}.relative.txt -o {tmp_dir}/Barplot-of-Group-Mean/ -p ${{category_1}}_${{n7}}_mean_ -b T;
            {R_path} {bayegy_home}/Function_DunnTest.r -m {mapping_file} -c $category_1 -i {tmp_dir}/Relative/otu_table.${{n7}}.relative.txt -o {tmp_dir}/DunnTest/ -p ${{n7}}_;
            {R_path} {bayegy_home}/write_data_for_lefse.R -i  {tmp_dir}/Relative/otu_table.${{n7}}.relative.txt -m  {mapping_file} -c  $category_1 -o  {tmp_dir}/Lefse/${{n7}}/${{category_1}}_${{n7}}_lefse.txt -u l;
        done;
    done;

    echo "############################################Additional R related plot "
    sed 's/taxonomy/Consensus Lineage/' < {abundance_table} | sed 's/ConsensusLineage/Consensus Lineage/' > {tmp_dir}/feature-table.ConsensusLineage.txt
    for category_1 in {category_set};
        do echo $category_1;
        {R_path} {bayegy_home}/clean_na_of_inputs.R -m {mapping_file} --group $category_1 -o {tmp_dir}
        map="{tmp_dir}/cleaned_map.txt"
        {R_path} {bayegy_home}/venn_and_flower_plot.R -s F  -i {abundance_table} -m {mapping_file} -c $category_1 -o {tmp_dir}/VennAndFlower;
        {R_path} {bayegy_home}/pcoa_and_nmds.R  -i {tmp_dir}/feature-table.ConsensusLineage.txt -m $map -c $category_1 -o {tmp_dir}/PCoA-NMDS;
        done;
        """, frequency=1000, category_set=self.categories.replace(',', ' '), exclude=exclude.replace(';', ' '))

        with OSEnv(pythonpath=self.lefse_pylib_home, path=self.lefse_py_home, r_libs=self.lefse_rlib_home):
            self.system("""
    echo "#####################Run Lefse for taxa abundance";
    for n7 in "Phylum" "Class" "Order" "Family" "Genus" "Species";
        do echo $n7;
            for category_1 in {category_set};
                do echo $category_1;
                    base="{tmp_dir}/Lefse/${{n7}}/${{category_1}}_${{n7}}_lefse_LDA2"; lefse-format_input.py {tmp_dir}/Lefse/${{n7}}/${{category_1}}_${{n7}}_lefse.txt ${{base}}.lefseinput.txt -c 2 -u 1 -o 1000000; run_lefse.py ${{base}}.lefseinput.txt ${{base}}.LDA.txt -l 2;
                    {bayegy_home}/mod_lefse-plot_res.py --category $category_1 --map {mapping_file} --max_feature_len 200 --orientation h --format pdf --left_space 0.3 --dpi 300 ${{base}}.LDA.txt ${{base}}.pdf; {bayegy_home}/mod_lefse-plot_cladogram.py ${{base}}.LDA.txt --category $category_1 --map {mapping_file} --dpi 300 ${{base}}.cladogram.pdf --format pdf;
                    base4="{tmp_dir}/Lefse/${{n7}}/${{category_1}}_${{n7}}_lefse_LDA4";
                    # lefse-format_input.py ${{category_1}}_${{n7}}_lefse.txt ${{base}}.lefseinput.txt -c 2 -u 1 -o 1000000; run_lefse.py ${{base}}.lefseinput.txt ${{base}}.LDA.txt -l 4;
                    {base_dir}/lda22ldamt.py ${{base}}.LDA.txt ${{base4}}.LDA.txt 4;
                    {bayegy_home}/mod_lefse-plot_res.py --category $category_1 --map {mapping_file}  --max_feature_len 200 --orientation h --format pdf --left_space 0.3 --dpi 300 ${{base4}}.LDA.txt ${{base4}}.pdf; {bayegy_home}/mod_lefse-plot_cladogram.py ${{base4}}.LDA.txt --category $category_1 --map {mapping_file} --dpi 300 ${{base4}}.cladogram.pdf --format pdf;
                done;
        done;

    mkdir -p {out_dir}/1-AbundanceSummary/1-RelativeAbundance \
        {out_dir}/1-AbundanceSummary/2-Barplots \
        {out_dir}/1-AbundanceSummary/3-Heatmaps \
        {out_dir}/2-AbundanceComparison/VennAndFlower \
        {out_dir}/2-AbundanceComparison/ANCOM \
        {out_dir}/2-AbundanceComparison \
        {out_dir}/2-AbundanceComparison/LEfSe \
        {out_dir}/3-DiversityAnalysis \

    mv {tmp_dir}/Relative/Classified_stat_relative.png {out_dir}/1-AbundanceSummary/
    mv {tmp_dir}/Relative/*relative.txt {out_dir}/1-AbundanceSummary/1-RelativeAbundance/
    mv {tmp_dir}/All.Taxa.OTU.taxa-bar-plots {tmp_dir}/All.Taxa.OTU.taxa-bar-plots.1000  {tmp_dir}/Barplot-of-Group-Mean {tmp_dir}/Taxa-bar-plots-top20  {out_dir}/1-AbundanceSummary/2-Barplots/
    mv {tmp_dir}/Heatmap/* {out_dir}/1-AbundanceSummary/3-Heatmaps/
    mv {tmp_dir}/ANCOM/*ANCOM* {out_dir}/2-AbundanceComparison/ANCOM/
    mv {tmp_dir}/Lefse/* {out_dir}/2-AbundanceComparison/LEfSe/
    mv {tmp_dir}/DunnTest {out_dir}/2-AbundanceComparison/
    mv {tmp_dir}/VennAndFlower/*  {out_dir}/2-AbundanceComparison/VennAndFlower/
    mv {tmp_dir}/core-metrics/bray_curtis_emperor {tmp_dir}/PCoA-NMDS/* {out_dir}/3-DiversityAnalysis/
    mv {tmp_dir}/CorrelationAnalysis {out_dir}/4-CorrelationAnalysis
                """, category_set=self.categories.replace(',', ' '))

    def __visualize_without_group__(self):
        os.system("bash {}visualize_otu_table_without_group.sh {} {} none {} {} {} {}".format(
            self._base_dir, self.abundance_table, self.mapping_file, self.prefix, "NONE", self.path[
                'bayegy_home'], self.tmp_dir
        ))
