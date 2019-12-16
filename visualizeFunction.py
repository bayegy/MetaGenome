#!usr/bin/python3
import os
from circos import Circos
from visualize import Visualize
import re
# import pdb


class VisualizeFunction(Visualize):
    """
    Sample usage:
        from visualizeFunction import VisualizeFunction
        v=VisualizeFunction('test/otu_table.Genus.absolute.txt','test/mapping_file.txt','Group1')
    """

    def __visualize_with_group__(self, exclude='all'):
        categories = [g.strip() for g in re.split(',', self.categories)]
        exclude_list = [] if exclude == "all" else exclude.split(';')
        for g in categories:
            print(g)
            Circos(self.abundance_table, mapping_file=self.mapping_file, category=g, by_group_mean=False,
                   prefix=self.prefix + g + '_', out_dir=self.out_dir + '3-Circos').visualize()
            # Circos(self.abundance_table, mapping_file=self.mapping_file, category=g, by_group_mean=True,
            #        prefix=self.prefix + g + '_groupMean_', out_dir=self.out_dir + 'Circos').visualize()
            os.system('''
Rscript {SCRIPTPATH}/abundance_barplot.R -n 20 -m {mapping_file} -c {category} -i {abundance_table} -o {outdir}1-Barplots/ -p {prefix}{category}_ -b F;
Rscript {SCRIPTPATH}/abundance_barplot.R -n 20 -m {mapping_file} -c {category} -i {abundance_table} -o {outdir}1-Barplots/ -p {prefix}{category}_groupMean_ -b T;
Rscript {SCRIPTPATH}/abundance_heatmap.R  -m {mapping_file} -c {category} -n 20 -i {abundance_table} -o {outdir}2-Heatmaps -p {prefix}{category}_nocluster_ -l F -t F;
Rscript {SCRIPTPATH}/abundance_heatmap.R -m {mapping_file} -c {category}  -n 20 -i {abundance_table} -o {outdir}2-Heatmaps -p {prefix}{category}_clustered_ -l F -t F -u T;
Rscript {SCRIPTPATH}/abundance_heatmap.R  -m {mapping_file} -c {category} -n 20 -i {abundance_table} -o {outdir}2-Heatmaps -p {prefix}{category}_groupMean_ -l F -t T -b T;
Rscript {SCRIPTPATH}/write_data_for_lefse.R -i  {abundance_table} -m  {mapping_file} -c  {category} -o  {outdir}4-SignificanceAnalysis/LEfSe/{prefix}{category}_lefse.txt -u f;
Rscript {SCRIPTPATH}/Function_DunnTest.r -i  {abundance_table} -m  {mapping_file} -c  {category} -o  {outdir}4-SignificanceAnalysis/DunnTest -p {prefix};
source lefse;
{SCRIPTPATH}/mod_format_input.py {outdir}4-SignificanceAnalysis/LEfSe/{prefix}{category}_lefse.txt {base}2.lefseinput.txt -c 2 -u 1 -o 1000000; run_lefse.py {base}2.lefseinput.txt {base}2.LDA.txt -l 2 -a 0.03 -w 0.03;
{SCRIPTPATH}/mod_lefse-plot_res.py --map {mapping_file} --category {category}  --max_feature_len 200 --orientation h --format pdf --left_space 0.3 --dpi 300 {base}2.LDA.txt {base}2.pdf;
#{SCRIPTPATH}/mod_lefse-plot_cladogram.py {base}2.LDA.txt  --map {mapping_file} --category {category} --dpi 300 {base}2.cladogram.pdf --clade_sep 1.8 --format pdf --right_space_prop 0.45 --label_font_size 10;
#{SCRIPTPATH}/mod_format_input.py {outdir}4-SignificanceAnalysis/LEfSe/{prefix}{category}_lefse.txt {base}4.lefseinput.txt -c 2 -u 1 -o 1000000; run_lefse.py {base}4.lefseinput.txt {base}4.LDA.txt -l 4 -a 0.03 -w 0.03;
lda22ldamt.py {base}2.LDA.txt {base}4.LDA.txt 4
{SCRIPTPATH}/mod_lefse-plot_res.py --map {mapping_file} --category {category}  --max_feature_len 200 --orientation h --format pdf --left_space 0.3 --dpi 300 {base}4.LDA.txt {base}4.pdf;
#{SCRIPTPATH}/mod_lefse-plot_cladogram.py {base}4.LDA.txt --map {mapping_file} --category {category} --dpi 300 {base}4.cladogram.pdf --clade_sep 1.8 --format pdf --right_space_prop 0.45 --label_font_size 10;
source delefse'''.format(base='{}4-SignificanceAnalysis/LEfSe/{}{}_lefse_LDA'.format(self.out_dir, self.prefix, g), SCRIPTPATH=self.path['bayegy_home'], abundance_table=self.abundance_table, mapping_file=self.mapping_file, category=g, outdir=self.out_dir, prefix=self.prefix))
            # if not exclude == 'all':
            for el in exclude_list:
                prefix = self.prefix if el == 'none' else self.prefix + \
                    re.sub('\(|\)|%|\\|\/', "", el.replace(',', '_')) + '_excluded_'
                os.system('''
if [ ! -f {outdir}5-CorrelationAnalysis/CorrelationHeatmap/{prefix}spearman_rank_correlation_matrix.txt ]; then Rscript {SCRIPTPATH}/cor_heatmap.R -i {abundance_table} -o {outdir}5-CorrelationAnalysis/CorrelationHeatmap -n 25 -m {mapping_file} -e {eld} -p "{prefix}";fi;
Rscript {SCRIPTPATH}/RDA.R -i {abundance_table} -m {mapping_file} -c {category} -o {outdir}5-CorrelationAnalysis/RDA -n 25 -e {eld} -p "{prefix}";
'''.format(SCRIPTPATH=self.path['bayegy_home'], abundance_table=self.abundance_table, mapping_file=self.mapping_file, category=g, outdir=self.out_dir, prefix=prefix, eld=el))
            # os.system("bash {}".format(self.running_bash))

    def __visualize_without_group__(self):
        Circos(self.abundance_table, prefix=self.prefix, out_dir=self.out_dir + 'Circos').visualize()
        os.system('''
Rscript {SCRIPTPATH}/abundance_barplot.R -n 20 -i {abundance_table} -o {outdir}1-Barplots/ -p {prefix} -b F;
Rscript {SCRIPTPATH}/abundance_heatmap.R -n 20 -i {abundance_table} -o {outdir}2-Heatmaps  -p {prefix} -l F -t F -u T;
'''.format(SCRIPTPATH=self.path['bayegy_home'], abundance_table=self.abundance_table, outdir=self.out_dir, prefix=self.prefix))
