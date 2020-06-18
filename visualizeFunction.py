#!usr/bin/python3
import os
from circos import Circos
from visualize import Visualize
import re
from osEnv import OSEnv
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
            self.system('''
{R_path} {bayegy_home}/abundance_barplot.R -n 20 -m {mapping_file} -c {category} -i {abundance_table} -o {out_dir}1-Barplots/ -p {prefix}{category}_ -b F;
{R_path} {bayegy_home}/abundance_barplot.R -n 20 -m {mapping_file} -c {category} -i {abundance_table} -o {out_dir}1-Barplots/ -p {prefix}{category}_groupMean_ -b T;
{R_path} {bayegy_home}/abundance_heatmap.R  -m {mapping_file} -c {category} -n 20 -i {abundance_table} -o {out_dir}2-Heatmaps -p {prefix}{category}_nocluster_ -l F -t F;
{R_path} {bayegy_home}/abundance_heatmap.R -m {mapping_file} -c {category}  -n 20 -i {abundance_table} -o {out_dir}2-Heatmaps -p {prefix}{category}_clustered_ -l F -t F -u T;
{R_path} {bayegy_home}/abundance_heatmap.R  -m {mapping_file} -c {category} -n 20 -i {abundance_table} -o {out_dir}2-Heatmaps -p {prefix}{category}_groupMean_ -l F -t T -b T;
{R_path} {bayegy_home}/write_data_for_lefse.R -i  {abundance_table} -m  {mapping_file} -c  {category} -o  {out_dir}4-SignificanceAnalysis/LEfSe/{prefix}{category}_lefse.txt -u f;
{R_path} {bayegy_home}/Function_DunnTest.r -i  {abundance_table} -m  {mapping_file} -c  {category} -o  {out_dir}4-SignificanceAnalysis/DunnTest -p {prefix};
            ''', category=g)

            with OSEnv(pythonpath=self.lefse_pylib_home, path=self.lefse_py_home, r_libs=self.lefse_rlib_home):
                lefse_base = '{}4-SignificanceAnalysis/LEfSe/{}{}_lefse_LDA'.format(
                    self.out_dir, self.prefix, g)
                self.system("""
{bayegy_home}/mod_format_input.py {out_dir}4-SignificanceAnalysis/LEfSe/{prefix}{category}_lefse.txt {base}2.lefseinput.txt -c 2 -u 1 -o 1000000 && \
{bayegy_home}/mod_run_lefse.py {base}2.lefseinput.txt {base}2.LDA.txt -l 2 -a 0.05 -w 0.05
                """, base=lefse_base, category=g)
            # if not exclude == 'all':
                lda2_result = "{base}2.LDA.txt".format(base=lefse_base)
                if os.path.exists(lda2_result):
                    self.system("""
{bayegy_home}/mod_lefse-plot_res.py --map {mapping_file} --category {category}  --max_feature_len 200 --orientation h --format pdf --left_space 0.3 --dpi 300 {base}2.LDA.txt {base}2.pdf;
{base_dir}/lda22ldamt.py {base}2.LDA.txt {base}4.LDA.txt 4
{bayegy_home}/mod_lefse-plot_res.py --map {mapping_file} --category {category}  --max_feature_len 200 --orientation h --format pdf --left_space 0.3 --dpi 300 {base}4.LDA.txt {base}4.pdf;
                """, base=lefse_base, category=g)
                else:
                    with open(lda2_result, 'w') as lefout, open(self.abundance_table, 'r') as abcin:
                        for num, line in enumerate(abcin):
                            if num > 0:
                                lefout.write(
                                    "{}\t-\t\t\t-\n".format(line.split()[0]))

            for el in exclude_list:
                el, *aprefix = el.split(":") + [""]
                # info = re.sub(r'\(|\)|%|\\|\/', "", el.replace(',', '_')) + '_excluded_'
                # info = info if info < 50 else hash(info)
                prefix = self.prefix + "".join(aprefix)
                self.system('''
if [ ! -f {out_dir}5-CorrelationAnalysis/CorrelationHeatmap/{prefix}spearman_rank_correlation_matrix.txt ]; then {R_path} {bayegy_home}/cor_heatmap.R -i {abundance_table} -o {out_dir}5-CorrelationAnalysis/CorrelationHeatmap -n 25 -m {mapping_file} -e {eld} -p "{prefix}";fi;
{R_path} {bayegy_home}/RDA.R -i {abundance_table} -m {mapping_file} -c {category} -o {out_dir}5-CorrelationAnalysis/RDA -n 25 -e {eld} -p "{prefix}";
''', category=g, prefix=prefix, eld=el)
            # os.system("bash {}".format(self.running_bash))

    def __visualize_without_group__(self):
        Circos(self.abundance_table, prefix=self.prefix,
               out_dir=self.out_dir + 'Circos').visualize()
        self.system('''
{R_path} {bayegy_home}/abundance_barplot.R -n 20 -i {abundance_table} -o {out_dir}1-Barplots/ -p {prefix} -b F;
{R_path} {bayegy_home}/abundance_heatmap.R -n 20 -i {abundance_table} -o {out_dir}2-Heatmaps  -p {prefix} -l F -t F -u T;
''')
