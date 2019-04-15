#!usr/bin/python3
import os
import re
from circos import Circos
import json


class VisualizeFunction(object):
    """docstring for VisualizaFunction"""

    def __init__(self, abundance_table, mapping_file=False, categories=False, prefix=False, out_dir=False):
        out_dir = out_dir if out_dir else os.path.dirname(abundance_table)
        self.out_dir = os.path.abspath(out_dir) + '/'
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        self._base_dir = os.path.dirname(__file__) + '/'
        with open(self._base_dir + "pipconfig/path.conf") as f:
            self.path = json.load(f)
        self.abundance_table = os.path.abspath(abundance_table)
        self.mapping_file = os.path.abspath(mapping_file)
        self.categories = [g.strip() for g in re.split(',', categories)]
        self.prefix = prefix if prefix else re.sub('\.[^\.]*$', '', os.path.basename(self.abundance_table)) + '_'
        self.running_bash = self.out_dir + '.visualize_function.sh'

    def visualiza_with_group(self):
        for g in self.categories:
            with open(self.running_bash, 'w') as shscript:
                print('''
SCRIPTPATH=%s
abundance_table=%s
mapping_file=%s
category=%s
outdir=%s
prefix=%s
Rscript ${SCRIPTPATH}/abundance_barplot.R -n 20 -m $mapping_file -c $category -i $abundance_table -o ${outdir}Barplots/ -p ${prefix}${category}_ordered_ -b F;
Rscript ${SCRIPTPATH}/abundance_barplot.R -n 20 -m $mapping_file -c $category -i $abundance_table -o ${outdir}Barplots/ -p ${prefix}${category}_mean_ -b T;
Rscript ${SCRIPTPATH}/abundance_heatmap.R  -m $mapping_file -c $category -n 20 -i $abundance_table -o ${outdir}Heatmaps -p ${prefix}${category}_nocluster_ -l F -t F;
Rscript ${SCRIPTPATH}/abundance_heatmap.R -m $mapping_file -c $category  -n 20 -i $abundance_table -o ${outdir}Heatmaps -p ${prefix}${category}_clustered_ -l F -t F -u T;
Rscript ${SCRIPTPATH}/abundance_heatmap.R  -m $mapping_file -c $category -n 20 -i $abundance_table -o ${outdir}Heatmaps -p ${prefix}${category}_groupMean_ -l F -t T -b T;
source lefse;
Rscript ${SCRIPTPATH}/write_data_for_lefse.R -i  $abundance_table -m  $mapping_file -c  $category -o  ${outdir}LEfSe/${prefix}${category}_lefse.txt -u f;
base="${outdir}LEfSe/${prefix}${category}_lefse_LDA2"; format_input.py ${outdir}LEfSe/${prefix}${category}_lefse.txt ${base}.lefseinput.txt -c 2 -u 1 -o 1000000; run_lefse.py ${base}.lefseinput.txt ${base}.LDA.txt -l 2;
plot_res.py  --max_feature_len 200 --orientation h --format pdf --left_space 0.3 --dpi 300 ${base}.LDA.txt ${base}.pdf; plot_cladogram.py ${base}.LDA.txt --dpi 300 ${base}.cladogram.pdf --clade_sep 1.8 --format pdf --right_space_prop 0.45 --label_font_size 10;
base="${outdir}LEfSe/${prefix}${category}_lefse_LDA4"; format_input.py ${outdir}LEfSe/${prefix}${category}_lefse.txt ${base}.lefseinput.txt -c 2 -u 1 -o 1000000; run_lefse.py ${base}.lefseinput.txt ${base}.LDA.txt -l 4;
plot_res.py  --max_feature_len 200 --orientation h --format pdf --left_space 0.3 --dpi 300 ${base}.LDA.txt ${base}.pdf; plot_cladogram.py ${base}.LDA.txt --dpi 300 ${base}.cladogram.pdf --clade_sep 1.8 --format pdf --right_space_prop 0.45 --label_font_size 10;
source delefse''' % (self.path['bayegy_home'], self.abundance_table, self.mapping_file, g, self.out_dir, self.prefix), file=shscript)
                os.system("bash {}".format(self.running_bash))
