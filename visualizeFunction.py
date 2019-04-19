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

    def __visualize_with_group__(self):
        categories = [g.strip() for g in re.split(',', self.categories)]
        for g in categories:
            # pdb.set_trace()
            Circos(self.abundance_table, mapping_file=self.mapping_file, category=g, by_group_mean=False,
                   prefix=self.prefix + g + '_', out_dir=self.out_dir + 'Circos').visualize()
            Circos(self.abundance_table, mapping_file=self.mapping_file, category=g, by_group_mean=True,
                   prefix=self.prefix + g + '_groupMean_', out_dir=self.out_dir + 'Circos').visualize()
            os.system('''
SCRIPTPATH=%s
abundance_table=%s
mapping_file=%s
category=%s
outdir=%s
prefix=%s
Rscript ${SCRIPTPATH}/abundance_barplot.R -n 20 -m $mapping_file -c $category -i $abundance_table -o ${outdir}Barplots/ -p ${prefix}${category}_ -b F;
Rscript ${SCRIPTPATH}/abundance_barplot.R -n 20 -m $mapping_file -c $category -i $abundance_table -o ${outdir}Barplots/ -p ${prefix}${category}_groupMean_ -b T;
Rscript ${SCRIPTPATH}/abundance_heatmap.R  -m $mapping_file -c $category -n 20 -i $abundance_table -o ${outdir}Heatmaps -p ${prefix}${category}_nocluster_ -l F -t F;
Rscript ${SCRIPTPATH}/abundance_heatmap.R -m $mapping_file -c $category  -n 20 -i $abundance_table -o ${outdir}Heatmaps -p ${prefix}${category}_clustered_ -l F -t F -u T;
Rscript ${SCRIPTPATH}/abundance_heatmap.R  -m $mapping_file -c $category -n 20 -i $abundance_table -o ${outdir}Heatmaps -p ${prefix}${category}_groupMean_ -l F -t T -b T;
Rscript ${SCRIPTPATH}/write_data_for_lefse.R -i  $abundance_table -m  $mapping_file -c  $category -o  ${outdir}LEfSe/${prefix}${category}_lefse.txt -u f;
source lefse;
base="${outdir}LEfSe/${prefix}${category}_lefse_LDA2"; format_input.py ${outdir}LEfSe/${prefix}${category}_lefse.txt ${base}.lefseinput.txt -c 2 -u 1 -o 1000000; run_lefse.py ${base}.lefseinput.txt ${base}.LDA.txt -l 2;
plot_res.py  --max_feature_len 200 --orientation h --format pdf --left_space 0.3 --dpi 300 ${base}.LDA.txt ${base}.pdf;
#plot_cladogram.py ${base}.LDA.txt --dpi 300 ${base}.cladogram.pdf --clade_sep 1.8 --format pdf --right_space_prop 0.45 --label_font_size 10;
base="${outdir}LEfSe/${prefix}${category}_lefse_LDA4"; format_input.py ${outdir}LEfSe/${prefix}${category}_lefse.txt ${base}.lefseinput.txt -c 2 -u 1 -o 1000000; run_lefse.py ${base}.lefseinput.txt ${base}.LDA.txt -l 4;
plot_res.py  --max_feature_len 200 --orientation h --format pdf --left_space 0.3 --dpi 300 ${base}.LDA.txt ${base}.pdf;
#plot_cladogram.py ${base}.LDA.txt --dpi 300 ${base}.cladogram.pdf --clade_sep 1.8 --format pdf --right_space_prop 0.45 --label_font_size 10;
source delefse''' % (self.path['bayegy_home'], self.abundance_table, self.mapping_file, g, self.out_dir, self.prefix))
            # os.system("bash {}".format(self.running_bash))

    def __visualize_without_group__(self):
        Circos(self.abundance_table, prefix=self.prefix, out_dir=self.out_dir + 'Circos').visualize()
        os.system('''
SCRIPTPATH=%s
abundance_table=%s
outdir=%s
prefix=%s
Rscript ${SCRIPTPATH}/abundance_barplot.R -n 20 -i $abundance_table -o ${outdir}Barplots/ -p ${prefix} -b F;
Rscript ${SCRIPTPATH}/abundance_heatmap.R -n 20 -i $abundance_table -o ${outdir}Heatmaps  -p ${prefix} -l F -t F -u T;
''' % (self.path['bayegy_home'], self.abundance_table, self.out_dir, self.prefix))
