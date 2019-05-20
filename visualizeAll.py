import os
import json
import re
from visualizeSpecies import VisualizeSpecies
from visualizeFunction import VisualizeFunction
from pyutils.read import read_to_html_table, get_kingdom_ratio, format_file
from colorMap import ColorMap


class VisualizeAll(VisualizeSpecies):
    """
    Sample usage:
        from visualizeAll import VisualizeAll
        VisualizeAll("/home/cheng/Projects/rll_testdir/mapping_file.txt","Group1").visualize()
    """

    def __init__(self, mapping_file, categories=False, out_dir=False, exclude_species="Environmentalsamples"):
        super(VisualizeAll, self).__init__("multi_table", mapping_file,
                                           categories, out_dir=out_dir, exclude_species=exclude_species)
        self.out_dir = out_dir if out_dir else os.path.dirname(self.mapping_file) + '/'
        self._base_dir = os.path.dirname(__file__) + '/'

    def visualize(self, exclude='none'):

        os.system("{}/piputils/write_colors_plan.py -i {} -c {} -p {}/piputils/group_color.list -o {}colors_plan.json".format(
            self.path['bayegy_home'], self.mapping_file, self.categories, self.path['bayegy_home'], self.out_dir))
        os.environ['COLORS_PLAN_PATH'] = self.out_dir + 'colors_plan.json'
        """
        for abundance_table in [self.out_dir + 'FMAP/' + f for f in ('All.Function.abundance.KeepID.KO.txt',
                                                                     'All.Function.abundance.KeepID.Module.txt',
                                                                     'All.Function.abundance.KeepID.Pathway.txt',
                                                                     'All.Function.abundance.KeepID.Pathway.Level1.txt',
                                                                     'All.Function.abundance.KeepID.Pathway.Level2.txt')]:
            VisualizeFunction(abundance_table, self.mapping_file, self.categories).visualize(exclude)

        VisualizeFunction(self.out_dir + 'AMR' + '/All.AMR.abundance.txt',
                          self.mapping_file, self.categories).visualize(exclude)

        VisualizeSpecies(self.out_dir + 'Kraken2/All.Taxa.OTU.txt', self.mapping_file,
                         self.categories, exclude_species=self.exclude_species).visualize(exclude)
        """
        for g in self.categories.split(','):
            ColorMap(ko_lefse_lda=self.out_dir + 'FMAP/LEfSe/KO_{}_lefse_LDA2.LDA.txt'.format(g),
                     ko_abundance_table=self.out_dir + 'FMAP/All.Function.abundance.KeepID.KO.txt', mapping_file=self.mapping_file, category=g, prefix=g + '_', out_dir=self.out_dir + 'FMAP/ColoredMaps/{}'.format(g)).plot_all()

        os.system("cd {}&&bash {}orgnize_dir_structure.sh {} {}".format(self.out_dir,
                                                                        self._base_dir, self.mapping_file, re.sub(',.*$', '', self.categories)))
        format_file(self.out_dir + 'Result_Metagenomics/FiguresTablesForReport/src/pages/main_cleaned.html', table1=read_to_html_table(
            self.out_dir + 'Report/reads_summary.txt'), species_ratio=get_kingdom_ratio(self.out_dir + 'Kraken2/All.Taxa.OTU.txt'))
