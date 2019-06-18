import os
import json
import re
from visualizeSpecies import VisualizeSpecies
from visualizeFunction import VisualizeFunction
from visualizeHumann import VisualizeHumann
from pyutils.read import read_to_html_table, get_kingdom_ratio, format_file
from colorMap import ColorMap
from mapInfo import MapInfo


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

    def map_ko_annotation(self):
        ko_file = self.out_dir + 'FMAP/All.KO.abundance_unstratified.tsv'
        out = os.path.dirname(ko_file)
        os.system('''
perl {SCRIPTPATH}/ConvergeKO2Module.pl {ko_file} > {out_dir}/All.KEGG.Module.txt
perl {SCRIPTPATH}/ConvergeKO2Pathway.pl {ko_file} > {out_dir}/All.KEGG.Pathway.txt
perl {SCRIPTPATH}/ConvergePathway2Level1.pl {out_dir}/All.KEGG.Pathway.txt > {out_dir}/All.KEGG.Pathway.Level1.txt
perl {SCRIPTPATH}/ConvergePathway2Level2.pl {out_dir}/All.KEGG.Pathway.txt > {out_dir}/All.KEGG.Pathway.Level2.txt
            '''.format(SCRIPTPATH=self.path['cii_home'], ko_file=ko_file, out_dir=out))

    def map_func_definition(self):
        FMAP_data = self.path['fmap_home'] + '/FMAP_data/'
        FMAP = self.out_dir + "FMAP/"
        mi = self.mi = MapInfo()
        mi.mapping(FMAP + 'All.KEGG.Pathway.txt', [FMAP_data + f for f in ['KEGG_Pathway2Level1.txt', 'KEGG_Pathway2Level2.txt', 'KEGG_pathway.txt']],
                   out_file=FMAP + 'All.KEGG.Pathway.mapping.txt', mapped_headers=["Level1", "Level2", "Level3"])
        mi.mapping(FMAP + 'All.KO.abundance_unstratified.tsv', [FMAP_data + 'KEGG_orthology.txt'],
                   out_file=FMAP + 'All.KO.mapping.txt', adjust_func=mi.ajust_ko_info, mapped_headers=['Gene_name\tEnzyme_number\tDefinition'])
        mi.mapping(FMAP + 'All.KO.mapping.txt', [FMAP_data + f for f in [
                   'KEGG_orthology2module.txt', 'KEGG_orthology2pathway.txt']], mapped_headers=["Module", "KEGG Pathway"], add=True)

        datas = [FMAP + f for f in ["All.KO.abundance_unstratified.tsv",
                                    "All.KEGG.Module.txt",
                                    "All.KEGG.Pathway.txt"]] + [self.out_dir + f for f in ["EC/All.LEVEL4EC.abundance_unstratified.tsv",
                                                                                           "EggNOG/All.EGGNOG.abundance_unstratified.tsv",
                                                                                           "GO/All.GO.abundance_unstratified.tsv"]]

        mapping_sources = [FMAP_data +
                           f for f in ["KEGG_orthology.txt",
                                       "KEGG_module.txt",
                                       "KEGG_pathway.txt"]] + [self.path['humann2_utility_mapping'] + '/' + f for f in ["map_level4ec_name.txt.gz",
                                                                                                                        "map_eggnog_name.txt.gz",
                                                                                                                        "map_go_name.txt.gz"]]

        for data, mapping_source in zip(datas, mapping_sources):
            mi.mapping(data, [mapping_source])
        mi.mapping(self.out_dir + 'AMR/All.AMR.abundance.txt', [FMAP_data + '/aro.csv'],
                   pattern="ARO[^\|]+", first_pattern="[^\|]+$", add_sid_to_info=True, map_column=-1)

    def visualize(self, exclude='none'):

        os.system("{}/piputils/write_colors_plan.py -i {} -c {} -p {}/piputils/group_color.list -o {}colors_plan.json".format(
            self.path['bayegy_home'], self.mapping_file, self.categories, self.path['bayegy_home'], self.out_dir))
        os.environ['COLORS_PLAN_PATH'] = self.out_dir + 'colors_plan.json'
        """
        os.system('''
humann2_renorm_table --input {0}/RPK.All.UniRef90.genefamilies.tsv --units cpm -o {0}/All.UniRef90.genefamilies.tsv -s n
humann2_renorm_table --input {0}/RPK.All.Metacyc.pathabundance.tsv --units cpm -o {0}/All.Metacyc.pathabundance.tsv -s n
'''.format(self.out_dir + 'Metagenome/Humann'))

        for group, out_dir in zip(['ko', 'level4ec', 'go', 'eggnog'], ['FMAP', 'EC', 'GO', 'EggNOG']):
            VisualizeHumann(self.out_dir + 'Metagenome/Humann/All.UniRef90.genefamilies.tsv', self.mapping_file,
                            self.categories, regroup='uniref90_' + group, out_dir=self.out_dir + out_dir).visualize(exclude)

        VisualizeHumann(self.out_dir + 'Metagenome/Humann/All.UniRef90.genefamilies.tsv', self.mapping_file,
                        self.categories, custom_to='cazy', out_dir=self.out_dir + 'CAZy').visualize(exclude)

        self.map_ko_annotation()
        self.map_func_definition()

        for abundance_table in [self.out_dir + 'FMAP/' + f for f in (
            # 'All.KO.abundance_unstratified.tsv',
            'All.KEGG.Module.txt',
            'All.KEGG.Pathway.txt',
            'All.KEGG.Pathway.Level1.txt',
            'All.KEGG.Pathway.Level2.txt'
        )]:
            VisualizeFunction(abundance_table, self.mapping_file, self.categories).visualize(exclude)

        VisualizeFunction(self.out_dir + 'AMR' + '/All.AMR.abundance.txt',
                          self.mapping_file, self.categories).visualize(exclude)

        VisualizeSpecies(self.out_dir + 'Kraken2/All.Taxa.OTU.txt', self.mapping_file,
                         self.categories, exclude_species=self.exclude_species).visualize(exclude)
        
        for abundance_table in [self.out_dir + 'Metagenome/Humann/' + f for f in [
            # 'All.UniRef90.genefamilies.tsv',
            'All.Metacyc.pathabundance.tsv'
        ]]:
            VisualizeHumann(abundance_table, self.mapping_file, self.categories, id_with_name=True).visualize(exclude)
        """
        categories_list = self.categories.split(',')
        """
        for g in categories_list:
            ColorMap(ko_lefse_lda=self.out_dir + 'FMAP/4-SignificanceAnalysis/LEfSe/KO_{}_lefse_LDA2.LDA.txt'.format(g),
                     ko_abundance_table=self.out_dir + 'FMAP/All.KO.abundance_unstratified.tsv', mapping_file=self.mapping_file, category=g, prefix=g + '_', out_dir=self.out_dir + 'FMAP/ColoredMaps/{}'.format(g)).plot_all()
        """
        os.system("cd {}&&bash {}orgnize_dir_structure.sh {} {}".format(self.out_dir,
                                                                        self._base_dir, self.mapping_file, categories_list[0]))
        page = self.out_dir + 'Result_Metagenomics/FiguresTablesForReport/src/pages/main_cleaned.html'
        format_file(page, page, table1=read_to_html_table(
            self.out_dir + 'Report/reads_summary.txt', table_class=['table', 'table-striped', 'table-sm'], thead_class=['thead-dark']),
            species_ratio=get_kingdom_ratio(self.out_dir + 'Kraken2/All.Taxa.OTU.txt'), report_category=categories_list[0])
