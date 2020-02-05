import os
import re
from visualizeSpecies import VisualizeSpecies
from visualizeFunction import VisualizeFunction
from visualizeHumann import VisualizeHumann
from visualizeAssembly import VisualizeAssembly
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
        root_dir = self.out_dir + 'Result/'
        self.set_path(force=True,
                      root_dir=root_dir,
                      qc_dir=root_dir + '00-QCStats',
                      )
        self.set_path(force=False,
                      FMAP_data=self.fmap_home + '/FMAP_data/',
                      )

        self.mi = MapInfo()

    def map_ko_annotation(self):
        ko_file = self.function_dir + '1-KEGG/All.KO.abundance_unstratified.tsv'
        # self.kegg_dir = os.path.dirname(ko_file) + '/'
        self.system('''
{perl_path} {base_dir}/ConvergeKO2Module.pl {ko_file} > {kegg_dir}/All.KEGG.Module.txt
{perl_path} {base_dir}/ConvergeKO2Pathway.pl {ko_file} > {kegg_dir}/All.KEGG.Pathway.txt
{perl_path} {base_dir}/ConvergePathway2Level1.pl {kegg_dir}/All.KEGG.Pathway.txt > {kegg_dir}/All.KEGG.Pathway.Level1.txt
{perl_path} {base_dir}/ConvergePathway2Level2.pl {kegg_dir}/All.KEGG.Pathway.txt > {kegg_dir}/All.KEGG.Pathway.Level2.txt
            ''', ko_file=ko_file)

    @staticmethod
    def extract_empper_cog(annotation):
        return ['ENOG41' + c if not c.startswith('COG') else c for c in re.findall('([^,]+)@NOG', annotation)]

    def map_func_definition(self):
        FMAP_data = self.FMAP_data
        FMAP = self.kegg_dir
        mi = self.mi
        asem = self.base_on_assembly
        mi.mapping(FMAP + 'All.KEGG.Pathway.txt', [FMAP_data + f for f in ['KEGG_Pathway2Level1.txt', 'KEGG_Pathway2Level2.txt', 'KEGG_pathway.txt']],
                   out_file=FMAP + 'All.KEGG.Pathway.mapping.txt', mapped_headers=["Level1", "Level2", "Level3"])
        mi.mapping(FMAP + 'All.KO.abundance_unstratified.tsv', [FMAP_data + 'KEGG_orthology.txt'],
                   out_file=FMAP + 'All.KO.mapping.txt', adjust_func=mi.ajust_ko_info, mapped_headers=['Gene_name\tEnzyme_number\tDefinition'])
        mi.mapping(FMAP + 'All.KO.mapping.txt', [FMAP_data + f for f in [
                   'KEGG_orthology2module.txt', 'KEGG_orthology2pathway.txt']], mapped_headers=["Module", "KEGG Pathway"], add=True)

        datas = [FMAP + f for f in ["All.KO.abundance_unstratified.tsv",
                                    "All.KEGG.Module.txt",
                                    "All.KEGG.Pathway.txt"]] + [self.function_dir + f for f in ["5-EC/All.LEVEL4EC.abundance_unstratified.tsv",
                                                                                                "{}/All.EGGNOG.abundance_unstratified.tsv".format(
                                                                                                    "2-EggNOG" if asem else"3-EggNOG"),
                                                                                                "{}/All.GO.abundance_unstratified.tsv".format("3-GO" if asem else "4-GO")]]

        mapping_sources = [FMAP_data +
                           f for f in ["KEGG_orthology.txt",
                                       "KEGG_module.txt",
                                       "KEGG_pathway.txt"]] + [self.path['humann2_utility_mapping'] + '/' + f for f in ["map_level4ec_name.txt.gz",
                                                                                                                        "map_eggnog_name.txt.gz",
                                                                                                                        "map_go_name.txt.gz"]]

        for data, mapping_source in zip(datas, mapping_sources):
            if os.path.exists(data):
                mi.mapping(data, [mapping_source])

        mi.mapping(self.cazy_dir + '/All.CAZY.abundance_unstratified.tsv',
                   [self.FMAP_data + '/fam_description.txt'])

    def init_out_dir(self, category):
        asem = self.base_on_assembly
        self.set_attr(
            category=category
        )
        categroy_dir = "{}/{}/".format(self.root_dir, self.category)
        self.set_path(force=True,
                      categroy_dir=categroy_dir,
                      taxa_dir=categroy_dir + '1-TaxaAundanceAnalysis/',
                      function_dir=categroy_dir + '2-FuctionAnalysis/',
                      amr_dir=categroy_dir + '3-AMRAnalysis/',
                      report_dir=categroy_dir + 'FiguresTablesForReport/',
                      )
        self.set_path(force=True,
                      kegg_dir=self.function_dir + '1-KEGG/',
                      go_dir=self.function_dir +
                      ('3-GO/' if asem else '4-GO/'),
                      eggnog_dir=self.function_dir +
                      ('2-EggNOG/' if asem else '3-EggNOG/'),
                      cazy_dir=self.function_dir +
                      ('4-CAZy/' if asem else '6-CAZy/'),
                      )

    def visualize(self, exclude='none', base_on_assembly=False):
        self.base_on_assembly = base_on_assembly
        self.system("{bayegy_home}/piputils/write_colors_plan.py -i {mapping_file} -c {categories} -p {bayegy_home}/piputils/group_color.list -o {out_dir}colors_plan.json")
        os.environ['COLORS_PLAN_PATH'] = self.out_dir + 'colors_plan.json'

        if base_on_assembly:
            self.set_path(force=True, asem_dir=self.root_dir + '01-Assembly')
            spec = """
mkdir -p {asem_dir}/2-ORFPrediction/ {asem_dir}/1-Quast/
cp {out_dir}Assembly_out/ORF* {asem_dir}/2-ORFPrediction/
cp -r {out_dir}Assembly_out/quast_results/* {asem_dir}/1-Quast/
cp {out_dir}/salmon_out/All.genes.abundance.txt {out_dir}Assembly_out/NR.nucleotide.fa {asem_dir}/
            """.format(**self.context)
        else:
            self.set_path(force=True, table_dir=self.root_dir + '01-BaseTables')
            self.system('''
{humann2_home}/humann2_renorm_table --input {b}/RPK.All.UniRef90.genefamilies.tsv --units cpm -o {b}/All.UniRef90.genefamilies.tsv -s n
{humann2_home}/humann2_renorm_table --input {b}/RPK.All.Metacyc.pathabundance.tsv --units cpm -o {b}/All.Metacyc.pathabundance.tsv -s n
            ''', b=self.out_dir + 'Metagenome/Humann')
            spec = """
cp {out_dir}Metagenome/Humann/All.UniRef90.genefamilies.tsv {table_dir}
cp {out_dir}Metagenome/Metaphlan/All.Metaphlan2.profile.txt {table_dir}/All.Metaphlan2.taxa.txt
cp {out_dir}Kraken2/All.Taxa.OTU.txt  {table_dir}/All.Kraken2.taxa.txt
{humann2_home}/humann2_split_stratified_table -i {table_dir}/All.UniRef90.genefamilies.tsv -o {table_dir}
            """.format(**self.context)

        self.system("""
mkdir -p {qc_dir}/1-QC_report_Rawfastq/ {qc_dir}/2-QC_report_Filtered/
cp {mapping_file} {root_dir}
cp {out_dir}/kneaddata_out/fastqc/*html {qc_dir}/1-QC_report_Rawfastq/
rm {qc_dir}/1-QC_report_Rawfastq/*unmatched*
mv {qc_dir}/1-QC_report_Rawfastq/*kneaddata* {qc_dir}/2-QC_report_Filtered/
cp {out_dir}Report/reads_summary.txt {qc_dir}/
{spec}
            """, spec=spec)
        for g in self.categories.split(','):
            self.init_out_dir(g)

            if base_on_assembly:
                VisualizeAssembly(self.out_dir + 'salmon_out/All.genes.abundance.txt', self.mapping_file, self.category, annotation_file=self.out_dir +
                                  'salmon_out/genes.emapper.annotations', prefix='KO_', annotation_column=6, out_dir=self.function_dir + '1-KEGG').visualize(exclude)
                VisualizeAssembly(self.out_dir + 'salmon_out/All.genes.abundance.txt', self.mapping_file, self.category, annotation_file=self.out_dir +
                                  'salmon_out/genes.emapper.annotations', prefix='GO_', annotation_column=5, out_dir=self.function_dir + '3-GO').visualize(exclude)

                VisualizeAssembly(self.out_dir + 'salmon_out/All.genes.abundance.txt', self.mapping_file, self.category, annotation_file=self.out_dir +
                                  'salmon_out/genes.emapper.annotations', prefix='EGGNOG_', annotation_column=9, out_dir=self.function_dir + '2-EggNOG', adjust_func=self.extract_empper_cog).visualize(exclude)
                VisualizeAssembly(self.out_dir + 'salmon_out/All.genes.abundance.txt', self.mapping_file, self.category, annotation_file=self.out_dir +
                                  'salmon_out/genes_cazy.f6', prefix='CAZY_', out_dir=self.function_dir + '4-CAZy', adjust_func=lambda x: re.search('\|([^\|_]+)', x).group(1)).visualize(exclude)
                vcard = VisualizeAssembly(self.out_dir + 'salmon_out/All.genes.abundance.txt', self.mapping_file, self.category, annotation_file=self.out_dir +
                                          'salmon_out/genes_card.f6', prefix='AMR_', out_dir=self.amr_dir, adjust_func=False)
            else:
                for group, out_dir in zip(['ko', 'level4ec', 'go', 'eggnog', 'cazy'], ['1-KEGG', '5-EC', '4-GO', '3-EggNOG', '6-CAZy']):
                    VisualizeHumann(self.out_dir + 'Metagenome/Humann/All.UniRef90.genefamilies.tsv', self.mapping_file,
                                    self.category, custom_to=group, out_dir=self.function_dir + out_dir).visualize(exclude)
                VisualizeHumann(self.out_dir + 'Metagenome/Humann/All.Metacyc.pathabundance.tsv',
                                self.mapping_file, self.category, id_with_name=True, out_dir=self.function_dir + '2-Metacyc').visualize(exclude)

            self.map_ko_annotation()

            if base_on_assembly:
                self.mi.mapping(self.amr_dir + 'All.AMR.abundance_unstratified.tsv', [self.FMAP_data + '/aro.csv'],
                                pattern="ARO[^\|]+",
                                # first_pattern="[^\|]+$",
                                first_pattern="ARO[^\|]+",
                                # add_sid_to_info=True,
                                add_sid_to_info=False,
                                map_column=-1)
                vcard.visualize(exclude)
            else:
                self.mi.mapping(self.out_dir + 'FMAP/All.AMR.abundance_unstratified.tsv', [self.FMAP_data + '/aro.csv'],
                                pattern="ARO[^\|]+",
                                # first_pattern="[^\|]+$",
                                first_pattern="ARO[^\|]+",
                                # add_sid_to_info=True,
                                add_sid_to_info=False,
                                map_column=-1)
                VisualizeFunction(self.out_dir + 'FMAP/All.AMR.abundance_unstratified.tsv',
                                  self.mapping_file, self.category, out_dir=self.amr_dir).visualize(exclude)

            for abundance_table in [self.kegg_dir + f for f in (
                # 'All.KO.abundance_unstratified.tsv',
                'All.KEGG.Module.txt',
                'All.KEGG.Pathway.txt',
                'All.KEGG.Pathway.Level1.txt',
                'All.KEGG.Pathway.Level2.txt'
            )]:
                VisualizeFunction(abundance_table, self.mapping_file, self.category,
                                  out_dir=self.function_dir + '1-KEGG').visualize(exclude)

            VisualizeSpecies(self.out_dir + 'Kraken2/All.Taxa.OTU.txt', self.mapping_file,
                             self.category, exclude_species=self.exclude_species, out_dir=self.taxa_dir, tmp_dir=self.out_dir + '/TMP_DIR').visualize(exclude)

            self.map_func_definition()

            ColorMap(ko_lefse_lda=self.function_dir + '1-KEGG/4-SignificanceAnalysis/LEfSe/KO_{}_lefse_LDA2.LDA.txt'.format(g),
                     ko_abundance_table=self.function_dir + '1-KEGG/All.KO.abundance_unstratified.tsv', mapping_file=self.mapping_file, category=g, out_dir=self.function_dir + '1-KEGG/5-ColoredMaps/').plot_all(map_level=self.function_dir + '1-KEGG/All.KEGG.Pathway.mapping.txt')

            reads_spec = """
cp -rp {base_dir}/Report/src {report_dir}
cp {base_dir}/Report/结题报告.html {categroy_dir}
cp {function_dir}/2-Metacyc/1-Barplots/Metacyc_{category}_barplot.pdf  Figure5-5.pdf
cpfirst "{function_dir}/5-EC/4-SignificanceAnalysis/LEfSe/SignificantFeatures/.pdf"   Figure5-8.pdf
            """.format(**self.context)

            asem_spec = """
cp -rp {base_dir}/Report_assembly/src {report_dir}
cp {base_dir}/Report_assembly/结题报告.html {categroy_dir}
            """.format(**self.context)

            self.system("""
cd {report_dir}
cp {taxa_dir}/1-AbundanceSummary/Classified_stat_relative.png Figure4-1.png
cp {taxa_dir}/1-AbundanceSummary/2-Barplots/Taxa-bar-plots-top20/Phylum_{category}_barplot.pdf Figure4-2.pdf
cp {taxa_dir}/1-AbundanceSummary/3-Heatmaps/Heatmap_top20_clustered/Phylum/{category}_heatmap.pdf  Figure4-3.pdf
cp {taxa_dir}/2-AbundanceComparison/LEfSe/Genus/{category}_Genus_lefse_LDA2.pdf  Figure4-4.pdf
cp {taxa_dir}/2-AbundanceComparison/VennAndFlower/{category}_Venn_plot.png Figure4-5.png
cp {function_dir}/1-KEGG/1-Barplots/KEGG.Pathway.Level1_{category}_barplot.pdf Figure5-1.pdf
cp {function_dir}/1-KEGG/4-SignificanceAnalysis/LEfSe/KEGG.Pathway_{category}_lefse_LDA2.pdf   Figure5-2.pdf
cpfirst "{function_dir}/1-KEGG/4-SignificanceAnalysis/LEfSe/SignificantFeatures/.pdf"   Figure5-4.pdf
cp {eggnog_dir}/2-Heatmaps/EGGNOG_{category}_clustered_heatmap.pdf   Figure5-6.pdf
cp {go_dir}/3-Circos/GO_{category}_circos.png   Figure5-7.png
cp {cazy_dir}/4-SignificanceAnalysis/LEfSe/CAZY_{category}_lefse_LDA2.pdf   Figure5-9.pdf
cp {amr_dir}/1-Barplots/AMR_{category}_barplot.pdf Figure6-1.pdf
cp {amr_dir}/2-Heatmaps/AMR_{category}_nocluster_heatmap.pdf Figure6-2.pdf
cp {amr_dir}/3-Circos/AMR_{category}_circos.png Figure6-3.png
cp {taxa_dir}/4-CorrelationAnalysis/RDA/Genus/{category}_RDA_features_location_plot.pdf Figure7-1.pdf
cp {amr_dir}/5-CorrelationAnalysis/CorrelationHeatmap/AMR_Correlation_heatmap.pdf Figure7-2.pdf
mv {kegg_dir}/5-CorrelationAnalysis {kegg_dir}6-CorrelationAnalysis
{spec}
if [ -f Figure4-2.pdf ];then echo "Converting pdf to png"; for pdfs in *.pdf; do echo $pdfs; base=$(basename $pdfs .pdf); convert  -density 300 -quality 80 $pdfs ${{base}}.png; rm $pdfs;done;fi;
                """, spec=asem_spec if base_on_assembly else reads_spec)

            page = self.report_dir + \
                'src/pages/main_cleaned.html' if base_on_assembly else self.categroy_dir + "结题报告.html"

            format_file(page, page, table1=read_to_html_table(
                self.out_dir + 'Report/reads_summary.txt', table_class=['table', 'table-striped', 'table-sm'], thead_class=['thead-dark']),
                species_ratio=get_kingdom_ratio(self.out_dir + 'Kraken2/All.Taxa.OTU.txt'), report_category=self.category)

        self.system(
            "{base_dir}/change_suffix.py {root_dir} -o txt -s 'All.Taxa.OTU.taxa-bar-plots,bray_curtis_emperor' ")
        self.system(
            "{base_dir}/change_suffix.py {root_dir} -o tsv -s 'All.Taxa.OTU.taxa-bar-plots,bray_curtis_emperor' ")
