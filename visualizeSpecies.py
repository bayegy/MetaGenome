import os
from visualize import Visualize


class VisualizeSpecies(Visualize):
    """docstring for VisualizeSpecies"""

    def __init__(self, abundance_table, mapping_file, categories=False, prefix=False, out_dir=False, exclude_species="UNREALTAX", tmp_dir='./'):
        super(VisualizeSpecies, self).__init__(abundance_table, mapping_file, categories, prefix, out_dir)
        self.exclude_species = exclude_species
        # self.tmp_dir = tmp_dir

    def __visualize_with_group__(self, exclude='none'):
        os.system("bash {}visualize_otu_table_with_group.sh {} {} {} {} {} {} {} '{}'".format(
            self._base_dir, self.abundance_table, self.mapping_file, self.categories, self.prefix, self.exclude_species, self.path[
                'bayegy_home'], self.TMP_DIR, exclude
        ))
        os.system("""
mkdir -p {out_dir}/1-AbundanceSummary/1-RelativeAbundance \
    {out_dir}/1-AbundanceSummary/2-Barplots \
    {out_dir}/1-AbundanceSummary/4-Abundance_Metaphlan2 \
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
            """.format(tmp_dir=self.TMP_DIR, out_dir=self.out_dir))

    def __visualize_without_group__(self):
        os.system("bash {}visualize_otu_table_without_group.sh {} {} none {} {} {} {}".format(
            self._base_dir, self.abundance_table, self.mapping_file, self.prefix, self.exclude_species, self.path[
                'bayegy_home'], self.tmp_dir
        ))
