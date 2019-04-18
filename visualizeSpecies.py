import os
from visualize import Visualize


class VisualizeSpecies(Visualize):
    """docstring for VisualizeSpecies"""

    def __init__(self, abundance_table, mapping_file, categories=False, prefix=False, out_dir=False, exclude_species="UNREALTAX"):
        super(VisualizeSpecies, self).__init__(abundance_table, mapping_file, categories, prefix, out_dir)
        self.exclude_species = exclude_species

    def __visualize_with_group__(self):
        os.system("bash {}visualize_otu_table_with_group.sh {} {} {} {} {} {} {}".format(
            self._base_dir, self.abundance_table, self.mapping_file, self.categories, self.prefix, self.exclude_species, self.path[
                'bayegy_home'], self.out_dir
        ))

    def __visualize_without_group__(self):
        os.system("bash {}visualize_otu_table_without_group.sh {} {} none {} {} {} {}".format(
            self._base_dir, self.abundance_table, self.mapping_file, self.prefix, self.exclude_species, self.path[
                'bayegy_home', self.out_dir]
        ))
