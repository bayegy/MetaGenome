import os
from visualize import Visualize


class VisualizeSpecies(Visualize):
    """docstring for VisualizeSpecies"""

    def visualize_with_group(self):
        os.system("{}visualize_otu_table_with_group.sh {} {} {} {} {}".format(
            self._base_dir, self.path['bayegy_home'], self.abundance_table, self.mapping_file, self.categories, self.prefix
        ))

    def visualize_without_group():
        pass
