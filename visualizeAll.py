import os
from visualizeSpecies import VisualizeSpecies
from visualizeFunction import VisualizeFunction


class VisualizeAll(VisualizeSpecies):
    """docstring for VisualizeAll"""

    def __init__(self, mapping_file, categories=False, out_dir=False, exclude_species="UNREALTAX"):
        super(VisualizeAll, self).__init__("multi_table", mapping_file,
                                           categories, out_dir=out_dir, exclude_species=exclude_species)
        self.out_dir = out_dir if out_dir else os.path.dirname(self.mapping_file) + '/'

    def visualize(self):
        for abundance_table in [self.out_dir + 'FMAP/' + f for f in ('All.Function.abundance.KeepID.KO.txt', 'All.Function.abundance.KeepID.Module.txt', 'All.Function.abundance.KeepID.Pathway.txt', 'All.Function.abundance.KeepID.Pathway.Level1.txt', 'All.Function.abundance.KeepID.Pathway.Level2.txt')] + [self.out_dir + 'AMR' + '/All.AMR.abundance.txt']:
            VisualizeFunction(abundance_table, self.mapping_file, self.categories).visualize()
        VisualizeSpecies(self.out_dir + 'Kraken2/All.Taxa.OTU.txt', self.mapping_file,
                         self.categories, exclude_species=self.exclude_species).visualize()
