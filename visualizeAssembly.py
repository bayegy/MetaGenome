from visualizeFunction import VisualizeFunction
import pandas as pd
from mapInfo import MapInfo
from numpy import nan
# import pdb


class VisualizeAssembly(VisualizeFunction):
    """"""

    def __init__(self, *args, annotation_file, prefix, annotation_column=1, adjust_func=lambda x: x.split(',')[0], **kwargs):
        super(VisualizeAssembly, self).__init__(*args, **kwargs, prefix=prefix)
        m = MapInfo()
        m.load_map(annotation_file, adjust_func=adjust_func,
                   map_column=annotation_column)
        gene_abdc = pd.read_csv(self.abundance_table, sep='\t')
        gene_abdc.iloc[:, 0] = [m.map[n][0] if n in m.map.keys(
        ) else nan for n in gene_abdc.iloc[:, 0]]
        gene_abdc.iloc[:, 0] = [n or nan for n in gene_abdc.iloc[:, 0]]
        columns = list(gene_abdc.columns)
        columns[0] = "Genes"
        gene_abdc.columns = columns
        # pdb.set_trace()
        gene_abdc = gene_abdc.groupby("Genes").sum()
        self.set_attr(abundance_table=self.out_dir +
                      "All.{}.abundance_unstratified.tsv".format(prefix.strip('_')))
        gene_abdc.to_csv(self.abundance_table, sep='\t', index=True)
        self.mapping_df = gene_abdc
