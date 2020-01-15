from visualizeFunction import VisualizeFunction
import pandas as pd
from numpy import nan
# import pdb


class VisualizeAssembly(VisualizeFunction):
    """"""

    def __init__(self, *args, annotation_file, prefix="Gene", annotation_column=1, adjust_func=lambda x: x.split(','), **kwargs):
        super(VisualizeAssembly, self).__init__(*args, **kwargs, prefix=prefix)
        self.set_attr(tmp_map=self.out_dir + "/map_{}_to_genes.tsv".format(prefix.strip("_")))
        ft_map = {}
        with open(annotation_file) as f:
            for line in f:
                li = line.strip().split('\t')
                annotation = li[annotation_column]
                if annotation:
                    if adjust_func:
                        annotation = adjust_func(annotation)
                    if not isinstance(annotation, list):
                        annotation = [annotation]
                    for anno in annotation:
                        ft_map.setdefault(anno, []).append(li[0])
        with open(self.tmp_map, 'w') as fout:
            for k, v in ft_map.items():
                fout.write("{}\t{}\n".format(k, '\t'.join(v)))

        pre_abundance_table = self.abundance_table
        self.set_attr(abundance_table=self.out_dir +
                      "All.{}.abundance_unstratified.tsv".format(prefix.strip('_')))
        self.system("{humann2_home}/humann2_regroup_table --input {pre_abundance_table} --custom {tmp_map} --ungrouped N --output {abundance_table}",
                    pre_abundance_table=pre_abundance_table)
