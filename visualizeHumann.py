from visualizeFunction import VisualizeFunction
import os
import pandas as pd
import re
from osEnv import OSEnv
# import pdb


class VisualizeHumann(VisualizeFunction):
    """"""

    def __init__(self, *args, regroup=False, prefix=False, id_with_name=False, custom_to=False, **kwargs):
        super(VisualizeHumann, self).__init__(*args, **kwargs, prefix=prefix)
        self.regroup = regroup
        self.set_attr(save_colon='T' if id_with_name else 'F')

        if regroup or custom_to:
            ob = custom_to.upper() if custom_to else regroup.replace('uniref90_', '').upper()
            self.set_attr(prefix=prefix or ob + '_',
                          pre_abundance_table="{}All.{}.abundance.tsv".format(
                              self.out_dir, ob)
                          )
            self.system("{humann2_home}/humann2_regroup_table --input {abundance_table} {par} --output {pre_abundance_table}",
                        par="-c {}/map_{}_uniref90.txt.gz".format(self.humann2_utility_mapping, custom_to) if custom_to else '--groups ' + regroup)

        else:
            self.set_attr(pre_abundance_table=self.abundance_table)

        self.system(
            "{humann2_home}/humann2_split_stratified_table -i {pre_abundance_table} -o {out_dir}")
        self.set_path(force=False, abundance_table="{}{}_unstratified.tsv".format(
            self.out_dir, self.get_file_name(self.pre_abundance_table)))
        df = pd.read_csv(self.abundance_table, sep='\t')
        if id_with_name:
            df['Description'] = df.iloc[:, 0]
        df = df.loc[(i not in ['UNMAPPED', 'UNINTEGRATED',
                               'UniRef90_unknown', 'UNGROUPED'] for i in df.iloc[:, 0]), :]
        if id_with_name:
            df.iloc[:, 0] = [re.search('^[^:]*', i).group()
                             for i in df.iloc[:, 0]]
        # pdb.set_trace()
        df.to_csv(self.abundance_table, sep='\t', index=False)

    def get_file_name(self, path):
        return os.path.splitext(os.path.basename(path))[0]

    def __visualize_with_group__(self, *args, **kwargs):
        super(VisualizeHumann, self).__visualize_with_group__(*args, **kwargs)
        categories = [g.strip() for g in re.split(',', self.categories)]

        for g in categories:
            print(g)
            bar_out = "{}4-SignificanceAnalysis/LEfSe/SignificantFeatures/".format(
                self.out_dir)
            self.set_path(force=True, bar_out=bar_out)

            humann2_ft_lefse_lda = '{}4-SignificanceAnalysis/LEfSe/{}{}_lefse_LDA2.LDA.txt'.format(
                self.out_dir, self.prefix, g)
            self.set_path(
                force=False, humann2_ft_lefse_lda=humann2_ft_lefse_lda)
            df = pd.read_csv(humann2_ft_lefse_lda, sep='\t',
                             header=None, index_col=0)[2]
            df = df[df.notna()]
            bar_table = "{}table_for_stratification_bar.txt".format(bar_out)
            self.set_attr(bar_table=bar_table)
            self.system("{R_path} {bayegy_home}/write_data_for_lefse.R -i  {pre_abundance_table} -m  {mapping_file} -c  {category} -o  {bar_table} -u f -j F -e {isenzyme} -n {save_colon}",
                        category=g, isenzyme=("T" if (self.regroup == 'uniref90_level4ec' or self.regroup == 'uniref50_level4ec') else "F"))
            features = [re.search('^[^:\|]*', f).group().strip()
                        for f in pd.read_csv(self.bar_table, sep='\t', index_col=0).index if not f.find('|') == -1]
            features = list(set(features))

            with OSEnv(path=self.lefse_py_home, pythonpath=self.lefse_pylib_home):
                # pdb.set_trace()
                for f in df.index:
                    print(f)
                    # f = f.replace('_', '-')
                    if f in features:
                        self.system("{base_dir}/humann2_barplot --input {bar_table} --focal-feature {f} --focal-metadatum {g} --last-metadatum {g} --output {bar_out}/{prefix}{f}_stratification_bar.pdf",
                                    f=f, g=g)
