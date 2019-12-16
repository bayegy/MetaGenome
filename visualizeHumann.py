from visualizeFunction import VisualizeFunction
import os
import pandas as pd
import re
# import pdb


class VisualizeHumann(VisualizeFunction):
    """"""

    def __init__(self, *args, regroup=False, prefix=False, id_with_name=False, custom_to=False, **kwargs):
        super(VisualizeHumann, self).__init__(*args, **kwargs, prefix=prefix)
        self.regroup = regroup
        self.save_colon = 'T' if id_with_name else 'F'

        if regroup or custom_to:
            ob = custom_to.upper() if custom_to else regroup.replace('uniref90_', '').upper()
            self.pre_abundance_table = "{}All.{}.abundance.tsv".format(self.out_dir, ob)
            os.system("humann2_regroup_table --input {} {} --output {}".format(self.abundance_table,
                                                                               "-c {}/map_{}_uniref90.txt.gz".format(self.path['humann2_utility_mapping'], custom_to) if custom_to else '--groups ' + regroup, self.pre_abundance_table))
            self.prefix = prefix or ob + '_'
        else:
            self.pre_abundance_table = self.abundance_table

        os.system("humann2_split_stratified_table -i {} -o {}".format(self.pre_abundance_table, self.out_dir))
        self.abundance_table = "{}{}_unstratified.tsv".format(
            self.out_dir, self.get_file_name(self.pre_abundance_table))
        df = pd.read_csv(self.abundance_table, sep='\t')
        if id_with_name:
            df['Description'] = df.iloc[:, 0]
        df = df.loc[(i not in ['UNMAPPED', 'UNINTEGRATED', 'UniRef90_unknown', 'UNGROUPED'] for i in df.iloc[:, 0]), :]
        if id_with_name:
            df.iloc[:, 0] = [re.search('^[^:]*', i).group() for i in df.iloc[:, 0]]
        # pdb.set_trace()
        df.to_csv(self.abundance_table, sep='\t', index=False)

    def get_file_name(self, path):
        return os.path.splitext(os.path.basename(path))[0]

    def __visualize_with_group__(self, *args, **kwargs):
        super(VisualizeHumann, self).__visualize_with_group__(*args, **kwargs)
        categories = [g.strip() for g in re.split(',', self.categories)]

        for g in categories:
            print(g)
            bar_out = "{}4-SignificanceAnalysis/LEfSe/SignificantFeatures/".format(self.out_dir)
            if not os.path.exists(bar_out):
                os.makedirs(bar_out)
            humann2_ft_lefse_lda = '{}4-SignificanceAnalysis/LEfSe/{}{}_lefse_LDA2.LDA.txt'.format(
                self.out_dir, self.prefix, g)
            df = pd.read_csv(humann2_ft_lefse_lda, sep='\t', header=None, index_col=0)[2]
            df = df[df.notna()]
            bar_table = "{}table_for_stratification_bar.txt".format(bar_out)
            os.system("Rscript {SCRIPTPATH}/write_data_for_lefse.R -i  {abundance_table} -m  {mapping_file} -c  {category} -o  {bar_table} -u f -j F -e {isenzyme} -n {save_colon}".format(
                SCRIPTPATH=self.path['bayegy_home'], abundance_table=self.pre_abundance_table, save_colon=self.save_colon, mapping_file=self.mapping_file, category=g, bar_table=bar_table, isenzyme=("T" if (self.regroup == 'uniref90_level4ec' or self.regroup == 'uniref50_level4ec') else "F")))
            features = [re.search('^[^:\|]*', f).group().strip()
                        for f in pd.read_csv(bar_table, sep='\t', index_col=0).index if not f.find('|') == -1]
            features = list(set(features))
            # pdb.set_trace()
            for f in df.index:
                print(f)
                # f = f.replace('_', '-')
                if f in features:
                    os.system("{b}humann2_barplot --input {i} --focal-feature {f} --focal-metadatum {g} --last-metadatum {g} --output {bar_out}/{p}{f}_stratification_bar.pdf".format(
                        i=bar_table, f=f, g=g, bar_out=bar_out, p=self.prefix, b=self._base_dir))
