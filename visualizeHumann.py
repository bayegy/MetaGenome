from visualizeFunction import VisualizeFunction
import os
import pandas as pd
import re
import pdb


class VisualizeHumann(VisualizeFunction):
    """"""

    def __init__(self, *args, **kwargs):
        super(VisualizeHumann, self).__init__(*args, **kwargs)
        self.pre_abundance_table = self.abundance_table
        os.system("humann2_split_stratified_table -i {} -o {}".format(self.pre_abundance_table, self.out_dir))
        self.abundance_table = "{}{}_unstratified.tsv".format(
            self.out_dir, os.path.splitext(os.path.basename(self.pre_abundance_table))[0])
        df = pd.read_csv(self.abundance_table, sep='\t')
        df['Description'] = df.iloc[:, 0]
        df = df.loc[(i not in ['UNMAPPED', 'UNINTEGRATED', 'UniRef90_unknown'] for i in df.iloc[:, 0]), :]
        df.iloc[:, 0] = [re.search('^[^:]*', i).group() for i in df.iloc[:, 0]]
        # pdb.set_trace()
        df.to_csv(self.abundance_table, sep='\t', index=False)

    def __visualize_with_group__(self, *args, **kwargs):
        super(VisualizeHumann, self).__visualize_with_group__(*args, **kwargs)
        categories = [g.strip() for g in re.split(',', self.categories)]

        for g in categories:
            print(g)
            bar_out = "{}LEfSe/SignificantFeatures/{}/".format(self.out_dir, g)
            if not os.path.exists(bar_out):
                os.makedirs(bar_out)
            humann2_ft_lefse_lda = '{}LEfSe/{}{}_lefse_LDA2.LDA.txt'.format(self.out_dir, self.prefix, g)
            df = pd.read_csv(humann2_ft_lefse_lda, sep='\t', header=None, index_col=0)[2]
            df = df[df.notna()]
            bar_table = "{}table_for_stratification_bar.txt".format(bar_out)
            os.system("Rscript {SCRIPTPATH}/write_data_for_lefse.R -i  {abundance_table} -m  {mapping_file} -c  {category} -o  {bar_table} -u f -j F".format(
                SCRIPTPATH=self.path['bayegy_home'], abundance_table=self.pre_abundance_table, mapping_file=self.mapping_file, category=g, bar_table=bar_table))
            features = [re.search('^[^:]*', f).group()
                        for f in pd.read_csv(bar_table, sep='\t', index_col=0).index if not f.find('|') == -1]
            features = list(set(features))

            for f in df.index:
                print(f)
                # f = f.replace('_', '-')
                if f in features:
                    os.system("humann2_barplot --input {i} --focal-feature {f} --focal-metadatum {g} --last-metadatum {g} --output {bar_out}/{p}{f}_stratification_bar.pdf".format(
                        i=bar_table, f=f, g=g, bar_out=bar_out, p=self.prefix))
