import os
import re
import json
from abc import ABCMeta, abstractmethod
import pandas as pd
import shutil
import pdb
from pipconfig import settings


class Visualize(metaclass=ABCMeta):
    """docstring for Visualize"""

    def __init__(self, abundance_table, mapping_file=False, categories=False, prefix=False, out_dir=False):
        out_dir = out_dir if out_dir else os.path.dirname(abundance_table)
        self.out_dir = os.path.abspath(out_dir) + '/'
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        self._base_dir = os.path.dirname(__file__) + '/'
        self.path = settings.path
        self.abundance_table = os.path.abspath(abundance_table)
        self.mapping_file = os.path.abspath(mapping_file) if mapping_file else ''
        self.categories = categories if categories else ""
        self.TMP_DIR = None
        if prefix:
            self.prefix = prefix
        else:
            p_prefix = os.path.splitext(os.path.basename(self.abundance_table))[0]
            for r in ['pathabundance', 'abundance', 'KeepID', 'Humann2', "unstratified", 'All']:
                p_prefix = p_prefix.replace(r, '')
            p_prefix = re.sub('\.+', '.', p_prefix.strip('.|_'))
            self.prefix = p_prefix + '_'
        # self.running_bash = self.out_dir + 'visualize_function.sh'
        if mapping_file and categories and not abundance_table == 'multi_table':
            TMP_DIR = os.path.dirname(self.abundance_table) + '/TMP_DIR/'
            if not os.path.exists(TMP_DIR):
                os.makedirs(TMP_DIR)
            self.TMP_DIR = TMP_DIR
            abd_df = pd.read_csv(self.abundance_table, sep='\t')
            map_df = pd.read_csv(self.mapping_file, sep='\t')
            not_sample = [c not in map_df.iloc[:, 0].values for c in abd_df.columns]
            map_df = map_df.loc[map_df[self.categories].notna(), :]
            has_group = [c in map_df.iloc[:, 0].values for c in abd_df.columns]
            sel = [a or b for a, b in zip(not_sample, has_group)]
            self.mapping_df = map_df
            self.abundance_df = abd_df.loc[:, sel]
            is_sample = [c in self.mapping_df.iloc[:, 0].values for c in self.abundance_df.columns]
            # pdb.set_trace()
            self.abundance_df = self.abundance_df.loc[self.abundance_df.loc[:, is_sample].T.sum() > 0, :]
            self.abundance_table = TMP_DIR + os.path.basename(self.abundance_table)
            self.mapping_file = TMP_DIR + os.path.basename(self.mapping_file)
            self.mapping_df.to_csv(self.mapping_file, sep='\t', index=False)
            self.abundance_df.to_csv(self.abundance_table, sep='\t', index=False)

    @abstractmethod
    def __visualize_with_group__():
        pass

    @abstractmethod
    def __visualize_without_group__():
        pass

    def visualize(self, exclude='all'):
        print("Visualizing abundance table: {}".format(self.abundance_table))
        if self.categories and self.mapping_file:
            print("Visualize using group info...")
            self.__visualize_with_group__(exclude)
        else:
            print("No group info detected, visualize without group info")
            self.__visualize_without_group__()
        result_abundance = self.out_dir + os.path.basename(self.abundance_table)
        if not os.path.exists(result_abundance):
            self.abundance_df.to_csv(result_abundance, sep='\t', index=False)
        if self.TMP_DIR:
            shutil.rmtree(self.TMP_DIR)
