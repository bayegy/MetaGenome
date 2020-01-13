import os
import re
import json
from abc import ABCMeta, abstractmethod
import pandas as pd
import shutil
import pdb
from pipconfig import settings
from systemMixin import SystemMixin


class Visualize(SystemMixin, metaclass=ABCMeta):
    """docstring for Visualize"""

    def __init__(self, abundance_table, mapping_file=False, categories=False, prefix=False, out_dir=False):
        out_dir = out_dir if out_dir else os.path.dirname(
            mapping_file or abundance_table)
        self.set_path(force=True,
                      out_dir=out_dir,
                      tmp_dir=out_dir + '/TMP_DIR/',
                      )

        self.path = settings.path

        if prefix:
            p_prefix = prefix
        else:
            p_prefix = os.path.splitext(os.path.basename(abundance_table))[0]
            for r in ['pathabundance', 'abundance', 'KeepID', 'Humann2', "unstratified", 'All']:
                p_prefix = p_prefix.replace(r, '')
            p_prefix = re.sub('\.+', '.', p_prefix.strip('.|_'))
            p_prefix = p_prefix + '_'
        # self.running_bash = self.out_dir + 'visualize_function.sh'

        if mapping_file and categories and len(categories.split(",")) == 1 and not abundance_table == 'multi_table':

            abd_df = pd.read_csv(abundance_table, sep='\t')
            map_df = pd.read_csv(mapping_file, sep='\t')
            # not_sample = [c not in map_df.iloc[:, 0].values for c in abd_df.columns]
            map_df = map_df.loc[map_df[categories].notna(), :]
            has_group = [c in map_df.iloc[:, 0].values for c in abd_df.columns]
            has_group[0] = True
            # sel = [a or b for a, b in zip(not_sample, has_group)]
            self.mapping_df = map_df
            self.abundance_df = abd_df.loc[:, has_group].groupby(
                abd_df.columns[0]).sum()
            # is_sample = [c in self.mapping_df.iloc[:, 0].values for c in self.abundance_df.columns]
            # pdb.set_trace()
            self.abundance_df = self.abundance_df.loc[self.abundance_df.T.sum(
            ) > 0, :]
            abundance_table = self.tmp_dir + os.path.basename(abundance_table)
            mapping_file = self.tmp_dir + os.path.basename(mapping_file)
            self.mapping_df.to_csv(mapping_file, sep='\t', index=False)
            self.abundance_df.to_csv(abundance_table, sep='\t', index=True)

        self.set_attr(
            prefix=p_prefix,
            categories=categories if categories else "",
            abundance_table=abundance_table,
            mapping_file=mapping_file,
        )
        self.set_path(
            force=False,
            _base_dir=os.path.dirname(__file__) + '/',
            base_dir=os.path.dirname(__file__) + '/',
            **settings.path,
        )

    @staticmethod
    def array_not(array):
        return [not e for e in array]

    @abstractmethod
    def __visualize_with_group__():
        pass

    @abstractmethod
    def __visualize_without_group__():
        pass

    def visualize(self, exclude='none'):
        print("Visualizing abundance table: {}".format(self.abundance_table))
        if self.categories and self.mapping_file:
            print("Visualize using group info...")
            self.__visualize_with_group__(exclude)
        else:
            print("No group info detected, visualize without group info")
            self.__visualize_without_group__()
        result_abundance = self.out_dir + \
            os.path.basename(self.abundance_table)
        if not os.path.exists(result_abundance):
            self.abundance_df.to_csv(result_abundance, sep='\t', index=False)
        if os.path.exists(self.tmp_dir):
            shutil.rmtree(self.tmp_dir)
