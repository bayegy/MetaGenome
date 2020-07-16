import os
import re
from abc import ABCMeta, abstractmethod
import pandas as pd
import shutil
from pipconfig import settings
from systemMixin import SystemMixin


class Visualize(SystemMixin, metaclass=ABCMeta):
    """docstring for Visualize"""

    def __init__(self, abundance_table, mapping_file=False, categories=False, prefix=False, out_dir=False, filter_abc_description=False, normalize_abc_description=False):
        self.valid = True
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
            p_prefix = re.sub(r'\.+', '.', p_prefix.strip('.|_'))
            p_prefix = p_prefix + '_'
        # self.running_bash = self.out_dir + 'visualize_function.sh'

        if mapping_file and categories and len(categories.split(",")) == 1 and not abundance_table == 'multi_table':

            abd_df = pd.read_csv(abundance_table, sep='\t')
            map_df = pd.read_csv(mapping_file, sep='\t')
            not_sample = [c not in map_df.iloc[:,
                                               0].values for c in abd_df.columns]
            abd_df_other = abd_df.iloc[:, not_sample]
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

            """
            # some sample do not contain any genes
            if sum(self.abundance_df.sum() == 0) > 0:
                self.abundance_df = self.abundance_df + 1
            """

            if sum(not_sample) > 1:
                abd_df_other = abd_df_other.iloc[self.array_not(
                    abd_df_other.iloc[:, 0].duplicated()), :]
                abd_df_other.index = abd_df_other.iloc[:, 0]
                abd_df_other.drop(
                    abd_df_other.columns[0], axis=1, inplace=True)
                self.abundance_df = self.abundance_df.join(abd_df_other)
                # filter abundance table description
                if filter_abc_description:
                    flag, ftaxa = filter_abc_description.split(':')
                    ftaxa = ftaxa.split(',')
                    if flag == "keep":
                        self.abundance_df = self.abundance_df.iloc[self.iterfind(
                            self.abundance_df.iloc[:, -1], ftaxa), :]
                    else:
                        self.abundance_df = self.abundance_df.iloc[[
                            not b for b in self.iterfind(self.abundance_df.iloc[:, -1], ftaxa)], :]

                if normalize_abc_description:
                    self.abundance_df.iloc[:, -1] = list(map(lambda x: x.replace("; ", ";").replace(" ", "_"),
                                                             self.abundance_df.iloc[:, -1]))
            nrow, ncol = self.abundance_df.shape
            if nrow < 2:
                print("Warning: features number is less than 2")
                self.valid = False

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

    def visualize(self, exclude='none', **kwargs):
        if not self.valid:
            print("Warning: visualize function was passed!")
            return
        print("Visualizing abundance table: {}".format(self.abundance_table))
        if self.categories and self.mapping_file:
            print("Visualize using group info...")
            self.__visualize_with_group__(exclude, **kwargs)
        else:
            print("No group info detected, visualize without group info")
            self.__visualize_without_group__(**kwargs)
        result_abundance = self.out_dir + \
            os.path.basename(self.abundance_table)
        if not os.path.exists(result_abundance) and hasattr(self, "abundance_df"):
            self.abundance_df.to_csv(result_abundance, sep='\t', index=True)
        if os.path.exists(self.tmp_dir):
            shutil.rmtree(self.tmp_dir)

    def set_colors(self, colors=False):
        if colors:
            colors_list = '\n'.join([colors[k] for k in sorted(colors.keys(), key=str.lower)])
            colors_list_file = os.path.join(self.out_dir, 'group_color.list')
            with open(colors_list_file, 'w') as f:
                f.write(colors_list)
        else:
            colors_list_file = "{bayegy_home}/piputils/group_color.list".format(**self.context)
        self.system(
            "{bayegy_home}/piputils/write_colors_plan.py -i {mapping_file} -c {categories} \
            -p {colors_list_file} -o {out_dir}/colors_plan.json", colors_list_file=colors_list_file)
        os.environ['COLORS_PLAN_PATH'] = os.path.join(self.out_dir, 'colors_plan.json')

    @staticmethod
    def iterfind(a, b):
        out = []
        for i in a:
            found = False
            for j in b:
                if re.search(j, i, flags=re.IGNORECASE):
                    found = True
                    break
            out.append(found)
        return out
