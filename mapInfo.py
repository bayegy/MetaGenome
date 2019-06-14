import re
import os
import pandas as pd
from numpy import dtype
import gzip


class MapInfo(object):
    """docstring for MapInfo"""

    def __init__(self):
        pass
        # data = data
        # self.map = {}
        # mapping_source_list = mapping_source_list

    def load_map(self, mapping_source, header=False, adjust_func=False, map_column=False):
        self.map = {}

        def add(fid, ldef):
            try:
                self.map[fid].append(ldef)
            except KeyError:
                self.map[fid] = [ldef]
        is_gz = os.path.splitext(mapping_source)[1] == '.gz'
        with (gzip.open(mapping_source, 'r') if is_gz else open(mapping_source, 'r')) as f:
            for number, line in enumerate(f):
                li = re.split('\t', line.decode().strip() if is_gz else line.strip())
                fid = li[0].strip()
                fdef = li[map_column] if map_column else '\t'.join(li[1:])
                ldef = adjust_func(fdef) if adjust_func else fdef
                if number == 0:
                    self.mapped_header = fdef if header else False
                    if not header:
                        add(fid, ldef)
                else:
                    add(fid, ldef)

    @staticmethod
    def ajust_ko_info(info):
        info = re.split('; *', info)
        # gene_name = re.sub('EC?[\.\d]+', '', info[0])
        # gene_name = re.sub('K\d{5}', '', gene_name).strip(' |,')
        gene_name = info[0].strip()
        ez = re.search('\[(E.*)\] *$', info[1])
        enzyme_number = ez.group(1) if ez else ""
        definition = re.sub(' *\[E.*\] *$', '', info[1]).strip()
        return '\t'.join([gene_name, enzyme_number, definition])

    def map_data(self, data, pattern=False, first_pattern=False, add_sid_to_info: bool=False, mapped_header: str=False, out_file=False):
        out = ""
        with open(data, 'r') as f:
            for number, line in enumerate(f):
                line = line.strip('\n')
                if number == 0:
                    out += '{}\t{}\n'.format(line, mapped_header
                                             ) if mapped_header else line + "\t{}\n".format(self.mapped_header or 'Description')
                else:
                    li = re.split('\t', line)
                    sid = re.search(pattern, li[0]).group() if pattern else li[0].strip()
                    li[0] = re.search(first_pattern, li[0]).group() if first_pattern else sid
                    try:
                        mdef = ', '.join(self.map[sid])
                        li.append("{}; {}\n".format(sid, mdef) if add_sid_to_info else mdef + '\n')
                    except Exception as e:
                        li.append("{}; {}\n".format(sid, "")
                                  if add_sid_to_info else "" + '\n')
                    out += '\t'.join(li)
        with open(out_file or data, 'w') as f:
            f.write(out)

    def mapping(self, data, mapping_source_list: list, pattern=False, first_pattern=False, add_sid_to_info=False, header: bool=False, adjust_func=False, mapped_headers: list=False, out_file=False, add: bool=False, map_column=False):
        """
        用来合并几个表格或者给一个表格匹配信息，默认所有表格第一列为匹配的id所在，匹配的表格的id可以不唯一

        参数：

            data: 原表格路径

            mapping_source_list: 一个包含要匹配的表格（可能有多个）路径的list

            pattern: 原表格的id是否需要正则提取（从原表格第一列），如需要，请提供正则表达式

            first_pattern: 结果表格的第一列是否需要正则提取（从原表格第一列），如需要，请提供正则表达式

            add_sid_to_info：是否需要把匹配的id加入匹配得到信息

            header: 要匹配的表格是否有表头

            adjust_func：是否需要用函数对匹配得到信息进行修正，如需要，请提供函数

            mapped_headers：是否给匹配得到的信息（可能有多个），提供表头，如需要，请提供一个包含表头的list

            out_file: 匹配结果文件的名字，如果为False，覆盖原表格

            add: 如果原表格最后一列已经有匹配信息，是否再追加

            map_column: 要匹配的表格的匹配信息所在列，如果为False，匹配所有列


        """

        df = pd.read_csv(data, sep='\t')
        if df.dtypes[-1] == dtype("O") and (not add) and (not out_file):
            print("Last column is string object, MapInfo will treat it as definition, and won't override.")
        else:
            mapped_headers = mapped_headers if mapped_headers else [False] * len(mapping_source_list)
            for number, map_data in enumerate(zip(mapped_headers, mapping_source_list)):
                mapped_header, mapping_source = map_data
                self.load_map(mapping_source, header, adjust_func, map_column)
                if number == 0:
                    self.map_data(data, pattern, first_pattern, add_sid_to_info, mapped_header, out_file)
                else:
                    self.map_data(out_file or data, pattern, first_pattern, add_sid_to_info, mapped_header)
