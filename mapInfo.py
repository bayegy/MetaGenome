import re
import pandas as pd
from numpy import dtype


class MapInfo(object):
    """docstring for MapInfo"""

    def __init__(self):
        pass
        # data = data
        # self.map = {}
        # mapping_source_list = mapping_source_list

    def load_map(self, mapping_source, header=False, adjust_func=False):
        self.map = {}

        def add(fid, ldef):
            try:
                self.map[fid].append(ldef)
            except KeyError:
                self.map[fid] = [ldef]

        with open(mapping_source, 'r') as f:
            for number, line in enumerate(f):
                li = re.split('\t', line.strip())
                fid = li[0].strip()
                fdef = '\t'.join(li[1:])
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
                        li.append("{}: {}\n".format(sid, mdef) if add_sid_to_info else mdef + '\n')
                    except Exception as e:
                        li.append("{}: {}\n".format(sid, "")
                                  if add_sid_to_info else "" + '\n')
                    out += '\t'.join(li)
        with open(out_file or data, 'w') as f:
            f.write(out)

    def mapping(self, data, mapping_source_list: list, pattern=False, first_pattern=False, add_sid_to_info=False, header: bool=False, adjust_func=False, mapped_headers: list=False, out_file=False, add: bool=False):
        df = pd.read_csv(data, sep='\t')
        if df.dtypes[-1] == dtype("O") and (not add) and (not out_file):
            print("Last column is string object, MapInfo will treat it as definition, and won't override.")
        else:
            mapped_headers = mapped_headers if mapped_headers else [False] * len(mapping_source_list)
            for number, map_data in enumerate(zip(mapped_headers, mapping_source_list)):
                mapped_header, mapping_source = map_data
                self.load_map(mapping_source, header, adjust_func)
                if number == 0:
                    self.map_data(data, pattern, first_pattern, add_sid_to_info, mapped_header, out_file)
                else:
                    self.map_data(out_file or data, pattern, first_pattern, add_sid_to_info, mapped_header)
