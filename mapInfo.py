import re
import pandas as pd
from numpy import dtype


class MapInfo(object):
    """docstring for MapInfo"""

    def __init__(self, data, mapping_source):
        self.data = data
        self.mapping_source = mapping_source
        self.map = {}

    def load_map(self):
        with open(self.mapping_source, 'r') as f:
            for line in f:
                try:
                    li = re.split('\t', line.strip())
                    self.map[li[0]] = li[-1]
                except Exception as e:
                    print("Can not map at: " + line)

    def map_data(self, pattern=False, first_pattern=False, add_sid_to_info=False):
        out = ""
        with open(self.data, 'r') as f:
            for number, line in enumerate(f):
                line = line.strip()
                if number == 0:
                    out += line + "\tDefinition\n"
                else:
                    li = re.split('\t', line)
                    sid = re.search(pattern, li[0]).group() if pattern else li[0]
                    li[0] = re.search(first_pattern, li[0]).group() if first_pattern else li[0]
                    try:
                        li.append("{}: {}\n".format(sid, self.map[sid]) if add_sid_to_info else self.map[sid] + '\n')
                    except Exception as e:
                        li.append("{}: {}\n".format(sid, "No definition")
                                  if add_sid_to_info else "No definition" + '\n')
                    out += '\t'.join(li)
        with open(self.data, 'w') as f:
            f.write(out)

    def mapping(self, pattern=False, first_pattern=False, add_sid_to_info=False):
        df = pd.read_csv(self.data, sep='\t')
        if df.dtypes[-1] == dtype("O"):
            print("Last column is string object, MapInfo will treat it as definition, and won't override.")
        else:
            self.load_map()
            self.map_data(pattern, first_pattern, add_sid_to_info)
