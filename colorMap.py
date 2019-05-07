import sys
import os
import json
import seaborn as sns
import pandas as pd
import numpy as np
import re
from PIL import Image
from pyutils.colors import rgb2hex, hex2color
import pdb


class ColorMap(object):
    """docstring for ColorMap"""

    def __init__(self, ko_lefse_lda, map_abundance_table=False, ko_abundance_table=False, mapping_file=False, category=False, out_dir='./'):
        self.out_dir = os.path.abspath(out_dir) + '/'
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        self._base_dir = os.path.dirname(__file__) + '/'
        with open(self._base_dir + "pipconfig/path.conf") as f:
            self.path = json.load(f)
        sys.path.append(self.path['bayegy_home'])
        from getColors import get_lefse_colors
        self.sig_gps = pd.read_csv(ko_lefse_lda, sep='\t', header=None, index_col=0)[2]
        self.gps_colors = get_lefse_colors(category, mapping_file, ko_lefse_lda, return_dict=True) if (
            mapping_file and category) else {gp: color for gp, color in zip(list(set(self.sig_gps[self.sig_gps.notna()])), sns.color_palette('Accent', 12).as_hex())}
        self.kos_colors = {ko: self.gps_colors[self.sig_gps[ko]]
                           if not self.sig_gps[ko] is np.nan else "#666666" for ko in self.sig_gps.index}

    def cac_map_colors(self, mapid):
        self.current_mapid = mapid
        self.color_data = []
        with open("{}/{}.conf".format(self.path['map_conf'], self.current_mapid), 'r') as f:
            for line in f:
                if line.startswith('rect'):
                    li = line.split('\t')
                    coordinate = re.findall('\d+', li[0])
                    kos = re.findall('K\d+', li[1])
                    colors = [self.kos_colors[ko] for ko in kos if ko in self.kos_colors.keys()]
                    if not len(colors) == 0:
                        colors = list(set([c for c in colors if not c == "#666666"]))
                        color = colors[0] if len(colors) == 1 else "#666666"
                        self.color_data.append([coordinate, color])

    @staticmethod
    def apply(array, func, **kargs):
        """This func will conside only the first and second dimension of array"""
        shape = array.shape
        return np.array([[func(array[i, j], **kargs) for j in range(shape[1])] for i in range(shape[0])])

    def color_map(self):
        plot = np.array(Image.open("{}/{}.png".format(self.path['map_conf'], self.current_mapid)).convert("RGB"))
        pdb.set_trace()
        plot = ColorMap.apply(plot, rgb2hex)
