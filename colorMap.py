import sys
import os
import json
import seaborn as sns
import pandas as pd
import numpy as np
import re
from PIL import Image, ImageFont, ImageDraw
from pyutils.colors import rgb2hex, hex2color
from pyutils.tools import dupply, time_counter
from pyutils.read import read_abundance


class ColorMap(object):
    """
    This class is used to color the KEGG map according to LEfSe analysis results.

    arguments:

        ko_lefse_lda: LEfSe analysis results of KOs

        ko_abundance_table: the abundance_table of KOs, Optional, if not passed, there will not be text in maps.

        map_abundance_table: map abundance table, optional, if not passed, method plot_all can not be used, still can use method plot_map.

        mapping_file: mapping file, optional, if not passed the color will be auto-asigned.

        category: name of category, optional, if not passed the color will be auto-asigned.

        prefix: prefix of out pngs, optional.

        out_dir: where to store the results, optional.

    Sample usage:
        First open the python3 command line.
        to plot all maps in map_abundance_table
	(base) cheng@ps-Super-Server [~/pipelines/MetaGenome]$python3 

            >>>from colorMap import ColorMap
            >>>c=ColorMap("KO_Group1_lefse_LDA2.LDA.txt","All.Function.abundance.KeepID.KO.txt","All.Function.abundance.KeepID.Pathway.txt",out_dir="map_test")
            >>>c.plot_all()

        to plot a single map:

            >>>from colorMap import ColorMap
            >>>c=ColorMap("KO_Group1_lefse_LDA2.LDA.txt")
            >>>c.plot_map("map00010")


    """

    def __init__(self, ko_lefse_lda, ko_abundance_table=False, map_abundance_table=False, mapping_file=False, category=False, prefix="", out_dir='./'):
        self.out_dir = os.path.abspath(out_dir) + '/'
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        self.out_dir += prefix
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
        self.ko_abundance_table = ko_abundance_table
        self.map_abundance_table = map_abundance_table

    def __cac_map_colors__(self, mapid):
        self.current_mapid = mapid
        self.color_data = []
        self.coord_data = []
        with open("{}/{}.conf".format(self.path['map_conf'], self.current_mapid), 'r') as f:
            for line in f:
                if line.startswith('rect'):
                    li = line.split('\t')
                    coordinate = [int(d) for d in re.findall('\d+', li[0])]
                    kos = re.findall('K\d+', li[1])
                    colors = [self.kos_colors[ko] for ko in kos if ko in self.kos_colors.keys()]
                    if not len(colors) == 0:
                        self.coord_data.append([coordinate, kos])
                        colors = list(set([c for c in colors if not c == "#666666"]))
                        color = colors[0] if len(colors) == 1 else "#666666"
                        self.color_data.append([coordinate, color])

    def __color_map__(self):
        plot = np.array(Image.open("{}/{}.png".format(self.path['map_conf'], self.current_mapid)).convert("RGB"))
        plot = dupply(plot, rgb2hex)
        for coordinate, color in self.color_data:
            x1, y1, x2, y2 = coordinate
            plot[y1:y2, x1:x2][plot[y1:y2, x1:x2] == '#ffffff'] = color
        plot = np.uint8(dupply(plot, hex2color))
        self.plot = Image.fromarray(plot)

    def __cac_map_text__(self):
        abundance = read_abundance(self.ko_abundance_table, return_sum=1)
        assert abundance.shape[0] == self.sig_gps.shape[0]
        # pdb.set_trace()
        self.text_data = [[coordinate, "({},{})".format(str(kos[0]), str(
            int(np.ceil(sum([abundance[ko] for ko in kos if ko in abundance.index])))))] for coordinate, kos in self.coord_data]

    def __text_map__(self):
        for coordinate, text in self.text_data:
            x1, y1, x2, y2 = coordinate
            font = ImageFont.truetype('arial.ttf', 9)
            draw = ImageDraw.Draw(self.plot)
            # pdb.set_trace()
            draw.text((x1 - 9, y2 + 1), text, font=font, fill="#FF0000")

    @time_counter
    def plot_map(self, mapid):
        self.__cac_map_colors__(mapid)
        self.__color_map__()
        if self.ko_abundance_table:
            self.__cac_map_text__()
            self.__text_map__()
        self.plot.save("{}{}.png".format(self.out_dir, mapid))

    def plot_all(self):
        if self.map_abundance_table:
            for mapid in read_abundance(self.map_abundance_table).index:
                print(mapid)
                self.plot_map(mapid)
        else:
            print("Please pass a map id table(map_abundance_table)")
