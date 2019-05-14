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
from mapInfo import MapInfo
import pdb


class ColorMap(object):
    """
    This class is used to color the KEGG map according to LEfSe analysis results.

    arguments:

        ko_lefse_lda: LEfSe analysis results of KOs

        ko_abundance_table: the abundance_table of KOs, Optional, if not passed, abundance of kos will not be showed in maps.

        mapping_file: mapping file, optional, if not passed the color will be auto-asigned.

        category: name of category, optional, if not passed the color will be auto-asigned.

        prefix: prefix of out output file name, optional.

        out_dir: where to store the results, optional.

    Sample usage:
        First open the python3 command line.
        (base) cheng@ps-Super-Server [~/pipelines/MetaGenome]$python3

        to plot all maps with abundance showed in map:

            >>>from colorMap import ColorMap
            >>>c=ColorMap("test/KO_Group1_lefse_LDA2.LDA.txt","test/All.Function.abundance.KeepID.KO.txt",out_dir="map_test")
            >>>c.plot_all()

        to plot a single map:

            >>>from colorMap import ColorMap
            >>>c=ColorMap("test/KO_Group1_lefse_LDA2.LDA.txt")
            >>>c.plot_map("map00010")


    """

    def __init__(self, ko_lefse_lda, ko_abundance_table=False, mapping_file=False, category=False, prefix="", out_dir='./'):
        self.out_dir = os.path.abspath(out_dir) + '/'
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        self.out_dir += prefix
        self._base_dir = os.path.dirname(__file__) + '/'
        with open(self._base_dir + "pipconfig/path.conf") as f:
            self.path = json.load(f)
        sys.path.append(self.path['bayegy_home'])
        from getColors import get_lefse_colors
        self.user_kos = pd.read_csv(ko_lefse_lda, sep='\t', header=None, index_col=0)[2]
        self.gps_colors = get_lefse_colors(category, mapping_file, ko_lefse_lda, return_dict=True) if (
            mapping_file and category) else {gp: color for gp, color in zip(list(set(self.user_kos[self.user_kos.notna()])), sns.color_palette('Accent', 12).as_hex())}
        # pdb.set_trace()
        self.kos_colors = {ko: self.gps_colors[self.user_kos[ko]]
                           if notna else "#999999" for ko, notna in zip(self.user_kos.index, self.user_kos.notna())}
        self.ko_abundance_table = ko_abundance_table
        maps = []
        mi = MapInfo()
        self.prefix = prefix
        mi.load_map(self.path['fmap_home'] + '/FMAP_data/KEGG_orthology2pathway.txt')
        for ko in self.user_kos.index:
            try:
                maps += mi.map[ko]
            except KeyError:
                pass
        self.maps = list(set(maps))

    def get_map_conf(self, mapid):
        self.current_mapid = mapid
        self.coord_kos, self.coord_gene_names, self.coord_enzyme, self.coord_reaction = [], [], [], []
        with open("{}/{}.conf".format(self.path['map_conf'], self.current_mapid), 'r') as f:
            for line in f:
                if line.startswith('rect'):
                    li = line.split('\t')
                    coordinate = [int(d) for d in re.findall('\d+', li[0])]
                    p_kos = re.findall('[^\(]*(K\d+)[^\)]*', li[2])
                    p_gene_names = re.findall('\(([^\(\)]+)\)', li[2])
                    kos, gene_names = [], []
                    for ko, name in zip(p_kos, p_gene_names):
                        if ko in self.user_kos.index:
                            kos.append(ko)
                            gene_names.append(name)
                    if not len(kos) == 0:
                        enzyme = [i.strip() for i in re.findall(',([^,]+),[^,]*$', li[2])]
                        reaction = re.findall('(R\d+) *$', li[2])
                        self.coord_enzyme.append([coordinate, enzyme])
                        self.coord_reaction.append([coordinate, reaction])
                        self.coord_kos.append([coordinate, kos])
                        self.coord_gene_names.append([coordinate, gene_names])

    def __cac_map_colors__(self):
        self.color_data = []
        for coordinate, kos in self.coord_kos:
            colors = [self.kos_colors[ko] for ko in kos]
            colors = list(set([c for c in colors if not c == "#999999"]))
            color = colors[0] if len(colors) == 1 else "#999999"
            self.color_data.append([coordinate, color])

    def __color_map__(self, override=True):
        """
            override: When Ture, override the enzyme name in map.png, this is faster.
        """
        plot = np.array(Image.open("{}/{}.png".format(self.path['map_conf'], self.current_mapid)).convert("RGB"))
        if override:
            for coordinate, color in self.color_data:
                x1, y1, x2, y2 = coordinate
                # pdb.set_trace()
                plot[y1:y2, x1:x2] = hex2color(color)
        else:
            plot = dupply(plot, rgb2hex)
            for coordinate, color in self.color_data:
                x1, y1, x2, y2 = coordinate
                plot[y1:y2, x1:x2][plot[y1:y2, x1:x2] == '#ffffff'] = color
            plot = np.uint8(dupply(plot, hex2color))
        self.plot = Image.fromarray(plot)

    def __cac_map_text__(self, use_text="gene", show_abundance=True):
        """
        aguments:
            use_text: 'gene' or 'enzyme' or 'ko' or 'reaction'
        """
        text_data = self.coord_gene_names if use_text == "gene" else (
            self.coord_enzyme if use_text == "enzyme" else (self.coord_kos if use_text == "ko" else self.coord_reaction))
        if self.ko_abundance_table:
            self.text_data = []
            abundance = read_abundance(self.ko_abundance_table, return_sum=1)
            # assert abundance.shape[0] == self.user_kos.shape[0]
            # pdb.set_trace()
            for coord_text, coord_kos in zip(text_data, self.coord_kos):
                coord, text = coord_text
                kos_abundance = [abundance[ko] if ko in abundance.index else 0 for ko in coord_kos[1]]
                if len(text) > 1 and len(text) == len(kos_abundance):
                    max_abundance, showtext = kos_abundance[0], text[0]
                    for abd, txt in zip(kos_abundance, text):
                        if abd > max_abundance:
                            max_abundance, showtext = abd, txt
                else:
                    showtext = text[0] if text else ""
                self.text_data.append([coord, "({},{})".format(showtext, str(
                    int(np.ceil(sum(kos_abundance))))) if show_abundance else showtext])
        else:
            self.text_data = [[coordinate, texts[0] if texts else ""] for coordinate, texts in text_data]

    def __text_map__(self, position="center", color="#000000", fontsize=9):
        """
        aguments:
            position: 'bottom' or 'center' or 'top'.
        """
        for coordinate, text in self.text_data:
            x1, y1, x2, y2 = coordinate
            font = ImageFont.truetype('arial.ttf', fontsize)
            width, height = font.getsize(text)
            offsetx, offsety = font.getoffset(text)
            w, h = width + offsetx, height + offsety

            draw = ImageDraw.Draw(self.plot)
            # pdb.set_trace()
            cord = ((x1 + x2 - w) / 2, y2) if position == "bottom" else (((x1 + x2 - w) / 2, (y1 + y2 - h) / 2)
                                                                         if position == "center" else ((x1 + x2 - w) / 2, y1 - h))
            # pdb.set_trace()
            draw.text(cord, text, font=font, fill=color)

    @time_counter
    def plot_map(self, mapid, use_text="gene", position="center", color="#000000", fontsize=9, show_abundance=False):
        """
        To plot a single map:

        arguments:

            mapid: map id

            see also plot_all
        """
        self.get_map_conf(mapid)
        self.__cac_map_colors__()
        self.__color_map__()
        self.__cac_map_text__(use_text, show_abundance)
        self.__text_map__(position, color, fontsize)
        self.plot.save("{}{}{}.png".format(self.out_dir, self.prefix, mapid))

    def plot_all(self, use_text="gene", position="center", color="#000000", fontsize=9, show_abundance=False):
        """
        To plot all maps in map_abundance_table

            argument:
                position: 'bottom' or 'center' or 'top'.

                use_text: 'gene' or 'enzyme' or 'ko' or 'reaction'

                show_abundance: if True, show abundance in plot

        """
        filter_maps = ["map01100", "map01110", "map01120", "map01130", "map00312"]
        for mapid in self.maps:
            if mapid not in filter_maps:
                print(mapid)
                self.plot_map(mapid, use_text, position, color, fontsize, show_abundance)
