#!/usr/bin/env python3
import sys
import os
import json
import seaborn as sns
import pandas as pd
import numpy as np
import re
import cv2 as cv
from PIL import Image, ImageFont, ImageDraw
from pyutils.colors import rgb2hex, hex2color
from pyutils.tools import dupply, time_counter
from pyutils.read import read_abundance
from mapInfo import MapInfo
from pyutils.read import update_html_properties
try:
    from pipconfig import settings
except ImportError:
    pass

# import pdb


class ColorMap(object):
    """
    This class is used to color the KEGG map according to LEfSe analysis results.

    arguments:

        feature_list_path: LEfSe analysis results of KOs

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
            >>>c.plot_all(show_abundance=True)

        to plot a single map:

            >>>from colorMap import ColorMap
            >>>c=ColorMap("test/KO_Group1_lefse_LDA2.LDA.txt")
            >>>c.plot_map("map00010")

    """

    def __init__(self, feature_list_path, map_conf_path=False, colors=False, column=0, ko_abundance_table=False, prefix="", out_dir='./'):
        self.out_dir = os.path.abspath(out_dir) + '/'
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        self.out_dir += prefix
        self._base_dir = os.path.dirname(__file__) + '/'
        user_kos = pd.read_csv(feature_list_path, sep=',', header=0, index_col=0)
        # pdb.set_trace()
        self.user_kos = user_kos.iloc[:, column] if isinstance(column, int) else user_kos.loc[:, column]
        is_compound = [i.startswith('C') for i in self.user_kos.index]
        self.is_compound = sum(is_compound) == len(is_compound)
        if self.is_compound:
            print("Inputs are compounds.")
        self.annoted_kos = self.user_kos[self.user_kos.notna()]
        if isinstance(colors, dict):
            self.gps_colors = colors
        elif isinstance(colors, str):
            gps, cols = colors.split(";")
            self.gps_colors = dict(((gp, col) for gp, col in zip(gps.split(","), cols.split(","))))
        else:
            self.gps_colors = dict(zip(list(set(self.annoted_kos)), sns.color_palette('Accent', 12).as_hex()))
        # pdb.set_trace()
        self.kos_colors = {ko: self.gps_colors[self.user_kos[ko]]
                           if notna else "#999999" for ko, notna in zip(self.user_kos.index, self.user_kos.notna())}
        self.ko_abundance_table = ko_abundance_table
        self.map_conf_path = map_conf_path if map_conf_path else settings.path['map_conf']
        self.prefix = prefix

    def get_map_conf(self, mapid, margin_right=60, clean_frame=True):
        self.current_mapid = mapid
        plot = np.array(Image.open("{}/{}.png".format(self.map_conf_path, self.current_mapid)).convert("RGB"))
        if margin_right:
            margin = np.zeros((plot.shape[0], margin_right, 3), dtype='uint8')
            margin[:] = 255
            if clean_frame:
                plot[:, (0, plot.shape[1] - 1)] = [255, 255, 255]
                plot[(0, plot.shape[0] - 1), :] = [255, 255, 255]
            self.plot = np.concatenate((plot, margin), axis=1)
        else:
            self.plot = plot
        with open("{}/{}.conf".format(self.map_conf_path, self.current_mapid), 'r') as f:
            self.coord_kos = []
            if self.is_compound:
                for line in f:
                    if line.startswith("circ"):
                        li = line.strip().split("\t")
                        ko = re.findall('^C\d+', li[2].strip())
                        if ko and ko[0] in self.user_kos.index:
                            coordinate = [int(d) for d in re.findall('\d+', li[0])]
                            coordinate = coordinate[:2]
                            self.coord_kos.append([coordinate, ko])
            else:
                self.coord_enzyme, self.coord_reaction, self.coord_gene_names = [], [], []
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
        color_used = set()
        self.legend_data = []
        for coordinate, kos in self.coord_kos:
            # pdb.set_trace()
            colors = [self.kos_colors[ko] for ko in kos]
            colors = list(set([c for c in colors if not c == "#999999"]))
            color = colors[0] if len(colors) == 1 else "#999999"
            self.color_data.append([coordinate, color])
            # for legend
            if not color == "#999999":
                color_used.add(color)
        for clr in list(color_used):
            for k, v in self.gps_colors.items():
                if clr == v:
                    self.legend_data.append([k, v])
        self.legend_data.sort()

    def __color_map__(self, color_data, override=True):
        """
            override: When Ture, override the enzyme name in map.png, this is faster.
        """
        plot = np.array(self.plot) if not isinstance(self.plot, np.ndarray) else self.plot.copy()
        if self.is_compound:
            for coordinate, color in color_data:
                cv.circle(plot, tuple(coordinate), 6, tuple(hex2color(color)), -1)
        else:
            # color_data = color_data or self.color_data
            if override:
                for coordinate, color in color_data:
                    x1, y1, x2, y2 = coordinate
                    # pdb.set_trace()
                    plot[y1:y2, x1:x2] = hex2color(color)
            else:
                plot = dupply(plot, rgb2hex)
                for coordinate, color in color_data:
                    x1, y1, x2, y2 = coordinate
                    plot[y1:y2, x1:x2][plot[y1:y2, x1:x2] == '#ffffff'] = color
                plot = np.uint8(dupply(plot, hex2color))
        self.plot = plot

    def __cac_map_text__(self, use_text="gene", show_abundance=True):
        """
        aguments:
            use_text: 'gene' or 'enzyme' or 'ko' or 'reaction'
        """
        if self.is_compound:
            print("Text is not supported for compounds")
        else:
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

    def __text_map__(self, text_data, position="center", color="#000000", fontsize=9):
        """
        aguments:
            position: 'bottom' or 'center' or 'top' or 'left'.

            text_data self-sapplied text data
        """
        plot = Image.fromarray(self.plot) if not isinstance(self.plot, Image.Image) else self.plot
        # text_data = text_data or self.text_data
        for coordinate, text in text_data:
            x1, y1, x2, y2 = coordinate
            font = ImageFont.truetype('arial.ttf', fontsize)
            width, height = font.getsize(text)
            offsetx, offsety = font.getoffset(text)
            w, h = width + offsetx, height + offsety

            draw = ImageDraw.Draw(plot)
            # pdb.set_trace()
            cord = ((x1 + x2 - w) / 2, y2) if position == "bottom" else (((x1 + x2 - w) / 2, (y1 + y2 - h) / 2)
                                                                         if position == "center" else (((x1 + x2 - w) / 2, y1 - h) if position == 'top' else (x1, (y1 + y2 - h) / 2)))
            # pdb.set_trace()
            draw.text(cord, text, font=font, fill=color)
        self.plot = plot

    def __cac_legend__(self, off_right=120, legend_data=False):
        """
        position: coordinate of the left top point of legend.

        """
        position = (np.array(self.plot).shape[1] - off_right, 20)
        self.legend_color, self.legend_text = [], []
        legend_data = legend_data or self.legend_data
        rect = np.array([position, [position[0] + 46, position[1] + 17]])
        for gp, color in legend_data:
            rect_text = rect.copy()
            rect_text[:, 0] += 56  # interval of text to rect
            if self.is_compound:
                self.legend_color.append([[int(rect[:, 0].mean()), int(rect[:, 1].mean())], color])
            else:
                self.legend_color.append([rect.flatten().tolist(), color])
            self.legend_text.append([rect_text.flatten().tolist(), gp])
            rect[:, 1] += 25   # interval of rect to rect

    def write_report(self, mapid, report_detail=True):
        # pdb.set_trace()
        link_data = {'img[name=pathwayimage]': {"src": '{}{}.png'.format(self.prefix, mapid)}}

        out_report = "{}{}.html".format(self.out_dir, mapid)
        in_report = "{}/{}.html".format(self.map_conf_path, mapid)

        # update_html_properties(in_report, link_data, out_report)
        if report_detail and not self.is_compound:
            tip_data = {
                'area[coords={}]'.format(','.join([str(co) for co in coord])): {'title': '''{value}

The following KOs were found in your samples[KO number(Group of feature)]:

    %s''' % (', '.join([("{}(biomarker of group {})".format(ko, self.annoted_kos[ko]) if ko in self.annoted_kos.index else ko) for ko in kos]))} for coord, kos in self.coord_kos
            }

            link_data.update(tip_data)

        update_html_properties(in_report, link_data, out_report)

    def show(self):
        plot = Image.fromarray(self.plot) if not isinstance(self.plot, Image.Image) else self.plot
        plot.show()

    def save(self, fp):
        plot = Image.fromarray(self.plot) if not isinstance(self.plot, Image.Image) else self.plot
        plot.save(fp)

    @time_counter
    def plot_map(self, mapid, use_text="gene", position="center", color="#000000", fontsize=9, show_abundance=False, legend_fontsize=12, margin_right=80, off_right=140, report_detail=True):
        """
        To plot a single map:

        arguments:

            mapid: map id

            position: 'bottom' or 'center' or 'top'.

            use_text: 'gene' or 'enzyme' or 'ko' or 'reaction'

            show_abundance: if True, show abundance in plot

            fontsize: font size in map

            legend_fontsize: font size of legend text

            margin_right: the margin at right

            off_right: adjust the postion of legend

            report_detail: report details of KOs in map if True
        """
        self.get_map_conf(mapid, margin_right)
        self.__cac_map_colors__()
        self.__cac_legend__(off_right=off_right)
        self.__color_map__(self.color_data)
        self.__color_map__(self.legend_color)
        self.__text_map__(self.legend_text, 'left', "#000000", legend_fontsize)
        if not self.is_compound:
            self.__cac_map_text__(use_text, show_abundance)
            self.__text_map__(self.text_data, position, color, fontsize)
        self.plot.save("{}{}.png".format(self.out_dir, mapid))
        self.write_report(mapid, report_detail)


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(
        description="This script is used to plot RDA of species. The numeric enviroment factors must be encluded in maping file. The categories will be filterd before RDA")
    p.add_argument('-i', '--input', dest='input', metavar='<path>',
                   help='Feature list file path')
    p.add_argument('-m', '--mapid', dest='mapid', metavar='<str>',
                   help='Map id')
    p.add_argument('-d', '--database', dest='database', metavar='<path>', default="/home/cheng/Databases/map",
                   help='Database path')
    p.add_argument('-o', '--output', dest='output', metavar='<directory>', default='./',
                   help='Output directory')
    p.add_argument('-c', '--colors', dest='colors', metavar='<str>', default=False,
                   help='Colors with the format of "group1,group2;color1,color2"')
    p.add_argument('-n', '--column', dest='column', metavar='<str>', default=0,
                   help='Which column specify the feature classes')
    p.add_argument('-p', '--prefix', dest='prefix', metavar='<int>', default="",
                   help='The prefix of output files, default if null')
    options = p.parse_args()
    c = ColorMap(feature_list_path=options.input, map_conf_path=options.database, column=options.column,
                 out_dir=options.output, prefix=options.prefix, colors=options.colors)
    c.plot_map(options.mapid)
