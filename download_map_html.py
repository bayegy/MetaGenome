import sys
import os
from myutils.net import Net
from myutils.progress_bar import ProgressBar

maplist, outdir = sys.argv[1], sys.argv[2]

n = Net()
if not os.path.exists(outdir):
    os.makedirs(outdir)
with open(maplist, 'r') as f:
    map_list = f.read().split('\n')
    bar = ProgressBar(len(map_list))
    for mapid in map_list:
        out_html = '{}/{}.html'.format(outdir, mapid)
        if not os.path.exists(out_html):
            print(mapid)
            result = n.lrequests('https://www.kegg.jp/kegg-bin/show_pathway?' + mapid, return_tree=False)
            with open(out_html, 'w') as o:
                o.write(result)
            bar.move()
        else:
            bar.move()
