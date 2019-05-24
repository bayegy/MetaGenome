import os
import sys
import re
from pyutils.read import format_html_properties


map_html_dir = sys.argv[1]

map_html_list = ['{}/{}'.format(map_html_dir, f) for f in os.listdir(map_html_dir) if f.endswith('.html')]

for map_html in map_html_list:
    print(map_html)
    with open(map_html) as f:
        f = f.read()
    tb = re.search('(<table[\S\s]*</table>)[\S\s]*(<table[\S\s]*</table>)[\S\s]*<table[\S\s]*</table>', f)
    if tb:
        d = tb.group(1)
        m = tb.group(2)
    else:
        tb = re.search('(<table[\S\s]*</table>)[\S\s]*<table[\S\s]*</table>', f)
        if tb:
            d = tb.group(1)
            m = ''
        else:
            print(map_html + ' is not a regular report. ignore it')
            continue

    fo = re.sub('<table[\S\s]*</table>', d + m, f)
    fo = re.sub('display: *none;', '', fo)
    fo = re.sub("window.open\('/", "window.open('https://www.kegg.jp/", fo)
    with open(map_html, 'w') as o:
        o.write(fo)

    link_data = {'*': {'href': 'https://www.kegg.jp{value}'}}
    format_html_properties(map_html, link_data, map_html,
                           filter_function=lambda x: True if x.startswith('/') else False, use_selector=True)
