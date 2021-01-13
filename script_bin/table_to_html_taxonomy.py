#!/usr/bin/env python3
import sys, re, os
ms, infile, outfile = sys.argv
with open(outfile, 'w') as o :
    with open(infile) as f:
        num = 1
        for line in f:
            cell = re.split(r'\t', line)
            if num == 1 :
                cell[0] = "BinID"
                cell[1] = "Taxonomy"
                o.write("<tr><th>{}</th><th>{}</th><tr>\n".format(cell[0], cell[1]))
                num += 1
            else:
                o.write("\t\t\t\t\t<tr><td>{}</td><td>{}</td><tr>\n".format(cell[0], cell[1]))
