#!/usr/bin/env python3
import sys, re, os
ms, infile, outfile = sys.argv
with open(outfile, 'w') as o :
    with open(infile) as f:
        num = 1
        for line in f:
            line = line.strip()
            cell = re.split(r'\t', line)
            if num == 1 :
                o.write("<tr><th>{}</th><th>{}</th><th>{}</th><th>{}</th><th>{}</th><th>{}</th><th>{}</th><th>{}</th><th>{}</th><tr>\n".format(cell[0], cell[1], cell[2], cell[3], cell[4], cell[5], cell[6], cell[7], cell[8]))
                num += 1
            else:
                o.write("\t\t\t\t\t<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><tr>\n".format(cell[0], cell[1], cell[2], cell[3], cell[4], cell[5], cell[6], cell[7], cell[8]))
