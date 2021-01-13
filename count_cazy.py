#!/usr/bin/env python3.8
# -*- coding: utf-8 -*-
import re
import sys

sp, *infiles = sys.argv

cazys = {}
queries = []
for file in infiles:
    with open(file) as inf:
        for line in inf:
            czs = line.split()
            cz = czs[1].split('|')[1].split('_')[0]
            if czs[0] not in queries:
                if cz in cazys:
                    cazys[cz] += 1
                else:
                    cazys[cz] = 1
            queries.append(czs[0])

print("entries\tgroup\tcount")
for k, v in cazys.items():
    g = re.search(r'^\D+', k).group()
    print("{}\t{}\t{}".format(k, g, v))
