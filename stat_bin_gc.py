#!/usr/bin/env python3.8

import os
import sys
import re

sp, indir, outfile = sys.argv

with open(outfile, 'w') as outfile:
    outfile.write("BinID\tSize\tGC_percent\n")
    bins = os.listdir(indir)
    for bi in bins:
        bi_path = "{}/{}".format(indir, bi)
        if os.path.isfile(bi_path):
            with open(bi_path) as bif:
                seq = bif.read()
            seq = re.sub('>[^\n]*\n', '', seq)
            seq = re.sub(r'\s', '', seq)
            nG = seq.count('G')
            nC = seq.count('C')
            total = len(seq)
            p_GC = (nC + nG) / total
            bin_n = re.sub(r'\.fa$', '', bi)
            outfile.write("{}\t{}\t{}\n".format(bin_n, total, p_GC))
