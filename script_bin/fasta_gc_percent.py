#!/usr/bin/env python3
import os, sys, re
ms, infile, outfile = sys.argv
with open(infile) as f:
    Dict = {}
    for line in f:
        if line[0] == ">":
            key = re.sub('>', '', line.strip())
            Dict[key] = []
        else:
            Dict[key].append(line.strip())
            
with open(outfile, 'w') as o:
    o.write("id\tGC_percent\n")
    for key, value in Dict.items():
        seq = ''.join(value)
        nG = seq.count("G") + seq.count("g")
        nC = seq.count("C") + seq.count("c")
        gc_percent = (nG + nC)/len(seq)
        o.write("{}\t{}\n".format(key, gc_percent))
