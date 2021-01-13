#!/usr/bin/env python3

from MetaGenome.pyutils.read import iter_fa
import argparse
import os
import re


p = argparse.ArgumentParser(
    description="rename the contigs name to contig_1, contig_2, and so on")


p.add_argument('--inputs', '-i', nargs='+', help="input fasta files.",
               default=False, metavar='<path>')
p.add_argument('--minlen', '-l', help="min contigs length to keep in output",
               default=False, metavar='<path>')
p.add_argument('--replaceinput', action='store_true', help="replace input files")

options = p.parse_args()

i = 1
for ipt in options.inputs:
    with open(ipt + '.renamed', 'w') as fout:
        for header, seq in iter_fa(ipt):
            if options.minlen:
                if len(seq) <= options.minlen:
                    continue

            header = re.sub(r'^\S+', ">contig_" + str(i), header)
            fout.write(header + seq)
            i += 1

if options.replaceinput:
    for ipt in options.inputs:
        os.remove(ipt)
        os.rename(ipt + '.renamed', ipt)
