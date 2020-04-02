#!/usr/bin/env python3

from pyutils.read import read_file_n_lines
import argparse
import os
import re


p = argparse.ArgumentParser(
    description="Change the file name suffix of all files under a specific directory. ")


p.add_argument('--input', '-i', help="input contigs fasta",
               default=False, metavar='<path>')
p.add_argument('--minlen', '-l', help="min contigs length to keep in output",
               default=False, metavar='<path>')
p.add_argument('--output', '-o', help="output file path",
               default=False, metavar='<path>')
p.add_argument('--replaceinput', action='store_true', help="replace input files")

options = p.parse_args()


with open(options.output, 'w') as fout:
    i = 1
    for header, seq in read_file_n_lines(options.input, 2):
        if options.minlen:
            if len(seq) <= options.minlen:
                continue

        header = re.sub('^\S+', ">contig_" + str(i), header)
        fout.write(header + seq)
        i += 1

if options.replaceinput:
    os.remove(options.input)
