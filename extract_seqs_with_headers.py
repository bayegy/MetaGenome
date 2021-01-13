#!/usr/bin/env python3.8
from MetaGenome.pyutils.read import iter_fa
import argparse
p = argparse.ArgumentParser(
    description="extract sequences from an fasta file with sequence headers")

p.add_argument('-s', '--headers', required=True,
               help='file containing sequence headers.')
p.add_argument('-f', '--fasta', required=True,
               help='fasta file containing sequences.')
p.add_argument('-o', '--out_file', required=True,
               help='file to save the extracted sequences. fasta format.')
options = p.parse_args()

headers_dict = {}
with open(options.headers) as file_handler:
    for line in file_handler:
        headers_dict[line.split()[0].split('_')[1]] = True


sequences_dict = {}

with open(options.out_file, 'w') as fh:
    for header, seq in iter_fa(options.fasta, trim_line_break=True):
        if sequences_dict.get(seq):
            continue
        if headers_dict.get(header.split()[0].split('_')[1]):
            fh.write("{}\n{}\n".format(header, seq))
            sequences_dict[seq] = True
