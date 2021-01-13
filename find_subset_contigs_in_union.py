#!/usr/bin/env python3.8
from MetaGenome.pyutils.read import iter_fa
import argparse
p = argparse.ArgumentParser(
    description="find subset contigs in union contigs.")

p.add_argument('-s', '--subset', required=True,
               help='fasta file containing subset contigs')
p.add_argument('-u', '--union', required=True,
               help='fasta file containing union contigs')
p.add_argument('-o', '--out_file', required=True,
               help='file to save the found match contig names')

options = p.parse_args()

subset_dict = {}
for header, seq in iter_fa(options.subset, trim_line_break=True):
    subset_dict[seq] = True

with open(options.out_file, 'w') as file_handler:
    for header, seq in iter_fa(options.union, trim_line_break=True):
        if subset_dict.get(seq):
            file_handler.write("{}\n".format(header))
