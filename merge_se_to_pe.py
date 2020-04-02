#!/usr/bin/env python3
import os
# import io
import argparse
from pyutils.read import read_file_n_lines

p = argparse.ArgumentParser(
    description="Change the file name suffix of all files under a specific directory. ")


p.add_argument('--se1', '-1', help="reads 1",
               default=False, metavar='<path>')
p.add_argument('--se2', '-2', help="reads 2",
               default=False, metavar='<path>')
p.add_argument('--out', '-o', help="The %s in path will be replaced by 1 or 2.",
               default=False, metavar='<path>')
p.add_argument('--replaceinput', action='store_true', help="replace input files")


options = p.parse_args()


se1 = {}

for lines in read_file_n_lines(options.se1, 4):
    k = lines[0].split("/")[0]
    se1[k] = lines


with open(options.out % (1), 'w') as pe1, open(options.out % (2), 'w') as pe2:
    # with io.TextIOWrapper(ppe1, encoding='ascii') as pe1, io.TextIOWrapper(ppe2, encoding='ascii') as pe2:
    for lines in read_file_n_lines(options.se2, 4):
        k = lines[0].split("/")[0]
        if k in se1:
            pe2.write("".join(lines))
            pe1.write("".join(se1[k]))


if options.replaceinput:
    os.remove(options.se1)
    os.remove(options.se2)

print("SE merged to PE")
