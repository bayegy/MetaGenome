import os
import argparse

p = argparse.ArgumentParser(
    description="Change the file name suffix of all files under a specific directory. ")


p.add_argument('--se1', '-1', help="reads 1",
               default=False, metavar='<path>')
p.add_argument('--se2', '-2', help="reads 2",
               default=False, metavar='<path>')
p.add_argument('--out', '-o', help="The %s in path will be replaced by 1 or 2.",
               default=False, metavar='<path>')
p.add_argument('--removeinput', action='store_true', help="remove input files")


options = p.parse_args()


def read_file_n_lines(file, n):
    """ Read a file n lines at a time """

    line_set = []
    with open(file) as file_handle:
        for line in file_handle:
            if len(line_set) == n:
                yield line_set
                line_set = []
            line_set.append(line)

    # yield the last set
    if len(line_set) == n:
        yield line_set


se1 = {}

for lines in read_file_n_lines(options.se1, 4):
    k = lines[0].split("/")[0]
    se1[k] = lines


with open(options.out % (1), 'w') as pe1, open(options.out % (2), 'w') as pe2:
    for lines in read_file_n_lines(options.se2, 4):
        k = lines[0].split("/")[0]
        if k in se1:
            pe2.write("".join(lines))
            pe1.write("".join(se1[k]))


if options.removeinput:
    os.remove(options.se1)
    os.remove(options.se2)