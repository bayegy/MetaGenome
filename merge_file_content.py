#!/usr/bin/env python3.6

import os
import sys
import argparse


p = argparse.ArgumentParser(
    description="This script is used for the analysis of metagenome sequencing.")
p.add_argument('-i', '--input', dest='input', metavar='<directory>', default=False,
               help='Directory of files')
p.add_argument('-f', '--filter', dest='filter', metavar='<str>', default="txt,pdf,pyc,tsv,xls,pm,git",
               help='the type of file excepted to exclude from merge.')
p.add_argument('-s', '--space', dest='space', metavar='<int>', default=2,
               help='The space size  between files content')
p.add_argument('-o', '--outfile', dest='outfile', default="merged_conten.txt", metavar='<file name>',
               help="Output file name")
options = p.parse_args()


# sp, indir, filter_type, outfile = sys.argv


def iswanted(file):
    for t in options.filter.split(','):
        if file.endswith(t):
            return False
    return True


with open(options.outfile, 'w', encoding='utf-8') as out:
    space = '\n' * int(options.space)
    for root, dirs, files in os.walk(options.input):
        for f in files:
            if iswanted(f):
                fp = f"{root}/{f}"
                print(fp)
                with open(fp) as fr:
                    try:
                        out.write(f"{space}Script {fp} file content:\n{space}")
                        out.write(fr.read())
                    except Exception as e:
                        pass
