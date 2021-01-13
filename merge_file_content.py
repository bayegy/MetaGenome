#!/usr/bin/env python3.8

import os
# import sys
import argparse


p = argparse.ArgumentParser(
    description="This script is used for the analysis of metagenome sequencing.")
p.add_argument('-i', '--input', dest='input', metavar='<directory>', default=False,
               help='Directory of files')
p.add_argument('-e', '--exclude', dest='exclude', metavar='<str>', default="txt,pdf,pyc,tsv,xls,pm,git",
               help='the type of file excepted to exclude from merge.')
p.add_argument('-I', '--include', dest='include', metavar='<str>', default=False,
               help='the type of file excepted to include to merge.')
p.add_argument('-s', '--space', dest='space', metavar='<int>', default=2,
               help='The space size  between files content')
p.add_argument('-o', '--outfile', dest='outfile', default="merged_conten.txt", metavar='<file name>',
               help="Output file name")
options = p.parse_args()


# sp, indir, filter_type, outfile = sys.argv


def filter_exclude(file):
    for t in options.exclude.split(','):
        if file.endswith(t):
            return False
    return True


def filter_include(file):
    for t in options.include.split(','):
        if file.endswith(t):
            return True
    return False


def filter(file):
    if options.include:
        return filter_include(file)
    return filter_exclude(file)


with open(options.outfile, 'w', encoding='utf-8') as out:
    space = '\n' * int(options.space)
    for root, dirs, files in os.walk(options.input):
        for f in files:
            if filter(f):
                fp = f"{root}/{f}"
                print(fp)
                with open(fp) as fr:
                    try:
                        out.write(f"{space}Script {fp} file content:\n{space}")
                        out.write(fr.read())
                    except Exception:
                        pass
