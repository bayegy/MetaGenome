#!/usr/bin/env python3

import sys
import json


current_path = []


def get_all_children(tree):
    for child in tree['children']:
        current_path.append(child['name'])
        if 'children' in child:
            get_all_children(child)
        else:
            entry, name = child['name'].strip().split(' ', 1)
            current_path[-1] = name
            print("{}\t{}".format(entry, "\t".join(current_path)))
        current_path.pop()


if __name__ == '__main__':
    input_json = sys.argv[1]
    with open(input_json) as infile:
        root = json.load(infile)
    get_all_children(root)
