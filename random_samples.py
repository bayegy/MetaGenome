#!/usr/bin/env python3
import random
import sys
import numpy as np
import pdb

abundance_table = sys.argv[1]
group_size = int(sys.argv[2])
out_file = sys.argv[3]

sep = ','


def get_u_sigma(scale):
    scale = int(np.ceil(scale * 0.2 * 10000))
    st = int(np.ceil(scale * 0.01))
    floor = int(np.ceil(scale * 0.1))
    # pdb.set_trace()
    return (random.randrange(-scale * 2, scale * 2, st) / 10000, random.randrange(floor, floor * 3, st) / 10000)


def random_samples(file):
    with open(file) as infile, open(out_file, 'w') as outfile:
        for num, line in enumerate(infile):
            li = line.strip().split(sep)
            line_header, line_body = [li[0]], li[1:]
            sample_number = len(line_body)
            wv = 0.2 if num % 2 == 1 else -0.2
            if num == 0:
                group_number = np.ceil(sample_number / group_size)
                outfile.write(line)
            else:
                seed = np.array(line_body[:group_size], dtype="float64")
                scale = seed.mean()
                new_line = seed.tolist()
                for i in range(int(group_number - 1)):
                    # pdb.set_trace()
                    new_line += (seed + np.random.normal(wv * scale, 0.06 * scale, group_size)).tolist()
                new_line = new_line[:sample_number]
                new_line = [str(i) for i in new_line]
                outfile.write(sep.join(line_header + new_line) + '\n')


random_samples(abundance_table)
