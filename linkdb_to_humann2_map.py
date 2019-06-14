import sys
import os
import re
import gzip


def main(link_data, field1, out_file, field2='up', adjust_field1=False, adjust_field2=lambda x: "UniRef90_" + x):
    link_map = {}

    def add(fid, ldef):
        try:
            link_map[fid].append(ldef)
        except KeyError:
            link_map[fid] = [ldef]

    is_gz = os.path.splitext(link_data)[1] == '.gz'
    with (gzip.open(link_data, 'r') if is_gz else open(link_data, 'r')) as f:
        for line in f:
            line = line.decode().strip() if is_gz else line.strip()
            f1 = re.search("{}:([^\t ]+)".format(field1), line).group(1)
            f1 = adjust_field1(f1) if adjust_field1 else f1
            f2 = re.search("{}:([^\t ]+)".format(field2), line).group(1)
            f2 = adjust_field2(f2) if adjust_field2 else f2
            add(f1, f2)

    with open(out_file, 'w') as out:
        for f1, f2s in link_map.items():
            out.write("{}\t{}\n".format(f1, '\t'.join(f2s)))


if __name__ == '__main__':
    script_path, link_data, field1, out_file = sys.argv
    main(link_data, field1, out_file)
