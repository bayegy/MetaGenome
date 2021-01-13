#!/usr/bin/env python3
import os, re, sys
import pandas as pd
ms, infile, outfile = sys.argv
#infile = "summary_rawdata.txt"
#outfile = "summary_rawdata_html.txt"
with open(infile, 'r') as f:
    f = pd.read_csv(f, index_col=0, header=0, sep="\t")
    f = f.sort_index(axis=0, ascending=True)
    f.to_csv("temp.txt", sep="\t", index=True)
    with open(outfile, 'w') as o:
        with open("./temp.txt") as temp:
            num = 1
            for line in temp:
                line = line.strip()
                cell = re.split(r'\t', line)
                if num == 1:
                    o.write("<tr><th>{}</th><th>{}</th><th>{}</th><th>{}</th><th>{}</th><th>{}</th><th>{}</th><th>{}</th></tr>\n".format("SampleID", "InsertSize(bp)", "SeqStrategy", "Reads(#)", "Base(GB)", "GC(%)", "MaxLength", "MinLength"))
                    num += 1
                else:
                    o.write("\t\t\t\t\t<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>\n".format(cell[0], cell[1], cell[2], cell[3], cell[4], cell[5], cell[6], cell[7]))
os.remove("temp.txt")
