#!/usr/bin/env python3
import sys, re, os
ms, temple, table_id, infile, outfile= sys.argv
with open(temple) as f_html:
    with open(infile) as f_table:
        f_html = f_html.read()
        f_table = f_table.read()
        judge = f_html.find(table_id)
        if judge != -1:
            result = re.sub(table_id, f_table, f_html)
            with open(outfile, 'w') as o:
                o.write(result)
        else:
            print("\033[31m error: there is no {} !\033[0m".format(table_id))

if judge != -1:
    os.remove(temple)
    os.rename(outfile, temple)
