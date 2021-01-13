import os
import re
import sys
db_path = sys.argv[1]
base_name = re.sub(r'\.[^\.]+$', '', os.path.basename(db_path))
out_path = './' + base_name + '_cleaned.fasta'

with open(db_path, 'r') as db, open(out_path, 'w') as odb:
    line_drop = []
    for num, line in enumerate(db):
        if re.search(r'^>', line):
            header = re.search('^ *[^ ]+', line).group()
            header = re.sub(r"[^\|]+$", '', header)
            tail = re.search(r'\[.+\]', line).group()
            body = line.replace(tail, '').replace(header, '')
            header = header.strip()
            tail = tail.strip()
            body = body.strip()
            if body.find(' ') == -1:
                odb.write(line)
            else:
                gene1 = re.search(r'[A-Za-z\-0-9]+[A-Z][A-Za-z\-0-9]*', body)
                gene2 = re.search('[A-Za-z]+[0-9]+[A-Za-z]*', body)

                gene_hint = ""
                if gene1 or gene2:
                    gene_hint = gene1.group() if gene1 else gene2.group()
                elif body.count(' ') == 1:
                    gene_hint = re.sub(' ', '_', body)

                print("Invlid gene name detected, the current gene definition is: '{}'; the bacteria name is {}".format(body, tail))
                ipt = input(
                    "Input 'y' or '' to use '{}' (drop sequence if gene name ='') as gene name, input 'n' to drop this sequence from database, or input a gene name defined by yourself\n: ".format(gene_hint))

                ipt == ipt.strip()
                gene = ''
                if any((ipt == 'y', ipt == '', ipt == '\n')):
                    gene = gene_hint
                elif not ipt == 'n':
                    gene = ipt

                gene = gene.strip()
                if gene:
                    odb.write(header + gene + ' ' + tail + '\n')
                else:
                    line_drop.append(num + 1)
        else:
            if num not in line_drop:
                odb.write(line)
