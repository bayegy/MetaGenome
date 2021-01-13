#!/usr/bin/env python3
import os, re, sys
ms, infold, outfile = sys.argv
#infold = "rawdata"
#outfile = "rawdata_summary.txt"

def count_fastq(fastq_1, fastq_2):
    with open(fastq_1) as f1, open(fastq_2) as f2:
        result = []
        lengths = []
        read_num = 0
        base_num = 0
        gc_num = 0
        
        line_num = 0
        for line in f1:
            line_num += 1
            if line_num % 4 == 2:
                read_num += 1
                base_num += len(line)
                gc_num += line.count("G") + line.count("g") + line.count("C") + line.count("c")
                lengths.append(len(line))
                
        line_num = 0
        for line in f2:
            line_num += 1
            if line_num % 4 == 2:
                read_num += 1
                base_num += len(line)
                gc_num += line.count("G") + line.count("g") + line.count("C") + line.count("c")
                lengths.append(len(line))
        gc_percent = format(gc_num/base_num*100, '0.2f')
        result.append(read_num)
        result.append(base_num)
        result.append(gc_percent)
        result.append(max(lengths))
        result.append(min(lengths))
        return(result)
        
file_name = []
for each in os.listdir(infold):
    each = re.findall(r'(^.*)_', each)
    file_name.append(''.join(each))

postfix_1 = re.findall(r'_(.*)', sorted(os.listdir(infold))[0])
postfix_2 = re.findall(r'_(.*)', sorted(os.listdir(infold))[1])

with open(outfile, 'w') as o:
    o.write("\tInsertSize(bp)\tSeqStrategy\tReads(#)\tBase(GB)\tGC(%)\tMaxLength\tMinLength\n")
    for each in list(set(file_name)):
        fastq_1 = "{}/{}_{}".format(infold, each, ''.join(postfix_1))
        fastq_2 = "{}/{}_{}".format(infold, each, ''.join(postfix_2))
        reads_info = count_fastq(fastq_1, fastq_2)
        read_num = reads_info[0]
        base_num = format(reads_info[1]/10**9, '0.2f')
        gc_percent = reads_info[2]
        max_length = reads_info[3]
        min_length = reads_info[4]
        o.write("{}\t350\t(150:150)\t{}\t{}\t{}\t{}\t{}\n".format(each, read_num, base_num, gc_percent, max_length, min_length))

