#!/usr/bin/env python3
import pandas as pd
import sys

sp, lda2, lda4, thresh = sys.argv


with open(lda2) as lda2, open(lda4, 'w') as lda4:
    for line in lda2:
        feature, raw_lda, sig_group, adj_lda, p_value = line.strip().split('\t')
        if adj_lda == "" or float(adj_lda) >= float(thresh):
            lda4.write(line)
        else:
            lda4.write('{}\t{}\t\t\t{}\n'.format(feature, raw_lda, p_value))
