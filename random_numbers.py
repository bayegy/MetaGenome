#!/usr/bin/env python3.8
import sys
import pandas as pd
import numpy as np
# import random
sp, input_file, output_file = sys.argv

df = pd.read_csv(input_file, sep="\t")

nrow, ncol = df.shape

for index, dtype in enumerate(df.dtypes):
    if dtype != np.dtype("O"):
        df.iloc[:, index] = df.iloc[:, index] + np.random.normal(10000, 1000, nrow)

df.to_csv(output_file, sep="\t", index=False)
