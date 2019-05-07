import numpy as np
import pandas as pd


def read_abundance(file_path, index_col=0, return_sum=-1, return_mean=-1):
    """
    Abundance table should be tab seprated.
    args:
        file_path: file path of abundance table
        index_col: which column is the indexes (row names), e.g: -1 for last column, 0 for the first column
        return_sum: -1(return df.DataFrame of abundance table), or 1(sums of each row), or 0(sums of each column)
        return_mean: -1(return df.DataFrame of abundance table), or 1(means of each row), or 0(means of each column),
                    if not return_sum==-1, this argument will be disabled.
    """
    df = pd.read_csv(file_path, sep="\t", index_col=index_col)
    df = df.loc[:, df.dtypes != np.object]
    return df.sum(axis=return_sum) if not return_sum == -1 else (df if return_mean == -1 else df.mean(axis=return_mean))
