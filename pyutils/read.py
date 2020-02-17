import numpy as np
import pandas as pd
import re
import pdb
from bs4 import BeautifulSoup
# import matplotlib.pyplot as plt


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


def read_to_html_table(file_path, table_class=False, thead_class=False):
    file_type = re.search('[^\.]+$', file_path).group()
    if file_type == 'txt' or file_type == 'csv':
        out_df = pd.read_csv(file_path, sep='\t')
    elif file_type == "html":
        out_df = pd.read_html(file_path)[0]
    # pdb.set_trace()
    table = out_df.to_html(index=False)
    if table_class or thead_class:
        table = BeautifulSoup(table, 'lxml')
        if table_class:
            table.table['class'] = table_class
        if thead_class:
            table.thead['class'] = thead_class
        table = table.prettify()
    return table


def get_kingdom_ratio(species_abundance_table):
    df = read_abundance(species_abundance_table, index_col=-1, return_sum=1)
    kingdom = [re.search('^[^;]+', s).group().replace('k__', '')
               for s in df.index]
    # pdb.set_trace()
    df = df.groupby(kingdom).sum()
    try:
        del df['Environmentalsamples']
    except Exception:
        pass

    def series_to_str(series, value_format=':.2%'):
        return ', '.join([('{}({%s})' % (value_format)).format(k, v) for k, v in series.items()])

    # max_abundance = max(df)
    # explode = [0.1 if i == max_abundance else 0 for i in df]
    # plt.pie(df, explode=explode, labels=df.index, labeldistance=1.1,
    #         autopct='%2.2f%%', shadow=False, startangle=90, pctdistance=0.6)
    # plt.show()
    return "您样本中，检测到的物种有：{}; 物种占比情况：{}。".format(series_to_str(df, value_format=':.0f'), series_to_str(df / df.sum()))


def format_file(in_fp, out_fp, **kwargs):
    with open(in_fp, 'r') as f:
        out = f.read()
        for k, v in kwargs.items():
            out = out.replace("{{%s}}" % (k), v)
    with open(out_fp, 'w') as f:
        f.write(out)


def update_html_properties(html_fp, format_dict, out_fp, filter_function=False, use_selector=False):
    """
    in a format_dict, you should never use '' or "" in a CSS selector
    """
    with open(html_fp, "r", encoding="utf-8") as file:
        fcontent = file.read()
    sp = BeautifulSoup(fcontent, 'lxml')
    for k, v in format_dict.items():
        # pdb.set_trace()
        if use_selector:
            label_list = sp.select(k)
        else:
            tag = re.search('^[^\[]+', k).group()
            attrs = re.search('\[(.+)\]', k)
            attrs = dict([attr.split('=')
                          for attr in attrs.group(1).split('&')]) if attrs else {}
            label_list = sp.findAll(tag, attrs=attrs)
        for label in label_list:
            for p, n in v.items():
                if label.has_attr(p) and (filter_function(label[p]) if filter_function else True):
                    label[p] = n.format(value=label[p])
    with open(out_fp, "w", encoding="utf-8") as out:
        out.write(sp.prettify())
