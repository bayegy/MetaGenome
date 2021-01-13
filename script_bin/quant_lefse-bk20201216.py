#!/usr/bin/env python3
import re, sys, os, shutil
import pandas as pd
from task.libs.groupCompare import lefse
ms, metadata, all_group, color_list = sys.argv

with open(metadata) as meta_file:
    with open("Bin_all/Bin_quant/bin_abundance_table.txt") as data_file:
        meta = pd.read_csv(meta_file, sep = "\t", header = 0, keep_default_na = False)
        data = pd.read_csv(data_file, sep = "\t", header = 0, keep_default_na = False)
        data.rename(columns = {'Genomic bins':'SampleID'}, inplace = True)
        data_new = data.set_index(data.SampleID).iloc[:, 1:data.shape[1]]
        data_new = data_new.sort_index(axis=1)
        med = meta.set_index(meta.SampleID).iloc[:, 1:(meta.shape[1]-1)]
        meta_new = pd.DataFrame(med.values.T, index=med.columns, columns=med.index)
        new_list = []
        new_list.append(meta_new)
        new_list.append(data_new)
        new_df = pd.concat(new_list)
        
cate_num = sum(map(lambda cell : 'Category' in cell, new_df.index))

if os.path.isdir("Bin_all/Bin_quant/lefse"):
    shutil.rmtree("Bin_all/Bin_quant/lefse", ignore_errors=True)
    print("\033[32mlefse already existed! make new lefse!\033[0m")
    os.makedirs("Bin_all/Bin_quant/lefse")
else:
    print("\033[32mNo lefse! make lefse!\033[0m")
    os.makedirs("Bin_all/Bin_quant/lefse")

new_df.to_csv('Bin_all/Bin_quant/lefse/bin_group_abundance.txt', sep='\t', index = True)

def colorDict(all_group, color_list):
    cate_group = list(set(new_df.loc[all_group].tolist()))
    group_num = len(cate_group)
    with open(color_list) as color:
        color = color.readlines()
        color_new = []
        for each in color[:group_num]:
            each = re.sub('\n', '', each)
            color_new.append(each)
    return(dict(zip(cate_group, color_new)))

color_dict = colorDict(all_group, color_list)

for i in range(cate_num):
    n = i + 1
    index_name = "Category{}".format(n)
    fold = "Bin_all/Bin_quant/lefse/{}".format(index_name)
    os.makedirs(fold)
    lefse_input = new_df.ix[[i] + list(range(cate_num, len(new_df))), (new_df.loc[index_name] != '')]
    lefse_input.insert(0, 'SampleID', list(lefse_input.index))
    lefse_input_name = "{}/{}_abundance.txt".format(fold, index_name)
    lefse_input.to_csv(lefse_input_name, sep='\t', index = False)
    lefse(lefse_input_name, colors=color_dict, lda_th=2, format="pdf", out_dir=fold)
    os.remove("{}/lefse.cladogram.pdf".format(fold))
    os.remove("{}/lefse.lefseinput.txt".format(fold))
    os.system("convert -density 400 -quality 200 {}/lefse.pdf {}/lefse.png".format(fold, fold))
    


