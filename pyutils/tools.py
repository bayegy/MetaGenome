import numpy as np
import re
import os
import pandas as pd
import time


def dupply(array, func, **kwargs):
    """This func will conside only the first and second dimension of array"""
    shape = array.shape
    return np.array([[func(array[i, j], **kwargs) for j in range(shape[1])] for i in range(shape[0])])


def time_counter(func):
    def wrapper(*args, **kwargs):
        t1 = time.time()
        out = func(*args, **kwargs)
        print("{} done, time used: {}".format(str(func), str(time.time() - t1)))
        return out
    return wrapper


def generate_span(number_list: [], start=0) -> []:
    """Please sapply iterable number list"""
    step_sum = []
    current_sum = start
    for e in number_list:
        current_sum += e
        step_sum.append(current_sum)
    flat = list(np.array(step_sum) - np.array(number_list))
    return [[m, n] for m, n in zip(flat, step_sum)]


def split_list(init_list, each=5):
    ll = len(init_list)
    rem = ll % each
    rem = [rem] if rem else []
    num_list = ([each] * int(np.floor(ll / each))) + rem
    index = generate_span(num_list)
    return [init_list[l[0]:l[1]] for l in index]


def parse_premap(raw_fqs_dir, pre_mapping_file, forward_regex, reverse_regex, sample_regex) ->dict:
    # get the mapping relation of id to description
    print("Assert no duplicated sample names: please make sure no duplicated sample names and did not use capital and small letter to distinguish samples if error happened\n")
    id_map = {}
    with open(pre_mapping_file, 'r') as mapping:
        for line_number, line in enumerate(mapping):
            li = re.split('\t', line.strip())
            if line_number > 0:
                first_element = li[0].strip().lower()
                assert first_element not in id_map.keys()
                id_map[first_element] = li[len(li) - 1].strip()

    # find the corresponding fastq file
    fq_succeeded = 0
    fq_failed = 0
    matched_id_set = []
    matched_fq_info = []
    print("\nThe following samples were found in fastq files, but not found in pre_mapping file:")
    for root, dirs, files in os.walk(raw_fqs_dir):
        if not len(files) == 0:
            root = os.path.abspath(root)
            for fl in files:
                if re.search(sample_regex, fl):
                    fq_path = "%s/%s" % (root, fl)
                    fq_id = re.search(sample_regex, fl).group(1)
                    pre_id = fq_id.lower()
                    try:
                        new_id = id_map[pre_id]
                        if re.search(reverse_regex, fl):
                            matched_fq_info.append([fq_path, new_id, "R2"])
                            fq_succeeded += 1
                        elif re.search(forward_regex, fl):
                            matched_fq_info.append([fq_path, new_id, "R1"])
                            matched_id_set.append(pre_id)
                            fq_succeeded += 1
                    except KeyError:
                        print('    ' + fq_id)
                        fq_failed += 1
    matched_fq_info = pd.DataFrame(matched_fq_info, columns=["Fastq_path", "New_SampleID", "Direction"])
    matched_fq_info = matched_fq_info.sort_values(by="New_SampleID")

    print("\nThe following samples were found in pre_mapping file, but not found in fastq files:")
    matched_map = []
    with open(pre_mapping_file, 'r') as mapping:
        map_failed_id = 0
        for line_number, line in enumerate(mapping):
            li = re.split('\t', line.strip())
            li = [l.strip() for l in li]
            fc = len(li) - 1
            if line_number == 0:
                li[0] = "#SampleID"
                li[fc] = "Description"
                columns = li
            else:
                if li[0].lower() in matched_id_set:
                    li[0] = li[fc]
                    matched_map.append(li)
                else:
                    print('    ' + li[0])
                    map_failed_id += 1

    matched_map = pd.DataFrame(matched_map, columns=columns)
    matched_map = matched_map.sort_values(by="#SampleID")
    print('''
    In summary:
        %s fastq files were matched;
        %s fastq files were filterd as the ids of samples were not found in pre_mapping file;
        %s samples were filterd from pre_mapping file, as the ids of samples were not found among forward fastq files\' names'.
        ''' % (fq_succeeded, fq_failed, map_failed_id))
    return {"fastq": matched_fq_info, "map": matched_map}
