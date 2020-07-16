import pickle
import os
import shutil
import uuid
from concurrent.futures import ProcessPoolExecutor


class GroupGene(object):
    """docstring for GroupGene"""

    def __init__(self, annotation_file, by=1, adjust_func=False, tmp_dir="/dev/shm"):
        """
        aguments:
            annotation_file: annotation file from diamond, bowtie2, emmaper(without header), and so no
            by: column contains annotaions in annotation file
            adjust_func: function used to transform annotation, returning a list of adjusted annotation
        """
        self.tmp_dir = os.path.join(tmp_dir, "groupGene_tmpdir_" + str(uuid.uuid1()))
        print("tmp_dir at {}".format(self.tmp_dir))
        self.__load_annotation(annotation_file, by, adjust_func)

    def __load_annotation(self, annotation_file, by, adjust_func):
        print("start to load annotation...")
        ft_map = {}
        with open(annotation_file) as f:
            for line in f:
                li = line.strip().split('\t')
                annotation = li[by].strip()
                if annotation:
                    if adjust_func:
                        annotation = adjust_func(annotation)
                    if not isinstance(annotation, list):
                        annotation = [annotation]
                    for anno in annotation:
                        # use set to avoid duplicates
                        ft_map.setdefault(anno, set()).add(li[0])
        # list eat less memery than set
        # for k, v in ft_map.items():
        #    ft_map[k] = list(v)
        # self.annotation = ft_map
        print("annotation loaded!")
        print("dumping annotation to tmp_dir...")
        if not os.path.exists(self.tmp_dir):
            os.makedirs(self.tmp_dir)
        self.dump_file = os.path.join(self.tmp_dir, str(uuid.uuid1()))
        with open(self.dump_file, 'wb') as f:
            pickle.dump(obj=ft_map, file=f)
        print("annotation dumped!")

    def group(self, gene_abundance, out_file, column=3, header=True):
        """
        aguments:
            gene_abundance: salmon quant.sf path
            out_file: output file path
            column: column in gene_abundance for group abundance calculation
            header: header in first line of gene_abundance file (True) or not (False)
        """
        print("Regrouping {} to {}...".format(gene_abundance, out_file))

        gene_abc_dict = self.load_abc(gene_abundance, column, header)

        print("reload annotation from dumped file...")
        with open(self.dump_file, 'rb') as f:
            annotation = pickle.load(file=f)
        print("annotation reloaded, start to group...")
        with open(out_file, "w") as fout:
            fout.write("# Ortholog Group\tAbundance\n")
            for anno, genes in annotation.items():
                sum_abc = sum([gene_abc_dict.get(gene, 0) for gene in genes])
                if sum_abc > 0:
                    fout.write("{}\t{}\n".format(anno, sum_abc))
        print("Regroup {} to {} done!".format(gene_abundance, out_file))

    def load_abc(self, gene_abundance, column=3, header=True):
        abc_tmp = self.get_abc_tmp(gene_abundance, column)
        if os.path.exists(abc_tmp):
            print("loading dumped abundance dict from {}".format(abc_tmp))
            with open(abc_tmp, 'rb') as f:
                return pickle.load(file=f)
        print("loading abundance dict from {}".format(gene_abundance))
        abc_tmp_dir = os.path.dirname(abc_tmp)
        if not os.path.exists(abc_tmp_dir):
            os.makedirs(abc_tmp_dir)
        gene_abc_dict = {}
        with open(gene_abundance, "r") as abc:
            for number, line in enumerate(abc):
                if number == 0 and header:
                    continue
                li = line.strip().split("\t")
                gene_abc_dict[li[0]] = float(li[column])
        with open(abc_tmp, 'wb') as f:
            pickle.dump(obj=gene_abc_dict, file=f)
        return gene_abc_dict

    def get_abc_tmp(self, gene_abundance, column=3):
        return "{}/groupGene_abctmp/{}_column_{}".format(
            os.path.dirname(gene_abundance) or ".",
            os.path.basename(gene_abundance),
            column
        )

    def map(self, in_out_dict, column=3, header=True, processors=3):
        try:
            executor = ProcessPoolExecutor(max_workers=processors)
            for gene_abundance, out_file in in_out_dict.items():
                executor.submit(self.group, gene_abundance, out_file, column, header)
            executor.shutdown(True)
        finally:
            print("removing the tmp_dir, a moment please")
            self.remove_tmpdir()

    def remove_tmpdir(self):
        if os.path.exists(self.tmp_dir):
            shutil.rmtree(self.tmp_dir)
