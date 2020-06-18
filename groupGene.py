
class GroupGene(object):
    """docstring for GroupGene"""

    def __init__(self, annotation_file, by=1, adjust_func=False):
        """
        aguments:
            annotation_file: annotation file from diamond, bowtie2, emmaper(without header), and so no
            by: column number
            adjust_func: function used to transform annotation, returning a list of adjusted annotation
        """
        self.__load_annotation(annotation_file, by, adjust_func)

    def __load_annotation(self, annotation_file, by, adjust_func):
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
        for k, v in ft_map.items():
            ft_map[k] = list(v)
        self.annotation = ft_map

    def group(self, gene_abundance, out_file, column=3, header=True):
        """
        aguments:
            gene_abundance: salmon quant.sf path
            out_file: output file path
            column: column in gene_abundance for group abundance calculation
            header: header in first line (True) or not (False)
        """
        print("Regroup {} to {}".format(gene_abundance, out_file))
        gene_abc_dict = {}
        with open(gene_abundance, "r") as abc:
            for number, line in enumerate(abc):
                if number == 0 and header:
                    continue
                li = line.strip().split("\t")
                gene_abc_dict[li[0]] = float(li[column])
        with open(out_file, "w") as fout:
            fout.write("# Ortholog Group\tAbundance\n")
            for anno, genes in self.annotation.items():
                sum_abc = sum([gene_abc_dict.get(gene, 0) for gene in genes])
                if sum_abc > 0:
                    fout.write("{}\t{}\n".format(anno, sum_abc))
