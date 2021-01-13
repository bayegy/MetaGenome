# from ..metagenomePipline import MetagenomePipline
from ..pyutils.controller import *
from ..pipconfig import settings
from .envShell import metabat2, prokka
import os
import re
import glob
from Bayegy.ampliconLibs.systemMixin import SystemMixin
from task.libs.__utils.script import script
import pandas as pd
import numpy as np
import shutil
from MetaGenome.pyutils.read import iter_fa
from MetaGenome.visualizeFunction import VisualizeFunction

class SummBinAnnos(SystemMixin):

    def __init__(self, inputs, out_dir, column=1, regroup=None,
                 adjust_func=None, skip=0, basename="gene",
                 id_pattern=r"bin\..+\.\d+", levels=None,
                 outputs=["summary"]):
        """
         regroup: path to file contain the gene map information which will be used to regroup genes
         levels: path to file contain the gene levels information which will be used to annotate genes
                 or a function take gene id as input and return a level tag
        """
        kwargs = locals()
        kwargs = {k: v for k, v in kwargs.items() if k not in ["self"]}
        self.set_attr(**kwargs)
        self.annotations = glob.glob(inputs)
        self.set_path(force=True,
                      out_dir=out_dir
                      )
        if regroup:
            self.regroup = self.load_map(regroup)
        if levels and isinstance(levels, str):
            levels = self.load_map(levels)
            self.levels = {k: '; '.join(v) for k, v in levels.items()}
        self.levels_paths = []

    def count_genes(self):
        regroup = self.regroup
        levels = self.levels
        bins_gene_count = {}
        index = set()
        for anno in self.annotations:
            bin_id = re.search(self.id_pattern, os.path.basename(anno))
            if not bin_id:
                continue
            bin_id = bin_id.group()
            orf_genes = {}
            gene_count = {}
            with open(anno) as fh:
                for num, line in enumerate(fh):
                    if num < self.skip or line.startswith("#"):
                        continue
                    li = line.strip().split('\t')
                    if len(li) <= self.column:
                        continue
                    genes = li[self.column]
                    if not genes:
                        continue
                    if self.adjust_func:
                        genes = self.adjust_func(genes)
                    if isinstance(genes, str):
                        genes = [genes]
                    if regroup:
                        genes = [regroup[gene] for gene in genes if regroup.get(gene)]
                        genes = self.flat_list(genes)
                    orf_id = li[0]
                    for gene in genes:
                        orf_genes.setdefault(orf_id, set()).add(gene)
            for value in orf_genes.values():
                for gene in value:
                    index.add(gene)
                    gene_count[gene] = gene_count.get(gene, 0) + 1
            bins_gene_count[bin_id] = gene_count
        index = list(index)
        bins_gene_count = {
            k: [v.get(g, 0) for g in index] for k, v in bins_gene_count.items()
        }
        df = pd.DataFrame(bins_gene_count)
        df.index = index
        df.index.set_names("#Gene Family", inplace=True)
        if levels:
            if isinstance(levels, dict):
                df['levels'] = [levels.get(k, "") for k in df.index]
            elif callable(levels):
                df['levels'] = [levels(k) for k in df.index]
        copy_path = os.path.join(self.out_dir, self.basename + 'GeneCount.xls')
        df.to_csv(copy_path, sep='\t')
        self.set_path(force=False, copy_path=copy_path)
        self.copy_df = df

    def load_map(self, path, skip=0):
        map_dict = {}
        with open(path) as fh:
            for num, line in enumerate(fh):
                if num < skip:
                    continue
                li = line.strip().split('\t')
                key, *vals = li
                map_dict.setdefault(key, []).extend(vals)
        return map_dict

    def flat_list(self, array: list):
        fl = set()
        for e in array:
            if not isinstance(e, list):
                fl.add(e)
                continue
            for se in e:
                fl.add(se)
        return list(fl)

    def summarize_levels(self):
        out_pattern = os.path.join(self.out_dir, self.basename + "GeneCount.L{}.xls")
        glob_pattern = out_pattern.replace('{}', '*')
        script.summarize_levels.py(
            input=self.copy_path,
            output_abs=out_pattern
        )
        self.levels_paths = sorted(glob.glob(glob_pattern))

    def plot_summary(self):
        out_dir = os.path.join(self.out_dir, self.basename + "summary")
        script.plot_table.R(
            input=self.copy_path,
            top=20,
            boxplot=out_dir,
            heatmap=out_dir,
            prefix=self.basename
        )
        script.plot_table.R(
            input=self.copy_path,
            top=20,
            stackbar=out_dir,
            other=True,
            prefix=self.basename
        )

    def plot_levels_summary(self):
        out_dir = os.path.join(self.out_dir, self.basename + "summary")
        for num, path in enumerate(self.levels_paths):
            level = num + 1
            prefix = self.basename + "L{}_".format(level)
            script.plot_table.R(
                input=path,
                top=20,
                boxplot=out_dir,
                heatmap=out_dir,
                prefix=prefix,
                annotate=level != 1
            )
            script.plot_table.R(
                input=self.copy_path,
                top=20,
                stackbar=out_dir,
                other=True,
                prefix=prefix
            )

    def detail(self, input, out_dir, annotate=True, transpose=False, percent=False):
        df = pd.read_csv(input, sep="\t", index_col=0)
        is_num = [t != np.dtype("O") for t in df.dtypes]
        df = df.loc[:, is_num]
        if transpose:
            df = df.T
        for column in df.columns:
            column_out = os.path.join(out_dir, column)
            script.plot_table.R(
                input=input,
                transpose=transpose,
                colname=column,
                top=50,
                bar=column_out,
                annotate=annotate
            )
            script.plot_table.R(
                input=input,
                transpose=transpose,
                colname=column,
                top=10,
                pie=column_out,
                other=True,
                percent=percent
            )

    def plot_bin_detail(self):
        self.detail(
            input=self.copy_path,
            out_dir=os.path.join(self.out_dir, self.basename + "detail_bin"),
            percent=True
        )

    def plot_gene_detail(self):
        self.detail(
            input=self.copy_path,
            out_dir=os.path.join(self.out_dir, self.basename + "detail_gene"),
            transpose=True
        )

    def plot_bin_detail_levels(self):
        for num, path in enumerate(self.levels_paths):
            level = num + 1
            self.detail(
                input=path,
                out_dir=os.path.join(
                    self.out_dir,
                    self.basename + "detail_bin_L{}".format(level)
                ),
                annotate=level != 1,
                percent=True
            )

    def plot_gene_detail_levels(self):
        for num, path in enumerate(self.levels_paths):
            level = num + 1
            self.detail(
                input=path,
                out_dir=os.path.join(
                    self.out_dir,
                    self.basename + "detail_gene_L{}".format(level)
                ),
                transpose=True
            )

    def pipeline(self):
        self.count_genes()
        if self.levels:
            self.summarize_levels()
        for output in self.outputs:
            handler = getattr(self, "plot_" + output)
            handler()
        if "bin_detail" in self.outputs:
            bin_detail = os.path.join(self.out_dir, self.basename + "detail_bin")
            for anno in self.annotations:
                bin_id = re.search(self.id_pattern, os.path.basename(anno))
                if not bin_id:
                    continue
                bin_id = bin_id.group()
                shutil.copyfile(
                    anno,
                    os.path.join(bin_detail, bin_id, bin_id + "_annotations.txt")
                )


synchronize = sge_decorator if settings.use_sge else pool_decorator


class VisualizeBinMixin(object):
    """docstring for VisualizeBin"""

    def visualize_bin_prepare(self):
        self.set_path(
            force=True,
            vb_result="{bin_result}/Results",
            vb_qc="{vb_result}/00-QCStats",
            vb_rawqc="{vb_qc}/1-QC_report_Rawfastq",
            vb_cleanqc="{vb_qc}/2-QC_report_Filtered",
            vb_bin="{vb_result}/01-Bin",
            vb_bin_all="{vb_bin}/1-Bin_raw",
            vb_bin_pick="{vb_bin}/2-Bin_pick",
            vb_plot="{vb_result}/02-Bin_Plot",
            vb_abc="{vb_result}/03-Bin_Abundance",
            vb_function="{vb_result}/04-Bin_Function",
            vb_tax="{vb_result}/05-Bin_Taxonomy",
            vb_circos="{vb_result}/06-Bin_Circos",
            vb_kegg="{vb_function}/KEGG",
            vb_go="{vb_function}/GO",
            vb_cazy="{vb_function}/CAZyme",
            vb_cog="{vb_function}/COG",
            vb_bin_sum="{bin_result}/Bin_summary",
            vb_bin_plot="{bin_result}/Bin_plot"
        )
        self.set_attr(
            bin_abc_path=os.path.join(self.vb_abc, "bin_abundance_table.xls"),
            bin_abc_tax_path=os.path.join(self.vb_tax, "bin_abundance_taxonomy.xls")
        )

    def gather_data_process(self):
        self.system("""
cp {out_dir}/fastqc_out/*/*html {vb_rawqc}
cp {out_dir}/kneaddata_out/*/fastqc/*kneaddata*html {vb_cleanqc}
cp {out_dir}Report/reads_summary.txt {vb_qc}/
cp -L {bin_all}/* {vb_bin_all}/
cp -L  {bin_pick}/* {vb_bin_pick}/

n90_file={vb_bin_sum}/N50.N90.txt
echo -e "\\033[32m统计N50和N90: \\033[0m"
echo -e "BinID\\tN50\\tN90" > $n90_file
for i in {bin_pick}/*.fa; do
    base=${{i##*/}}
    name=${{base%.fa}}
    echo -ne "$name\\t" >> $n90_file
    {perl_path} {metagenome_home}/script_bin/N50.N90.pl $i | \
    sed 's/N50: //g' | \
    sed 's/N90: //g' | \
    tr -s '\\n' '\\t' | \
    awk '{{printf"%s\\t%s\\n",$1,$2}}' >> $n90_file
    echo -e "\\033[32m$i Done...\\033[0m"
done

cat {bin_drep}/data/checkM/checkM_outdir/results.tsv | \
 sed 's/\\.fa\\t/\\t/g' | awk -F '\\t' 'BEGIN{{OFS="\\t"}}{{print $1,$12,$13}}' \
 > {vb_bin}/checkm_raw.txt

{python3_path} {metagenome_home}/stat_bin_gc.py {bin_pick}/ \
 {vb_bin_sum}/bin.gc.txt

{python3_path} {bayegy_home}/merge_tables.py  \
 {vb_bin_sum}/N50.N90.txt \
 {vb_bin_sum}/bin.gc.txt \
 {vb_bin}/checkm_raw.txt \
 - | awk -F "\\t" 'BEGIN{{OFS="\\t"}}{{print $1,$6,$7,$2,$3,$4,$5}}' \
 > {vb_bin}/bin_picked_summary.txt""")

    def bins_contig_gc(self, inputs, out_path, bin_id_pattern=r"bin\..+\.\d+"):
        array = []
        for bp in glob.glob(inputs):
            bin_id = re.search(bin_id_pattern, os.path.basename(bp)).group()
            for header, seq in iter_fa(bp, trim_line_break=True):
                contig_id = header.split()[0].lstrip(">")
                seq.upper()
                gc = (seq.count('G') + seq.count('C')) / len(seq)
                array.append([contig_id, bin_id, gc])
        df = pd.DataFrame(array, columns=["ContigID", "BinID", "GC"])
        df.to_csv(out_path, sep="\t", index=False)

    def plot_bin(self):
        self.bins_contig_gc(
            inputs=os.path.join(self.bin_pick, "*.fa"),
            out_path=os.path.join(self.vb_bin_plot, "bins_contig_gc.txt")
        )
        self.system("""
cat {bin_out}/*/depth.txt | sed -e '/totalAvgDepth/d' | \
 awk -F"\\t" 'BEGIN{{OFS="\\t"}}{{print $1, $3}}'  | \
 sed '1i ContigID\\ttotalAvgDepth' > {vb_bin_plot}/bins_contig_depth.txt
{python3_path} {bayegy_home}/merge_tables.py \
 {vb_bin_plot}/bins_contig_gc.txt \
 {vb_bin_plot}/bins_contig_depth.txt \
 {vb_plot}/bins_contig_summary.txt""")
        script.scatter.R(
            input=os.path.join(self.vb_plot, "bins_contig_summary.txt"),
            axisX="GC",
            axisY="totalAvgDepth",
            category="BinID",
            outpath=os.path.join(self.vb_plot, "bins_contig_summary.pdf"),
            labelX="GC Percent",
            labelY="Average Depth"
        )

    def analyze_bin_abc(self):
        self.system("cp {bin_quant}/bin_abundance_table.tab {bin_abc_path}")
        VisualizeFunction(
            self.bin_abc_path,
            self.mapping_file,
            self.categories,
            prefix="bin_",
            out_dir=self.vb_abc
        ).visualize(
            exclude=self.exclude,
            orders=self.orders
        )

    def create_bin_consensus(self):
        abc_df = pd.read_csv(self.bin_abc_path, sep="\t", index_col=0)
        tax_dict = {}
        with open(os.path.join(self.bin_phylo, "Bin_phylo.tsv")) as fh:
            for num, line in enumerate(fh):
                if num == 0:
                    continue
                li = line.strip().split("\t")
                tax_dict[li[0]] = li[1].split(":")[2].replace("|", ";")
        abc_df["Consensus Lineage"] = [tax_dict.get(b, "") for b in abc_df.index]
        abc_df.index.set_names("#OTU ID", inplace=True)
        abc_df.to_csv(self.bin_abc_tax_path, sep="\t")

    def plot_bin_phylo(self):
        self.create_bin_consensus()
        self.system("cp {bin_phylo}/bin_tree/RAxML_bestTree.bin_faas_refined.tre \
            {vb_tax}/bin_tree.tre")
        script.phylotree_and_heatmap.R(
            input=self.bin_abc_tax_path,
            tree=os.path.join(self.vb_tax, "bin_tree.tre"),
            output=self.vb_tax
        )

    def function_vb(self):
        # KEGG
        SummBinAnnos(
            inputs=os.path.join(self.bin_kegg, "*/*kegg_raw.txt"),
            out_dir=self.vb_kegg,
            basename="KO_",
            outputs=["summary", "bin_detail"]
        ).pipeline()
        SummBinAnnos(
            inputs=os.path.join(self.bin_kegg, "*/*kegg_raw.txt"),
            out_dir=self.vb_kegg,
            basename="Pathway_",
            regroup=os.path.join(self.fmap_home, "FMAP_data/KEGG_orthology2pathway.txt"),
            levels=os.path.join(self.mapfiles_home, "kegg/kegg_pathways_merge_levels.tsv"),
            outputs=["summary", "levels_summary", "bin_detail_levels"]
        ).pipeline()
        # CAZyme
        SummBinAnnos(
            inputs=os.path.join(self.bin_cazy, "*.raw"),
            out_dir=self.vb_cazy,
            adjust_func=lambda x: [c.split("_")[0] for c in x.split("|")[1:] if re.search("^[A-Z]", c)],
            levels=lambda x: re.search(r'^\D*', x).group(),
            basename="CAZy_",
            outputs=["summary", "bin_detail", "gene_detail", "levels_summary", "bin_detail_levels"]
        ).pipeline()
        # GO
        SummBinAnnos(
            inputs=os.path.join(self.bin_emapper, "*/*.annotations"),
            out_dir=self.vb_go,
            column=5,
            adjust_func=lambda x: x.split(','),
            levels=os.path.join(self.mapfiles_home, "go/go.annotation.txt"),
            basename="GO_",
            outputs=["summary", "bin_detail", "levels_summary"]
        ).pipeline()
        # COG
        SummBinAnnos(
            inputs=os.path.join(self.bin_prokka, "*/*.tsv"),
            out_dir=self.vb_cog,
            skip=1,
            column=5,
            basename="COG_",
            outputs=["summary"]
        ).pipeline()
        SummBinAnnos(
            inputs=os.path.join(self.bin_prokka, "*/*.tsv"),
            out_dir=self.vb_cog,
            skip=1,
            column=5,
            regroup=os.path.join(self.mapfiles_home, "cog/COG.funccat.txt"),
            levels=os.path.join(self.mapfiles_home, "cog/funccat_levels.txt"),
            basename="functional_category_",
            outputs=["summary", "bin_detail", "levels_summary"]
        ).pipeline()

    def bin_visualize(self):
        # self.function_vb()
        # self.gather_data_process()
        # self.plot_bin()
        # self.analyze_bin_abc()
        self.plot_bin_phylo()


class BinningMixin(VisualizeBinMixin):
    """docstring for Binning"""

    def bin_prepare(self):
        self.set_path(
            force=True,
            bin_result="{out_dir}/Bin_all",
            bin_out="{bin_result}/bin_out",
            bin_sum="{bin_result}/bin_summary",
            bin_all="{bin_sum}/bin_all",
            bin_drep="{bin_sum}/bin_drep",
            bin_quant="{bin_result}/Bin_quant",
            bin_index="{bin_quant}/bin_index",
            bin_alignment="{bin_quant}/alignment_files",
            bin_pick="{bin_result}/Bin_pick",
            bin_prokka="{bin_result}/Bin_prokka",
            bin_cazy="{bin_result}/Bin_cazy",
            bin_phylo="{bin_result}/Bin_phylo",
            bin_kegg="{bin_result}/kegg",
            bin_emapper="{bin_result}/emapper",
            bin_faas="{bin_result}/bin_faas",
            bin_tmp="{out_dir}/tmp"
        )
        self.bin_faa_pattern = "{bin_prokka}/{bin}/{bin}.faa"
        self.bin_ffn_pattern = "{bin_prokka}/{bin}/{bin}.ffn"
        self.visualize_bin_prepare()

    def resub_file_names(pattern, replace, directory):
        for file in os.listdir(directory):
            new_file = re.sub(pattern, replace, file)
            file = os.path.join(directory, file)
            new_file = os.path.join(directory, new_file)
            os.rename(file, new_file)

    @synchronize
    def binning(self, fq_list, threads=2, mem=40, refinem="yes", mix_contigs=False):
        parsed_fqs = self.parse_fq_list(fq_list)
        sample_dir = os.path.join(self.bin_out, parsed_fqs['sample'])
        if mix_contigs:
            contigs = os.path.join(self.assembly_out, "final.contigs.fa")
        else:
            contigs = os.path.join(
                self.assembly_out,
                parsed_fqs['sample'],
                "final.contigs.fa"
            )
        with metabat2():
            self.system("""
echo 'mkdir -p {sample_dir}/bowtie2_db {sample_dir}/Bin {sample_dir}/stats \
 {sample_dir}/outliers {sample_dir}/filtered && \
{bowtie2_home}/bowtie2-build {contigs} {sample_dir}/bowtie2_db/{sample} \
 --threads {threads} --quiet && \
{bowtie2_home}/bowtie2 --mm -1 {r1} -2 {r2} -p {threads} \
 -x {sample_dir}/bowtie2_db/{sample} | \
 samtools sort --threads {threads} -o {sample_dir}/align.sorted.bam - && \
jgi_summarize_bam_contig_depths --outputDepth {sample_dir}/depth.txt \
 {sample_dir}/align.sorted.bam && \
metabat2 -m 1500 -t {threads} -i {contigs} \
 -a {sample_dir}/depth.txt -o {sample_dir}/Bin/bin.{sample} -v && \
if [ yes = {refinem} ]
then
 samtools index {sample_dir}/align.sorted.bam && \
 {refinem_path} scaffold_stats --genome_ext fa -c {threads} {contigs} \
  {sample_dir}/Bin/ {sample_dir}/stats {sample_dir}/align.sorted.bam && \
 {refinem_path} outliers {sample_dir}/stats/scaffold_stats.tsv {sample_dir}/outliers && \
 {refinem_path} filter_bins --genome_ext fa {sample_dir}/Bin/ \
  {sample_dir}/outliers/outliers.tsv {sample_dir}/filtered && \
 rm {sample_dir}/filtered/refinem.log;
else
 mv {sample_dir}/Bin/* {sample_dir}/filtered/;
fi && \
touch {bin_out}/{sample}/done;
 ' | qsub -l h_vmem={mem}G -pe {sge_pe} {threads} -q {sge_queue} -V -N {sample} \
 -o {bin_out} -e {bin_out}
            """, **parsed_fqs, **locals())

    def run_binning(self):
        self.binning(
            self.clean_paired_list,
            clean_before=["{bin_out}/{sample}"],
            pass_if_exists=["{bin_out}/{sample}/done"],
            **self.alloc_src("metabat2")
        )

    def drep_bins(self, threads=45):
        for dr in [self.bin_all, self.bin_pick]:
            if os.listdir(dr):
                self.system("rm {}/*".format(dr))
        self.system("""
for file in {bin_out}/*/filtered/*.fa;
    do base=$(basename $file) && \
     base=${{base//.filtered.fa/.fa}} && \
     ln -s $file {bin_all}/$base;
done && \
{drep_path} dereplicate {bin_drep} -g {bin_all}/*.fa -p {threads} --debug && \
ln -s {bin_drep}/dereplicated_genomes/* {bin_pick}/
            """, **locals())

    def index_bin_contigs(self, threads=45, entire_contigs="no"):
        self.system("""
if [ {entire_contigs} = no ]; then
    cat {bin_pick}/* > {bin_quant}/binned_contigs.fa
else
    ln -s {assembly_out}/final.contigs.fa {bin_quant}/binned_contigs.fa
fi && \
{salmon_path} index -p {threads} -t {bin_quant}/binned_contigs.fa \
 -i {bin_index}""", **locals())

    @synchronize
    def quant_bin_contigs(self, fq_list, threads=2, mem=40):
        self.system("""
echo '{salmon_path} quant -i {bin_index} --libType IU -1 {r1} -2 {r2} \
 -o {bin_alignment}/{sample}.quant --meta -p {threads} && \
touch {bin_alignment}/{sample}.quant/done' \
 | qsub -l h_vmem={mem}G -pe {sge_pe} {threads} -q {sge_queue} -V -N {sample} \
 -o {bin_quant} -e {bin_quant}""", **self.parse_fq_list(fq_list), **locals())

    def sum_bin_quants(self):
        self.system("bash {base_dir}/sum_bin_quants.sh {metawrap_scripts_home} \
 {bin_quant} {bin_pick} {bin_quant}/binned_contigs.fa")

    def phylophlan_tax_bin(self, threads=45):
        self.system("{phylophlan_bin}/phylophlan_metagenomic \
 -i {bin_pick} -o {bin_phylo} --nproc {threads} -n 1 -d SGB.Sep20 \
 --database_folder {phylophlan_databases}", **locals())

    def phylophlan_tree_bin(self, threads=20):
        if os.listdir(self.bin_faas):
            self.system("rm {bin_faas}/*")
        self.system("ln -s {bin_prokka}/*/*.faa {bin_faas}/ && \
{phylophlan_bin}/phylophlan --diversity medium -d phylophlan \
 --accurate -f {phylophlan_config} -i {bin_faas}/ --output_folder {bin_phylo} -o bin_tree \
 -t a --nproc {threads} --verbose --databases_folder {phylophlan_databases}", **locals())

    def bin_tax_level(self, level="k"):
        pattern = r"{}__[^|]+".format(level)
        tax_file = os.path.join(self.bin_phylo, "Bin_phylo.tsv")
        tax_dict = {}
        with open(tax_file) as fh:
            for num, line in enumerate(fh):
                if num > 0:
                    li = line.split()
                    tax_dict[li[0]] = re.search(pattern, li[1]).group().lstrip(
                        "{}__".format(level)
                    )
        return tax_dict

    @property
    def bin_tax_kingdom(self):
        return self.bin_tax_level('k')

    def quant_bins(self):
        self.index_bin_contigs(threads=self.this_threads, entire_contigs="no")
        self.quant_bin_contigs(
            fq_list=self.clean_paired_list,
            clean_before=["{bin_alignment}/{sample}.quant"],
            pass_if_exists=["{bin_alignment}/{sample}.quant/done"],
            **self.alloc_src("salmon")
        )
        self.sum_bin_quants()

    @property
    def bin_ids(self):
        ids = []
        for f in os.listdir(self.bin_pick):
            if f.endswith('.fa'):
                ids.append(f.rstrip('.fa'))
        return ids

    def map_bin_list(self, pattern):
        return [pattern.format(bin=b, **self.context) for b in self.bin_ids]

    @synchronize
    def prokka_bins(self, fq_list, threads=2, mem=40):
        parsed_fqs = self.parse_fq_list(fq_list)
        kingdoms = ["Archaea", "Bacteria", "Mitochondria", "Viruses"]
        bin_id = parsed_fqs['bin']
        kingdom = self.bin_tax_kingdom.get(bin_id)
        print("Found kingdom for {}: {}".format(bin_id, kingdom))
        if kingdom not in kingdoms:
            print("Unspported kingdom: {}".format(kingdom))
            kingdom = "Bacteria"
        with prokka():
            self.system("""
echo 'prokka {bin_path} --outdir {bin_prokka}/{bin} --prefix {bin} \
 --metagenome --cpus {threads} --kingdom {kingdom} && \
touch {bin_prokka}/{bin}/done' \
 | qsub -l h_vmem={mem}G -pe {sge_pe} {threads} -q {sge_queue} -V -N {bin} \
 -o {bin_prokka} -e {bin_prokka}""", **parsed_fqs, **locals())

    def run_prokka_bins(self):
        self.prokka_bins(
            fq_list=self.map_bin_list("{bin_pick}/{bin}.fa"),
            pass_if_exists=["{bin_prokka}/{bin}/done"],
            clean_before=["{bin_prokka}/{bin}"],
            **self.alloc_src("prokka")
        )

    @synchronize
    def diamond_bin_cazy(self, fq_list, threads=2, mem=40):
        parsed_fqs = self.parse_fq_list(fq_list)
        bin_id = parsed_fqs['bin']
        self.diamond(
            fa=fq_list,
            database=self.cazy_database,
            out_path=os.path.join(self.bin_cazy, bin_id + '.raw'),
            threads=threads, mem=mem, top=False,
            done_file=bin_id + '.done'
        )

    def run_diamond_bin_cazy(self):
        self.diamond_bin_cazy(
            fq_list=self.map_bin_list(self.bin_faa_pattern),
            pass_if_exists=["{bin_cazy}/{bin}.done"],
            clean_before=["{bin_cazy}/{bin}.raw"],
            **self.alloc_src("diamond")
        )

    @synchronize
    def kofamscan_bin(self, fq_list, threads=2, mem=40):
        self.system("""
echo 'mkdir -p {bin_kegg}/{bin}/_tmp && \
{kofamscan_path} -f mapper -o {bin_kegg}/{bin}/{bin}_kegg_raw.txt {fq_list} \
 --profile={kofamscan_database}/profiles --ko-list={kofamscan_database}/ko_list \
 --cpu={threads} --tmp-dir={bin_kegg}/{bin}/_tmp && \
rm -r {bin_kegg}/{bin}/_tmp && \
touch {bin_kegg}/{bin}/done' \
 | qsub -l h_vmem={mem}G -pe {sge_pe} {threads} -q {sge_queue} -V -N {bin} \
 -o {bin_kegg} -e {bin_kegg}""", **self.parse_fq_list(fq_list), **locals())

    def run_kofamscan_bin(self):
        self.kofamscan_bin(
            fq_list=self.map_bin_list(self.bin_faa_pattern),
            pass_if_exists=["{bin_kegg}/{bin}/done"],
            clean_before=["{bin_kegg}/{bin}"],
            **self.alloc_src("kofamscan")
        )

    @synchronize
    def emapper_bin(self, fq_list, threads=2, mem=40):
        self.system("""
echo 'mkdir -p {bin_emapper}/{bin}/_tmp && \
{emapper_path} -m diamond --seed_ortholog_evalue 0.00001 \
  --data_dir {emapper_database} --cpu {threads} --temp_dir {bin_emapper}/{bin}/_tmp \
  -i {fq_list} -o {bin_emapper}/{bin}/{bin} --usemem --override && \
rm -r {bin_emapper}/{bin}/_tmp && \
touch {bin_emapper}/{bin}/done' \
 | qsub -l h_vmem={mem}G -pe {sge_pe} {threads} -q {sge_queue} -V -N {bin} \
 -o {bin_emapper} -e {bin_emapper}""", **self.parse_fq_list(fq_list), **locals())

    def run_emapper_bin(self):
        self.emapper_bin(
            fq_list=self.map_bin_list(self.bin_faa_pattern),
            pass_if_exists=["{bin_emapper}/{bin}/done"],
            clean_before=["{bin_emapper}/{bin}"],
            **self.alloc_src("emapper")
        )

    def binning_pipeline(self):
        """
        self.run_binning()
        self.drep_bins(threads=self.this_threads)
        self.quant_bins()
        self.phylophlan_tree_bin(threads=self.this_threads)
        self.phylophlan_tax_bin(threads=self.this_threads)
        self.run_prokka_bins()
        self.run_diamond_bin_cazy()
        self.run_kofamscan_bin()
        self.run_emapper_bin()
        """
        self.bin_visualize()
