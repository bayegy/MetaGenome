import os
from .pyutils.tools import split_list, parse_premap
import re
from .visualizeAll import VisualizeAll
import pandas as pd
from Bayegy.ampliconLibs.systemMixin import SystemMixin
import numpy as np
from .pipconfig import settings
from lxml import html
from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor
from .groupGene import GroupGene
# import uuid
from multiprocessing import cpu_count
from .pyutils.controller import *
from .binning import BinningMixin
import logging
logger = logging.getLogger(__name__)


synchronize = sge_decorator if settings.use_sge else pool_decorator
parallel = sge_parallel if settings.use_sge else default_parallel


class RawDataNotPairedError(Exception):
    pass


class PathError(Exception):
    pass


class MetagenomePipline(SystemMixin, BinningMixin):
    """
    arguments:
        pre_mapping_file: The first column of pre_mapping_file should be smaple id in raw fastq files, the last column of pre_mapping_file should be revised id (or the same) of samples.

        categories: Categories seprated by ',' , optional, if not passed, the categories names should have the pattern of 'Category.*'

        host_type: Will control the de_host step.

        run_size: Control the max number of jobs submitted to sge each time

        raw_fqs_dir: Directory where the raw fastq file weres stored

        sample_regex: Regular expression to match sample id (contained by brackets)

        forward_regex: Regular expression to match forward fastq files

        reverse_regex: Regular expression to match reverse fastq files

        out_dir: Where to store the results

    sample usage:

    from metagnomePipline import MetagenomePipline
    m =MetagenomePipline('/home/cheng/Projects/rll_testdir/1.rawdata/','/home/cheng/Projects/rll_testdir/mapping_file.txt',out_dir="/home/cheng/Projects/rll_testdir/")
    m.run()
    """

    def __init__(self, raw_fqs_dir, pre_mapping_file, categories=False,
                 host_db="/home/bayegy/Databases/host_genome/hg38/hg38",
                 sample_regex=r"(.+)_.*_[12]\.fq\.?g?z?", forward_regex=r"_1\.fq\.?g?z?$",
                 reverse_regex=r"_2\.fq\.?g?z?$", out_dir='./', base_on_assembly=False,
                 exclude='none', filter_species=False, zip_kneaddata=True,
                 mix_asem=False, transcriptome=False, orders=False,
                 if_binning="no"):
        self.context = dict()

        self.set_path(force=True,
                      out_dir=out_dir,
                      raw_dir=out_dir + '/raw_dir',
                      salmon_out=out_dir + '/salmon_out',
                      fastqc_out=out_dir + '/fastqc_out',
                      sortmerna_out=out_dir + '/sortmerna_out',
                      kneaddata_out=out_dir + '/kneaddata_out',
                      kraken2_out=out_dir + '/Kraken2',
                      amr_out=out_dir + '/AMR',
                      humann2_out=out_dir + '/Metagenome/Humann/',
                      fmap_out=out_dir + '/FMAP/',
                      assembly_out=out_dir + "/Assembly_out",
                      metaphlan_out=out_dir + "/Metagenome/Metaphlan",
                      report_out=out_dir + "/Report",
                      assembled_db="{assembly_out}/assembled_db",
                      unassembled_out="{out_dir}/unassembled_out"
                      )

        self.set_attr(base_on_assembly=base_on_assembly,
                      # assembly_out=out_dir + "/Assembly_out",
                      mapping_file=self.out_dir + '/mapping_file.txt',
                      host_db=' -db '.join(host_db.split(',')),
                      exclude=exclude,
                      base_dir=os.path.dirname(__file__) + '/',
                      contigs_path=self.assembly_out + 'final.contigs.fa',
                      orf_path=self.assembly_out + 'final.contigs.fa.fna',
                      nucleotide_path=self.assembly_out + 'NR.nucleotide.fa',
                      protein_path=self.assembly_out + 'NR.protein.fa',
                      if_binning=if_binning,
                      # cdhit_out=self.assembly_out + "NR.nucleotide.fa",
                      **self.get_attrs(settings))

        # not nessesary for:
        # tmp_dir=os.path.join('/dev/shm/', str(uuid.uuid1()) + '_metagenomepipeline_tmp_dir')

        if not hasattr(self, 'tmp_dir'):
            self.set_path(force=True, tmp_dir="/dev/shm/")

        self.set_path(force=False, **self.path)
        self.parsed_map = parse_premap(
            raw_fqs_dir, pre_mapping_file, forward_regex, reverse_regex, sample_regex)
        self.fq_info = self.parsed_map['fastq'].values
        self.new_ids = list(self.parsed_map['map']['#SampleID'].values)
        self.new_ids.sort()
        self.parsed_map['map'].to_csv(self.mapping_file, sep='\t', index=False)
        self.set_attr(categories=categories if categories else ','.join(
            [g for g in self.parsed_map['map'].columns if not g.lower().find('category') == -1]))
        logger.info("The detected categories are: \n\n{}\n".format(self.categories))
        if orders:
            orders = orders.split(',')
        else:
            orders = [g for g in self.parsed_map['map'].columns if not g.find('Order') == -1]
        logger.info("The Orders are: {}".format(str(orders)))
        self.orders = orders

        self.raw_pattern = "{raw_dir}/{sample}_{direction}.fq"
        self.raw_list = self.paired_data(self.raw_pattern)
        self.sortmerna_pattern = "{sortmerna_out}/{sample}/{sample}_{direction}.fq"
        self.sortmerna_list = self.paired_data(self.sortmerna_pattern)
        self.unzip_clean_paired_pattern = "{kneaddata_out}/{sample}/{sample}_R1_kneaddata_paired_{direction_num}.fastq"
        self.clean_paired_pattern = self.unzip_clean_paired_pattern + \
            '.gz' if zip_kneaddata else self.unzip_clean_paired_pattern
        self.clean_r1_list = self.map_list(
            self.clean_paired_pattern, use_direction='R1')
        self.clean_paired_list = self.paired_data(self.clean_paired_pattern)
        self.kracken2_reports_pattern = "{kraken2_out}/{sample}.report"
        self.unassembled_pattern = "{unassembled_out}/{sample}/{sample}_R{direction_num}_unassembled_merged.fq.gz"
        self.kracken2_reports_list = self.map_list(
            self.kracken2_reports_pattern,
            use_direction='R1'
        )
        self.salmon_quant_pattern = "{salmon_out}/{sample}.quant/quant.sf"
        # self.salmon_quant_list = self.map_list(self.salmon_quant_pattern, use_direction="R1")
        self.running_list = self.out_dir + '.running_list'
        self.escape_sge = not settings.use_sge  # global setting of escape_sge
        self.this_threads = int(cpu_count() * 0.9)
        self.filter_species = filter_species
        self.zip_kneaddata = zip_kneaddata
        self.mix_asem = mix_asem
        self.transcriptome = transcriptome
        if if_binning != "no":
            self.bin_prepare()

    def format_raw(self, processors=3):
        executor = ThreadPoolExecutor(max_workers=processors)
        for fq_path, new_id, direction in self.fq_info:
            if fq_path.endswith('.gz'):
                executor.submit(self.unzip, fq_path)
        executor.shutdown(True)
        for fq_path, new_id, direction in self.fq_info:
            formated_fq = self.raw_pattern.format(
                sample=new_id, direction=direction, **self.context)
            if not os.path.exists(formated_fq):
                fq_path = fq_path.rstrip('.gz')
                os.system("ln -s {} {}".format(fq_path, formated_fq))

    def unzip(self, fq_path):
        print("Unzip file {}...".format(fq_path))
        os.system("gunzip {}".format(fq_path))

    def map_list(self, pattern=False, each=False, use_direction="both"):
        out = []
        for fq_path, new_id, direction in self.fq_info:
            ele = pattern.format(sample=new_id, direction=direction,
                                 direction_num=re.search(r"\d$", direction).group(),
                                 **self.context) if pattern else [new_id, direction]
            if direction == "R2" and (use_direction == "both" or use_direction == "R2"):
                out.append(ele)
            elif direction == "R1" and (use_direction == "both" or use_direction == "R1"):
                out.append(ele)
        out.sort()
        return split_list(out, each) if each else out

    def paired_data(self, pattern, each=False):
        r1 = self.map_list(pattern=pattern, each=False, use_direction="R1")
        r2 = self.map_list(pattern=pattern, each=False, use_direction="R2")
        if not len(r1) == len(r2):
            raise RawDataNotPairedError(
                "The file number of Reads1 and Reads2 is not equal")
        out = np.array([r1, r2]).T.tolist()
        return split_list(out, each) if each else out

    @staticmethod
    def mirror_move(source, target):
        for root, dirs, files in os.walk(source):
            target_dir = root.replace(source.rstrip('/'), target, 1)
            if target_dir.startswith('~'):
                raise PathError("'~' is illegal in path!")
            if not os.path.exists(target_dir):
                print("Making directory: {}".format(target_dir))
                os.system("mkdir -p '{}'".format(target_dir))
            for fl in files:
                source_file = os.path.join(root, fl)
                os.system("mv '{}' '{}'".format(source_file, target_dir))

    def homized_cmd(self, cmd, home=False):
        home = home or self.base_dir
        return "cd {}&&{}".format(home, cmd)

    def parse_fq_list(self, fq_list):
        if isinstance(fq_list, list):
            return {
                "r1": fq_list[0],
                "r2": fq_list[1],
                "sample": os.path.basename(fq_list[0]).split('_')[0]
            }

        basename = os.path.basename(fq_list)
        if basename.startswith('bin.'):
            return dict(
                bin_path=fq_list,
                bin=os.path.splitext(basename)[0]
            )

        return {
            "r1": fq_list,
            "sample": basename.split('_')[0]
        }

    @synchronize
    def run_fastqc(self, fq_list, threads=10, mem=20):
        self.system("""
echo 'mkdir -p {fastqc_out}/{sample}/temp_dir&&\
{fastqc_home}/fastqc {r1} -o {fastqc_out}/{sample} -t {threads} \
 -d {fastqc_out}/{sample}/temp_dir&&\
{fastqc_home}/fastqc {r2} -o {fastqc_out}/{sample} -t {threads} \
 -d {fastqc_out}/{sample}/temp_dir&&\
rm -r {fastqc_out}/{sample}/temp_dir {fastqc_out}/{sample}/*zip&&\
touch {fastqc_out}/{sample}/done' \
  | qsub -l h_vmem={mem}G -pe {sge_pe} {threads} -q {sge_queue} -V -N {sample} \
 -o {fastqc_out} -e {fastqc_out}""", **self.parse_fq_list(fq_list), **locals())

    @synchronize
    def run_sortmerna(self, fq_list, threads=10, mem=20):
        self.system("""
echo 'mkdir -p {sortmerna_out}/{sample}/temp_dir&&{sortmerna_path} \
--ref {sortmerna_databases}/*.fasta \
--reads {r1} \
--reads {r2} \
--fastx \
--paired_out \
--threads {threads} \
-v \
--out2 \
--workdir {sortmerna_out}/{sample}/temp_dir \
--other {sortmerna_out}/{sample}/{sample} && \
mv {sortmerna_out}/{sample}/{sample}_fwd.fq {sortmerna_out}/{sample}/{sample}_R1.fq&& \
mv {sortmerna_out}/{sample}/{sample}_rev.fq {sortmerna_out}/{sample}/{sample}_R2.fq&& \
rm -r {sortmerna_out}/{sample}/temp_dir&& touch {sortmerna_out}/{sample}/done' \
  | qsub -l h_vmem={mem}G -pe {sge_pe} {threads} -q {sge_queue} -V -N {sample} \
 -o {sortmerna_out} -e {sortmerna_out}""", **self.parse_fq_list(fq_list), **locals())

    @synchronize
    def run_kneaddata(self, fq_list, threads=10, mem=20, gzip=True, remove_input=False):
        inputs = self.parse_fq_list(fq_list)
        r1o = self.unzip_clean_paired_pattern.format(
            **inputs,
            direction_num=1,
            **self.context
        )
        r2o = self.unzip_clean_paired_pattern.format(
            **inputs,
            direction_num=2,
            **self.context
        )
        cmds = []
        if gzip:
            cmds.append("gzip {}&&gzip {}".format(r1o, r2o))
        cmds.append('touch {kneaddata_out}/{sample}/done'.format(
            **inputs,
            kneaddata_out=self.kneaddata_out
        ))
        if remove_input:
            cmds.append("rm {r1} {r2}".format(**inputs))
        a_cmd = "&&".join(cmds)
        self.system("""
echo 'mkdir -p {kneaddata_out}/{sample}&&{kneaddata_path} \
 -i {r1} -i {r2} -o {kneaddata_out}/{sample} -v -db {host_db} \
 --remove-intermediate-output  -t {threads} --run-fastqc-end \
 --trimmomatic {trimmomatic_home} --trimmomatic-options \
  "ILLUMINACLIP:{adapters_path}:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50" --max-memory {mem}g \
 --bowtie2 {bowtie2_home} \
 --fastqc  {fastqc_home} && {a_cmd}' \
 | qsub -l h_vmem={mem}G -pe {sge_pe} {threads} -q {sge_queue} -V -N {sample} \
 -o {kneaddata_out} -e {kneaddata_out}""", **self.parse_fq_list(fq_list), **locals())

    @pool_decorator
    def gzip_kneaddata(self, fq_list):
        self.system("gzip {r1}&&gzip {r2}", **self.parse_fq_list(fq_list))

    @synchronize
    def run_kraken2(self, fq_list, threads=2, mem=70):
        self.system("""
echo '{kraken2_path} --db {kraken2_database} --threads {threads} --confidence 0.2 \
 --report {kraken2_out}/{sample}.report --paired {r1} {r2} --output - && \
touch {kraken2_out}/{sample}_done' \
  | qsub -l h_vmem={mem}G -pe {sge_pe} {threads} -q {sge_queue} -V -N {sample} \
 -o {kraken2_out} -e {kraken2_out}""", **self.parse_fq_list(fq_list), **locals())

    def run_bracken(self):
        for file in [
            "{kraken2_out}/All.Taxa.OTU.txt",
            "{kraken2_out}/All.Taxa.txt"
        ]:
            file = file.format(**self.context)
            if os.path.exists(file):
                os.remove(file)
        bracken_list = [e + '.bracken' for e in self.kracken2_reports_list]
        for report in self.kracken2_reports_list:
            self.system(
                "python2 {bracken_path} -i {report} -k {bracken_database} \
                 -l S -o {report}.bracken", report=report)
        with open(self.running_list, 'w') as f:
            f.write('\n'.join(bracken_list))
        os.system(
            self.homized_cmd("perl Braken_to_OTUtable.pl {} {}".format(
                self.ncbi_taxaID_path,
                self.running_list)
            )
        )

    def clean_header(self, tb_name, pattern="_.*$", skip=[]):
        df = pd.read_csv(tb_name, sep='\t')
        df.columns = [re.sub(pattern, '', c) for c in df.columns]
        df.drop(index=skip).to_csv(tb_name, sep="\t", index=False)

    @synchronize
    def run_metaphlan2(self, fq_list, threads, mem):
        self.system("""
echo '{metaphlan2_home}/metaphlan2.py --input_type multifastq --bowtie2_exe {bowtie2_home}/bowtie2 \
 --bowtie2_build {bowtie2_home}/bowtie2-build  --tmp_dir {tmp_dir} \
 --bowtie2db {metaphlan2_database} \
 --nproc {threads} --bowtie2out {tmp_dir}/{sample}.bt2.out.txt {r1} \
 -o {metaphlan_out}/{sample}.metaphlan.profile.txt --sample_id {sample}' \
  | qsub -l h_vmem={mem}G -pe {sge_pe} {threads} -q {sge_queue} -V -N {sample} \
 -o {metaphlan_out} -e {metaphlan_out}""", **self.parse_fq_list(fq_list), **locals())

    def join_metaphlan(self):
        tb_name = self.metaphlan_out + '/All.Metaphlan2.profile.txt'
        self.system(
            "python2 {base_dir}/merge_metaphlan_tables.py {metaphlan_out}/*profile.txt > {tb_name}",
            tb_name=tb_name
        )
        self.clean_header(tb_name, pattern=r'\.metaphlan\.profile$', skip=[0])

    @synchronize
    def run_humann2(self, fq_list, threads=2, mem=40):
        self.system("""
echo 'mkdir -p {humann2_out}/{sample}&&{humann2_home}/humann2 --input {r1} \
  --bowtie2 {bowtie2_home} --metaphlan {metaphlan2_home} --diamond {diamond_home} \
  --protein-database {humann2_protein_database} --nucleotide-database {humann2_nucleotide_database} \
  --threads {threads} --memory-use  maximum \
  --output-basename {sample} \
  --output {humann2_out}/{sample} \
  --metaphlan-options "-t rel_ab --bowtie2db {metaphlan2_database} --sample_id {sample}" && \
mv {humann2_out}/{sample}/*/*bugs_list.tsv {humann2_out}/{sample}/*genefamilies.tsv \
 {humann2_out}/{sample}/*pathabundance.tsv {humann2_out}/&& \
rm -r {humann2_out}/{sample} && touch {humann2_out}/{sample}_done' \
  | qsub -l h_vmem={mem}G -pe {sge_pe} {threads} -q {sge_queue} -V \
 -N {sample} -o {humann2_out} -e {humann2_out}""", **self.parse_fq_list(fq_list), **locals())

    def join_humann(self):
        out_files = [f.format(**self.context) for f in [
            "{humann2_out}/RPK.All.UniRef90.genefamilies.tsv",
            "{humann2_out}/RPK.All.Metacyc.pathabundance.tsv",
            "{metaphlan_out}/All.Metaphlan2.profile.txt"
        ]]
        for file in out_files:
            if os.path.exists(file):
                os.remove(file)
            norm_file = file.replace('RPK.', '')
            if os.path.exists(norm_file):
                os.remove(norm_file)
        self.system('''
{humann2_home}/humann2_join_tables -i {humann2_out}/ \
 -o {humann2_out}/RPK.All.UniRef90.genefamilies.tsv --file_name genefamilies.tsv
{humann2_home}/humann2_join_tables -i {humann2_out}/ \
 -o {humann2_out}/RPK.All.Metacyc.pathabundance.tsv --file_name pathabundance.tsv
{humann2_home}/humann2_join_tables -i {humann2_out}/ \
 -o {metaphlan_out}/All.Metaphlan2.profile.txt --file_name bugs_list.tsv''')
        for file in out_files:
            self.clean_header(file)

    def set_fmap_db(self, database):
        with open(self.fmap_home + '/FMAP_data/database', 'w') as f:
            f.write(database)

    @synchronize
    def run_fmap(self, fq_list, threads=4, mem=24, evalue=0.00001):
        self.system("""
echo 'perl {fmap_home}/FMAP_mapping.pl -e {evalue}  -p {threads} {r1} \
 > {fmap_out}/{sample}.mapping.txt' | \
 qsub -l h_vmem={mem}G -pe {sge_pe} {threads} -q {sge_queue} -V \
  -N {sample} -o {fmap_out} -e {fmap_out}""", **self.parse_fq_list(fq_list), **locals())

    def quantify_fmap(self, all_name="all.txt", print_definition=False):
        p = "" if print_definition else "-n"
        for new_id in self.new_ids:
            self.system(
                "perl {fmap_home}/FMAP_quantification.pl {fmap_out}/{sample}.mapping.txt \
            > {fmap_out}/{sample}.abundance.txt", sample=new_id)
        args = [
            "{sample}={fmap_out}/{sample}.abundance.txt".format(
                sample=new_id,
                fmap_out=self.fmap_out
            ) for new_id in self.new_ids]
        self.system("perl {fmap_home}/FMAP_table.pl {p} {args} > {fmap_out}/{all_name}",
                    p=p, args=" ".join(args), all_name=all_name)

    def fmap_wrapper(self, fq_list, run_type="KEGG", threads=6, mem=10, evalue=0.00001):
        fmap_db = {
            "KEGG": "orthology_uniref90_2_2157_4751.20190412161853",
            "AMR": "protein_fasta_protein_homolog_model_cleaned",
            "ARDB": "ARDB.20180725064354",
            "MGE": "MGEs_FINAL_99perc_trim"
        }
        self.set_fmap_db(fmap_db[run_type])
        self.run_fmap(fq_list=fq_list, threads=threads,
                      mem=mem, evalue=evalue)
        self.quantify_fmap(
            all_name="All.{}.abundance_unstratified.tsv".format(run_type))

    @synchronize
    def assembly(self, fq_list, threads=4, mem=24):
        # 单样品组装
        # parsed_fqs = self.parse_fq_list(fq_list)
        # sample_out = os.path.join(self.assembly_out, parsed_fqs['sample'])
        # if os.path.exists(os.path.join(sample_out, 'all_done')):
        #     return
        # if os.path.exists(sample_out):
        #     self.system("rm -r {}".format(sample_out))
        mem_p = mem * 1000000000
        self.system("""
echo 'mkdir -p {tmp_dir}/{sample}&&{megahit_path} --k-list 21,29,39,59,79,99,119,141 \
 --tmp-dir {tmp_dir}/{sample} -m {mem_p} --mem-flag 1  \
 -1 {r1} -2 {r2} --min-contig-len 500 -t {threads} -o {assembly_out}/{sample} && \
rm -r {assembly_out}/{sample}/intermediate_contigs {tmp_dir}/{sample} && \
touch {assembly_out}/{sample}/done' | \
qsub -l h_vmem={mem}G -pe {sge_pe} {threads} -q {sge_queue} -V -N {sample} \
 -o {assembly_out} -e {assembly_out}""", **self.parse_fq_list(fq_list), **locals())

    def build_assembled(self, threads=4):
        self.system("""
cat {assembly_out}/*/final.contigs.fa > {assembly_out}/all.samples.contigs.fa && \
{python3_path} {base_dir}/rename_contigs.py -i {assembly_out}/all.samples.contigs.fa \
  --replaceinput && \
{bowtie2_home}/bowtie2-build --threads \
 {threads} {assembly_out}/all.samples.contigs.fa {assembled_db}/assembled""", **locals())

    @synchronize
    def unassembled(self, fq_list, threads=4, mem=24):
        self.system("""
echo 'mkdir -p {unassembled_out}/{sample}&&{bowtie2_home}/bowtie2 \
 --mm --threads {threads} -x {assembled_db}/assembled \
 -U {r1} --end-to-end --sensitive -S {sam_out} \
 --un {unassembled_out}/{sample}/{sample}_R1_unassembled.fastq && \
{bowtie2_home}/bowtie2 --mm --threads {threads} -x {assembled_db}/assembled \
 -U {r2} --end-to-end --sensitive -S {sam_out} \
 --un {unassembled_out}/{sample}/{sample}_R2_unassembled.fastq && \
{python3_path} {base_dir}/merge_se_to_pe.py \
 -1 {unassembled_out}/{sample}/{sample}_R1_unassembled.fastq \
 -2 {unassembled_out}/{sample}/{sample}_R2_unassembled.fastq \
 -o {unassembled_out}/{sample}/{sample}_R%s_unassembled_merged.fq  --replaceinput && \
gzip {unassembled_out}/{sample}/{sample}_R*_unassembled_merged.fq && \
touch {unassembled_out}/{sample}/done' | \
 qsub -l h_vmem={mem}G -pe {sge_pe} {threads} -q {sge_queue} -V -N {sample} \
  -o {unassembled_out} -e {unassembled_out}""",
                    **self.parse_fq_list(fq_list), **locals(), sam_out=os.devnull)

    def mix_assembly(self, threads=50):
        # 混合组装
        self.system(
            "{megahit_path} --continue  --kmin-1pass \
             --k-list 21,29,39,59,79,99,119,141 -m 0.94 --mem-flag 0 \
             -1 {r1_list} -2 {r2_list} --min-contig-len 500 -t {threads} \
             -o {assembly_out}/mixed_assembly",
            r1_list=','.join(self.map_list(
                self.unassembled_pattern, use_direction="R1")),
            r2_list=','.join(self.map_list(
                self.unassembled_pattern, use_direction="R2")),
            threads=threads
        )

    def cat_assembly(self, threads=50):
        # 基因预测，去冗余
        self.system(
            """
{python3_path} {base_dir}/rename_contigs.py -i {assembly_out}/*/final.contigs.fa \
  --replaceinput && \
cat {assembly_out}/*/final.contigs.fa > {assembly_out}/final.contigs.fa
""")

    @parallel
    def prodigal(self, fa, out_dir, threads=4, mem=24, escape_sge=False):
        self.system("""
echo '{prodigal_path}  -q -i {fa} -a {out_dir}/final.contigs.fa.faa \
 -d {out_dir}/final.contigs.fa.fna \
 -f gff -p meta -o {out_dir}/final.contigs.fa.gff && touch {out_dir}/done' \
| qsub -l h_vmem={mem}G -pe {sge_pe} {threads} -q {sge_queue} -V -N prodigal \
 -o {out_dir} -e {out_dir}""", **locals())

    def run_prodigal(self):
        self.prodigal(fa=self.contigs_path,
                      out_dir=self.assembly_out,
                      cat_file="final.contigs.fa.fna",
                      **self.alloc_src_para("prodigal")
                      )

    @parallel
    def cdhit(self, fa, out_dir, threads=4, mem=24, escape_sge=False):
        mem_mb = mem * 1000
        # name = out_dir.strip('/').split('/')[-1]
        self.system("""
echo '{cdhit_path} -T {threads} -i {fa} -d 0 -M {mem_mb} \
 -o {out_dir}/NR.nucleotide.fa && touch {out_dir}/done' \
| qsub -l h_vmem={mem}G -pe {sge_pe} {threads} -q {sge_queue} -V -N cdhit \
 -o {out_dir} -e {out_dir}""", **locals())

    def run_cdhit(self):
        self.cdhit(
            fa=self.orf_path,
            out_dir=self.assembly_out,
            cat_file="NR.nucleotide.fa",
            **self.alloc_src_para("cdhit")
        )

    def plot_orf(self):
        self.system(
            """cd {assembly_out} && \
grep '>' NR.nucleotide.fa > NR.nucleotide.fa.header && \
{acmd} && \
{R_path} {base_dir}/ORF_header_summary.R -i {assembly_out}""",
            acmd=self.homized_cmd(
                "{perl_path} ORF_generate_input_stats_file.pl \
                         {assembly_out}/NR.nucleotide.fa.header".format(
                    **self.context
                )
            )
        )

    def run_quast(self):
        self.system(
            "python2 {quast_path} -o {assembly_out}/quast_results/  {contigs_path}")

    def create_gene_db(self, threads=7):
        self.system('''
{salmon_path} index -t {assembly_out}/NR.nucleotide.fa -p {threads} \
 -k 31 -i {assembly_out}/salmon_index --keepFixedFasta && \
{transeq_path} -sequence {assembly_out}/salmon_index/ref_k31_fixed.fa \
 -outseq {assembly_out}/NR.protein.fa -trim Y && \
sed -i 's/_1\\b//' {assembly_out}/NR.protein.fa
        ''', threads=threads)

    @synchronize
    def quant_gene(self, fq_list, threads=7, mem=80):
        self.system("""
echo '{salmon_path} quant --validateMappings -i {assembly_out}/salmon_index \
 -l A -p {threads} --meta -1 {r1} -2 {r2} \
 -o {salmon_out}/{sample}.quant && touch {salmon_out}/{sample}.quant/done' \
| qsub -l h_vmem={mem}G -pe {sge_pe} {threads} -q {sge_queue} -V -N {sample} \
 -o {salmon_out} -e {salmon_out}""", **self.parse_fq_list(fq_list), **locals())

    def join_gene(self):
        self.system(
            "{salmon_path} quantmerge --quants {salmon_out}/*.quant \
             -o {salmon_out}/All.genes.abundance.txt")
        self.clean_header(
            self.out_dir + "salmon_out/All.genes.abundance.txt", ".quant$")

    def diamond(self, fa, database, out_path, threads=4, mem=24, max_target_seqs=1,
                evalue=0.00001, identity=80, top=3, escape_sge=False, done_file="done"):
        out_dir = os.path.dirname(out_path)
        fa_name = os.path.basename(fa)
        if top:
            hits_control = "--top " + str(top)
        else:
            hits_control = "--max-target-seqs " + str(max_target_seqs)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        self.system("""
echo '{diamond_home}/diamond blastp --db {database} \
 --query {fa} --outfmt 6 --threads {threads} \
 -e {evalue} --id {identity} {hits_control} --block-size 200 --index-chunks 1 --quiet \
  --out {out_path} && touch {out_dir}/{done_file}' \
| qsub -l h_vmem={mem}G -pe {sge_pe} {threads} -q {sge_queue} \
 -V -N {fa_name} -o {out_dir} -e {out_dir}""", **locals())

    @parallel
    def diamond_cazy(self, fa, out_dir, threads=4, mem=24, escape_sge=False):
        out_path = os.path.join(out_dir, "genes_cazy.f6")
        database = self.cazy_database
        kwargs = locals()
        kwargs = {k: v for k, v in kwargs.items() if k not in ["self", "out_dir"]}
        self.diamond(**kwargs)

    def run_diamond_cazy(self):
        self.diamond_cazy(
            fa=self.protein_path,
            out_dir=self.salmon_out,
            cat_file="genes_cazy.f6",
            **self.alloc_src_para("diamond")
        )

    @parallel
    def diamond_card(self, fa, out_dir, threads=4, mem=24, escape_sge=False):
        out_path = os.path.join(out_dir, "genes_card.f6")
        database = os.path.join(
            self.fmap_home,
            "FMAP_data/protein_fasta_protein_homolog_model_cleaned.dmnd"
        )
        kwargs = locals()
        kwargs = {k: v for k, v in kwargs.items() if k not in ["self", "out_dir"]}
        self.diamond(**kwargs)

    def run_diamond_card(self):
        self.diamond_card(
            fa=self.protein_path,
            out_dir=self.salmon_out,
            cat_file="genes_card.f6",
            **self.alloc_src_para("diamond")
        )

    @parallel
    def emapper(self, fa, out_dir, threads=4, mem=24, escape_sge=False):
        self.system("""
echo '{emapper_path} --no_file_comments -m diamond --seed_ortholog_evalue 0.00001 \
  --data_dir {emapper_database} --cpu {threads} --temp_dir {out_dir} \
  -i {fa} -o {out_dir}/genes --usemem --override && touch {out_dir}/done' \
| qsub -l h_vmem={mem}G -pe {sge_pe} {threads} -q {sge_queue} -V \
 -N emapper -o {out_dir} -e {out_dir}""", **locals())

    def run_emapper(self):
        self.emapper(
            fa=self.protein_path,
            out_dir=self.salmon_out,
            cat_file="genes.emapper.annotations",
            **self.alloc_src_para("emapper")
        )

    def alloc_src(self, proc):
        return {
            "threads": self.threads[proc],
            "mem": self.memorys[proc]
        }

    def alloc_src_para(self, proc):
        return {
            **self.alloc_src(proc),
            "splits": self.splits[proc]
        }

    def find_file(self, directory, pattern):
        for file in os.listdir(directory):
            if re.search(pattern, file):
                return os.path.join(directory, file)
        return None

    def search_fqc(self, fqc_file):
        with open(fqc_file) as f:
            fqc = f.read()
        tree = html.fromstring(fqc)
        num_reads = tree.xpath(
            '//td[text()="Total Sequences"]/following-sibling::td[1]/text()')
        gc_content = tree.xpath(
            '//td[text()="%GC"]/following-sibling::td[1]/text()')
        return int(num_reads[0]), int(gc_content[0])

    def q20_q30(self, read_file):
        qr = list(os.popen("{base_dir}/sim_fqstat {r1}".format(base_dir=self.base_dir, r1=read_file)))
        bases, q20_cnt, q30_cnt = [int(i) for i in qr[0].split()]
        return np.round(q20_cnt / bases * 100, 2), np.round(q30_cnt / bases * 100, 2)

    def get_qc_stats(self, sample_id):
        results = OrderedDict([
            ("Sample ID", sample_id),
            ("InsertSize(bp)", "350"),
            ("SeqStrategy", "(150:150)"),
        ])
        raw_fqc = "{fastqc_out}/{sample}/{sample}_R1_fastqc.html".format(**self.context, sample=sample_id)
        print("The raw fastqc file of {} is {}".format(sample_id, raw_fqc))
        raw_num_reads, raw_gc_content = self.search_fqc(raw_fqc)

        clean_fqc = "{kneaddata_out}/{sample}/fastqc/{sample}_R1_kneaddata_paired_1_fastqc.html".format(
            **self.context, sample=sample_id)
        print("The clean fastqc file of {} is {}".format(sample_id, clean_fqc))
        clean_num_reads, clean_gc_content = self.search_fqc(clean_fqc)
        raw_q20, raw_q30 = self.q20_q30(
            read_file=self.raw_pattern.format(sample=sample_id, direction="R1", **self.context))
        clean_q20, clean_q30 = self.q20_q30(
            read_file=self.clean_paired_pattern.format(sample=sample_id, direction_num=1, **self.context))
        results.update([
            ("RawReads(#)", raw_num_reads),
            ("Raw Base(GB)", np.round(raw_num_reads * 300 / 10**9, 2)),
            ("%GC", raw_gc_content),
            ("Raw Q20(%)", raw_q20),
            ("Raw Q30(%)", raw_q30),
            ("Clean Reads(#)", clean_num_reads),
            ("Cleaned(%)", np.round(clean_num_reads / raw_num_reads * 100, 2)),
            ("Clean Q20(%)", clean_q20),
            ("Clean Q30(%)", clean_q30),
        ])
        return results

    def generate_qc_report(self, processors=2):
        executor = ThreadPoolExecutor(max_workers=processors)
        futures = []
        for sample in self.new_ids:
            future = executor.submit(self.get_qc_stats, sample)
            futures.append(future)
        executor.shutdown(True)
        with open(os.path.join(self.report_out, "reads_summary.txt"), 'w') as f:
            for num, future in enumerate(futures):
                stat = future.result()
                if num == 0:
                    f.write('\t'.join(stat.keys()) + '\n')
                vals = [str(v) for v in stat.values()]
                f.write('\t'.join(vals) + '\n')

    def group_gene(self, db, processors=10):
        annotation_files = {
            "AMR": os.path.join(self.salmon_out, "genes_card.f6"),
            "CAZY": os.path.join(self.salmon_out, "genes_cazy.f6"),
            "KO": os.path.join(self.salmon_out, "genes.emapper.annotations"),
            "EGGNOG": os.path.join(self.salmon_out, "genes.emapper.annotations"),
            "GO": os.path.join(self.salmon_out, "genes.emapper.annotations")
        }
        adjust_funcs = {
            "KO": lambda x: x.split(','),
            "GO": lambda x: x.split(','),
            # "EGGNOG": lambda x: ['ENOG41' + c if not c.startswith('COG') else c for c in re.findall('([^,]+)@NOG', x)],
            "EGGNOG": lambda x: [c for c in re.findall('([^,]+)@NOG', x) if c.startswith('COG')],
            "CAZY": lambda x: [c.split("_")[0] for c in x.split("|")[1:] if re.search("^[A-Z]", c)],
            "AMR": False
        }
        columns = {
            "KO": 6,
            "GO": 5,
            "EGGNOG": 9,
            "CAZY": 1,
            "AMR": 1
        }
        gg = GroupGene(annotation_file=annotation_files[db],
                       by=columns[db],
                       adjust_func=adjust_funcs[db])

        in_out_dict = {}
        for sample_id in self.new_ids:
            gene_abundance = self.salmon_quant_pattern.format(sample=sample_id, **self.context)
            out_file = os.path.join(self.salmon_out, "{}_{}_abundance.txt".format(sample_id, db))
            in_out_dict[gene_abundance] = out_file
            # executor.submit(gg.group, gene_abundance, out_file)
            # futures.append(future)
        # executor.shutdown(True)
        gg.map(in_out_dict, processors=processors)
        self.system("""
{humann2_home}/humann2_join_tables -i {salmon_out}/ -o \
 {salmon_out}/All.{db}.abundance_unstratified.tsv --file_name _{db}_abundance.txt
            """, db=db)
        self.clean_header(os.path.join(self.salmon_out,
                                       "All.{}.abundance_unstratified.tsv".format(db)))

    def visualize(self, **kwargs):
        VisualizeAll(
            self.mapping_file,
            self.categories,
            out_dir=self.out_dir,
            filter_species=self.filter_species
        ).visualize(
            self.exclude, self.base_on_assembly,
            transcriptome=self.transcriptome,
            orders=self.orders, **kwargs
        )

    # def clean(self):
    #     self.system("rm -r {tmp_dir}")

    def run(self, **kwargs):
        """
        电脑内存：250G

        电脑逻辑CPU个数：72

        目前:
            kraken2每个样本需要内存：数据库大小（70G）

            Metaphlan2每个样本需要内存：数据库大小（1.5G）

            FMAP每个样本需要内存：数据库大小（uniref90, 2.5G; ARDB, 100M）× threads 个数
        """
        """
        self.format_raw(processors=3)
        self.run_fastqc(self.raw_list,
                        pass_if_exists=["{fastqc_out}/{sample}/done"],
                        clean_before=["{fastqc_out}/{sample}"],
                        first_check=1, **self.alloc_src("fastqc"))

        if self.transcriptome:
            self.run_sortmerna(self.raw_list,
                               pass_if_exists=["{sortmerna_out}/{sample}/done"],
                               clean_before=["{sortmerna_out}/{sample}"],
                               first_check=1, **self.alloc_src("sortmerna"))
        self.run_kneaddata(
            self.sortmerna_list if self.transcriptome else self.raw_list,
            pass_if_exists=["{kneaddata_out}/{sample}/done"],
            clean_before=["{kneaddata_out}/{sample}"],
            gzip=self.zip_kneaddata,
            remove_input=True if self.transcriptome else False,
            **self.alloc_src("kneaddata"))

        self.generate_qc_report(processors=3)

        if self.base_on_assembly or self.if_binning != "no":
            self.assembly(self.clean_paired_list,
                          pass_if_exists=["{assembly_out}/{sample}/done"],
                          clean_before=["{assembly_out}/{sample}"],
                          **self.alloc_src("megahit"))
            if self.mix_asem:
                self.build_assembled(threads=self.this_threads)
                self.unassembled(self.clean_paired_list,
                                 pass_if_exists=["{unassembled_out}/{sample}/done"],
                                 clean_before=["{unassembled_out}/{sample}"],
                                 **self.alloc_src("bowtie2"))
                self.mix_assembly(threads=self.this_threads)
            self.cat_assembly()
        """
        if self.if_binning != "no":
            self.binning_pipeline()
        if self.if_binning != "only":
            self.run_kraken2(self.clean_paired_list,
                             pass_if_exists=["{kraken2_out}/{sample}_done"],
                             clean_before=[
                                 "{kraken2_out}/{sample}.report",
                                 "{kraken2_out}/{sample}.report.bracken",
                                 "{kraken2_out}/{sample}.txt",
                                 "{kraken2_out}/{sample}_bracken.report"
                             ],
                             **self.alloc_src("kraken2"))
            self.run_bracken()
            if self.base_on_assembly:
                self.run_prodigal()
                self.run_cdhit()
                self.plot_orf()
                self.run_quast()
                self.create_gene_db(threads=self.this_threads)
                self.quant_gene(self.clean_paired_list,
                                pass_if_exists=["{salmon_out}/{sample}.quant/done"],
                                clean_before=["{salmon_out}/{sample}.quant"],
                                first_check=10, **self.alloc_src("salmon"))
                self.run_diamond_card()
                self.run_diamond_cazy()
                self.run_emapper()

                for num, db in enumerate([
                        "AMR",
                        "CAZY",
                        "EGGNOG",
                        "KO",
                        "GO"
                ]):
                    if num == 0:
                        # first run will require lots of IO, less processors are appropriate
                        self.group_gene(db, processors=3)
                    else:
                        self.group_gene(db, processors=10)
            else:
                self.run_humann2(
                    self.clean_r1_list,
                    pass_if_exists=["{humann2_out}/{sample}_done"],
                    clean_before=["{humann2_out}/{sample}"],
                    **self.alloc_src("humann2"))
                self.join_humann()

                self.fmap_wrapper(self.clean_r1_list,
                                  run_type="AMR", **self.alloc_src("fmap"))

            self.visualize(**kwargs)
