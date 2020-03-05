import os
import time
from pyutils.tools import split_list, parse_premap
import re
from visualizeAll import VisualizeAll
import pandas as pd
from systemMixin import SystemMixin
import numpy as np
from pipconfig import settings
from lxml import html
from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor
# import uuid


def wait_sge(first_check=1):
    time.sleep(first_check * 60)
    while True:
        if not list(os.popen("qstat")):
            break
        time.sleep(60)


def sge_decorator(func):
    """
    对于被修饰的函数：fq_list参数值应是单个fq文件的路径（或路径对）
    对于最终函数：fq_list参数是多个fq文件路径的二(至少)维数组(list)
    """

    def wfunc(self, fq_list, first_check=10, callback=False, callback_kwargs={}, max_workers=10, **kwargs):
        fq_list = split_list(fq_list, max_workers)
        print("######################Running " + str(func))
        start_time = time.time()
        for fqs in fq_list:
            for fq in fqs:
                func(self, fq_list=fq, **kwargs)
            if callback:
                wait_sge(first_check)
                print("This run done! checking callback...")
                callback(**callback_kwargs)
        if not callback:
            wait_sge(first_check)
        end_time = time.time()
        time_used = (end_time - start_time) / 60
        print("######################" + str(func) +
              " done; time used: {} min".format(time_used))
    return wfunc


def pool_decorator(func):
    """
    对于被修饰的函数：fq_list参数值应是单个fq文件的路径（或路径对）
    对于最终函数：fq_list参数是多个fq文件路径的二(至少)维数组(list)
    """

    def wfunc(self, fq_list, first_check=1, callback=False, callback_kwargs={}, max_workers=10, **kwargs):
        fq_list = split_list(fq_list, max_workers)
        print("######################Running " + str(func))
        start_time = time.time()
        if not callback:
            executor = ThreadPoolExecutor(max_workers=max_workers)
        for fqs in fq_list:
            if callback:
                executor = ThreadPoolExecutor(max_workers=max_workers)
            for fq in fqs:
                executor.submit(func, self, fq_list=fq, **kwargs)
            if callback:
                executor.shutdown(True)
                print("This run done! checking callback...")
                callback(**callback_kwargs)
        if not callback:
            executor.shutdown(True)
        end_time = time.time()
        time_used = (end_time - start_time) / 60
        print("######################" + str(func) +
              " done; time used: {} min".format(time_used))
    return wfunc


synchronize = sge_decorator if settings.use_sge else pool_decorator


class RawDataNotPairedError(Exception):
    pass


class PathError(Exception):
    pass


class MetagenomePipline(SystemMixin):
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

    def __init__(self, raw_fqs_dir, pre_mapping_file, categories=False, host_db="/home/cheng/Databases/hg38/hg38", sample_regex="(.+)_.*_[12]\.fq\.?g?z?", forward_regex="_1\.fq\.?g?z?$", reverse_regex="_2\.fq\.?g?z?$", out_dir='./', base_on_assembly=False, exclude='none', filter_species=False, zip_kneaddata=True):
        self.context = dict()

        self.set_path(force=True,
                      out_dir=out_dir,
                      raw_dir=out_dir + '/raw_dir',
                      salmon_out=out_dir + '/salmon_out',
                      kneaddata_out=out_dir + '/kneaddata_out',
                      kraken2_out=out_dir + '/Kraken2',
                      amr_out=out_dir + '/AMR',
                      humann2_out=out_dir + '/Metagenome/Humann/',
                      fmap_out=out_dir + '/FMAP/',
                      assembly_out=out_dir + "/Assembly_out",
                      metaphlan_out=out_dir + "/Metagenome/Metaphlan",
                      report_out=out_dir + "/Report",
                      )

        self.set_attr(base_on_assembly=base_on_assembly,
                      # assembly_out=out_dir + "/Assembly_out",
                      mapping_file=self.out_dir + '/mapping_file.txt',
                      host_db=' -db '.join(host_db.split(',')),
                      exclude=exclude,
                      base_dir=os.path.dirname(__file__) + '/',
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
        self.categories = categories if categories else ','.join(
            [g for g in self.parsed_map['map'].columns if not g.find('Category') == -1])
        print("The detected categories are: \n\n{}\n".format(self.categories))

        self.raw_pattern = "%s/{sample_id}_{direction}.fq" % (self.raw_dir)
        self.raw_list = self.paired_data(self.raw_pattern)
        self.unzip_clean_paired_pattern = "%s/{sample_id}_R1_kneaddata_paired_{direction_num}.fastq" % (
            self.kneaddata_out)
        self.clean_paired_pattern = self.unzip_clean_paired_pattern + \
            '.gz' if zip_kneaddata else self.unzip_clean_paired_pattern
        self.clean_r1_list = self.map_list(
            self.clean_paired_pattern, use_direction='R1')
        self.clean_paired_list = self.paired_data(self.clean_paired_pattern)
        self.unzip_clean_paired_list = self.paired_data(
            self.unzip_clean_paired_pattern)
        self.kracken2_reports_pattern = "%s/{sample_id}.report" % (
            self.kraken2_out)
        self.unassembled_pattern = "%s/{sample_id}/{sample_id}_R{direction_num}_unassembled_merged.fq.gz" % (
            self.assembly_out)
        self.kracken2_reports_list = self.map_list(
            self.kracken2_reports_pattern, use_direction='R1')
        self.running_list = self.out_dir + '.running_list'
        self.escape_sge = not settings.use_sge
        self.threads_single = int(self.threads / self.hosts)
        self.filter_species = filter_species
        self.zip_kneaddata = zip_kneaddata

    def format_raw(self, processors=3):
        executor = ThreadPoolExecutor(max_workers=processors)
        for fq_path, new_id, direction in self.fq_info:
            if fq_path.endswith('.gz'):
                executor.submit(self.unzip, fq_path)
        executor.shutdown(True)
        for fq_path, new_id, direction in self.fq_info:
            formated_fq = self.raw_pattern.format(
                sample_id=new_id, direction=direction)
            if not os.path.exists(formated_fq):
                fq_path = fq_path.rstrip('.gz')
                os.system("ln -s {} {}".format(fq_path, formated_fq))

    def unzip(self, fq_path):
        print("Unzip file {}...".format(fq_path))
        os.system("gunzip {}".format(fq_path))

    def map_list(self, pattern=False, each=False, use_direction="both"):
        out = []
        for fq_path, new_id, direction in self.fq_info:
            ele = pattern.format(sample_id=new_id, direction=direction, direction_num=re.search(
                "\d$", direction).group()) if pattern else [new_id, direction]
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
        else:
            return {
                "r1": fq_list,
                "sample": os.path.basename(fq_list).split('_')[0]
            }

    @synchronize
    def run_kneaddata(self, fq_list, threads=10, mem=20):
        self.system("""
echo '{kneaddata_path} -i {r1} -i {r2} -o {kneaddata_out} -v -db {host_db} \
 --remove-intermediate-output  -t {threads} --run-fastqc-start --run-fastqc-end \
 --trimmomatic {trimmomatic_home} --trimmomatic-options "ILLUMINACLIP:{adapters_path}:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50" --max-memory {mem}g \
 --bowtie2 {bowtie2_home} \
 --fastqc  {fastqc_home}' \
 | qsub -l h_vmem={mem}G -pe {sge_pe} {threads} -q {sge_queue} -V -N {sample} -o {kneaddata_out} -e {kneaddata_out}
            """, **self.parse_fq_list(fq_list), kneaddata_out=self.kneaddata_out, threads=threads, mem=mem, escape_sge=self.escape_sge)

    def kneaddata_callback(self):
        self.mirror_move(self.tmp_dir, self.kneaddata_out)

    @pool_decorator
    def gzip_kneaddata(self, fq_list):
        self.system("gzip {r1}&&gzip {r2}", **self.parse_fq_list(fq_list))

    @synchronize
    def run_kraken2(self, fq_list, threads=2, mem=70):
        self.system("""
echo '{kraken2_path} --db {kraken2_database} --threads {threads} --confidence 0.2 \
 --report {kraken2_out}/{sample}.report --paired {r1} {r2} --output -' \
  | qsub -l h_vmem={mem}G -pe {sge_pe} {threads} -q {sge_queue} -V -N {sample} -o {kraken2_out} -e {kraken2_out}
            """, **self.parse_fq_list(fq_list), threads=threads, mem=mem, escape_sge=self.escape_sge)

    def run_bracken(self):
        bracken_list = [l + '.bracken' for l in self.kracken2_reports_list]
        for report in self.kracken2_reports_list:
            self.system(
                "python2 {bracken_path} -i {report} -k {bracken_database} -l S -o {report}.bracken", report=report)
        with open(self.running_list, 'w') as f:
            f.write('\n'.join(bracken_list))
        os.system(
            self.homized_cmd("perl Braken_to_OTUtable.pl {} {}".format(self.ncbi_taxaID_path, self.running_list)))

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
  | qsub -l h_vmem={mem}G -pe {sge_pe} {threads} -q {sge_queue} -V -N {sample} -o {metaphlan_out} -e {metaphlan_out}
            """, **self.parse_fq_list(fq_list), threads=threads, mem=mem, escape_sge=self.escape_sge)

    def join_metaphlan(self):
        tb_name = self.metaphlan_out + '/All.Metaphlan2.profile.txt'
        self.system(
            "python2 {base_dir}/merge_metaphlan_tables.py {metaphlan_out}/*profile.txt > {tb_name}", tb_name=tb_name)
        self.clean_header(tb_name, pattern='\.metaphlan\.profile$', skip=[0])

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
mv {humann2_out}/{sample}/*/*bugs_list.tsv {humann2_out}/{sample}/*genefamilies.tsv {humann2_out}/{sample}/*pathabundance.tsv {humann2_out}/&& \
rm -r {humann2_out}/{sample}' \
  | qsub -l h_vmem={mem}G -pe {sge_pe} {threads} -q {sge_queue} -V -N {sample} -o {humann2_out} -e {humann2_out}
            """, **self.parse_fq_list(fq_list), threads=threads, mem=mem, escape_sge=self.escape_sge)

    def humann2_callback(self):
        self.system(
            "mv {humann2_out}/*/*/*bugs_list.tsv {humann2_out}/")
        for e in os.listdir(self.humann2_out):
            path = os.path.join(self.humann2_out, e)
            if os.path.isdir(path):
                os.system("rm -r {}".format(path))

    def join_humann(self):
        self.system('''
{humann2_home}/humann2_join_tables -i {humann2_out}/ -o {humann2_out}/RPK.All.UniRef90.genefamilies.tsv --file_name genefamilies.tsv
{humann2_home}/humann2_join_tables -i {humann2_out}/ -o {humann2_out}/RPK.All.Metacyc.pathabundance.tsv --file_name pathabundance.tsv
{humann2_home}/humann2_join_tables -i {humann2_out}/ -o {metaphlan_out}/All.Metaphlan2.profile.txt --file_name bugs_list.tsv''')
        self.clean_header(self.humann2_out +
                          "/RPK.All.UniRef90.genefamilies.tsv")
        self.clean_header(self.humann2_out +
                          "/RPK.All.Metacyc.pathabundance.tsv")
        self.clean_header(self.metaphlan_out + "/All.Metaphlan2.profile.txt")

    def set_fmap_db(self, database):
        with open(self.fmap_home + '/FMAP_data/database', 'w') as f:
            f.write(database)

    @synchronize
    def run_fmap(self, fq_list, threads=4, mem=24):
        self.system("""
echo 'perl {fmap_home}/FMAP_mapping.pl -e 0.00001 -p {threads} {r1} > {fmap_out}/{sample}.mapping.txt' | \
 qsub -l h_vmem={mem}G -pe {sge_pe} {threads} -q {sge_queue} -V -N {sample} -o {fmap_out} -e {fmap_out}
            """, **self.parse_fq_list(fq_list), threads=threads, mem=mem, escape_sge=self.escape_sge)

    def quantify_fmap(self, all_name="all.txt", print_definition=False):
        p = "" if print_definition else "-n"
        for new_id in self.new_ids:
            self.system(
                "perl {fmap_home}/FMAP_quantification.pl {fmap_out}/{sample}.mapping.txt > {fmap_out}/{sample}.abundance.txt", sample=new_id)
        args = [
            "{sample}={fmap_out}/{sample}.abundance.txt".format(sample=new_id, fmap_out=self.fmap_out) for new_id in self.new_ids]
        self.system("perl {fmap_home}/FMAP_table.pl {p} {args} > {fmap_out}/{all_name}",
                    p=p, args=" ".join(args), all_name=all_name)

    def fmap_wrapper(self, fq_list, run_type="KEGG", threads=6, mem=10, max_workers=10):
        fmap_db = {
            "KEGG": "orthology_uniref90_2_2157_4751.20190412161853",
            "AMR": "protein_fasta_protein_homolog_model_cleaned",
            "ARDB": "ARDB.20180725064354",
            "MGE": "MGEs_FINAL_99perc_trim"
        }
        self.set_fmap_db(fmap_db[run_type])
        self.run_fmap(fq_list=fq_list, threads=threads,
                      mem=mem, max_workers=max_workers)
        self.quantify_fmap(
            all_name="All.{}.abundance_unstratified.tsv".format(run_type))

    @synchronize
    def assembly(self, fq_list, threads=4, mem=24):
        # 单样品组装
        parsed_fqs = self.parse_fq_list(fq_list)
        sample_out = os.path.join(self.assembly_out, parsed_fqs['sample'])
        if os.path.exists(os.path.join(sample_out, 'all_done')):
            return
        if os.path.exists(sample_out):
            self.system("rm -r {}".format(sample_out))
        self.system("""
echo 'mkdir -p {tmp_dir}/{sample}/bowtie2_db&&{megahit_path} --k-list 21,29,39,59,79,99,119,141 \
 --tmp-dir {tmp_dir}/{sample} -m {mem_p} --mem-flag 1  \
 -1 {r1} -2 {r2} --min-contig-len 500 -t {threads} -o {assembly_out}/{sample} && \
{bowtie2_home}/bowtie2-build --threads {threads} {assembly_out}/{sample}/final.contigs.fa {tmp_dir}/{sample}/bowtie2_db/{sample} && \
{bowtie2_home}/bowtie2 --threads {threads} -x {tmp_dir}/{sample}/bowtie2_db/{sample}  -U {r1} --end-to-end --sensitive -S {sam_out} \
 --un {assembly_out}/{sample}/{sample}_R1_unassembled.fastq && \
{bowtie2_home}/bowtie2 --threads {threads} -x {tmp_dir}/{sample}/bowtie2_db/{sample}  -U {r2} --end-to-end --sensitive -S {sam_out} \
 --un {assembly_out}/{sample}/{sample}_R2_unassembled.fastq && \
rm -r {assembly_out}/{sample}/intermediate_contigs {tmp_dir}/{sample} && \
{base_dir}/merge_se_to_pe.py -1 {assembly_out}/{sample}/{sample}_R1_unassembled.fastq -2 {assembly_out}/{sample}/{sample}_R2_unassembled.fastq \
-o {assembly_out}/{sample}/{sample}_R%s_unassembled_merged.fq  --removeinput && \
gzip {assembly_out}/{sample}/{sample}_R*_unassembled_merged.fq && \
echo {sample}_done > {assembly_out}/{sample}/all_done' | \
 qsub -l h_vmem={mem}G -pe {sge_pe} {threads} -q {sge_queue} -V -N {sample} -o {assembly_out} -e {assembly_out}
            """, **parsed_fqs, threads=threads, mem_p=mem * 1000000000, mem=mem,
                    escape_sge=self.escape_sge,
                    sam_out=os.devnull)

    def sum_assembly(self, threads=50):
        # 混合组装
        self.system(
            "{megahit_path} --continue  --kmin-1pass --presets meta-large -m 0.94 --mem-flag 0 \
             -1 {r1_list} -2 {r2_list} --min-contig-len 500 -t {threads} -o {assembly_out}/mixed_assembly",
            r1_list=','.join(self.map_list(
                self.unassembled_pattern, use_direction="R1")),
            r2_list=','.join(self.map_list(
                self.unassembled_pattern, use_direction="R2")),
            threads=threads
        )

        # 基因预测，去冗余
        self.system(
            """
cd {assembly_out}
cat */final.contigs.fa > final.contigs.fa
{prodigal_path} -i final.contigs.fa -a final.contigs.fa.faa -d final.contigs.fa.fna  -f gff -p meta -o final.contigs.fa.gff
{cdhit_path} -i final.contigs.fa.fna -d 0 -M 0 -o NR.nucleotide.fa -T 0
grep '>' NR.nucleotide.fa > NR.nucleotide.fa.header
{acmd}
{R_path} {base_dir}/ORF_header_summary.R -i {assembly_out}""",
            acmd=self.homized_cmd(
                "{perl_path} ORF_generate_input_stats_file.pl {assembly_out}/NR.nucleotide.fa.header".format(**self.context)),
        )

    def run_quast(self):
        self.system(
            "python2 {quast_path} -o {assembly_out}/quast_results/  {assembly_out}/final.contigs.fa")

    def create_gene_db(self, threads=7):
        self.system('''
{transeq_path} -sequence {assembly_out}/NR.nucleotide.fa -outseq {assembly_out}/NR.protein.fa -trim Y
sed -i 's/_1 / /' {assembly_out}/NR.protein.fa
{salmon_path} index -t {assembly_out}/NR.nucleotide.fa -p {threads} -k 31 -i {assembly_out}/salmon_index''', threads=threads)

    @synchronize
    def quant_gene(self, fq_list, threads=7, mem=80):
        self.system("""
echo '{salmon_path} quant --validateMappings -i {assembly_out}/salmon_index -l A -p {threads} --meta -1 {r1} -2 {r2} -o {salmon_out}/{sample}.quant' \
| qsub -l h_vmem={mem}G -pe {sge_pe} {threads} -q {sge_queue} -V -N {sample} -o {salmon_out} -e {salmon_out}
            """, **self.parse_fq_list(fq_list), threads=threads, mem=mem, escape_sge=self.escape_sge)

    def join_gene(self):
        self.system(
            "{salmon_path} quantmerge --quants {salmon_out}/*.quant -o {salmon_out}/All.genes.abundance.txt")
        self.clean_header(
            self.out_dir + "salmon_out/All.genes.abundance.txt", ".quant$")

    def diamond_gene(self, database, out_file, threads=66):
        self.system("""
{diamond_home}/diamond blastp --db {database} --query {assembly_out}/NR.protein.fa --outfmt 6 --threads {threads} \
 -e 0.00001 --id 80 --top 3 --block-size 200 --index-chunks 1 --quiet --out {salmon_out}/{out_file}
            """, out_file=out_file, threads=threads, database=database)

    def map_gene(self, threads=70):
        self.system('''
{emapper_path} --no_file_comments -m diamond --seed_ortholog_evalue 0.00001 \
  --data_dir {emapper_database} --cpu {threads} \
  -i {assembly_out}/NR.protein.fa -o {salmon_out}/genes --usemem --override 

# {emapper_path} --annotate_hits_table {salmon_out}/genes.emapper.seed_orthologs --no_file_comments \
  -o {salmon_out}/genes --cpu {threads} --data_dir {emapper_database} --usemem --override
sed -i '1 i Name\teggNOG\tEvalue\tScore\tGeneName\tGO\tKO\tBiGG\tTax\tOG\tBestOG\tCOG\tAnnotation' {salmon_out}/genes.emapper.annotations
            ''', threads=threads)
        self.diamond_gene(self.cazy_database, 'genes_cazy.f6', threads)
        self.diamond_gene(self.fmap_home + '/FMAP_data/protein_fasta_protein_homolog_model_cleaned.dmnd',
                          'genes_card.f6', threads)

    def alloc_src(self, proc, threads=False, sam_num=False):
        memery_needs = self.memery_needs[proc]
        total_sample = len(self.new_ids)
        if not sam_num:
            sample_number = self.memery // memery_needs
            sample_number = sample_number if sample_number < settings.max_workers[
                proc] else settings.max_workers[proc]
            sample_number = sample_number if sample_number > 1 else 1
            sample_number = sample_number if sample_number < total_sample else total_sample
            runs = np.ceil(total_sample / sample_number)
            sample_number = int(np.ceil(total_sample / runs))
        else:
            sample_number = sam_num
        threads = threads if threads else (self.threads // sample_number)
        print("For {}:\nSample number per run is: {}\n Threads number per sample is {}\n Memery size per sample is {}".format(
            proc, sample_number, threads, memery_needs))
        return {
            "max_workers": sample_number,
            "threads": threads,
            "mem": memery_needs

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
        qr = list(os.popen(
            "{base_dir}/sim_fqstat {r1}".format(base_dir=self.base_dir, r1=read_file)))
        bases, q20_cnt, q30_cnt = [int(i) for i in qr[0].split()]
        return np.round(q20_cnt / bases * 100, 2), np.round(q30_cnt / bases * 100, 2)

    def get_qc_stats(self, sample_id):
        target_dir = os.path.join(self.kneaddata_out, "fastqc")
        results = OrderedDict([
            ("Sample ID", sample_id),
            ("InsertSize(bp)", "350"),
            ("SeqStrategy", "(150:150)"),
        ])
        raw_fqc = self.find_file(
            target_dir, "^reformatted_identifier.+_{}_R1_fastqc.html$".format(sample_id))
        print("The raw fastqc file of {} is {}".format(sample_id, raw_fqc))
        raw_num_reads, raw_gc_content = self.search_fqc(raw_fqc)
        clean_fqc = self.find_file(
            target_dir, "^{}_R1_kneaddata_paired_1_fastqc.html$".format(sample_id))
        print("The clean fastqc file of {} is {}".format(sample_id, clean_fqc))
        clean_num_reads, clean_gc_content = self.search_fqc(clean_fqc)
        raw_q20, raw_q30 = self.q20_q30(
            self.raw_pattern.format(sample_id=sample_id, direction="R1"))
        clean_q20, clean_q30 = self.q20_q30(
            self.clean_paired_pattern.format(sample_id=sample_id, direction_num=1))

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

    def visualize(self):
        VisualizeAll(self.mapping_file, self.categories, out_dir=self.out_dir, filter_species=self.filter_species).visualize(
            self.exclude, self.base_on_assembly)

    # def clean(self):
    #     self.system("rm -r {tmp_dir}")

    def run(self):
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
        self.run_kneaddata(
            self.raw_list, callback=False, **self.alloc_src("kneaddata"))
        self.generate_qc_report(processors=3)
        """
        if self.zip_kneaddata:
            self.gzip_kneaddata(self.unzip_clean_paired_list, max_workers=6)
        """
        self.run_kraken2(self.clean_paired_list, **self.alloc_src("kraken2"))
        self.run_bracken()
        """
        if self.base_on_assembly:
            """
            self.assembly(self.clean_paired_list, **self.alloc_src("megahit"))
            self.sum_assembly(threads=self.threads_single)
            self.run_quast()
            self.create_gene_db(threads=self.threads_single)
            """
            self.quant_gene(self.clean_paired_list,
                            first_check=5, **self.alloc_src("salmon"))
            self.join_gene()
            self.map_gene(threads=self.threads_single)
        else:

            self.run_humann2(
                self.clean_r1_list, callback=False, **self.alloc_src("humann2"))
            self.join_humann()

            self.fmap_wrapper(self.clean_r1_list,
                              run_type="AMR", **self.alloc_src("fmap"))

        self.visualize()
