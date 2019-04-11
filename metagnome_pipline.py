import os
import time
from pyutils.tools import split_list
from pyutils.tools import parse_premap
import json
import re


class MetagenomePipline(object):
    """
    arguments:
        run_size: control the max number of jobs submitted to sge each time
        raw_fqs_dir: directory where the raw fastq file were stored
        sample_regex: regular expression to match sample id (contained by brackets)
        forward_regex: regular expression to match forward fastq files
        reverse_regex: regular expression to match reverse fastq files
        out_results_dir: where to store the results
    sample usage:
    from metagnome_pipline import MetagenomePipline
    m =MetagenomePipline('/home/cheng/Projects/rll_testdir/1.rawdata/','/home/cheng/Projects/rll_testdir/mapping_file.txt',out_results_dir="/home/cheng/Projects/rll_testdir/")
    m.fmap_wrapper("AMR",8)
    """

    def __init__(self, raw_fqs_dir, pre_mapping_file, run_size=4, host_type="hg38", sample_regex="(.+)_.*_[12]\.fq\.gz", forward_regex="_1\.fq\.gz$", reverse_regex="_2\.fq\.gz$", out_results_dir='/home/cheng/Projects/rll_testdir/test/'):
        self.out_dir = os.path.abspath(out_results_dir) + '/'
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        self._base_dir = os.path.dirname(__file__) + '/'
        with open(self._base_dir + "pipconfig/path.conf") as f:
            self.path = json.load(f)
        self.parsed_map = parse_premap(raw_fqs_dir, pre_mapping_file, forward_regex, reverse_regex, sample_regex)
        self.fq_info = self.parsed_map['fastq'].values
        self._init_outdir_()
        self.parsed_map['map'].to_csv(self.out_dir + 'mapping_file.txt', sep='\t', index=False)
        raw_pattern = self.out_dir + "Raw_fastq/{}_{}.fq.gz"
        trimmed_pattern = self.out_dir + "Primer_trimmed/{}_{}.fq.gz"
        filtered_pattern = self.out_dir + "Filtered/{}_{}.good.fastq.gz"
        de_host_pattern = self.out_dir + "Host_subtracted/bowtie/%s/{}_{}.%s.unmapped.fastq.gz" % (host_type, host_type)
        kracken2_reports_pattern = self.out_dir + "Kraken2/{}_{}.report"
        self.merged_pe_pattern = self.out_dir + \
            "Host_subtracted/bowtie/%s/{}_R1.%s.unmapped.{}.fastq.gz" % (host_type, host_type)
        self.raw_list = self.map_list(raw_pattern, run_size)
        self.trimmed_list = self.map_list(trimmed_pattern, run_size)
        self.filtered_list = self.map_list(filtered_pattern, run_size)
        self.de_host_r1_list = self.map_list(de_host_pattern, run_size, only_r1=True)
        self.merged_pe_r1_list = self.map_list(self.merged_pe_pattern, run_size, only_r1=True)
        self.kracken2_reports_list = self.map_list(kracken2_reports_pattern, only_r1=True)
        global running_list
        running_list = self.out_dir + '.running_list'

    def _init_outdir_(self):
        os.system('perl ' + self.path['cii_home'] + 'create_dir_structure.pl ' + self.out_dir)
        for fq_path, new_id, direction in self.fq_info:
            os.system(("ln -s {} " + self.out_dir + "Raw_fastq/{}_{}.fq.gz").format(fq_path, new_id, direction))

    def map_list(self, pattern=False, each=False, only_r1=False):
        out = []
        for fq_path, new_id, direction in self.fq_info:
            ele = pattern.format(new_id, direction) if pattern else [new_id, direction]
            if direction == "R2":
                if not only_r1:
                    out.append(ele)
            else:
                out.append(ele)
        out.sort()
        return split_list(out, each) if each else out

    def synchronize(func):

        def wfunc(*args, fq_list, first_check=10, **kwargs):
            print("######################Running " + str(func))
            global running_list
            # running_list = MetagenomePipline.out_dir + '.running_list'
            start_time = time.time()
            for fq in fq_list:
                with open(running_list, 'w') as f:
                    f.write('\n'.join(fq))
                    f.flush()
                func(*args, fq_list=running_list, **kwargs)
                time.sleep(first_check * 60)
                while True:
                    if not list(os.popen("qstat")):
                        break
                    time.sleep(60)
            end_time = time.time()
            time_used = (end_time - start_time) / 60
            print("######################" + str(func) + " done; time used: {} min".format(time_used))
        return wfunc

    def synchronize2(func):

        def wfunc(*args, fq_list, first_check=10, **kwargs):
            print("######################Running " + str(func))
            # running_list = MetagenomePipline.out_dir + '.running_list'
            start_time = time.time()
            for fq in fq_list:
                func(*args, fq_list=fq, **kwargs)
                time.sleep(first_check * 60)
                while True:
                    if not list(os.popen("qstat")):
                        break
                    time.sleep(60)
            end_time = time.time()
            time_used = (end_time - start_time) / 60
            print("######################" + str(func) + " done; time used: {} min".format(time_used))
        return wfunc

    def homized_cmd(self, cmd, home="/home/cheng/pipelines/CII_meta/"):
        return "cd {}&&".format(home) + cmd

    @synchronize
    def run_fastqc(self, fq_list: list, processor=2):
        """
        fq_list: 2 dimension list of fastq files
        processor: number of processor
        """
        os.system(self.homized_cmd("perl fastqc_wrapper.pl {} {}".format(fq_list, str(processor))))

    @synchronize
    def run_trim(self, fq_list: list):
        os.system(self.homized_cmd("perl cutadapt_wrapper.pl -f {} -b adaptor_Illumina.list".format(fq_list)))

    @synchronize
    def run_filter(self, fq_list: list):
        os.system(self.homized_cmd("perl prinseq_wrapper.pl --aim Filter --single-end-fq-list {}".format(fq_list)))

    @synchronize
    def run_de_host(self, fq_list: list, processor=2):
        os.system("perl " + self.path['de_host_path'] +
                  ' --aim Host --config human --single-end-fq-list {} --threads {} -m N'.format(fq_list, str(processor)))

    @synchronize
    def run_merge_se_to_pe(self, fq_list: list):
        os.system("perl " + self.path['merge_se_path'] + ' --paired-end-fq-list {}'.format(fq_list))

    @synchronize
    def run_kraken2(self, fq_list: list, processor=2):
        os.system(self.homized_cmd("perl run_kraken2.pl {} {} {} PE".format(
            fq_list, str(processor), self.path['kraken2_database'])))

    def run_bracken(self):
        global running_list
        bracken_list = [l + '.bracken' for l in self.kracken2_reports_list]
        for report in self.kracken2_reports_list:
            os.system("python2 {} -i {} -k {} -l S -o {}.bracken".format(
                self.path['bracken_path'], report, self.path['bracken_database'], report))
        with open(running_list, 'w') as f:
            f.write('\n'.join(bracken_list))
        os.system(
            "perl {}/Braken_to_OTUtable.pl {} {}".format(self.path['cii_home'], self.path['ncbi_taxaID_path'], running_list))

    @synchronize
    def run_metaphlan2(self, fq_list: list, processor=2):
        os.system(self.homized_cmd("perl run_metaphlan2.pl {} PE {}".format(
            fq_list, str(processor))))

    @synchronize
    def run_humann(self, fq_list: list):
        os.system(self.homized_cmd('perl run_humann.pl {} SE'.format(fq_list)))

    @synchronize2
    def run_fmap(self, fq_list, database="protein_fasta_protein_homolog_model", processor=4, out_dir='FMAP'):
        """
        arguments:
            database: database at FMAP_data
        """
        with open(self.path['fmap_home'] + '/FMAP_data/database', 'w') as f:
            f.write(database)
        out = self.out_dir + out_dir + '/'
        if not os.path.exists(out):
            os.makedirs(out)
        for fq in fq_list:
            cmd = "echo 'perl {}/FMAP_mapping.pl -p {} {} > {}.mapping.txt' | qsub -V -N {} -cwd -l h_vmem=24G -o {} -e {} -pe smp 4".format(
                self.path['fmap_home'], str(processor), fq, out + os.path.basename(fq), os.path.basename(fq), out, out)
            print("submit:\n {}\n\n".format(cmd))
            os.system(cmd)

    def quantify_fmap(self, out_dir='FMAP', all_name="all.txt", print_definition=False):
        out = self.out_dir + out_dir + '/'

        p = "" if print_definition else "-n"
        file_pattern = out + re.search("[^/]+$", self.merged_pe_pattern).group()
        info_list = self.map_list(only_r1=True)
        for new_id, direction in info_list:
            cmd = "perl {}/FMAP_quantification.pl {}.mapping.txt > {}.abundance.txt".format(
                self.path['fmap_home'], file_pattern.format(new_id, direction), out + new_id)
            print(cmd)
            os.system(cmd)
        args = ["{}={}.abundance.txt".format(new_id, out + new_id) for new_id, direction in info_list]
        os.system("perl {}/FMAP_table.pl {} {} > {}".format(
            self.path['fmap_home'], p, " ".join(args), out + all_name))

    def fmap_wrapper(self, run_type="KEGG", processor=4):
        if run_type == "KEGG":
            self.run_fmap(fq_list=self.merged_pe_r1_list,
                          database="orthology_uniref90_2_2157_4751.20180725040837", processor=processor, out_dir="FMAP")
            self.quantify_fmap(out_dir="FMAP", all_name="All.Function.abundance.KeepID.KO.txt")
        elif run_type == "AMR":
            self.run_fmap(fq_list=self.merged_pe_r1_list, database="protein_fasta_protein_homolog_model_cleaned",
                          processor=processor, out_dir="AMR")
            self.quantify_fmap(out_dir="AMR", all_name="All.AMR.abundance.txt")
        elif run_type == "ARDB":
            self.run_fmap(fq_list=self.merged_pe_r1_list, database="ARDB.20180725064354",
                          processor=processor, out_dir="ARDB")
            self.quantify_fmap(out_dir="ARDB", all_name="All.ARDB.abundance.txt", print_definition=True)

    def run_pipline(self, processor=3):
        self.run_fastqc(fq_list=self.raw_list, processor=processor)
        self.run_trim(fq_list=self.raw_list)
        self.run_filter(fq_list=self.trimmed_list)
        self.run_fastqc(fq_list=self.filtered_list, processor=processor)
        self.run_de_host(fq_list=self.filtered_list, processor=processor)
        self.run_merge_se_to_pe(fq_list=self.de_host_r1_list)
        self.run_kraken2(fq_list=self.merged_pe_r1_list, processor=processor)
        self.run_bracken()
        self.run_metaphlan2(fq_list=self.merged_pe_r1_list, processor=processor)
        self.run_humann(fq_list=self.merged_pe_r1_list)
        self.run_fmap(fq_list=self.merged_pe_r1_list, processor=processor)
        self.fmap_wrapper(run_type="KEGG", processor=processor * 2)
        self.fmap_wrapper(run_type="AMR", processor=processor * 2)
