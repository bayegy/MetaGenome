import os
import time
from pyutils.tools import split_list
from pyutils.tools import parse_premap
import json


class MetagenomePipline(object):
    """
    out_dir: the out results directory
    """

    def __init__(self, raw_fqs_dir, pre_mapping_file, run_size=10, sample_regex="(.+)_.*_[12]\.fq\.gz", forward_regex="_1\.fq\.gz$", reverse_regex="_2\.fq\.gz$", out_results_dir='./'):
        self._base_dir = os.path.dirname(__file__) + '/'
        self.out_dir = os.path.abspath(out_results_dir) + '/'
        self.parsed_map = parse_premap(raw_fqs_dir, pre_mapping_file, forward_regex, reverse_regex, sample_regex)
        self.parsed_map['map'].to_csv(self.out_dir + 'mapping_file.txt', sep='\t', index=False)
        self.fq_info = self.parsed_map['fastq'].values
        raw_pattern = self.out_dir + "Raw_fastq/{}_{}.fq.gz"
        trimmed_pattern = self.out_dir + "Primer_trimmed/{}_{}.fq.gz"
        filtered_pattern = self.out_dir + "Filtered/{}_{}.good.fastq.gz"
        de_host_pattern = self.out_dir + "Host_subtracted/bowtie/hg38/{}_{}.hg38.unmapped.fastq.gz"
        merged_pe_pattern = self.out_dir + "Host_subtracted/bowtie/hg38/{}_R1.hg38.unmapped.{}.fastq.gz"
        self.raw_list = self.map_list(self.fq_info, raw_pattern, run_size)
        self.trimmed_list = self.map_list(self.fq_info, trimmed_pattern, run_size)
        self.filtered_list = self.map_list(self.fq_info, filtered_pattern, run_size)
        self.de_host_r1_list = self.map_list(self.fq_info, de_host_pattern, run_size, only_r1=True)
        self.merged_pe_r1_list = self.map_list(self.fq_info, merged_pe_pattern, run_size, only_r1=True)
        with open(self._base_dir + "pipconfig/path.conf") as f:
            self.path = json.load(f)
        self.running_list = self.out_dir + '.running_list'

    def map_list(self, matched_fq_info, pattern, each, only_r1=False):
        out = []
        for fq_path, new_id, direction in matched_fq_info:
            if direction == "R2":
                if not only_r1:
                    out.append(pattern.format(new_id, direction))
            else:
                out.append(pattern.format(new_id, direction))
        out.sort()
        return split_list(out, each)

    @staticmethod
    def synchronize(func):

        def wfunc(*agrs, fq_list, **kwargs):
            print("######################Running " + str(func))
            running_list = MetagenomePipline.out_dir + '.running_list'
            start_time = time.time()
            for fq in fq_list:
                with open(running_list, 'w') as f:
                    f.write('\n'.join(fq))
                    f.flush()
                func(*agrs, fq_list=running_list, **kwargs)
                time.sleep(10 * 60)
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
    def run_fastq(self, fq_list: list, processor=1):
        """
        fq_list: 2 dimension list of fastq files
        processor: number of processor
        """
        os.system(self.homized_cmd("perl fastqc_wrapper.pl {} {}".format(fq_list, str(processor))))

    @synchronize
    def run_trim(self, fq_list: list):
        os.system(self.homized_cmd("perl cutadapt_wrapper.pl -f {} -b adaptor_Illumina.list".format(fq_list)))
