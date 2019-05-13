import os
import time
from pyutils.tools import split_list, parse_premap
import json
import re
from visualizeAll import VisualizeAll
from mapInfo import MapInfo
import pandas as pd


class MetagenomePipline(object):
    """
    arguments:
        pre_mapping_file: The first column of pre_mapping_file should be smaple id in raw fastq files, the last column of pre_mapping_file should be new id (or the same) of samples.

        categories: Categories seprated by ',' , optional, if not passed, the categories names should have the pattern of 'Group.*'

        host_type: Will control the de_host step.

        run_size: Control the max number of jobs submitted to sge each time

        raw_fqs_dir: Directory where the raw fastq file were stored

        sample_regex: Regular expression to match sample id (contained by brackets)

        forward_regex: Regular expression to match forward fastq files

        reverse_regex: Regular expression to match reverse fastq files

        out_dir: Where to store the results

    sample usage:

    from metagnomePipline import MetagenomePipline
    m =MetagenomePipline('/home/cheng/Projects/rll_testdir/1.rawdata/','/home/cheng/Projects/rll_testdir/mapping_file.txt',out_dir="/home/cheng/Projects/rll_testdir/")
    m.run()
    """

    def __init__(self, raw_fqs_dir, pre_mapping_file, categories=False, run_size=10, host_type="hg38", sample_regex="(.+)_.*_[12]\.fq\.gz", forward_regex="_1\.fq\.gz$", reverse_regex="_2\.fq\.gz$", out_dir='/home/cheng/Projects/rll_testdir/test/'):
        self.out_dir = os.path.abspath(out_dir) + '/'
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        self._base_dir = os.path.dirname(__file__) + '/'
        with open(self._base_dir + "pipconfig/path.conf") as f:
            self.path = json.load(f)
        self.parsed_map = parse_premap(raw_fqs_dir, pre_mapping_file, forward_regex, reverse_regex, sample_regex)
        self.fq_info = self.parsed_map['fastq'].values
        self._init_outdir_()
        self.mapping_file = self.out_dir + 'mapping_file.txt'
        self.parsed_map['map'].to_csv(self.mapping_file, sep='\t', index=False)
        self.categories = categories if categories else ','.join(
            [g for g in self.parsed_map['map'].columns if not g.find('Group') == -1])
        print("The detected categories are: \n    {}\n".format(self.categories))
        self.raw_pattern = self.out_dir + "Raw_fastq/{}_{}.fq.gz"
        trimmed_pattern = self.out_dir + "Primer_trimmed/{}_{}.fq.gz"
        filtered_pattern = self.out_dir + "Filtered/{}_{}.good.fastq.gz"
        de_host_pattern = self.out_dir + \
            "Host_subtracted/bowtie/%s/{}_{}.%s.unmapped.fastq.gz" % (host_type, host_type)
        kracken2_reports_pattern = self.out_dir + "Kraken2/{}_{}.report"
        self.merged_pe_pattern = self.out_dir + \
            "Host_subtracted/bowtie/%s/{}_R1.%s.unmapped.{}.fastq.gz" % (host_type, host_type)
        self.raw_list = self.map_list(self.raw_pattern, run_size)
        self.trimmed_list = self.map_list(trimmed_pattern, run_size)
        self.filtered_list = self.map_list(filtered_pattern, run_size)
        self.de_host_r1_list = self.map_list(de_host_pattern, run_size, use_direction='R1')
        self.merged_pe_r1_list = self.map_list(self.merged_pe_pattern, run_size, use_direction='R1')
        self.kracken2_reports_list = self.map_list(kracken2_reports_pattern, use_direction='R1')
        global running_list
        running_list = self.out_dir + '.running_list'

    def _init_outdir_(self):
        os.system('perl ' + self.path['cii_home'] + 'create_dir_structure.pl ' + self.out_dir)
        for fq_path, new_id, direction in self.fq_info:
            formated_fq = self.out_dir + "Raw_fastq/{}_{}.fq.gz".format(new_id, direction)
            if not os.path.exists(formated_fq):
                os.system("ln -s {} {}".format(fq_path, formated_fq))

    def map_list(self, pattern=False, each=False, use_direction="both"):
        out = []
        for fq_path, new_id, direction in self.fq_info:
            ele = pattern.format(new_id, direction) if pattern else [new_id, direction]
            if direction == "R2" and (use_direction == "both" or use_direction == "R2"):
                out.append(ele)
            elif direction == "R1" and (use_direction == "both" or use_direction == "R1"):
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

    def homized_cmd(self, cmd, home=False):
        home = home or self.path['cii_home']
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
        os.system(self.homized_cmd('perl host_subtraction_wrapper.pl --aim Host --config human --single-end-fq-list {} --threads {} -m N'.format(
            fq_list, str(processor)), home=self.path['de_host_home']))

    def generate_summary(self):
        os.system("mv {}Host_subtracted/bowtie/QC_report/* {}QC_report/".format(self.out_dir, self.out_dir))
        global running_list
        with open(running_list, 'w') as f:
            f.write('\n'.join(self.map_list(self.raw_pattern, use_direction="R1")))
        os.system(self.homized_cmd(
            "perl generate_summary_wrapper.pl --aim Host --file-list {}".format(running_list)))

    @synchronize
    def run_merge_se_to_pe(self, fq_list: list):
        os.system(self.homized_cmd(
            'perl merge_se_into_pe_fastq_wrapper.pl --paired-end-fq-list {}'.format(fq_list), home=self.path['merge_se_home']))

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
            self.homized_cmd("perl Braken_to_OTUtable.pl {} {}".format(self.path['ncbi_taxaID_path'], running_list)))

    @synchronize
    def run_metaphlan2(self, fq_list: list, processor=2):
        os.system(self.homized_cmd("perl run_metaphlan2.pl {} PE {}".format(
            fq_list, str(processor))))

    def clean_header(self, tb_name, pattern, skip=[]):
        df = pd.read_csv(tb_name, sep='\t')
        df.columns = [re.sub(pattern, '', c) for c in df.columns]
        df.drop(index=skip).to_csv(tb_name, sep="\t", index=False)

    def join_metaphlan(self):
        tb_name = self.out_dir + 'Metagenome/Metaphlan/All.Metaphlan2.profile.txt'
        os.system("python2 {}/merge_metaphlan_tables.py {}Metagenome/Metaphlan/*profile.txt > {}".format(
            self.path['cii_home'], self.out_dir, tb_name))
        self.clean_header(tb_name, pattern='_R1.metaphlan.profile$', skip=[0])

    @synchronize
    def run_humann(self, fq_list: list):
        os.system(self.homized_cmd('perl run_humann.pl {} SE'.format(fq_list)))

    def join_humann(self):
        os.system('''
out_dir=%s
ln -s ${out_dir}Metagenome/Humann/*/*genefamilies.tsv ${out_dir}Metagenome/Humann/*/*pathabundance.tsv ${out_dir}Metagenome/Humann/
humann2_join_tables -i ${out_dir}Metagenome/Humann/ -o ${out_dir}Metagenome/Humann/All.Humann2.genefamilies.tsv --file_name genefamilies.tsv
humann2_join_tables -i ${out_dir}Metagenome/Humann/ -o ${out_dir}Metagenome/Humann/All.Humann2.pathabundance.tsv --file_name pathabundance.tsv''' % (self.out_dir))
        self.clean_header(self.out_dir + "Metagenome/Humann/All.Humann2.genefamilies.tsv",
                          pattern='_R1\..*\.unmapped\.R1_Abundance-RPKs$')
        self.clean_header(self.out_dir + "Metagenome/Humann/All.Humann2.pathabundance.tsv",
                          pattern='_R1\..*\.unmapped.R1_Abundance$')

    @synchronize2
    def run_fmap(self, fq_list, database="protein_fasta_protein_homolog_model", processor=4, out_dir='AMR'):
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
        info_list = self.map_list(use_direction='R1')
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
                          database="orthology_uniref90_2_2157_4751.20190412161853", processor=processor, out_dir="FMAP")
            self.quantify_fmap(out_dir="FMAP", all_name="All.Function.abundance.KeepID.KO.txt")
        elif run_type == "AMR":
            self.run_fmap(fq_list=self.merged_pe_r1_list, database="protein_fasta_protein_homolog_model_cleaned",
                          processor=processor, out_dir="AMR")
            self.quantify_fmap(out_dir="AMR", all_name="All.AMR.abundance.txt")
        elif run_type == "ARDB":
            self.run_fmap(fq_list=self.merged_pe_r1_list, database="ARDB.20180725064354",
                          processor=processor, out_dir="ARDB")
            self.quantify_fmap(out_dir="ARDB", all_name="All.ARDB.abundance.txt", print_definition=True)

    def run_assembly(self, processor=50):
        os.system(
            "{} -1 {} -2 {} --min-contig-len 1000 -t {} -o {}Assembly/Assembly".format(
                self.path['megahit_path'], ','.join(self.map_list(self.merged_pe_pattern, use_direction="R1")), ','.join(self.map_list(self.merged_pe_pattern, use_direction="R2")), processor, self.out_dir)
        )
        os.system(
            "cd {}Assembly/Assembly&&prodigal -i final.contigs.fa -a final.contigs.fa.faa -d final.contigs.fa.fna  -f gff -p meta -o final.contigs.fa.gff&&\
            {} -i final.contigs.fa.fna -c 0.95 -G 0 -aS 0.9 -g 1 -d 0 -M 80000 -o final.contigs.fa.fna.out -T 0&&\
            grep '>' final.contigs.fa.fna.out > final.contigs.fa.fna.out.header&&\
            {}&&\
            Rscript ORF_header_summary.R -i {}Assembly/Assembly/".format(
                self.out_dir, self.path['cdhit_path'],
                self.homized_cmd(
                    "perl ORF_generate_input_stats_file.pl {}Assembly/Assembly/final.contigs.fa.fna.out.header".format(self.out_dir)),
                self.out_dir
            )
        )

    def map_ko_annotation(self):
        ko_file = self.out_dir + 'FMAP/All.Function.abundance.KeepID.KO.txt'
        out = os.path.dirname(ko_file)
        os.system('''
SCRIPTPATH=%s
ko_file=%s
out_dir=%s
perl ${SCRIPTPATH}/ConvergeKO2Module.pl $ko_file > ${out_dir}/All.Function.abundance.KeepID.Module.txt
perl ${SCRIPTPATH}/ConvergeKO2Pathway.pl $ko_file > ${out_dir}/All.Function.abundance.KeepID.Pathway.txt
perl ${SCRIPTPATH}/ConvergePathway2Level1.pl ${out_dir}/All.Function.abundance.KeepID.Pathway.txt > ${out_dir}/All.Function.abundance.KeepID.Pathway.Level1.txt
perl ${SCRIPTPATH}/ConvergePathway2Level2.pl ${out_dir}/All.Function.abundance.KeepID.Pathway.txt > ${out_dir}/All.Function.abundance.KeepID.Pathway.Level2.txt
            ''' % (self.path['cii_home'], ko_file, out))

    def map_func_definition(self):
        FMAP_data = self.path['fmap_home'] + '/FMAP_data/'
        FMAP = self.out_dir + "FMAP/"
        mi = MapInfo()
        mi.mapping(FMAP + 'All.Function.abundance.KeepID.Pathway.txt', [FMAP_data + f for f in ['KEGG_Pathway2Level1.txt', 'KEGG_Pathway2Level2.txt', 'KEGG_pathway.txt']],
                   out_file=FMAP + 'All.Function.abundance.Pathway.full_info.txt', mapped_headers=["Level1", "Level2", "Level3"])
        mi.mapping(FMAP + 'All.Function.abundance.KeepID.KO.txt', ['/home/cheng/softwares/FMAP/FMAP_data/KEGG_orthology.txt'],
                   out_file=FMAP + 'All.Function.abundance.KO.full_info.txt', adjust_func=mi.ajust_ko_info, mapped_headers=['Gene_name\tEnzyme_number\tDefinition'])
        mi.mapping(FMAP + 'All.Function.abundance.KO.full_info.txt', [FMAP_data + f for f in [
                   'KEGG_orthology2module.txt', 'KEGG_orthology2pathway.txt']], mapped_headers=["Module", "KEGG Pathway"], add=True)
        datas = [FMAP + f for f in ["All.Function.abundance.KeepID.KO.txt",
                                    "All.Function.abundance.KeepID.Module.txt", "All.Function.abundance.KeepID.Pathway.txt"]]
        mapping_sources = [FMAP_data +
                           f for f in ["KEGG_orthology.txt", "KEGG_module.txt", "KEGG_pathway.txt"]]
        for data, mapping_source in zip(datas, mapping_sources):
            mi.mapping(data, [mapping_source])
        mi.mapping(self.out_dir + 'AMR/All.AMR.abundance.txt', [FMAP_data + '/aro.csv'],
                   pattern="ARO[^\|]+", first_pattern="[^\|]+$", add_sid_to_info=True)

    def run_quast(self):
        os.system("python2 {} -o {}Assembly/Assembly/quast_results/quast_results/  {}Assembly/Assembly/final.contigs.fa".format(
            self.path['quast_path'], self.out_dir, self.out_dir))

    def run(self, processor=2, base_on_assembly=False):
        self.run_fastqc(fq_list=self.raw_list, processor=processor, first_check=5)
        self.run_trim(fq_list=self.raw_list)
        self.run_filter(fq_list=self.trimmed_list)
        self.run_fastqc(fq_list=self.filtered_list, processor=processor, first_check=5)
        self.run_de_host(fq_list=self.filtered_list, processor=processor)
        self.run_merge_se_to_pe(fq_list=self.de_host_r1_list)
        self.run_fastqc(fq_list=self.de_host_r1_list, processor=processor, first_check=5)
        self.generate_summary()
        self.run_kraken2(fq_list=self.merged_pe_r1_list, processor=processor)
        self.run_bracken()
        self.run_metaphlan2(fq_list=self.merged_pe_r1_list, processor=processor)
        self.join_metaphlan()
        self.run_humann(fq_list=self.merged_pe_r1_list)
        self.join_humann()
        if base_on_assembly:
            self.run_assembly(processor=processor * 20)
            self.run_quast()
        else:
            self.fmap_wrapper(run_type="KEGG", processor=processor * 4)
            self.fmap_wrapper(run_type="AMR", processor=processor * 4)
        self.map_ko_annotation()
        self.map_func_definition()
        VisualizeAll(self.mapping_file, self.categories).visualize()
