path = {
    "bracken_path": "/home/bayegy/pipelines/metagenome/softwares/bracken/src/est_abundance.py",
    "bracken_database": "/home/bayegy/Databases/kraken_bracken/database150mers.kmer_distrib",
    "fmap_home": "/home/bayegy/pipelines/metagenome/softwares/FMAP/",
    "bayegy_home": "/home/bayegy/pipelines/metagenome/Bayegy/",
    "megahit_path": "/home/bayegy/pipelines/metagenome/miniconda2/bin/megahit",
    # conda install -c bioconda cd-hit # cd-hit installed by conda does not suport multi-threading
    "cdhit_path": "/home/bayegy/pipelines/metagenome/softwares/cd-hit-v4.8.1-2019-0228/cd-hit-est",
    "circos_path": "/home/bayegy/pipelines/metagenome/miniconda2/bin/circos",
    "circos_etc": "/home/bayegy/pipelines/metagenome/miniconda2/etc",
    "ncbi_taxaID_path": "/home/bayegy/Databases/kraken_bracken/taxid2OTU_ranks.txt",
    "quast_path": "/home/bayegy/pipelines/metagenome/miniconda2/bin/quast.py",
    "kraken2_database": "/home/bayegy/Databases/kraken_bracken/",
    "emapper_database": "/home/bayegy/Databases/emapper/eggnog",
    "cazy_database": "/home/bayegy/Databases/diamond/CAZyDB.07312018.dmnd",
    "humann2_home": "/home/bayegy/pipelines/metagenome/miniconda2/bin/",
    # conda install humann2 -c biobackery # do not use bioconda channel
    "humann2_protein_database": "/home/bayegy/Databases/humann2/uniref90_full",
    "humann2_nucleotide_database": "/home/bayegy/Databases/humann2/chocophlan_full",
    "humann2_utility_mapping": "/home/bayegy/Databases/humann2/utility_mapping",
    "map_conf": "/home/bayegy/Databases/colormap",
    "kneaddata_path": "/home/bayegy/pipelines/metagenome/miniconda2/bin/kneaddata",
    "trimmomatic_home": "/home/bayegy/pipelines/metagenome/softwares/Trimmomatic-0.39",
    "bowtie2_home": "/home/bayegy/pipelines/metagenome/miniconda2/bin/",
    "metaphlan2_home": "/home/bayegy/pipelines/metagenome/miniconda2/bin/",
    "metaphlan2_database": "/home/bayegy/pipelines/metagenome/miniconda2/bin/databases",
    "diamond_home": "/home/bayegy/pipelines/metagenome/miniconda2/bin/",
    "salmon_path": "/home/bayegy/pipelines/metagenome/softwares/salmon-latest_linux_x86_64/bin/salmon",
    "kraken2_path": "/home/bayegy/pipelines/metagenome/miniconda2/bin/kraken2",
    "adapters_path": "/home/bayegy/pipelines/metagenome/MetaGenome/data/adaptor_Illumina.fa",
    # conda install emboss
    "transeq_path": "/home/bayegy/pipelines/metagenome/miniconda2/bin/transeq",
    # conda install -c bioconda eggnog-mapper
    "emapper_path": "/home/bayegy/pipelines/metagenome/miniconda2/bin/emapper.py",
    "fastqc_home": "/home/bayegy/pipelines/metagenome/miniconda2/bin/",
    "python3_path": "/home/bayegy/pipelines/metagenome/miniconda2/bin/python3.8",
    "R_path": "/home/bayegy/pipelines/metagenome/miniconda2/bin/Rscript",
    "perl_path": "/home/bayegy/pipelines/metagenome/miniconda2/bin/perl",
    "qiime2_home": "/home/bayegy/pipelines/metagenome/miniconda2/envs/qiime2-2019.10",
    "lefse_pylib_home": "/home/bayegy/pipelines/metagenome/miniconda2/share/lefse-1.0.8.post1-1",
    "lefse_rlib_home": "/home/bayegy/pipelines/metagenome/miniconda2/lib/R/library",
    "lefse_py_home": "/home/bayegy/pipelines/metagenome/miniconda2/bin/",
    "prodigal_path": "/home/bayegy/pipelines/metagenome/miniconda2/bin/prodigal",
    "sortmerna_path": "/home/bayegy/pipelines/metagenome/miniconda2/envs/sortmerna/bin/sortmerna",
    "sortmerna_databases": "/home/bayegy/Databases/rRNA_databases"
    # sudo apt-get install ttf-mscorefonts-installer #安装字体
    # /etc/ImageMagick-6/policy.xml replace the value of PDF rights to "read|write"
    # do not use conda to install salmon, it's too old
    # sudo apt-get install mpich #cd-hit nedd this
    # humann2 installed by bioconda is too old!
    # sudo apt-get install libopenblas-dev # to solve 'libopenblas.so.0: cannot open shared object file'

}

memery_needs = {
    "kneaddata": 160,
    "fastqc": 160,
    "sortmerna": 160,
    "kraken2": 100,
    "humann2": 100,
    "fmap": 30,
    "salmon": 1000,
    "megahit": 200,
    "prodigal": 50,
    "cdhit": 100,
    "diamond_cazy": 100,
    "diamond_card": 100,
    "emapper": 100,
}

memery = 3000

threads = 270

hosts = 3

use_sge = True

sge_pe = "smp"

sge_queue = "metagqueue"

max_workers = {
    "kneaddata": 9,
    "fastqc": 9,
    "sortmerna": 9,
    "kraken2": 9,
    "humann2": 9,
    "fmap": 18,
    "salmon": 3,
    "megahit": 9,
    "prodigal": 30,
    "cdhit": 1,
    "diamond_cazy": 12,
    "diamond_card": 12,
    "emapper": 21,
}


each = {
    "prodigal": 30,
    "cdhit": 9,
    "diamond_cazy": 3,
    "diamond_card": 3,
    "emapper": 3,
}
