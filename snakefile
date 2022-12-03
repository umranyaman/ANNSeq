import os
import Bio
import shutil
from os import path
from re import search
from pathlib import Path

from snakemake.utils import validate
from snakemake.utils import min_version
min_version("5.3")

# ----------------------------------------------------------------

configfile: "config.yml"
validate(config, schema="schema/config_schema.yaml")
workdir: config["workdir"]

WORKDIR = config["workdir"]
SNAKEDIR = path.dirname(workflow.snakefile)

shutil.copy2(SNAKEDIR + "/config.yml", WORKDIR)

sample = config["sample_name"]


# ----------------------------------------------------------------

target_list = [
    "Nanostat/stat_out.txt",
    "Talon/" + sample + "_talon_abundance_filtered.tsv",
    "Talon/" + sample + "_talon.gtf"
]

if config.get("reads_fastq") == "":
    target_list.remove("Nanostat/stat_out.txt")



rule all:
    input:
        target_list


#########################################################################

rule concenate:
    input:
        fq = config["reads_fastq"]
        
    output:
        outfq = "processed_reads/input_reads.fq"
    
    shell:
        """
        mkdir -p processed_reads;
        if [[ {params.concat} == "True" ]];
        then
            find {input.fq}  -regextype posix-extended -regex '.*\.(fastq|fq)$' -exec cat {{}} \\; > processed_reads/input_reads.fq
        else
            ln -s `realpath {input.fq}` processed_reads/input_reads.fq
        fi
        """
        
rule nanostat:
    input:
        fq = rules.concenate.output.outfq

    output:
        ns = "Nanostat/stat_out.txt"

    threads: config["threads"]

    shell:
        """
        NanoStat -n {output.ns} -t {threads} --tsv --fastq {input.fq}
        """


rule pychopper:
    input:
        fq = rules.concenate.output.outfq

    output:
        pyfq = "Pychopper/" + sample + "_pychopped.fastq"

    params:
        outpath = "Pychopper",
        prefix = sample + "_pychopped.fastq"

    log: WORKDIR + "/Pychopper/" + sample + "_pychopped.log"

    threads: config["threads"]

    shell:
        """
        cd {params.outpath};
        (cdna_classifier.py -t {threads} -r report.pdf -u unclassified.fq -w rescued.fq {input.fq} {params.prefix}) 2> {log}
        """


rule minimap_mapping:
    input:
        genome = config["genome"],
        fq = rules.pychopper.output.pyfq

    output:
        sam = "Mapping/" + sample + "_minimap.sam"

    log: "Mapping/" + sample + "_minimap.log"

    threads: config["threads"]

    shell:
        """
        (minimap2 -ax splice -uf -k14 --MD --secondary=no -t {threads} -o {output.sam} {input.genome} {input.fq}) 2> {log}
        """



rule sort_sam:
    input:
        sam = rules.minimap_mapping.output.sam

    output:
        sortedSam = "Mapping/" + sample + "_minimap.sorted.sam"

    threads: config["threads"]

    shell:
        """
        samtools sort -O SAM -o {output.sortedSam} -@ {threads} {input.sam}
        """


rule transcriptclean:
    input:
        genome = config["genome"],
        sam2 = rules.sort_sam.output.sortedSam
         
    output:
        clsam = "Mapping/" + sample + "_minimap.sortedcleaned.sam"
    
    shell:
        """
        python /home/DSUmranYaman/TranscriptClean/TranscriptClean.py --sam {input.sam2} --genome {input.genome}
        """
        
        

rule talon_initialize_database:
    input:
        gtf = config["gtf"]

    output:
        db = "Talon/talon_" + Path(config["gtf"]).stem + ".db"

    params:
        gtf_name = Path(config["gtf"]).stem,
        gbuild = "mm10",
        prefix = "Talon/talon_" + Path(config["gtf"]).stem

    log: "Talon/talon_db.log"

    shell:
        """
        (talon_initialize_database --f {input.gtf} --a {params.gtf_name} --g {params.gbuild} --o {params.prefix}) &> {log}
        """



rule talon_label_reads:
    input:
        sam = rules.sort_sam.output.sortedSam,
        fa = config["genome"]

    output:
        sam = "Talon/" + sample + "_labeled.sam"

    params:
        prefix = "Talon/" + sample

    log: "Talon/talon_labelReads.log"

    threads: config["threads"]

    shell:
        """
        (talon_label_reads --f {input.sam} --g {input.fa} --t {threads} --deleteTmp --o {params.prefix}) &> {log}
        """



rule talon_annotate:
    input:
        sam = rules.talon_label_reads.output.sam,
        talon_db = rules.talon_initialize_database.output.db

    output:
        qc = "Talon/" + sample + "_QC.log",
        tsv = "Talon/" + sample + "_talon_read_annot.tsv"

    params:
        prefix1 = sample,
        prefix2 = "Talon/" + sample,
        description = sample + "_sam",
        gbuild = "mm10"

    log: "Talon/talon_annotate.log"

    threads: config["threads"]

    shell:
        """
        echo '{params.prefix1},{params.description},ONT-Promethian,{input.sam}' > Talon/config_file.txt &&
        (talon --f Talon/config_file.txt --db {input.talon_db} --build {params.gbuild} --threads {threads} --o {params.prefix2}) &> {log}
        """


rule talon_filter_transcripts:
    input:
        talon_db = rules.talon_initialize_database.output.db,
        job_hold = rules.talon_annotate.output.tsv

    output:
        csv = "Talon/" + sample + "_whitelist.csv"

    params:
        gtf_name = Path(config["gtf"]).stem

    log: "Talon/talon_filter.log"

    shell:
        """
        (talon_filter_transcripts --db {input.talon_db} -a {params.gtf_name} --maxFracA 0.5 --minCount 5 --o {output.csv}) &> {log}
        """



rule talon_abundance:
    input:
        talon_db = rules.talon_initialize_database.output.db,
        whitelist = rules.talon_filter_transcripts.output.csv

    output:
        tsv = "Talon/" + sample + "_talon_abundance_filtered.tsv"

    params:
        gtf_name = Path(config["gtf"]).stem,
        gbuild = "mm10",
        prefix = "Talon/" + sample

    log: "Talon/talon_abundance.log"

    shell:
        """
        (talon_abundance --db {input.talon_db} -a {params.gtf_name} --build {params.gbuild} --whitelist {input.whitelist} --o {params.prefix}) &> {log}
        """



rule talon_create_GTF:
    input:
        talon_db = rules.talon_initialize_database.output.db,
        whitelist = rules.talon_filter_transcripts.output.csv

    output:
        talon_gtf = "Talon/" + sample + "_talon.gtf"

    params:
        gtf_name = Path(config["gtf"]).stem,
        gbuild = "mm10",
        prefix = "Talon/" + sample

    log: "Talon/talon_gtf.log"

    shell:
        """
        (talon_create_GTF --db {input.talon_db} -a {params.gtf_name} --build {params.gbuild} --whitelist {input.whitelist} --o {params.prefix}) &> {log}
        """
