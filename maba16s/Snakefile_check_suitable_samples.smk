import time
import os
from shutil import copy2
import pathlib
import pandas as pd

configfile: "config/config.yaml"
SAMPLES = config['SAMPLES']

OUTDIR = config['parameters']['outdir'] + "/"

onstart:
    print("This is MABA16S")
    time.sleep(1)

    # copy the config file to output dir
    pathlib.Path(OUTDIR).mkdir(parents=True, exist_ok=True)
    copy2('config/config.yaml', OUTDIR)

    for i in SAMPLES.items():
        print(i[0], '\t', i[1])
    print(f'output directory is: {OUTDIR}')

onerror:
    print(f"these were the samples: {SAMPLES}")
    print("error has occured")


localrules: all, download_kraken2_db, combinereads

rule all:
    input:
        OUTDIR + "config/config_good_samples.yaml"


rule download_kraken2_db:
    output:
        directory("db/silva")
    threads: 1
    conda:
        "envs/kraken2.yaml"
    params:
        "db/silva"
    shell:
        "kraken2-build --db {params} --special silva"


rule build_blast_db:
    input:
        "db/silva/data/SILVA_138.1_SSURef_NR99_tax_silva.fasta"
    output:
        directory("db/blastDB")
    threads: 1
    conda:
        "envs/blast.yaml"
    log:
        "log/blastbuild_db.log"
    params:
        "db/blastDB"
    shell:
        "makeblastdb -in {input} -dbtype nucl -out {output}/blastDB 2> {log}"


rule combinereads:
    input:
        lambda wildcards: SAMPLES[wildcards.sample]
    output:
        temp(OUTDIR + "reads/{sample}.fastq.gz")
    threads: 1
    shell:
        "cat {input}/*fastq* > {output}"


rule filter_read_length:
    input:
        rules.combinereads.output
    output:
        OUTDIR + "filtered_reads/{sample}.fastq.gz"
    conda:
        "envs/filter_reads.yaml"
    threads: 1
    params:
        min_length = 1200
    log: 
        OUTDIR + "log/filter_read_length/{sample}.log"
    shell: 
        "filtlong --min_length {params.min_length} {input} | gzip > {output}  2> {log}"

rule kraken2:
    input:
        reads = rules.filter_read_length.output,
        db = rules.download_kraken2_db.output
    output:
        report = OUTDIR + "kraken2/{sample}/krakenreport.txt",
        out = OUTDIR + "kraken2/{sample}/output.txt"
    threads: 4
    conda:
        "envs/kraken2.yaml"
    log:
        OUTDIR + "log/{sample}/kraken2.txt"
    shell:
        "kraken2 --db {input.db} -t {threads} "
        "--report {output.report} {input.reads} > {output.out} "


rule strip_genera:
    input:
        rules.kraken2.output.report
    output:
        OUTDIR + "kraken2/{sample}/krakenreport_filtered.txt"
    params:
        reads_cutoff_genus = 50
    shell:
        '''awk '$4 == "G" {{print $0}}' {input} | awk '$2 > {params.reads_cutoff_genus} {{print $0}}' > {output}'''


rule write_good_samples:
    input:
        expand(OUTDIR + "kraken2/{sample}/krakenreport_filtered.txt", sample = SAMPLES)
    output:
        OUTDIR + "config/config_good_samples.yaml"
    run:
        with open(output[0], "w") as outfile:
            outfile.write("GOODSAMPLES: \n")
            for f in input:
                if os.path.getsize(f) > 0:
                    samplename = f.split("/")[-2]
                    outfile.write("  "+ samplename + ": " + SAMPLES[samplename] + "\n")
