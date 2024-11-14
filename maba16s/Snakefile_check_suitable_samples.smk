"""
MABA16S Preprocessing Pipeline
==============================

Description:
------------
This Snakemake pipeline is designed for the preprocessing and taxonomic 
classification of metagenomic 16S data, using tools like Kraken2 and BLAST.
The workflow includes steps for downloading databases, filtering reads, 
concatenating files, and generating a final sample configuration.

Main Steps:
-----------
1. Download Kraken2 SILVA database
2. Concatenate ONT read files
3. Filter reads by minimum length
4. Perform taxonomic classification with Kraken2
5. Select and filter Kraken2 output for genera
6. Calculate basic QC parameters - read counts and number of genera
7. Write a configuration file with samples for further analysis

"""

import time
import os
from shutil import copy2
import pathlib
import pandas as pd

configfile: "config/config.yaml"
SAMPLES = config['SAMPLES']

OUTDIR = config['parameters']['outdir'] + "/"

# workflow start: set up output directories & save the config file
onstart:
    print("This is MABA16S")
    time.sleep(1)

    # copy the config file to output dir
    pathlib.Path(OUTDIR).mkdir(parents=True, exist_ok=True)
    copy2('config/config.yaml', OUTDIR)

    for i in SAMPLES.items():
        print(i[0], '\t', i[1])
    print(f'The output directory is: {OUTDIR}')

# Error handling
onerror:
    print(f"Attempted to analyse the following samples: {SAMPLES}")
    print("An error has occured")


# Define local rules
localrules: all, download_kraken2_db, combinereads, strip_genera, qc_preprocessing, write_good_samples

# Master rule
rule all:
    input:
        OUTDIR + "config/config_good_samples.yaml",
        expand(OUTDIR + "QC/{sample}/{sample}_qcPreprocessing.txt", sample=SAMPLES)

# Download the database required for Kraken2
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

# Concetante ONT reads from multiple fastq files
rule combinereads:
    input:
        lambda wildcards: SAMPLES[wildcards.sample]
    output:
        OUTDIR + "reads/{sample}.fastq.gz"
    threads: 1
    shell:
        "cat {input}/*fastq* > {output}"

# Apply a minimum read length filter
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
        OUTDIR + "log/{sample}/filter_read_length.log"
    shell: 
        "filtlong --min_length {params.min_length} {input} | gzip > {output}  2> {log}"

# Taxonomic classification with kraken2
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

# Select and filter kraken2 output 
rule strip_genera:
    input:
        rules.kraken2.output.report
    output:
        OUTDIR + "kraken2/{sample}/krakenreport_filtered.txt"
    params:
        reads_cutoff_genus = 50
    shell:
        '''awk '$4 == "G" {{print $0}}' {input} | awk '$2 > {params.reads_cutoff_genus} {{print $0}}' > {output}'''

# QC pre-processing step
rule qc_preprocessing:
    input:
        raw = rules.combinereads.output,
        filtered = rules.filter_read_length.output,
        genera = rules.strip_genera.output
    output:
        OUTDIR + "QC/{sample}/{sample}_qcPreprocessing.txt"
    log:
        OUTDIR + "log/{sample}/qcPreprocessing.txt"
    shell:
        """       
        # Get raw reads count (line count divided by 4 for FASTQ format)
        rawReads=$(( $(zcat {input.raw} | wc -l) / 4 ))

        # Get filtered reads count (line count divided by 4 for FASTQ format)
        filteredReads=$(( $(zcat {input.filtered} | wc -l) / 4 ))

        # Count the number of genera in the filtered Kraken2 report
        nGenera=$(wc -l < {input.genera})

        # Write results to the output file
        echo -e "sampleID\trawReads\tfilteredReads\tnGenera" > {output}
        echo -e "{wildcards.sample}\t${{rawReads}}\t${{filteredReads}}\t${{nGenera}}" >> {output} 2> {log}
        """       

# Generate a new config file: samples for further analysis only
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
