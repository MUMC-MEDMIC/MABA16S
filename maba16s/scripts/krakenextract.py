"""
Kraken Extract Script
=====================

Description:
------------
This script processes a Kraken2 taxonomic classification report to extract 
reads associated with specific genera. It generates FASTA files for each 
genus by calling the `extract_kraken_reads.py` tool from the kraken-tools package.

Workflow:
---------
1. Read the Kraken2 report file.
2. Extract genus names and IDs.
3. Create output directories for each genus.
4. Extract reads for each genus using `extract_kraken_reads.py`.

Arguments:
----------
1. Input Kraken2 report file (path).
2. Input Kraken2 classification file (path).
3. Input FASTQ file (path).
4. Output directory (path).

"""

import pandas as pd
import os
import sys

def extract_reads(speciesID, name, krakenfile, genusID, fastqin, outdir):
    os.system(f'mkdir -p {outdir}')
    os.system(f'extract_kraken_reads.py -k {krakenfile} -o "{outdir + "/" + name + ".fasta"}" -s {fastqin} -t {genusID}')


def read_krakenreport(file, krakenfile, fastqin, outdir):
    df = pd.read_csv(file, sep = '\t', header = None)
    for i in df.iterrows():
        info = i[1]
        name = info[5].strip()
        genusID = info[4]
        extract_reads(genusID, name, krakenfile, genusID, fastqin, outdir)



def main():
    read_krakenreport(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])


if __name__ == "__main__":
    main()
