import pandas as pd
import os
import sys
import re

def extract_reads(speciesID, name, krakenfile, genusID, fastqin, outdir):
    os.system(f'mkdir -p {outdir}')
    os.system(f'extract_kraken_reads.py -k {krakenfile} -o {outdir + "/" + name + ".fasta"} -s {fastqin} -t {genusID}')


def read_krakenreport(file, krakenfile, fastqin, outdir):
    df = pd.read_csv(file, sep = '\t', header = None)
    for i in df.iterrows():
        info = i[1]
        name = re.sub(r'\W+', '', info[5])
        genusID = info[4]
        extract_reads(genusID, name, krakenfile, genusID, fastqin, outdir)



def main():
    read_krakenreport(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])


if __name__ == "__main__":
    main()
