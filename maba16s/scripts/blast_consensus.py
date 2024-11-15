"""
BLASTn Query and Output Script
==============================

Description:
------------
This script performs BLASTn searches for genus-specific FASTA files against a 
pre-built BLAST database and outputs the results in a tabular format.

Workflow:
---------
1. Reads input FASTA files from a specified directory.
2. Executes BLASTn for each FASTA file using the specified database and number 
   of threads.
3. Writes the BLASTn output in tabular format (`pident`, `length`, `bitscore`, 
   `stitle`) to the specified output directory.

Arguments:
----------
1. Input directory containing genus-specific FASTA files.
2. Number of threads to use for BLASTn.
3. Path to the BLAST database directory.
4. Output directory for BLASTn results.

Output:
-------
BLASTn result files for each input FASTA file in the specified output directory.
Each result file is named `<genusname>.txt` and contains the BLASTn output in 
tabular format.
"""

import sys
import os
import glob

def main():
    indir = sys.argv[1]
    threads = sys.argv[2]
    db = sys.argv[3]
    outdir = sys.argv[4]

    print(f'Blasting files from {indir}')

    os.system(f'mkdir -p {outdir}')
    files = glob.glob(os.path.join(indir, "*.fasta"))
    for fasta in files:
        genusname = os.path.basename(fasta).replace("_consensus.fasta", "")
        print(f'Blasting {genusname} fasta')
        os.system(f'blastn -db {db}/blastDB -query "{fasta}" -num_threads {threads} -out "{outdir}/{genusname}_BLASTn.txt" -outfmt "6 pident length bitscore stitle"')

if __name__ == "__main__":
    main()
