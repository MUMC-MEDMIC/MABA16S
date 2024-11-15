import sys
import os
import glob

def main():
    indir = sys.argv[1]
    threads = sys.argv[2]
    db = sys.argv[3]
    outdir = sys.argv[4]

    print(f'Blasting files from {indir}')

    files = glob.glob(os.path.join(indir, "*.fasta"))
    for fasta in files:
        os.system(f'mkdir -p {outdir}')
        genusname = os.path.basename(fasta).replace("_consensus.fasta", "")
        print(f'Blasting {genusname} fasta')
        os.system(f'blastn -db {db}/blastDB -query "{fasta}" -num_threads {threads} -out "{outdir}/{genusname}_BLASTn.txt" -outfmt "6 pident length bitscore stitle"')

if __name__ == "__main__":
    main()
