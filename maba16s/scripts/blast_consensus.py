import sys
import os
import glob

def main():
    inputfasta = sys.argv[1]
    threads = sys.argv[2]
    db = sys.argv[3]
    output = sys.argv[4]

    files = glob.glob(inputfasta + "*fasta")
    for fasta in files:
        genusname = os.path.basename(fasta).split(".")[0]
        os.system(f'mkdir -p {output}')
        os.system(f'blastn -db {db}/blastDB -query {fasta} -num_threads {threads} -out {output}/{genusname}.txt -outfmt "7 stitle"')


if __name__ == "__main__":
    main()
