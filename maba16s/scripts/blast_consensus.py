import sys
import os
import glob

def main():
    inputfasta = sys.argv[1]
    threads = [2]
    db = sys.argv[3]
    output = sys.argv[4]

    files = glob.glob(inputfasta + "*fasta")
    print(files)
    for fasta in files:
        genusname = os.path.basename(fasta).split(".")[0]
        print(genusname)
        os.system(f'mkdir -p {output}')
        os.system(f'blastn -db {db}/blastDB -query {fasta} -out {output}/{genusname}.txt -outfmt 11')


if __name__ == "__main__":
    main()
