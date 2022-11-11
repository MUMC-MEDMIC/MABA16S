import sys
import glob
import os

def catch_right_ref(genusname, ref, dumpdir):
    ID = os.popen(f'grep {genusname} {ref} | head -n 1 | cut -b 2-').read().split()[0]
    os.system(f'samtools faidx {ref} {ID} >  {dumpdir}/{genusname}_reference.fasta ')

def run_minimap2(reads, reference, outdir):
    os.system('minimap2 -ax')

def align_reads(indir, ref, outdir):
    os.system(f'mkdir -p {outdir}')
    for file in glob.glob(indir + "/*.fastq.gz"):
        genusname = os.path.basename(file).split(".")[0]
        catch_right_ref(genusname, ref, outdir)


def main():
    print(str(sys.argv[0]))
    align_reads(indir = sys.argv[1], ref = sys.argv[2], outdir = sys.argv[3])

if __name__ == "__main__":
    main()
