import sys
import glob
import os

def catch_right_ref(genusname, ref, dumpdir):
    ID = os.popen(f'grep {genusname} {ref} | head -n 1 | cut -b 2-').read().split()[0]
    os.system(f'samtools faidx {ref} {ID} >  {dumpdir}/{genusname}_reference.fasta ')

def run_minimap2(reads, reference, outdir, genusname):
    print(f'mapping reads {reads} to {reference}')
    os.system(f'minimap2 -ax map-ont {reference} {reads} > {outdir}/{genusname}_consensus.sam')
    os.system(f'samtools view -b -S {outdir}/{genusname}_consensus.sam | samtools sort > {outdir}/{genusname}_consensus_sort.bam')
    os.system(f'rm {outdir}/*consensus.sam')
    os.system(f'samtools consensus {outdir}/{genusname}_consensus_sort.bam > {outdir}/{genusname}_consensus.fasta')


def align_reads(indir, ref, outdir):
    os.system(f'mkdir -p {outdir}')
    for file in glob.glob(indir + "/*.fastq.gz"):
        genusname = os.path.basename(file).split(".")[0]
        catch_right_ref(genusname, ref, outdir)
        run_minimap2(
            reads=file,
            reference=outdir + genusname + '_reference.fasta',
            outdir=outdir,
            genusname=genusname)

def main():
    print(str(sys.argv[0]))
    align_reads(indir = sys.argv[1], ref = sys.argv[2], outdir = sys.argv[3])

if __name__ == "__main__":
    main()
