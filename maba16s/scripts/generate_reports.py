import sys
import pandas as pd
import glob

percentage_from_max_bitscore = 0.90
percentage_from_max_seqid = 1

def read_blastn(file):
    df = pd.read_csv(file, sep="\t")
    df.columns = ['percentage', 'length', 'bitscore', 'taxonomy']
    df = score_top_hits(df)
    return df


def score_top_hits(df):
    # cut off by taking top bitscore hits
    cutoff = max(df.bitscore) * percentage_from_max_bitscore
    df = df[df.bitscore > cutoff]

    # cut off from max seq alignment
    cutoff = max(df.percentage) - percentage_from_max_seqid
    df = df[df.percentage > cutoff]

    df.taxonomy = df.taxonomy.apply(get_species)
    print(str(df.taxonomy.value_counts()))
    return df

def get_species(value):
    species = value.split(";")[-1]
    return species


def main():
    blastindir = sys.argv[1]
    outdir = sys.argv[2]

    for file in glob.glob(blastindir + "/*consensus.txt"):
        read_blastn(file)



if __name__ == "__main__":
    main()
