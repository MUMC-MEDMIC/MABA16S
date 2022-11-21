import sys
import pandas as pd
import glob
import os

percentage_from_max_bitscore = 0.90
percentage_from_max_seqid = 0.35

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
    data = {
        "percentage":max(df.percentage),
        "length":max(df.length),
        "bitscore":max(df.bitscore),
        "blast_hit":df.taxonomy.value_counts()
    }
    return data

def get_species(value):
    species = value.split(";")[-1]
    return species

def get_read_count(file):
    data = pd.read_csv(file, sep="\t")
    data.columns = ["percentage", "reads_in","reads_out", "taxlvl", "taxid", "name"]
    data["name"] = data.name.apply(lambda x: x.strip())
    data = data.set_index("name")
    read_dict = data.to_dict()['reads_out']

    return read_dict


def main():
    blastindir = sys.argv[1]
    outxls = sys.argv[2]
    readfile = sys.argv[3]

    results = {}
    for file in glob.glob(blastindir + "/*consensus.txt"):
        data = read_blastn(file)
        genus = os.path.basename(file).split("_consensus.txt")[0]
        results[genus] = data

    results = pd.DataFrame(results).T
    readcount = get_read_count(readfile)
    results['num_reads'] = [readcount[x] for x in results.index]
    results = results.sort_values('num_reads', ascending = False)
    results.to_excel(outxls)



if __name__ == "__main__":
    main()
