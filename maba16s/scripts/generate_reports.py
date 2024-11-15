import sys
import pandas as pd
import glob
import os

# functionality to process the BLASTn output
percentage_from_max_bitscore = 0.99
percentage_from_max_seqid = 0.35

def read_blastn(blastfile):
    # Handle missing BLASTn information: no match due to empty consensus or NNNNNNN as consensus
    if os.path.getsize(blastfile) == 0:
        print(f'No BLASTn hits found, returning NAs')
        return pd.DataFrame([{
            "percentage": np.nan,
            "length": np.nan,
            "bitscore": np.nan,
            "taxonomy": "No hits or no good consensus"}])

    # When BLASTn has hits        
    print(f'Processing BLASTn hits from {file}')
    blast_df = pd.read_csv(blastfile, sep="\t")
    blast_df.columns = ['percentage', 'length', 'bitscore', 'taxonomy']
    blast_df = score_top_hits(blast_df)
    return blast_df


def score_top_hits(blast_df):
    # filter using bistscore value cut-off
    bitscore_cutoff = max(blast_df.bitscore) * percentage_from_max_bitscore
    blast_df = blast_df[blast_df.bitscore > bitscore_cutoff]

    # filter based on sequence identity score
    seqid_cutoff = max(blast_df.percentage) - percentage_from_max_seqid
    blast_df = blast_df[blast_df.percentage > seqid_cutoff]

    # format the taxonomy information
    blast_df.taxonomy = blast_df.taxonomy.apply(get_species)
    blast_hits = blast_df.taxonomy.value_counts()
    data = {
        "percentage":max(blast_df.percentage),
        "length":max(blast_df.length),
        "bitscore":max(blast_df.bitscore),
        "blast_hits":blast_hits,
        "blast_hit":blast_hits.index[0]
    }
    return data

def get_species(value):
    # retrieve the last entry of the taxonomy information
    species = value.split(";")[-1]
    # keep only the first two words
    species = " ".join((species.split(" ")[:2]))
    return species

# Functionality to process the filtered Kraken2 report to retrieve read counts
def get_read_count(readfile):
    kraken_data = pd.read_csv(readfile, sep="\t", header = None)
    kraken_data.columns = ["percentage", "reads_clade","reads_taxon", "taxon_level", "taxon_id", "name"]
    kraken_data["name"] = kraken_data.name.apply(lambda x: x.strip())
    kraken_data = kraken_data.set_index("name")
    kraken_read_dict = kraken_data[["reads_taxon", "reads_clade"]].to_dict(orient="index")

    return kraken_read_dict

# Main function
def main():
    blastindir = sys.argv[1]
    outxls = sys.argv[2]
    readfile = sys.argv[3]

    results = {}

    # Process the BLASTn output
    print(f'Collecting BLASTn reports from {blastindir}')
    
    for blastfile in glob.glob(os.path.join(blastindir, "*_BLASTn.txt")):
        blast_data = read_blastn(blastfile)
        genus = os.path.basename(blastfile).replace("_BLASTn.txt", "")
        print(f'Processing genus {genus}')
        results[genus] = blast_data

    results = pd.DataFrame(results).T

    # Retrieve the read counts from the filtered Kraken2 report
    print(f'Collecting read counts from {readfile}')
    readcount = get_read_count(readfile)
    results['num_reads'] = [readcount[x]['reads_taxon'] for x in results.index]
    results['reads_clade'] = [readcount[x]['reads_clade'] for x in results.index]

    # Sort and save the resulting table
    results = results.sort_values('num_reads', ascending = False)
    results.to_excel(outxls)


if __name__ == "__main__":
    main()
