"""
BLASTn and Kraken2 Processing Script
====================================

Description:
------------
This script processes BLASTn output files and filtered Kraken2 reports to 
generate a summary table with key metrics for each taxon, including the number 
of reads assigned directly to the taxon, the total reads within its clade, and 
the top BLASTn hit details.

Workflow:
---------
1. Parse BLASTn output files to extract top hits based on bitscore and 
   sequence identity thresholds.
2. Parse the Kraken2 filtered report to retrieve read counts for each taxon.
3. Combine BLASTn data and Kraken2 read counts into a final summary table.
4. Sort taxa by the number of reads assigned directly and export the table as 
   an Excel file.

Arguments:
----------
1. Input directory containing BLASTn output files (`*_BLASTn.txt`).
2. Path to the Kraken2 filtered report file.
3. Path to the output Excel file.

Output:
-------
An Excel file containing a summary table with the following columns:
- percentage: Maximum sequence identity percentage from BLASTn hits.
- length: Maximum alignment length from BLASTn hits.
- bitscore: Maximum bitscore from BLASTn hits.
- blast_hit: Top BLASTn hit taxonomy.
- reads_taxon: Number of reads assigned directly to the taxon.
- reads_clade: Total reads within the taxon's clade.
- num_reads: Same as reads_taxon (for compatibility).
"""

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
    print(f'Processing BLASTn hits from {blastfile}')
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
        results[genus] = blast_data

    results = pd.DataFrame(results).T

    # Retrieve the read counts from the filtered Kraken2 report
    print(f'Collecting read counts from {readfile}')
    readcount = get_read_count(readfile)
    results['num_reads'] = [readcount[x]['reads_taxon'] for x in results.index]

    # Sort and save the resulting table
    results = results.sort_values('num_reads', ascending = False)
    results.to_excel(outxls)

    print(f'Report generated')
    print(results)

if __name__ == "__main__":
    main()
