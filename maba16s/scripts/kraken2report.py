import pandas as pd
import sys


def read_kraken(reportfile, cutoff = 100):
    df = pd.read_csv(reportfile, sep = '\t', header = None)
    df.columns = ['percentage_hit', 'Num_reads', 'Num_reads_noHit', 'taxlvl', 'taxID', 'taxname']
    genus = df[df.taxlvl == 'G']
    abundantgenera = genus[genus.Num_reads > 100] 

    return abundantgenera


def write_kraken_report(df, output):
    df.to_csv(output, sep = '\t')

    
    






def main():
    report = read_kraken(sys.argv[1])
    write_kraken_report(report, sys.argv[2])


if __name__ == "__main__":
    main()


