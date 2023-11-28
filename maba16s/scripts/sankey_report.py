import pandas as pd
import sys
from plotly import graph_objects as go 

def make_sankey(excelfile, outputfile):
    df = pd.read_excel(excelfile, index_col = 0)
    
    taxa = list(df.index) + list(df['blast_hit'])
    taxa = pd.Series(taxa).unique()
    taxa = pd.Series(taxa)


    taxonomy = {}
    for num, tax in taxa.items():
        taxonomy[tax] = num
    

    source = [taxonomy[x] for x in df.index]
    target = [taxonomy[x] for x in df['blast_hit']]
    value = df['num_reads']


    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=list(taxonomy.keys())
            ),
        
        link = dict(
            source = source, # indices correspond to labels, eg A1, A2, A1, B1, ...
            target = target,
            value = value
    ))
    ])
    fig.write_html(outputfile)

    


def main():
    excelfile = sys.argv[1]
    outputfile = sys.argv[2]
    make_sankey(excelfile, outputfile)
    
if __name__ == "__main__":
    main()