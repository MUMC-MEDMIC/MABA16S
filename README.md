# MABA16S

A pipeline for analyzing Oxford Nanopore 16S rRNA sequencing data from clinical samples.

* Free software: MIT license

### Overview
----------

MABA16S processes 16S rRNA sequencing data to classify reads at the genus and species levels, 
providing detailed taxonomic identification for clinical microbiology. 

### Tools and Settings
----------

1. Filtlong:
   - Minimum read length: 1200 bp
2. Kraken2:
    - Genera with a minimum of 50 reads are processed further
3. extract_kraken_reads.py from kraken-tools
4. Minimap2: alignment to SILVA reference sequences for specific genera
5. Samtools consensus
6. BLASTn

### How does it work?
-----------------
1. reads are classified on genus level using kraken2 and SILVA database
2. reads for each genus are extracted
3. each genus readset is mapped to the first species in the SILVA database of this genus
4. consensus sequence is extracted and BLASTed to the SILVA database to obtain a species ID
5. results are compiled and written to a spreadsheet

### Quickstart
----------
As a quickstart to use this pipeline you need Python 3.6 or higher, conda environment manager  and snakemake.

__Usage__

    git clone https://github.com/MUMC-MEDMIC/MABA16S 
    cd MABA16S/maba16s
    python cli.py snakemake -i folders_containing_nanopore16s_reads -o my_output_directory --cores 1 

    # input are directories which hold your nanopore reads. Naming of the output will be done based on the names of these directories


Output File Structure
-----------------
The output directory contains the following structure:  
```
my_output_directory/ 
├── kraken2/  
│   ├── {sample}/  
│   │   ├── krakenreport_filtered.txt  # Filtered Kraken2 report  
│   │   ├── output.txt                # Full Kraken2 classification output  
│   │   ├── reads/                    # Genus-specific reads (FASTQ files)  
├── kraken2consensus/  
│   ├── {sample}/  
│   │   ├── reference_fastas/         # Reference FASTA files used for alignment  
│   │   ├── aligned_reads/            # BAM files for aligned reads  
│   │   ├── consensus_fastas/         # Consensus FASTA files  
├── BLAST/  
│   ├── {sample}/  
│   │   ├── *_BLASTn.txt              # BLASTn results for consensus sequences  
├── QC/  
│   ├── {sample}/  
│   │   ├── {sample}_qcPreprocessing.txt  # Preprocessing QC metrics  
│   │   ├── {sample}_qcPostAnalysis.txt   # Extended QC metrics  
├── reports/  
│   ├── {sample}.xlsx                 # Comprehensive report for each sample  
├── sankeys/  
    ├── {sample}_sankey.html          # Interactive Sankey diagrams
```



Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
