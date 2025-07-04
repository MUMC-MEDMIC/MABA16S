#!/usr/bin/env python3

import sys
from argparse import ArgumentParser
import yaml
import os
import subprocess

locationrepo = os.path.dirname(os.path.abspath(__file__))


def get_absolute_path(path):
    return os.path.abspath(path)


def file_name_generator(filepath):
    return os.path.splitext(os.path.basename(filepath))[0]


def snakemake_in(samples,  outdir, ):
    samplesdic = {}
    samplesdic['parameters'] = {}
    samplesdic['parameters']["outdir"] = get_absolute_path(outdir)
    samplesdic["SAMPLES"] = {}

    # generate the samples dictionary as input for snakemake
    for i in samples:
        samplename = file_name_generator(i)
        samplesdic["SAMPLES"][samplename] = get_absolute_path(i)
    data = yaml.dump(samplesdic, default_flow_style=False)

    # make and write config file location
    os.system(f"mkdir -p {locationrepo}/config")
    with open(f"{locationrepo}/config/config.yaml", 'w') as f:
        f.write(data)


def main(command_line=None):
    """Console script for MABA16S."""
    # add main parser object
    parser = ArgumentParser(description="Maastricht Bacterial 16S workflow")

    # add sub parser object
    subparsers = parser.add_subparsers(dest="mode")

    # add module rename barcodes to samples
    rename = subparsers.add_parser("rename",
                                   help='rename barcode.fasta to samplenames supplied in a spreadsheet')

    rename.add_argument("-i", 
                        required=True, 
                        dest='input_directory',
                        help='supply a raw data directory')
    
    rename.add_argument("--samplesheet",
                        required=True,
                        dest='samplesheet',
                        help='supply an excel sample sheet to rename samples')
    
    rename.add_argument("--technician",
                        required=False,
                        help='First email address contact to notify on failure (optional)')
    
    rename.add_argument("--mmb",
                        required=False,
                        help='Second email address contact to notify on failure (optional)')

    # add snakemake pipeline to completely run fasta to 16S report
    snakemake = subparsers.add_parser("snakemake",
                                      help='run the entire workflow on 16S sequenced samples')
    snakemake.add_argument("-i",
                           required=True,
                           dest="input_files",
                           nargs="+")
    
    snakemake.add_argument("--cores",
                           dest='cores',
                           required=True,
                           type=int,
                           help='Number of CPU cores to use')

    snakemake.add_argument("-o", required=True, dest="outdir")

    snakemake.add_argument("--snakemake-params",
                           required=False,
                           dest="smkparams")

####################
# parsing part
####################

    args = parser.parse_args(command_line)
    
    if args.mode == "rename":
        cmd = ["python3", "scripts/renamer.py",
               "--input-dir", args.input_directory,
               "--samplesheet", args.samplesheet]
        
        if args.technician:
            cmd += ["--technician", args.technician]
        if args.mmb:
            cmd += ["--mmb", args.mmb]
        
        subprocess.run(cmd, check=True)

    elif args.mode == "snakemake":
        snakemake_in(
                samples=args.input_files,
                outdir=args.outdir,
                )

        os.chdir(f"{locationrepo}")

        if args.smkparams == None:
            args.smkparams=""
        exitstatus = os.system(f"snakemake --cores {args.cores} --use-conda {args.smkparams} --snakefile Snakefile_check_suitable_samples.smk")
        if exitstatus > 0:
            sys.exit("pre-workflow crashed")
        os.system(f"snakemake -p --cores {args.cores} --use-conda {args.smkparams} --notemp --keep-going")

    else:
        parser.print_usage()
    return 0


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
