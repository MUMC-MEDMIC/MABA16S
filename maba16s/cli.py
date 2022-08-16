#!/usr/bin/env python3

import sys
from argparse import ArgumentParser
import yaml
import os
from scripts.renamer import renamer



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
    rename = subparsers.add_parser(
            "rename",
            help="""rename barcode.fasta to
            samplenames supplied in a spreadsheet""")

    rename.add_argument("-i",
                        required=True,
                        dest='input_directory',

                        )
    rename.add_argument("--spreadsheet",
                        required=True,
                        dest='spreadsheet',
                        help='supply a spreadsheet to rename samples'
                        )
                            
    # add snakemake pipeline to completely run fasta to 16S report
    snakemake = subparsers.add_parser("snakemake",
                                      help='''run the entire workflow on 16S
                                      sequenced samples''')
    snakemake.add_argument(
                "-i",
                required=True,
                dest="input_files",
                nargs="+"
                )
    snakemake.add_argument(
                            "--cores",
                            dest='cores',
                            required=True,
                            type=int,
                            help='Number of CPU cores to use'
                            )

    snakemake.add_argument("-o", required=True, dest="outdir")

####################
# parsing part
####################

    args = parser.parse_args(command_line)
    if args.mode == "rename":
        renamer(
                input_file=args.input_directory,
                spreadsheet=args.spreadsheet                
                )


    elif args.mode == "snakemake":
        snakemake_in(
                samples=args.input_files,
                outdir=args.outdir,
                )

        os.chdir(f"{locationrepo}")
        os.system(f"snakemake --cores {args.cores}")

    else:
        parser.print_usage()
    return 0


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
