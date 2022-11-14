FROM snakemake/snakemake:v7.18.2


LABEL base.image="snakemake/snakemake:v7.18.2"
LABEL dockerfile.version="2"
LABEL software="MABA16S"
LABEL description="16S sequencing of clinical samples: The Pipeline"
LABEL website="https://github.com/MUMC-MEDMIC/MABA16S"
LABEL license="https://github.com/MUMC-MEDMIC/MABA16S/blob/main/LICENSE"
LABEL maintainer="Casper Jamin"
LABEL maintainer.email="casperjamin@gmail.com"


# get Pipeline
RUN git clone https://github.com/MUMC-MEDMIC/MABA16S
WORKDIR MABA16S/maba16s

# install dependencies via conda
RUN snakemake --use-conda --conda-create-envs-only --cores 1 
