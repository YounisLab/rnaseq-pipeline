# rnaseq-pipeline

## Prerequisites

If using Docker, the only prerequisite is to install [docker CE](https://docs.docker.com/install/linux/docker-ce/ubuntu/).

If using outside of Docker, the requirements can be installed using

```
conda install -c bioconda --yes --file requirements.txt
```

## Building

```
docker image build -t rnaseq-pipe .
```

## Running

`rnaseq-pipeline` can be run with and without Docker, though the latter is recommended.

### Without Docker:

```
run-pipeline.sh <REF_DIR> <FASTQ_DIR> <STAR_INDEX_DIR> <GENOME_VERSION> <NUM_CORES> <OUTPUT_DIR>

        REF_DIR         Directory containing reference files
        FASTQ_FILE      Path to .fastq files to run through the pipeline
        STAR_INDEX_DIR  Directory containing STAR indices
        GENOME_VERSION  Human Genome version prefix used in REF_DIR files
        NUM_CORES       Number of CPU cores to use in pipeline.
        OUTPUT_DIR      Directory to save output from pipeline.

```


### With Docker:

## Development using Docker

```
docker run -it --rm -v $PWD:/home/ rnaseq-pipe
```

This maps the current directory to a folder inside the container. This makes it
so that there is no need to re-build the image upon making source code changes.
