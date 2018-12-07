# rnaseq-pipeline

A nextflow pipeline for aligning, mapping & estimating abundances of RNASeq data (.bam files). The following are the components of the pipeline:

1. STAR
2. regtools
3. Cufflinks

To get started, download the repo using:

`git clone --recursive https://github.com/YounisLab/rnaseq-pipeline`

See the 'Prerequisites' & 'Running' section for install & usage instructions.

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
Usage: ./run-pipeline.sh -R <REF_DIR> -F <FASTQ_FILE> -I <STAR_INDEX_DIR> -G <GENOME_VERSION> -P <NUM_CORES> -O <OUTPUT_DIR> -S <SAMPLE_NAME>

                -R <REF_DIR>         Directory containing reference files
                -F <FASTQ_FILE>      Path to .fastq file(s) to run through the pipeline.
                                     Separate paired-ends with spaces & replicates by comma.
                                     Eg: -F \"pair1_sample1,pair1_sample2...pair1_sampleN pair2_sample1,pair2_sample2...pair2_sampleN\"
                                     NOTE: MUST BE ENCLOSED IN QUOTES FOR PAIRED-END/REPLICATES.
                -I <STAR_INDEX_DIR>  Directory containing STAR indices
                -G <GENOME_VERSION>  Human Genome version prefix used in REF_DIR files
                -P <NUM_CORES>       Number of CPU cores to use in pipeline.
                -O <OUTPUT_DIR>      Directory to save output from pipeline.
                -S <SAMPLE_NAME>     Prefix name for output files.

```


### With Docker:

## Development using Docker

```
docker run -it --rm -v $PWD:/home/ rnaseq-pipe
```

This maps the current directory to a folder inside the container. This makes it
so that there is no need to re-build the image upon making source code changes.
