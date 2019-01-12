# rnaseq-pipeline

A nextflow pipeline for aligning, mapping & estimating abundances of RNASeq data (.bam files). The following are the components of the pipeline:

1. STAR
2. regtools
3. stringtie
4. [intron_analysis](https://github.com/YounisLab/splicing-analysis)

To get started, download the repo using:

`git clone --recursive https://github.com/YounisLab/rnaseq-pipeline`

See the 'Prerequisites' & 'Running' section for install & usage instructions.

## Prerequisites

Install [docker CE](https://docs.docker.com/install/linux/docker-ce/ubuntu/).

## Building

```
git clone https://github.com/YounisLab/splicing-analysis
cd splicing-analysis
docker image build -t splicing-analysis .
cd ../
git clone https://github.com/YounisLab/rnaseq-pipeline
cd rnaseq-pipeline
docker image build -t rnaseq-pipeline .
```

## Running

```
Usage: nextflow run rnaseq-pipeline.nf [OPTIONS] --ref_dir <REF_DIR> --fastq_dir <FASTQ_DIR> --star_index <STAR_INDEX_DIR> --genome <GENOME_VERSION> --cores <NUM_CORES> --output_dir <OUTPUT_DIR>
        Optional args:
            --single_end      Specifies that the --fastq_dir input contains single end reads only.
                              If enabled, files in the --fastq_files directory must be of the form
                              '<SAMPLENAME>_[A-Z].fastq', where [A-Z] refers to statistical replicates.

                              Example:
                                    SAMPLE_A.fastq SAMPLE_B.fastq SAMPLE_C.fastq SAMPLE_D.fastq

            --no_replicates   Specifies that each .fastq file in the input directory
                              should be run through the pipeline individually.

        Mandatory args:
            <REF_DIR>         Directory containing reference files
            <FASTQ_DIR>       Directory containing .fastq file(s) to run through the pipeline.
                              Files in this directory must be of the form '<SAMPLENAME>_[A-Z]_R{1,2}.fastq',
                              where [A-Z] refers to statistical replicates and R{1,2} refers to paired-ends.

                              Example:
                                    SAMPLE_A_R1.fastq SAMPLE_A_R2.fastq SAMPLE_B_R1.fastq SAMPLE_B_R2.fastq

                              Here 'SAMPLE_A_R1.fastq' and 'SAMPLE_A_R2.fastq' are paired end reads, while
                              'SAMPLE_A_R1.fastq' and 'SAMPLE_B_R1.fastq' are statistical replicates
                              belonging to the same paired end read.

            <STAR_INDEX_DIR>  Directory containing STAR indices
            <GENOME_VERSION>  Human Genome version prefix used in REF_DIR files
            <NUM_CORES>       Number of CPU cores to use in pipeline.
            <OUTPUT_DIR>      Directory to save output from pipeline.
```
