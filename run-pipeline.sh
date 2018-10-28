#!/usr/bin/env bash

set -e -u

if [ "$#" -ne "6" ]
then
    echo "Usage: run-pipeline.sh <REF_DIR> <FASTQ_DIR> <STAR_INDEX_DIR> <GENOME_VERSION> <NUM_CORES> <OUTPUT_DIR>

        REF_DIR         Directory containing reference files
        FASTQ_FILE      Path to .fastq file to run through the pipeline
        STAR_INDEX_DIR  Directory containing STAR indices
        GENOME_VERSION  Human Genome version prefix used in REF_DIR files
        NUM_CORES       Number of CPU cores to use in pipeline.
        OUTPUT_DIR      Directory to save output from pipeline.
    "
    exit
fi

ref_gene="$1/$4.refGene_gene_longest.gtf"
ref_fasta="$1/$4.fa"

cmd="nextflow run ./src/rnaseq-pipeline.nf --fastq $2 --STAR_index $3 --ref_dir $1 --ref_gene $ref_gene --ref_fasta $ref_fasta --output_dir $6"
echo $cmd
$cmd
