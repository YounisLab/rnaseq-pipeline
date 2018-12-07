#!/usr/bin/env bash

set -e -u

usage()
{
    echo "Usage: ./run-pipeline.sh -R <REF_DIR> -F <FASTQ_FILE> -I <STAR_INDEX_DIR> -G <GENOME_VERSION> -P <NUM_CORES> -O <OUTPUT_DIR>

                -R <REF_DIR>         Directory containing reference files
                -F <FASTQ_FILE>      Path to .fastq file to run through the pipeline
                -I <STAR_INDEX_DIR>  Directory containing STAR indices
                -G <GENOME_VERSION>  Human Genome version prefix used in REF_DIR files
                -P <NUM_CORES>       Number of CPU cores to use in pipeline.
                -O <OUTPUT_DIR>      Directory to save output from pipeline.
                -S <SAMPLE_NAME>     Prefix name for output files.
            "
}

OPTIONS=":R:F:I:G:P:O:S:h"

REF_DIR=''
FASTQ_FILE=''
STAR_INDEX_DIR=''
GENOME_VERSION=''
NUM_CORES=''
OUTPUT_DIR=''
SAMPLE_NAME=''

while getopts $OPTIONS opt; do
    case $opt in
    R)
        REF_DIR=$OPTARG ;;
    F)
        FASTQ_FILE=$OPTARG ;;
    I)
        STAR_INDEX_DIR=$OPTARG ;;
    G)
        GENOME_VERSION=$OPTARG ;;
    P)
        NUM_CORES=$OPTARG ;;
    O)
        OUTPUT_DIR=$OPTARG ;;
    S)
        SAMPLE_NAME=$OPTARG ;;
    h)
        usage
        exit 0 ;;
    \?)
        echo "Unrecognized option: -$OPTARG"
        usage
        exit 1 ;;
    :)
        echo "Error: -$OPTARG requires an argument."
        usage
        exit 1 ;;
    esac
done

shift $(expr $OPTIND - 1 )
if [[ $# -gt 0 ]]
then
    echo "Unrecognized options: $*"
    usage
    exit 1
fi

if [[ -z $REF_DIR ]] || [[ -z $FASTQ_FILE ]] || [[ -z $STAR_INDEX_DIR ]] || \
[[ -z $GENOME_VERSION ]] || [[ -z $NUM_CORES ]] || [[ -z $OUTPUT_DIR ]] ||
[[ -z $SAMPLE_NAME ]]
then
    usage
    exit 1
fi

echo "REF_DIR = $REF_DIR"
echo "FASTQ_FILE = $FASTQ_FILE"
echo "STAR_INDEX_DIR = $STAR_INDEX_DIR"
echo "GENOME_VERSION = $GENOME_VERSION"
echo "NUM_CORES = $NUM_CORES"
echo "OUTPUT_DIR = $OUTPUT_DIR"
echo "SAMPLE_NAME = $SAMPLE_NAME"

ref_gene="$REF_DIR/$GENOME_VERSION.refGene_gene_longest.gtf"
ref_fasta="$REF_DIR/$GENOME_VERSION.fa"

cmd="nextflow run ./src/rnaseq-pipeline.nf --fastq '$FASTQ_FILE' --STAR_index $STAR_INDEX_DIR --ref_dir $REF_DIR --ref_gene $ref_gene --ref_fasta $ref_fasta --output_dir $OUTPUT_DIR --cores $NUM_CORES --sample_name $SAMPLE_NAME"
$cmd
