#!/usr/bin/env bash

if [ "$#" -ne "3" ]
then
    echo "Not enough parameters..."
    echo "Usage: bam2bigwig.sh <accepted_hits.bam> <chrom.sizes.txt> <path_to_output_bigwig>"
    echo "Example: bam2bigwig.sh sample1.bam hg38.chrom.sizes ~/sample1/bigwig/sample1.bw"
    exit
fi

accepted_hits=$1
genome_file=$2

output_file_name=$(basename -- "$3")
output_file_name="${output_file_name%.*}"
output_path_without_ext=$(dirname -- "$3")/${output_file_name}

echo "Running genomeCoverageBed..."
genomeCoverageBed -split -bg -ibam $accepted_hits -g $genome_file > $output_file_name.interim.bedgraph

echo "Running sorting on bedgraph file..."
LC_COLLATE=C sort -k1,1 -k2,2n $output_file_name.interim.bedgraph > $output_path_without_ext.bedgraph

echo "Running bedGraphToBigWig"
bedGraphToBigWig $output_path_without_ext.bedgraph $genome_file $output_path_without_ext.bw

echo "Cleaning up"
rm $output_file_name.interim.bedgraph
rm $output_path_without_ext.bedgraph
