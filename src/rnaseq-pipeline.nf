#!/usr/bin/env nextflow

/* This script aligns a fastq file to a reference genome with annotations using the STAR aligner
 * and then performs quantification using Cufflinks.
 *
 * The script is NOT to be executed directly, instead the provided
 * 'run-pipeline.sh' wrapper script calls this script with the correct parameters.
 * If using docker, use 'run-pipeline-docker.sh'.
 *
 */

// fastq file directory
params.fastq             = ""

// path to STAR indexes
params.STAR_index        = ""

// path to hg_ref directory
params.ref_dir           = ""

// path to GTF refgene file for cufflinks.
// The corresponding .bed file must also be present in ref_dir
params.ref_gene          = ""

// path to reference fasta file for cufflinks
params.ref_fasta         = ""

params.output_dir        = ""

params.gene_version = "hg38"
params.cores    = 32

/*
 * Programmer's note:
 * To publish the output of a pipeline component to an output directory, we use publishDir.
 * But publishDir only publishes files pushed into the 'output' channel. Hence the last line of
 * every output block. Any channel named *_DIR is used for this purpose.
 */

process STAR {
        input:
        file fastq_file from Channel.fromPath( params.fastq)

        publishDir "${params.output_dir}/${fastq_file.baseName}/STAR_OUT", mode: 'link'

        output:
        set val(fastq_file), file("${fastq_file.baseName}_Aligned.sortedByCoord.out.bam") into STAR_out_1
        set val(fastq_file), file("${fastq_file.baseName}_Aligned.sortedByCoord.out.bam") into STAR_out_2
        set val(fastq_file), file("${fastq_file.baseName}_Aligned.sortedByCoord.out.bam") into STAR_out_3
        file '*' into STAR_DIR

        """
        STAR --genomeDir $params.STAR_index --runThreadN $params.cores --readFilesIn $fastq_file --outFileNamePrefix ${fastq_file.baseName}_ \
        --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif
        """
}

