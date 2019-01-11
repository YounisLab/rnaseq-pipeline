#!/usr/bin/env nextflow

// optional args
params.help = null
params.single_end = null
params.no_replicates = null

// manadatory args
params.ref_dir = null
params.fastq_dir = null
params.star_index = null
params.genome = null
params.cores = null
params.output_dir = null


def helpMessage() {
    log.info """
    Usage: ./run-pipeline.sh [OPTION] --ref_dir <REF_DIR> --fastq_dir <FASTQ_DIR> --star_index <STAR_INDEX_DIR> --genome <GENOME_VERSION> --cores <NUM_CORES> --output_dir <OUTPUT_DIR>
            Optional args:
                --single_end      Specifies that the --fastq_files input contained single end reads only.
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

                                    Here 'BL9_A_R1.fastq' and 'BL9_A_R2.fastq' are paired end reads, while
                                    'BL9_A_R1.fastq' and 'BL9_B_R1.fastq' are statistical replicates
                                    belonging to the same paired end read.

                <STAR_INDEX_DIR>  Directory containing STAR indices
                <GENOME_VERSION>  Human Genome version prefix used in REF_DIR files
                <NUM_CORES>       Number of CPU cores to use in pipeline.
                <OUTPUT_DIR>      Directory to save output from pipeline.
    """
}

Channel
    .fromFilePairs(params.fastq_dir +'/*[A-Z]_R{1,2}.fastq', size: -1)
    .map { key, files -> [key, files.findAll { it.toString().endsWith('R1.fastq') }, files.findAll { it.toString().endsWith('R2.fastq') }] }
    .ifEmpty { exit 1, "Cannot find any reads matching the glob!"}
    .into {raw_reads_fastq}

 process STAR {
        input:
        set val(sample), file(reads1), file(reads2) from raw_reads_fastq

        publishDir "${params.output_dir}/$sample/STAR", mode: 'copy'

        output:
        file '*' into STAR_DIR // For publishDir to publish everything

        script:
        """
        STR1="$reads1"
        STR2="$reads2"
        READS1=\$(echo \${STR1// /,})
        READS2=\$(echo \${STR2// /,})

        STAR --genomeDir $params.star_index --runThreadN $params.cores --readFilesIn \$READS1 \$READS2 --outFileNamePrefix ${sample}_ \
        --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif
        """
}
