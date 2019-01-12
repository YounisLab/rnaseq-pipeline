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
Usage: nextflow run rnaseq-pipeline.nf [OPTIONS] --ref_dir <REF_DIR> --fastq_dir <FASTQ_DIR> --star_index <STAR_INDEX_DIR> --genome <GENOME_VERSION> --cores <NUM_CORES> --output_dir <OUTPUT_DIR>
    Optional args:
        --single_end      Specifies that the --fastq_dir input contains single end reads only.
                            If enabled, files in the --fastq_files directory must be of the form
                            '<SAMPLENAME>_[A-Z].fastq', where [A-Z] refers to statistical replicates.

                            Example:
                                SAMPLE_A.fastq SAMPLE_B.fastq SAMPLE_C.fastq SAMPLE_D.fastq

                                This will be run in STAR as:
                                --readFilesIn SAMPLE_A.fastq,SAMPLE_B.fastq,SAMPLE_C.fastq,SAMPLE_D.fastq

        --no_replicates   Specifies that each .fastq file in the input directory
                            should be run through the pipeline individually.

    Mandatory args:
        <REF_DIR>         Directory containing reference files
        <FASTQ_DIR>       Directory containing .fastq file(s) to run through the pipeline.
                            Files in this directory must be of the form '<SAMPLENAME>_[A-Z]_R{1,2}.fastq',
                            where [A-Z] refers to statistical replicates and R{1,2} refers to paired-ends.

                            Example:
                                SAMPLE_A_R1.fastq SAMPLE_A_R2.fastq SAMPLE_B_R1.fastq SAMPLE_B_R2.fastq

                                This will be run in STAR as:
                                --readFilesIn SAMPLE_A_R1,SAMPLE_B_R1.fastq SAMPLE_A_R2.fastq,SAMPLE_B_R2.fastq

                            Here 'SAMPLE_A_R1.fastq' and 'SAMPLE_A_R2.fastq' are paired end reads, while
                            'SAMPLE_A_R1.fastq' and 'SAMPLE_B_R1.fastq' are statistical replicates
                            belonging to the same paired end read.

        <STAR_INDEX_DIR>  Directory containing STAR indices
        <GENOME_VERSION>  Human Genome version prefix used in REF_DIR files
        <NUM_CORES>       Number of CPU cores to use in pipeline.
        <OUTPUT_DIR>      Directory to save output from pipeline.
    """
}

if (params.help) {
    helpMessage()
    exit 0
}

if (!params.ref_dir) {
    helpMessage()
    exit 1, "REF_DIR not specified."
}

if (!params.fastq_dir) {
    helpMessage()
    exit 1, "FASTQ_DIR not specified."
}

if (!params.star_index) {
    helpMessage()
    exit 1, "STAR_INDEX_DIR not specified."
}

if (!params.genome) {
    helpMessage()
    exit 1, "GENOME_VERSION not specified."
}

if (!params.cores) {
    helpMessage()
    exit 1, "NUM_CORES not specified."
}

if (!params.output_dir) {
    helpMessage()
    exit 1, "OUTPUT_DIR not specified."
}

Channel
    .fromFilePairs(params.fastq_dir +'/*_[A-Z]_R{1,2}.fastq', size: -1)
    .map { key, files -> [key,
                          files.findAll { it.toString().endsWith('R1.fastq') },
                          files.findAll { it.toString().endsWith('R2.fastq') }] }
    .ifEmpty { exit 1, "Cannot find any reads matching the glob!"}
    .set {raw_reads_fastq}

 process STAR {
        input:
        set val(sample), file(reads1), file(reads2) from raw_reads_fastq
        file star_index from Channel.fromPath(params.star_index)

        publishDir "${params.output_dir}/$sample/STAR", mode: 'copy'

        output:
        set val(sample), file("${sample}_Aligned.sortedByCoord.out.bam") into \
            bam_for_regtools, bam_for_stringtie, bam_for_intron_analysis
        file '*' into STAR_DIR // Publish all files

        script:
        """
        STR1="$reads1"
        STR2="$reads2"
        READS1=\$(echo \${STR1// /,})
        READS2=\$(echo \${STR2// /,})

        STAR --genomeDir $star_index --runThreadN $params.cores --readFilesIn \$READS1 \$READS2 --outFileNamePrefix ${sample}_ \
        --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif
        """
}

process regtools {
    input:
    file ref_gene_bed from Channel.fromPath(params.ref_dir + "/" + params.genome + ".refGene_gene_longest.bed")
    set val(sample), file(bam_file) from bam_for_regtools

    publishDir "${params.output_dir}/$sample/regtools", mode: 'copy'

    output:
    file "${sample}_clean.bed" into bed_for_intron_analysis

    script:
    """
    samtools index $bam_file
    regtools junctions extract $bam_file -o ${sample}.bed
    remove_transgene.py $ref_gene_bed ${sample}.bed ${sample}_clean.bed
    """
}

process stringtie {
    input:
    set val(sample), file(bam_file) from bam_for_stringtie
    file ref_gene_gtf from Channel.fromPath(params.ref_dir + "/" + params.genome + ".refGene_gene_longest.gtf")

    publishDir "${params.output_dir}/$sample/stringtie", mode: 'copy'

    output:
    file "${sample}_gene_abund.tab" into fpkm_for_intron_analysis

    script:
    """
    stringtie -e -p 32 -A ${sample}_gene_abund.tab -G $ref_gene_gtf -o ${sample}_assembly.gtf $bam_file
    """
}

process intron_analysis {
    input:
    set val(sample), file(bam_file) from bam_for_intron_analysis
    file fpkm from fpkm_for_intron_analysis
    file junc_bed from bed_for_intron_analysis
    file ref_dir from Channel.fromPath(params.ref_dir)

    publishDir "${params.output_dir}/$sample/intron_analysis", mode: 'copy'

    output:
    set file("${sample}_intron_analysis.txt"), \
        file("${sample}_total_cvg.txt") into INTRON_ANALYSIS_DIR

    script:
    """
    echo "===> Computing coverage..."
    compute_coverage.sh $ref_dir $bam_file $junc_bed $params.genome $sample
    echo "===> Performing analysis..."
    analyze.py $ref_dir $fpkm $params.genome $sample
    """
}
