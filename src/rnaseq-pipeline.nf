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

params.cores             = ""

params.gene_version = "hg38"

// shortcuts to external commands
remove_transgene = PWD + "/src/scripts/remove_transgene.py"
splicing_analysis = PWD + "/src/splicing-analysis/splicing-analysis.sh"

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

process regtools {
        input:
        file ref_gene_file from Channel.fromPath(params.ref_gene)
        set val(fastq_file), file(bam_file) from STAR_out_1

        publishDir "${params.output_dir}/${fastq_file.baseName}/", mode: 'link'

        output:
        file "${fastq_file.baseName}_clean.bed" into REGTOOLS_out_1

        """
        samtools index $bam_file
        regtools junctions extract $bam_file -o ${fastq_file.baseName}.bed
        python $remove_transgene $params.ref_dir/${ref_gene_file.baseName}.bed ${fastq_file.baseName}.bed ${fastq_file.baseName}_clean.bed
        """
}

process cufflinks {
        input:
        set val(fastq_file), file(bam_file) from STAR_out_2

        publishDir "${params.output_dir}/${fastq_file.baseName}/CUFFLINKS_OUT/", mode: 'link'

        output:
        file 'genes.fpkm_tracking' into CUFFLINKS_out_1
        file '*' into CUFFLINKS_DIR

        """
        cufflinks -p $params.cores -G $params.ref_gene -b $params.ref_fasta -L experiment_descriptor -u $bam_file
        """
        
}

process intron_analysis {
        input:
        set val(fastq_file), file(bam_file) from STAR_out_3
        file(fpkm_file) from CUFFLINKS_out_1 
        file(junctions_file) from REGTOOLS_out_1

        publishDir "${params.output_dir}/${fastq_file.baseName}/"

        output:
        file("${fastq_file.baseName}_intron_analysis.txt") into ANALYSIS_DIR_1
        file("${fastq_file.baseName}_total_cvg.txt") into ANALYSIS_DIR_2

        """
        $splicing_analysis $params.ref_dir $bam_file $junctions_file $fpkm_file $params.gene_version $fastq_file.baseName
        """
}
