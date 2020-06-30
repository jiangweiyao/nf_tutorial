#!/usr/bin/env nextflow

params.in = "$baseDir/fastq/*.fastq.gz"
params.out = "$HOME/test"
params.reference = "$baseDir/Americas2.fasta"
params.model = "r941_min_high_g303"

ref = file(params.reference)
model = params.model
fastq_files = Channel.fromPath(params.in, type: 'file')
fastq_files2 = Channel.fromPath(params.in, type: 'file')

process fastqc {
    
    conda 'fastqc'

    //publishDir params.out, mode: 'copy', overwrite: true

    //Note to self: specifying the file name literally coerces the input file into that name. It doesn't select files matching pattern of the literal.
    input:
    file fastq from fastq_files
 
    output:
    file "*_fastqc.{zip,html}" into qc_files
    file "*_fastqc.{zip,html}" into qc_files1
    """
    fastqc ${fastq}
    """
}

//qc_files1.collect().print()

process multiqc {

    conda 'multiqc'

    publishDir params.out, mode: 'copy', overwrite: true

    input:
    file reports  from qc_files.collect().ifEmpty([])

    output:
    path "multiqc_report.html" into final1

    """
    multiqc .
    """
}

process minimap2 {

    conda 'minimap2'

//    publishDir params.out, mode: 'copy', overwrite: true

    input:
    file fastq from fastq_files2

    output:
    tuple file(fastq), path("${fastq.simpleName}.sam") into sam_files

    """
    minimap2 -ax map-ont $ref ${fastq} > ${fastq.simpleName}.sam
    """
}

process samtools {

    conda 'samtools'

//    publishDir params.out, mode: 'copy', overwrite: true

    input:
    tuple file(fastq), file(sam_file)from sam_files

    output:
    tuple file(fastq), path("${sam_file.simpleName}.coverage") into coverage_files

    """
    samtools view -S -b ${sam_file.simpleName}.sam > ${sam_file.simpleName}.bam
    samtools sort ${sam_file.simpleName}.bam -o ${sam_file.simpleName}.sorted.bam
    samtools index ${sam_file.simpleName}.sorted.bam
    samtools depth ${sam_file.simpleName}.sorted.bam > ${sam_file.simpleName}.coverage
    """
}

process scaffold {

    conda 'bioconductor-decipher'

    publishDir params.out, pattern: "${coverage_file.simpleName}_cov.txt", mode: 'copy', overwrite: true

    input:
    tuple file(fastq), file(coverage_file) from coverage_files

    output:
    tuple file(fastq), path("${coverage_file.simpleName}_scaffold.fasta") into scaffold_files
    path "${coverage_file.simpleName}_cov.txt" into cov_stat_files

    """
    $baseDir/scaffold_cutter.R 1 10 $ref ${coverage_file.simpleName}.coverage ${coverage_file.simpleName}_scaffold.fasta ${coverage_file.simpleName}_cov.txt   
    """
}

process medaka {

    conda 'medaka'

//    publishDir params.out, mode: 'copy', overwrite: true

    input:
    tuple file(fastq), file(scaffold_file) from scaffold_files

    output:
    tuple file(scaffold_file), path("${fastq.simpleName}_medaka/consensus.fasta") into medaka_files

    """
    medaka_consensus -i $fastq -d $scaffold_file -o ${fastq.simpleName}_medaka -m $model 
    """
}

process gapfixer {

    conda 'bioconductor-decipher'

    publishDir params.out, mode: 'copy', overwrite: true

    input:
    tuple file(scaffold_file), file(medaka_file) from medaka_files

    output:
    path("${scaffold_file.simpleName}_final.fasta") into final_files

    """
    $baseDir/gapfixer.R $scaffold_file $medaka_file ${scaffold_file.simpleName}_final.fasta
    """
}
