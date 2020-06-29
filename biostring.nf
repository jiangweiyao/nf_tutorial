#!/usr/bin/env nextflow

params.in = "$baseDir/fastq/*.fastq"
params.out = "$HOME/test"

fastq_files = Channel.fromPath(params.in, type: 'file')

process fastqc {
    
    conda 'bioconductor-biostrings'

    //publishDir params.out, mode: 'copy', overwrite: true

    input:
    file(fastq) from fastq_files
 
    output:
    stdout result

    """
    #!/usr/bin/env R --no-save
    
    library(Biostring)
    library(ShortRead)
    
    fq <- readFastq(${fastq})
    fq[1:5]
    """
}


/*
 * print the channel content
 */
result.subscribe { println it }

/*
process multiqc {

    conda 'multiqc'

    publishDir params.out, mode: 'copy', overwrite: true

    input:
    file '*' from qc_files.collect().ifEmpty([])

    output:
    path "multiqc_report.html" into records

    """
    multiqc .
    """
}
*/
