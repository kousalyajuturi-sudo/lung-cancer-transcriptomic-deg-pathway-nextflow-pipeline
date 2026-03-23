nextflow.enable.dsl=2

params.reads = "raw_data/*.fastq.gz"
params.index = "/home/dell/scRNA_project/reference/grch38_index"
params.gtf   = "/home/dell/scRNA_project/reference/gencode.v45.annotation.gtf"
params.meta  = "samples.tsv"

workflow {

    reads_ch = Channel.fromPath(params.reads)
     
    PRE_QC(reads_ch)
    trimmed_ch = TRIM(reads_ch)
    POST_QC(trimmed_ch)
    
    aligned_ch = ALIGN(trimmed_ch)
    sorted_ch  = SORT_BAM(aligned_ch)
    FLAGSTAT(sorted_ch)

    sorted_ch.collect().set{all_bams} 
    counts_ch= COUNT_MATRIX(all_bams)

    DESEQ2_OUT= DESEQ2(
        counts_ch,
        file(params.meta),
        file("scripts/deseq_scRNA.R")
    )

   // Add KEGG links process
   KEGG_OUT = KEGG_LINKS(
    DESEQ2_OUT,
    file("scripts/generate_kegg_links.R")
)

   FINAL_REPORT(
    file("scripts/final_report_scRNA.py"),
    DESEQ2_OUT
)
}


process PRE_QC {

    publishDir "results/pre_qc", mode: 'copy'

    input:
    path reads

    output:
    path "*_fastqc.*"

    script:
    """
    fastqc $reads
    """
}


process TRIM {

    publishDir "results/trimmed", mode: 'copy'

    input:
    path reads

    output:
    path "${reads.simpleName}_trimmed.fastq.gz"

    script:
    """
    fastp -i $reads -o ${reads.simpleName}_trimmed.fastq.gz --thread 2
    """
}


process POST_QC {

    publishDir "results/post_qc", mode: 'copy'

    input:
    path reads

    output:
    path "*_fastqc.*"

    script:
    """
    fastqc $reads
    """
}


process ALIGN {

    publishDir "results/aligned", mode: 'copy'

    input:
    path reads

    output:
    path "${reads.simpleName}.sam"

    script:
    """
    hisat2 -x ${params.index} -U $reads -S ${reads.simpleName}.sam -p 2
    """
}


process SORT_BAM {

    publishDir "results/sorted_bam", mode: 'copy'

    input:
    path sam

    output:
    path "${sam.simpleName}.bam"

    script:
    """
    samtools view -bS $sam | samtools sort -o ${sam.simpleName}.bam
    samtools index ${sam.simpleName}.bam
    """
}


process FLAGSTAT {

    publishDir "results/flagstat", mode: 'copy'

    input:
    path bam

    output:
    path "${bam.simpleName}_flagstat.txt"

    script:
    """
    samtools flagstat $bam > ${bam.simpleName}_flagstat.txt
    """
}


process COUNT_MATRIX {

    publishDir "results/counts", mode: 'copy'

    input:
    path bam_files

    output:
    path "counts.txt"

    script:
    """
    featureCounts -a ${params.gtf} -o counts.txt *.bam
    """
}
process DESEQ2 {

    publishDir "results/deseq", mode: 'copy'

    input:
    path count_file
    path metadata
    path rscript

    output:
    path "results/plots"

    script:
    """
    Rscript $rscript $count_file $metadata
    """
}
process KEGG_LINKS {
    publishDir "results/kegg", mode: 'copy'

    input:
    path deseq_results
    path rscript

    output:
    path "*.png"
    path "*.xml"
    path "kegg_links.csv"

    script:
    """
    Rscript $rscript $deseq_results
    """
}
process FINAL_REPORT {

    publishDir "results/final_report", mode: 'copy'

    input:
    path report_script
    path plot_dir

    output:
    path "scRNAseq_Report.html"

    script:
    """
    python3 ${report_script} ${plot_dir}
    """
}
