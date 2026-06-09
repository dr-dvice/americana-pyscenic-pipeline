process GET_FASTA {
    label 'pyscenic_conda'
    cpus 1
    memory '4.GB'
    publishDir "${params.outdir}/step2_fasta", mode: 'copy'

    input:
    path bed
    path genome_fasta
    path genome_fai

    output:
    path 'regulatory_regions.fa', emit: fasta

    script:
    """
    bedtools getfasta -fi ${genome_fasta} -bed ${bed} -s -name -fo raw.fa
    2_clean_fasta_headers.py raw.fa regulatory_regions.fa
    """
}
