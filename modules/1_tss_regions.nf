process TSS_REGIONS {
    label 'pyscenic_conda'
    cpus 1
    memory '2.GB'
    publishDir "${params.outdir}/step1_tss", mode: 'copy'

    input:
    path gff

    output:
    path 'regulatory_regions.bed', emit: bed

    script:
    """
    1_extract_tss_regions.py ${gff} regulatory_regions.bed \
        --upstream ${params.upstream_bp} \
        --downstream ${params.downstream_bp}
    """
}
