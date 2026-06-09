process H5AD_TO_LOOM {
    label 'pyscenic_conda'
    cpus 1
    memory '16.GB'
    publishDir "${params.outdir}/step7_loom", mode: 'copy'

    input:
    path h5ad

    output:
    path 'expression.loom', emit: loom

    script:
    """
    7_convert_h5ad_to_loom.py --input ${h5ad} --output expression.loom
    """
}
