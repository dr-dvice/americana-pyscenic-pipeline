process EXPORT_COUNTS {
    label 'pyscenic_conda'
    cpus 1
    memory '16.GB'
    publishDir "${params.outdir}/step0_export_counts", mode: 'copy'

    input:
    path h5ad

    output:
    path 'raw_counts.csv', emit: counts_csv

    script:
    def use_raw = params.use_raw_counts ? '--use_raw' : ''
    """
    0_export_raw_counts.py --input ${h5ad} --output raw_counts.csv ${use_raw}
    """
}
