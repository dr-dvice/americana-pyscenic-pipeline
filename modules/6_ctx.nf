process CTX {
    label 'singularity_pyscenic'
    cpus params.max_threads
    memory '32.GB'
    publishDir "${params.outdir}/step6_ctx", mode: 'copy'

    input:
    path grn_tsv
    path rankings_feather
    path motif_tbl
    path expression_loom

    output:
    path 'regulons.csv', emit: regulons

    script:
    """
    pyscenic ctx \
        ${grn_tsv} \
        ${rankings_feather} \
        --annotations_fname ${motif_tbl} \
        --expression_mtx_fname ${expression_loom} \
        --output regulons.csv \
        --num_workers ${task.cpus} \
        --mask_dropouts
    """
}
