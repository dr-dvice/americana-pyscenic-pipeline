process AUCELL {
    label 'singularity_pyscenic'
    cpus params.max_threads
    memory '32.GB'
    publishDir "${params.outdir}/step8_aucell", mode: 'copy'

    input:
    path expression_loom
    path regulons_csv

    output:
    path 'auc_matrix.loom', emit: auc_loom

    script:
    """
    pyscenic aucell \
        ${expression_loom} \
        ${regulons_csv} \
        --output auc_matrix.loom \
        --num_workers ${task.cpus}
    """
}
