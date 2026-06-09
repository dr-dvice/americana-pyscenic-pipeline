process GRN_LOCAL {
    label 'singularity_pyscenic'
    cpus params.max_threads
    memory '64.GB'
    time '48.h'
    publishDir "${params.outdir}/step4_grn", mode: 'copy'

    input:
    path expression_csv
    path tf_list

    output:
    path 'grn_network.tsv', emit: grn

    script:
    """
    # Extract gene IDs from TF list (handles both plain text and CSV with header)
    head -1 ${tf_list} | grep -q ',' && cut -d',' -f1 ${tf_list} | tail -n+2 > tf_genes.txt || cp ${tf_list} tf_genes.txt

    pyscenic grn \
        --method grnboost2 \
        --num_workers ${task.cpus} \
        --output grn_network.tsv \
        ${expression_csv} \
        tf_genes.txt
    """
}
