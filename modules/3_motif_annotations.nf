process MOTIF_ANNOTATIONS {
    label 'pyscenic_conda'
    cpus 1
    memory '4.GB'
    publishDir "${params.outdir}/step3_motif_annot", mode: 'copy'

    input:
    path fly_tbl
    path ortholog_tsv

    output:
    path 'motif_annotations.tbl', emit: tbl

    script:
    """
    3_create_motif_annotations.py \
        --fly_tbl ${fly_tbl} \
        --orthologs ${ortholog_tsv} \
        --output motif_annotations.tbl \
        --min_pident ${params.min_pident}
    """
}
