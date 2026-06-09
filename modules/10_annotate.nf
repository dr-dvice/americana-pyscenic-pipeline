process ANNOTATE {
    label 'pyscenic_conda'
    cpus 1
    memory '4.GB'
    publishDir "${params.outdir}/step10_annotate", mode: 'copy'

    input:
    path gff
    path ortholog_file
    path cluster_activity_csv
    path top_regulons_txt

    output:
    path 'regulon_annotations.csv',                emit: annotations
    path 'cluster_regulon_activity_annotated.csv',  emit: annotated_activity
    path 'top_regulons_per_cluster_annotated.txt',  emit: annotated_top
    path 'top20_regulons_overall.csv',              emit: top20

    script:
    """
    10_annotate_regulons.py \
        --gff ${gff} \
        --orthologs ${ortholog_file} \
        --cluster_activity ${cluster_activity_csv} \
        --top_regulons ${top_regulons_txt} \
        --output_dir .
    """
}
