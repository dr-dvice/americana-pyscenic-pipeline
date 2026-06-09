process CLUSTER_ANALYSIS {
    label 'pyscenic_conda'
    cpus 1
    memory '8.GB'
    publishDir "${params.outdir}/step9_cluster_analysis", mode: 'copy'

    input:
    path auc_loom
    path h5ad
    path cluster_names
    path ortholog_file

    output:
    path 'cluster_regulon_activity.csv', emit: activity_csv
    path 'top_regulons_per_cluster.txt', emit: top_regulons
    path '*.png', emit: heatmaps

    script:
    """
    9_analyze_cluster_regulons.py \
        --auc_loom ${auc_loom} \
        --h5ad ${h5ad} \
        --cluster_names ${cluster_names} \
        --orthologs ${ortholog_file} \
        --species_name "${params.species_name}" \
        --cluster_col ${params.cluster_col} \
        --cluster_name_col ${params.cluster_name_col} \
        --output_dir .
    """
}
