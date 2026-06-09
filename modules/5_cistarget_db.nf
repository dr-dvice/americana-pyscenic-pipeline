process CISTARGET_SINGLE {
    label 'cistarget_conda'
    cpus params.max_threads
    memory '32.GB'
    publishDir "${params.outdir}/step5_cistarget", mode: 'copy'

    input:
    path fasta
    path motifs_dir
    path motif_list
    path cbust_binary
    path cistarget_script
    val db_prefix
    val seed

    output:
    path '*.regions_vs_motifs.rankings.feather', emit: rankings

    script:
    """
    python ${cistarget_script} \
        -f ${fasta} \
        -M ${motifs_dir} \
        -m ${motif_list} \
        -o ${db_prefix} \
        -c ${cbust_binary} \
        -t ${task.cpus} \
        -s ${seed}
    """
}

process CISTARGET_PART {
    label 'cistarget_conda'
    cpus params.max_threads
    memory '32.GB'
    time '4.h'
    tag "part_${part_num}"

    input:
    path fasta
    path motifs_dir
    path motif_list
    path cbust_binary
    path cistarget_script
    val db_prefix
    val seed
    val total_parts
    val part_num

    output:
    path '*.scores.feather', emit: partial_scores

    script:
    """
    python ${cistarget_script} \
        -f ${fasta} \
        -M ${motifs_dir} \
        -m ${motif_list} \
        -o ${db_prefix} \
        -c ${cbust_binary} \
        -t ${task.cpus} \
        -s ${seed} \
        -p ${part_num} ${total_parts}
    """
}

process CISTARGET_COMBINE {
    label 'pyscenic_conda'
    cpus 1
    memory '16.GB'
    publishDir "${params.outdir}/step5_cistarget", mode: 'copy'

    input:
    path partial_scores
    val db_prefix
    val total_parts
    val seed

    output:
    path '*.regions_vs_motifs.rankings.feather', emit: rankings

    script:
    """
    5_combine_partial_dbs.py ${db_prefix} ${total_parts} ${seed}
    """
}
