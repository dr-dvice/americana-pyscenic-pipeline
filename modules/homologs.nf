process CREATE_HOMOLOGS {
    label 'pyscenic_conda'
    cpus 1
    memory '4.GB'
    publishDir "${params.outdir}/homologs", mode: 'copy'

    input:
    path blast_map
    path sp_headers
    path fly_headers
    path gff_key

    output:
    path 'species_to_fly_homologs.tsv', emit: ortholog_tsv

    script:
    """
    create_species_to_fly_homologs.py \
        --blast_map ${blast_map} \
        --sp_headers ${sp_headers} \
        --fly_headers ${fly_headers} \
        --gff_key ${gff_key} \
        --output species_to_fly_homologs.tsv
    """
}
