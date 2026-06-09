process SPREADSHEET {
    label 'pyscenic_conda'
    cpus 1
    memory '4.GB'
    publishDir "${params.outdir}/step11_spreadsheet", mode: 'copy'

    input:
    path regulons_csv
    path gff_key_tsv
    path ortholog_file

    output:
    path 'regulon_summary.xlsx', emit: xlsx

    script:
    """
    11_create_regulon_spreadsheet.py \
        --regulons ${regulons_csv} \
        --gff_key ${gff_key_tsv} \
        --orthologs ${ortholog_file} \
        --output regulon_summary.xlsx
    """
}
