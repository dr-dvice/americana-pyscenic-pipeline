#!/usr/bin/env nextflow

/*
 * pySCENIC Cross-Species Regulon Analysis Pipeline
 *
 * 11-step pipeline: TSS extraction → cisTarget DB → GRNBoost2 → ctx → AUCell → annotation
 * Supports importing pre-computed GRN and cisTarget DB from remote machines.
 *
 * Usage:
 *   nextflow run main.nf -profile americana                    # import mode (typical)
 *   nextflow run main.nf -profile americana --grn_import null  # recompute GRN locally
 *   nextflow run main.nf -profile americana -resume            # resume after failure
 */

nextflow.enable.dsl = 2

// ── Parameter validation ────────────────────────────────────────────────────

def required_params = ['gff', 'genome_fasta', 'atlas_h5ad',
                       'cluster_names', 'gff_key', 'species_name']
required_params.each { p ->
    if (!params[p]) {
        error "Missing required parameter: --${p}. Use a species profile (-profile americana) or provide via CLI."
    }
}

if (!params.ortholog_file) {
    ['blast_map', 'sp_headers', 'fly_headers'].each { p ->
        if (!params[p]) {
            error "Ortholog generation requires --${p}. Provide BLAST inputs or use --ortholog_file to import a pre-computed file."
        }
    }
}

if (!params.grn_import && !params.tf_list) {
    error "Local GRN computation requires --tf_list. Provide a TF list or use --grn_import to skip."
}

if (!params.cistarget_import) {
    ['cbust_binary', 'cistarget_script', 'motif_list'].each { p ->
        if (!params[p]) {
            error "Local cisTarget DB build requires --${p}. Provide it or use --cistarget_import to skip."
        }
    }
}

// ── Module includes ─────────────────────────────────────────────────────────

include { CREATE_HOMOLOGS }  from './modules/homologs'
include { EXPORT_COUNTS }    from './modules/0_export_counts'
include { TSS_REGIONS }      from './modules/1_tss_regions'
include { GET_FASTA }        from './modules/2_get_fasta'
include { MOTIF_ANNOTATIONS } from './modules/3_motif_annotations'
include { GRN_LOCAL }        from './modules/4_grn'
include { CISTARGET_SINGLE } from './modules/5_cistarget_db'
include { CISTARGET_PART }   from './modules/5_cistarget_db'
include { CISTARGET_COMBINE } from './modules/5_cistarget_db'
include { CTX }              from './modules/6_ctx'
include { H5AD_TO_LOOM }    from './modules/7_h5ad_to_loom'
include { AUCELL }           from './modules/8_aucell'
include { CLUSTER_ANALYSIS } from './modules/9_cluster_analysis'
include { ANNOTATE }         from './modules/10_annotate'
include { SPREADSHEET }      from './modules/11_spreadsheet'

// ── Main workflow ───────────────────────────────────────────────────────────

workflow {

    // ── Input channels (always needed) ──────────────────────────────────────

    // .first() converts queue channels to value channels so they can be reused
    // by multiple processes without being consumed on first use
    gff_ch        = channel.fromPath(params.gff, checkIfExists: true).first()
    genome_ch     = channel.fromPath(params.genome_fasta, checkIfExists: true).first()
    genome_fai_ch = channel.fromPath("${params.genome_fasta}.fai", checkIfExists: true).first()
    fly_tbl_ch    = channel.fromPath(params.fly_motif_tbl, checkIfExists: true).first()
    h5ad_ch       = channel.fromPath(params.atlas_h5ad, checkIfExists: true).first()
    cluster_ch    = channel.fromPath(params.cluster_names, checkIfExists: true).first()
    gff_key_ch    = channel.fromPath(params.gff_key, checkIfExists: true).first()

    // ── Ortholog file — import or generate from BLAST maps ──────────────────

    if (params.ortholog_file) {
        ortholog_ch = channel.fromPath(params.ortholog_file, checkIfExists: true).first()
    } else {
        blast_map_ch  = channel.fromPath(params.blast_map, checkIfExists: true).first()
        sp_headers_ch = channel.fromPath(params.sp_headers, checkIfExists: true).first()
        fly_headers_ch = channel.fromPath(params.fly_headers, checkIfExists: true).first()

        CREATE_HOMOLOGS(blast_map_ch, sp_headers_ch, fly_headers_ch, gff_key_ch)
        ortholog_ch = CREATE_HOMOLOGS.out.ortholog_tsv.first()
    }

    // ── Steps 1-3: TSS → FASTA → Motif annotations (parallel) ──────────────

    TSS_REGIONS(gff_ch)
    GET_FASTA(TSS_REGIONS.out.bed, genome_ch, genome_fai_ch)
    MOTIF_ANNOTATIONS(fly_tbl_ch, ortholog_ch)

    // ── Step 7: h5ad → loom (parallel with steps 1-3) ──────────────────────

    H5AD_TO_LOOM(h5ad_ch)

    // ── Step 4: GRN — import or compute locally ─────────────────────────────

    if (params.grn_import) {
        grn_ch = channel.fromPath(params.grn_import, checkIfExists: true)
    } else {
        // Export counts from h5ad, then run GRNBoost2
        EXPORT_COUNTS(h5ad_ch)
        tf_list_ch = channel.fromPath(params.tf_list, checkIfExists: true)
        GRN_LOCAL(EXPORT_COUNTS.out.counts_csv, tf_list_ch)
        grn_ch = GRN_LOCAL.out.grn
    }

    // ── Step 5: cisTarget DB — import or compute locally ────────────────────

    if (params.cistarget_import) {
        rankings_ch = channel.fromPath(params.cistarget_import, checkIfExists: true)
    } else {
        motifs_dir_ch       = channel.fromPath(params.motifs_dir, checkIfExists: true, type: 'dir')
        motif_list_ch       = channel.fromPath(params.motif_list, checkIfExists: true)
        cbust_ch            = channel.fromPath(params.cbust_binary, checkIfExists: true)
        cistarget_script_ch = channel.fromPath(params.cistarget_script, checkIfExists: true)

        if (params.cistarget_total_parts > 1) {
            // Partitioned mode — scatter across parts, then combine
            part_indices = channel.of(1..params.cistarget_total_parts)

            CISTARGET_PART(
                GET_FASTA.out.fasta,
                motifs_dir_ch,
                motif_list_ch,
                cbust_ch,
                cistarget_script_ch,
                params.cistarget_prefix,
                params.cistarget_seed,
                params.cistarget_total_parts,
                part_indices
            )

            CISTARGET_COMBINE(
                CISTARGET_PART.out.partial_scores.collect(),
                params.cistarget_prefix,
                params.cistarget_total_parts,
                params.cistarget_seed
            )
            rankings_ch = CISTARGET_COMBINE.out.rankings
        } else {
            // Single-run mode — produces rankings directly
            CISTARGET_SINGLE(
                GET_FASTA.out.fasta,
                motifs_dir_ch,
                motif_list_ch,
                cbust_ch,
                cistarget_script_ch,
                params.cistarget_prefix,
                params.cistarget_seed
            )
            rankings_ch = CISTARGET_SINGLE.out.rankings
        }
    }

    // ── Step 6: pySCENIC ctx ────────────────────────────────────────────────

    CTX(grn_ch, rankings_ch, MOTIF_ANNOTATIONS.out.tbl, H5AD_TO_LOOM.out.loom)

    // ── Step 8: AUCell ──────────────────────────────────────────────────────

    AUCELL(H5AD_TO_LOOM.out.loom, CTX.out.regulons)

    // ── Step 9: Cluster analysis ────────────────────────────────────────────

    CLUSTER_ANALYSIS(AUCELL.out.auc_loom, h5ad_ch, cluster_ch, ortholog_ch)

    // ── Step 10: Annotate regulons ──────────────────────────────────────────

    ANNOTATE(gff_ch, ortholog_ch, CLUSTER_ANALYSIS.out.activity_csv, CLUSTER_ANALYSIS.out.top_regulons)

    // ── Step 11: Excel spreadsheet ──────────────────────────────────────────

    SPREADSHEET(CTX.out.regulons, gff_key_ch, ortholog_ch)
}
