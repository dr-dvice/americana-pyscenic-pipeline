# pySCENIC cross-species regulon pipeline

Nextflow (DSL2) pipeline for gene regulatory network and regulon inference with
[pySCENIC](https://github.com/aertslab/pySCENIC), built to run across the species
in the paper using a fly-hub ortholog mapping. Companion code for the paper.

Part of the project hub: https://dr-dvice.github.io/americana/

## Overview

An 11-step pipeline that goes from a single-cell atlas (`.h5ad`) + genome to
per-cell regulon activity and annotated regulon tables:

```
homologs → export counts → TSS regions → get FASTA → motif annotations →
GRNBoost2 (GRN) → cisTarget DB → ctx (prune) → h5ad→loom → AUCell →
cluster analysis → annotate → spreadsheet
```

The heavy steps (GRNBoost2 GRN inference and the cisTarget DB build) can be
**imported** from pre-computed files via `--grn_import` / `--cistarget_import`,
so a repeated run skips them.

```
main.nf            Workflow + parameter validation
modules/*.nf       One Nextflow module per step (0–11 + homologs)
bin/*.py           Python helpers invoked by the modules
conf/*.config      Species profiles: americana, gregaria, test
envs/*.yml         Conda environments (see Requirements)
test_data/         Small 480-cell subset for `-profile test`
```

## Requirements

- **Nextflow** (DSL2)
- **Singularity** + the **pySCENIC 0.12.1** container (see below) — used for the
  GRN/ctx/AUCell steps (`singularity_pyscenic` label)
- **Conda** environments (provided in `envs/`):
  - `envs/pyscenic.yml` — pySCENIC tooling for the lighter Python steps
  - `envs/create_cistarget_databases.yml` — aertslab cisTarget DB builder
- aertslab **cisTarget motif collection** (v10nr_clust) + the fly motif
  annotation table, both from the aertslab resources site:
  [resources.aertslab.org/cistarget](https://resources.aertslab.org/cistarget/)
  (motif collections under `motif_collections/`, the fly motif-to-TF table under
  `motif2tf/`). Paths are set in `nextflow.config` / `conf/americana.config`.

### Obtaining the container
Download the pySCENIC 0.12.1 image, then point `params.singularity_image` at it:

```bash
singularity build pyscenic-0.12.1.sif docker://aertslab/pyscenic:0.12.1
```

## Usage

```bash
# Typical run (imports pre-computed GRN + cisTarget DB) for a species profile
nextflow run main.nf -profile americana

# Recompute the GRN locally instead of importing
nextflow run main.nf -profile americana --grn_import null

# Resume after a failure
nextflow run main.nf -profile americana -resume
```


```bash
nextflow run main.nf -profile test
```

Note that the test profile (`conf/test.config`) still references some external
inputs (genome, motif list, TF list, BLAST maps) by absolute path — see below.

## Inputs and paths

The profiles in `conf/` and `nextflow.config` contain **absolute paths specific
to the authors' system** (genome FASTA/GFF, atlas `.h5ad`, BLAST ortholog maps,
TF lists, the cisTarget motif collection, the fly motif table, conda env
locations, and the container path). They are kept as a record of the exact runs;
to reproduce, edit the relevant profile to point at your own inputs and
environments.

FUTURE UPDATE: SRA reference to be provided when available

## License

GNU General Public License v3.0 — see [`LICENSE`](LICENSE).
