# Non-Vaccinium generalization: field-grown rice

## Selected dataset

Edwards JA et al. *Compositional shifts in root-associated bacterial and
archaeal microbiota track the plant life cycle in field-grown rice*. PLOS
Biology (2018), DOI: `10.1371/journal.pbio.2003862`.

Data: Dryad DOI `10.5061/dryad.7q7k1`.

This dataset was selected because it provides:

- a complete 97% GreenGenes OTU count table;
- sample-level metadata;
- taxonomy and organelle reference objects;
- dense field sampling across three seasons and two geographic sites;
- root compartments and plant age, supporting membership, response, and
  prediction claims without reprocessing FASTQ files.

The primary analysis uses `Site × Season` as the independent cohort and
restricts samples to Arbuckle and Jonesboro, rhizosphere and endosphere, and
at least 5,000 reads.

## Required files

Download these four files from
<https://datadryad.org/dataset/doi:10.5061/dryad.7q7k1>:

- `lc_study_otu_table.tsv.gz`
- `lc_study_mapping_file.tsv`
- `gg_otus_tax.rds`
- `organelle.rds`

Expected MD5 checksums are recorded in
`inst/extdata/edwards_rice_dryad_manifest.tsv`.

## Server run

```bash
cd ~/projects/blueberry_meta

R CMD INSTALL microTransfer

export RICE_DATA_DIR="$PWD/external_generalization/edwards_rice"
export RICE_OUT_DIR="$PWD/external_generalization/edwards_rice/results"
mkdir -p "$RICE_OUT_DIR"

Rscript microTransfer/inst/scripts/run_edwards_rice_generalization.R \
  "$RICE_DATA_DIR/lc_study_otu_table.tsv.gz" \
  "$RICE_DATA_DIR/lc_study_mapping_file.tsv" \
  "$RICE_OUT_DIR" \
  "$RICE_DATA_DIR/gg_otus_tax.rds" \
  "$RICE_DATA_DIR/organelle.rds" \
  2>&1 | tee "$RICE_OUT_DIR/rice_generalization.log"
```

## Prespecified claims

1. **Membership:** features discovered without each site-year cohort are tested
   for prevalence in that held-out cohort and compared with a
   prevalence/abundance-matched random-subset null.
2. **Response:** age-associated features are discovered by study-specific
   models and random-effects pooling, then projected into a reserved cohort
   with direction-concordant nominal replication.
3. **Prediction:** the study-aware stable feature set is used to discriminate
   rhizosphere from endosphere samples in held-out site-year cohorts.

The OTU-level result is the prespecified primary generalization because all
site-year cohorts share a single reference-OTU table. Genus aggregation is a
sensitivity analysis and should only be enabled after confirming the taxonomy
object's feature identifier and genus column.
