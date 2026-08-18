# External extension runbook

## Objective

Expand the harmonized analysis from four studies and 211 primary paired
samples to as many as seven studies and 277 paired samples before quality
control. The new cohorts are deliberately assigned different roles:

- PRJNA835240: third pH-enabled cohort for study-level random-effects
  meta-analysis and external pH validation.
- PRJNA577971/PRJNA578171: 34-pair external validation cohort for core
  microbiome transportability.
- PRJNA1253324: small, recent 8-pair cultivar cohort used as a difficult
  external sensitivity test, not as a discovery cohort.

## Phase 1: reproduce the ENA audit

Run from `~/projects/blueberry_meta`:

```bash
python audit_external_extension.py
column -t external_extension/external_candidate_summary.tsv
```

Expected output:

- PRJNA835240: 48 runs, 24 paired biological samples.
- PRJNA1253324: 16 runs, 8 paired biological samples.
- PRJNA577971/PRJNA578171: 68 runs, 34 paired biological samples.
- Total download: approximately 1.49 GiB in 216 FASTQ files.

The audit fails immediately if marker multiplicity or sample pairing changes.

## Phase 2: download without mixing cohorts

```bash
for study in \
  PRJNA835240 \
  PRJNA1253324 \
  PRJNA577971_PRJNA578171
do
  mkdir -p "external_extension/$study/fastq" \
           "external_extension/$study/logs"
  nohup bash -c \
    "xargs -P 4 -n 1 wget -c -nv -P external_extension/$study/fastq \
     < external_extension/fastq_urls/${study}_fastq_urls.txt" \
    > "external_extension/$study/logs/download.log" 2>&1 &
  echo "$!" > "external_extension/$study/download.pid"
done
```

Monitor:

```bash
watch -n 60 '
for study in PRJNA835240 PRJNA1253324 PRJNA577971_PRJNA578171
do
  echo "===== $study ====="
  P=$(cat external_extension/$study/download.pid 2>/dev/null)
  ps -p "$P" -o pid,etime,%cpu,%mem 2>/dev/null
  find external_extension/$study/fastq -name "*.fastq.gz" | wc -l
  du -sh external_extension/$study/fastq 2>/dev/null
  tail -n 2 external_extension/$study/logs/download.log 2>/dev/null
done
'
```

Validate checksums from inside each cohort FASTQ directory:

```bash
for study in PRJNA835240 PRJNA1253324 PRJNA577971_PRJNA578171
do
  (
    cd "external_extension/$study/fastq" &&
    md5sum -c "../../fastq_urls/${study}_fastq_md5.txt"
  ) > "external_extension/$study/logs/md5.log" 2>&1
done
grep -R "FAILED" external_extension/*/logs/md5.log
```

## Phase 3: marker-specific processing

Do not combine raw reads across projects during error learning or denoising.
Process every study-marker combination independently, then harmonize at the
taxonomic genus level.

Important constraints:

- PRJNA835240 is single-end in ENA. Run single-end DADA2 independently for
  16S and ITS. Do not attempt paired-read merging.
- PRJNA1253324 and PRJNA577971/PRJNA578171 are paired-end.
- Confirm every primer marked `confirm_from_supplement` in
  `external_protocol_audit.tsv` before trimming.
- For ITS, retain the established two-table sensitivity design:
  UNITE-phylum-confident primary table and strict EUKARYOME fungal table.
- Use SILVA 138.1 and UNITE v10/19.02.2025 to match the current manuscript.
- Save filter, denoising, merging, chimera, taxonomic, and sample-depth tracks.

Minimum primary QC:

- Bacterial clean depth at least 5,000 reads.
- UNITE phylum-confident fungal depth at least 5,000 reads.
- Both markers present for a biological sample.

## Phase 4: rebuild the seven-study harmonized objects

Create:

```text
global_harmonized_extended/
  global_metadata_primary.tsv
  analysis_ready/
    bacteria_genus_model_counts.rds
    fungi_genus_model_counts.rds
```

Preserve the current `sample_uid` convention:

```text
study__biological_sample
```

Rerun the existing harmonization, prevalence filtering, study-effect,
cross-kingdom, core, and strict-fungal sensitivity scripts before adding new
models.

## Phase 5: new high-value analyses

Install the one new package if needed:

```bash
Rscript -e 'if(!requireNamespace("metafor",quietly=TRUE)) install.packages("metafor",repos="https://cloud.r-project.org")'
```

Run:

```bash
export IMETA_EXTENDED_ROOT=global_harmonized_extended

Rscript run_pH_random_effects_meta.R \
  > global_harmonized_extended/pH_random_effects_meta/run.log 2>&1

Rscript run_domain_shift_core_transport.R \
  > global_harmonized_extended/domain_shift_core_transport/run.log 2>&1
```

The pH meta-analysis deliberately requires at least three studies with at
least eight rhizosphere samples, three unique pH values, and a pH range of
at least 0.3. It estimates standardized study-specific linear and quadratic
pH effects and combines taxon effects with REML random-effects models,
Knapp-Hartung inference, heterogeneity, and prediction intervals.

The domain-shift analysis performs leakage-free leave-one-study-out feature
filtering and reports centroid shift, Jensen-Shannon divergence, prevalence
drift, and held-out validation of the training-study core.

## Decision rule

Promote the seven-study extension into the main manuscript only if:

- at least 50 of 66 new pairs pass both-marker QC;
- at least three studies remain eligible for pH meta-analysis;
- bacterial core validation remains at least 70% in most held-out studies;
- the principal management and pH conclusions do not reverse;
- protocol/domain-shift effects are reported rather than hidden.

Otherwise, retain the new cohorts as external sensitivity analyses.
