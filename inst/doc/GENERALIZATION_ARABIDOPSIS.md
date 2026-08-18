# Completed external demonstration: Arabidopsis

## Dataset and scope

The first completed non-*Vaccinium* demonstration uses the public
`CoreMicrobiomeR` 0.1.0 CRAN object derived from Lundberg et al.,
*Defining the core Arabidopsis thaliana root microbiome*
(DOI: `10.1038/nature11237`). The object contains 188 representative OTUs
across 1,439 samples. It is therefore a framework-level demonstration, not a
full biological reanalysis of the original 18,783-OTU table.

Sample identifiers encode soil, genotype, fraction, developmental stage,
experiment, and plate. The parser reproduced all five corresponding design
fields for all 103 rows in the packaged metadata subset.

## Prespecified analysis

- Primary blocks: C1, C2, M1, and M2.
- Compartments: endosphere and rhizosphere.
- QC: at least 100 counts across the 188 retained OTUs.
- Membership: prevalence at least 50% in every training block and total
  training abundance at least 20.
- Validation: leave-one-experiment-out with 999 prevalence- and
  abundance-matched null subsets.
- Response: partial old-versus-young correlations adjusted for compartment
  and genotype, random-effects discovery in training experiments, and
  direction-concordant nominal replication in the held-out experiment.
- Prediction: nearest-centroid endosphere-versus-rhizosphere discrimination
  using training-defined core features and 999 within-training-block label
  permutations.

## Results

- 978 samples passed the primary design and retained-feature-depth filters.
- At 50% prevalence, pooled analysis identified 154 core OTUs, whereas the
  study-aware definition retained 98.
- LOSO membership validation ranged from 86.0% to 98.0%; every fold exceeded
  its matched-null 95th percentile (all empirical P = 0.001).
- 138 of 203 fold-specific developmental-stage candidates replicated
  nominally with concordant direction (68.0%).
- Mean held-out compartment balanced accuracy was 0.978 and mean macro-F1 was
  0.980; the overall permutation P value was 0.001.
- Across retained-count thresholds of 50–500, mean core validation remained
  0.908–0.969 and mean balanced accuracy remained 0.943–0.981.

Fold C1 alone had a prediction permutation P value of 0.052; the fold-level
table is retained to prevent selective reporting.

## Packaged result tables

The compact summary tables are installed under `inst/extdata`:

- `arabidopsis_pooled_vs_study_aware_core.tsv`
- `arabidopsis_core_membership_loso.tsv`
- `arabidopsis_developmental_response_loso.tsv`
- `arabidopsis_compartment_prediction_loso.tsv`
- `arabidopsis_retained_depth_sensitivity.tsv`

The complete reproducibility bundle, including feature-level tables, SVG/PNG
figure, source manifest, Python audit implementation, and manuscript-ready
text, is distributed separately as
`arabidopsis_external_generalization_bundle_20260722.zip`.
