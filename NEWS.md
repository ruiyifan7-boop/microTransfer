# microTransfer 0.1.0

* Optimized prevalence-and-abundance-matched null sampling with vectorized
  stratum-level draws, avoiding repeated full-universe scans while preserving
  sampling without replacement and the specified matching strata.

- Added a claim declaration object for membership, response, and prediction.
- Added leave-one-study-out membership validation with simple and
  prevalence/abundance-matched random-subset nulls.
- Added discovery-projection response validation using study-specific
  feature models and inverse-variance random-effects pooling.
- Added held-out-study classification with balanced accuracy, macro-F1, and
  class-preserving label-noise sensitivity.
- Added an adapter and end-to-end script for the Edwards et al. field-grown
  rice lifecycle dataset (Dryad DOI: 10.5061/dryad.7q7k1).
