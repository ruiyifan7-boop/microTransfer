# microTransfer 0.1.0

`microTransfer` turns the study-aware framework used in the Vaccinium
benchmark into a reusable four-step workflow:

1. declare whether the claim concerns **membership**, **response**, or
   **prediction**;
2. discover candidate signals without using the held-out cohort;
3. match validation to the property;
4. audit matched nulls, filtering/depth choices, and label uncertainty before
   assigning a transferability label.

The package accepts a sample-by-feature count matrix and a metadata data frame.
It does not require `phyloseq`.

## Install from the source directory

```bash
R CMD INSTALL microTransfer
```

```r
library(microTransfer)
```

## Minimal membership example

```r
claim <- declare_transferability(
  property = "membership",
  study_col = "study",
  feature_level = "genus"
)

fit <- run_transferability_workflow(
  claim,
  counts = counts,
  metadata = metadata,
  sample_col = "sample_id",
  prevalence = 0.50,
  n_null = 999,
  seed = 20260722
)

fit$summary
classify_transferability(fit)
```

See `vignette("study-aware-workflow", package = "microTransfer")` for all
three properties. A completed non-*Vaccinium* demonstration in 978
Arabidopsis endosphere/rhizosphere samples is documented in
`inst/doc/GENERALIZATION_ARABIDOPSIS.md`; the prespecified full-OTU rice
generalization is in `inst/scripts/run_edwards_rice_generalization.R`.

## Reproducing the study analyses

The `analysis/` directory contains the scripts used for the *Vaccinium*, rice,
Arabidopsis, and FAPROTAX analyses. Large source data are not redistributed.
The scripts retrieve or read them from their persistent public records:

- raw *Vaccinium* amplicon reads: ENA/NCBI BioProjects listed in the article;
- rice lifecycle data: Dryad DOI `10.5061/dryad.7q7k1`;
- Arabidopsis demonstration data: CoreMicrobiomeR 0.1.0.

Derived data, analysis-ready tables, and exact study outputs are distributed
with the companion data record described in `DATA_AVAILABILITY.md`.

## Package checks

From the parent directory, run:

```sh
R CMD build microTransfer
R CMD check --as-cran microTransfer_0.1.0.tar.gz
```

The release archive includes the check summary and SHA-256 manifest generated
for version 0.1.0.

## Citation

Please cite the associated mSystems Research Article and the archived software
release. The release-specific DOI should be used to identify the exact code
version; it will be added to the article after Zenodo registration.

## License

The package and analysis code are released under the MIT License. External
datasets remain subject to the terms of their source repositories.
