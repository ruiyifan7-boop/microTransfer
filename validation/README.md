# Release validation

Validation was completed on 18 August 2026 with R 4.6.1 on Windows 11.

## Package check

The source package was built with its vignette and checked using
`R CMD check --as-cran --no-manual`. The final status was zero errors, zero
warnings, and one informational note (`New submission`). All testthat tests,
documentation checks, examples, and vignette rebuilds passed. The complete
check log is in `check/00check.log`.

## Rice source integrity

The four public Dryad files from DOI `10.5061/dryad.7q7k1` matched the archived
manifest:

| File | MD5 |
|---|---|
| `lc_study_mapping_file.tsv` | `70c776ab613800347a86a1b78f9153a5` |
| `lc_study_otu_table.tsv.gz` | `85bde2c7a81140907d2a9cceae5f55c4` |
| `gg_otus_tax.rds` | `cfb234a8be3f3f195b6b5b903d16782f` |
| `organelle.rds` | `0965d49f14366ed053e37e4943c61e27` |

These third-party source files are not redistributed in this software
release. The download and verification scripts are in
`analysis/rice_arabidopsis/`.

## End-to-end results

The official rice workflow retained 965 samples and 30,609 OTUs across four
site-year cohorts.

- Membership: 806 leave-one-site-year-out-stable features; 4/4 held-out
  cohorts passed the matched-null criterion.
- Response: held-out replication fractions were 0.7841, 0.7367, 0.7273, and
  0.5411.
- Prediction: mean held-out balanced accuracy was 0.9845; permutation
  `P = 0.0010` with 999 permutations.
- ComBat control: study R2 decreased from 0.1539 to 0.0120, while membership
  remained 1,540 pooled versus 806 study-aware features.

Exact TSV outputs and session information are retained in `rice_official/`
and `combat_official/`.
