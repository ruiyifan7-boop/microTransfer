# Analysis scripts

This directory contains study-analysis scripts that accompany the reusable
`microTransfer` package. They are retained outside the built R package so the
software API remains compact while the article workflow stays auditable.

## Directory map

- `vaccinium/`: harmonization, transferability, pH, management, sensitivity,
  external-cohort, figure, and FAPROTAX workflow scripts.
- `rice_arabidopsis/`: external-generalization adapters and reproduction
  scripts for the public rice and Arabidopsis datasets.
- `faprotax/`: sample-level FAPROTAX input, LOSO, and summary scripts.

The scripts expect the source data and directory layouts described in their
headers. Public source identifiers are recorded in `DATA_AVAILABILITY.md` and
in the relevant manifests. Exact derived inputs and outputs are supplied in
the companion data archive rather than duplicated in this software release.
