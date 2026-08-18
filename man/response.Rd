\name{response}
\alias{discover_response}
\alias{validate_response}
\title{Study-Aware Response Analysis}
\usage{
discover_response(
  counts, metadata, sample_col, study_col, outcome_col,
  discovery_studies = NULL, covariates = NULL, transform = "clr",
  pseudocount = 0.5, min_prevalence = 0.1, min_total = 20,
  alpha = 0.05
)

validate_response(
  counts, metadata, sample_col, study_col, outcome_col,
  discovery_studies, validation_studies, covariates = NULL,
  transform = "clr", pseudocount = 0.5, min_prevalence = 0.1,
  min_total = 20, alpha = 0.05, nominal_alpha = 0.05
)
}
\arguments{
  \item{counts}{Sample-by-feature counts.}
  \item{metadata}{Sample metadata.}
  \item{sample_col}{Sample identifier column.}
  \item{study_col}{Cohort column.}
  \item{outcome_col}{Numeric environmental response column.}
  \item{discovery_studies}{Cohorts used for discovery.}
  \item{validation_studies}{Cohorts reserved for projection.}
  \item{covariates}{Optional model covariates.}
  \item{transform}{Feature transformation.}
  \item{pseudocount}{CLR pseudocount.}
  \item{min_prevalence}{Within-study prevalence filter.}
  \item{min_total}{Within-study total-count filter.}
  \item{alpha}{Discovery FDR threshold.}
  \item{nominal_alpha}{Independent validation threshold.}
}
\value{A response discovery or validation object.}
\description{
Fits study-specific feature models, combines discovery effects, and tests
direction-concordant replication in reserved cohorts.
}
