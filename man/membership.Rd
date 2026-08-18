\name{membership}
\alias{discover_membership}
\alias{validate_membership}
\alias{pooled_membership_contrast}
\alias{membership_sensitivity}
\title{Study-Aware Membership Analysis}
\usage{
discover_membership(
  counts, metadata, sample_col, study_col, heldout = NULL,
  prevalence = 0.5, min_total = 0, min_studies = NULL
)

validate_membership(
  counts, metadata, sample_col, study_col, prevalence = 0.5,
  min_total = 0, n_null = 999, matched = TRUE, seed = 1L
)

pooled_membership_contrast(
  counts, metadata, sample_col, study_col,
  thresholds = c(0.1, 0.3, 0.5)
)

membership_sensitivity(
  counts, metadata, sample_col, study_col,
  prevalence = c(0.3, 0.5, 0.7), min_depth = c(0, 5000),
  min_total = c(0, 20), n_null = 0, matched = TRUE, seed = 1L
)
}
\arguments{
  \item{counts}{Sample-by-feature counts.}
  \item{metadata}{Sample metadata.}
  \item{sample_col}{Sample identifier column.}
  \item{study_col}{Cohort column.}
  \item{heldout}{Optional cohort omitted from discovery.}
  \item{prevalence}{Within-cohort prevalence threshold.}
  \item{min_total}{Minimum discovery total count.}
  \item{min_studies}{Number of discovery cohorts required to pass.}
  \item{n_null}{Number of random-subset draws.}
  \item{matched}{Use prevalence/abundance-matched nulls.}
  \item{seed}{Random seed.}
  \item{thresholds}{Thresholds for pooled versus study-aware comparison.}
  \item{min_depth}{Minimum sample-depth thresholds.}
}
\value{Discovery or validation objects, or a framework-value data frame.}
\description{
Discovers features without held-out cohorts and evaluates prevalence
transportability against cohort-specific random-subset nulls.
}
