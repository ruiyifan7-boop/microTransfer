\name{declare_transferability}
\alias{declare_transferability}
\title{Declare a Transferability Claim}
\usage{
declare_transferability(
  property = c("membership", "response", "prediction"),
  study_col,
  outcome_col = NULL,
  feature_level = "feature",
  validation = NULL,
  provisional = FALSE,
  notes = NULL
)
}
\arguments{
  \item{property}{Membership, response, or prediction.}
  \item{study_col}{Metadata column identifying independent cohorts.}
  \item{outcome_col}{Outcome column for response or prediction.}
  \item{feature_level}{Human-readable feature level.}
  \item{validation}{Optional validation design.}
  \item{provisional}{Whether labels are provisional.}
  \item{notes}{Optional notes.}
}
\value{An object of class \code{mt_claim}.}
\description{
Records the property being claimed before candidate discovery or validation.
}
