\name{run_transferability_workflow}
\alias{run_transferability_workflow}
\title{Run a Declared Transferability Workflow}
\usage{
run_transferability_workflow(claim, counts, metadata, sample_col, ...)
}
\arguments{
  \item{claim}{An object from \code{declare_transferability()}.}
  \item{counts}{Sample-by-feature counts.}
  \item{metadata}{Sample metadata.}
  \item{sample_col}{Sample identifier column.}
  \item{...}{Property-specific validation arguments.}
}
\value{A membership, response, or prediction validation result.}
\description{Dispatches to the validation design matched to the declared claim.}
