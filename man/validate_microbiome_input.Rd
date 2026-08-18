\name{validate_microbiome_input}
\alias{validate_microbiome_input}
\alias{clr_transform}
\title{Validate and Transform Microbiome Inputs}
\usage{
validate_microbiome_input(
  counts, metadata, sample_col, study_col,
  outcome_col = NULL, allow_zero_depth = FALSE
)

clr_transform(counts, pseudocount = 0.5)
}
\arguments{
  \item{counts}{Numeric sample-by-feature matrix.}
  \item{metadata}{Sample metadata data frame.}
  \item{sample_col}{Metadata sample identifier column.}
  \item{study_col}{Metadata cohort column.}
  \item{outcome_col}{Optional outcome column.}
  \item{allow_zero_depth}{Whether zero-depth samples may remain.}
  \item{pseudocount}{Positive CLR pseudocount.}
}
\value{
\code{validate_microbiome_input()} returns aligned counts and metadata.
\code{clr_transform()} returns a numeric matrix.
}
\description{Input validation, alignment, and centered log-ratio transformation.}
