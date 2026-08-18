\name{collapse_counts_by_taxonomy}
\alias{collapse_counts_by_taxonomy}
\title{Collapse Counts to a Taxonomic Rank}
\usage{
collapse_counts_by_taxonomy(
  counts, taxonomy, rank, feature_col = NULL,
  unclassified_label = "Unclassified"
)
}
\arguments{
  \item{counts}{Sample-by-feature count matrix.}
  \item{taxonomy}{Feature taxonomy data frame.}
  \item{rank}{Taxonomic rank column.}
  \item{feature_col}{Optional feature identifier column.}
  \item{unclassified_label}{Missing-rank label.}
}
\value{A sample-by-taxon count matrix.}
\description{Aggregates features sharing the same taxonomic assignment.}
