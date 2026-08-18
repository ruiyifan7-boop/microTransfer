\name{read_edwards_rice}
\alias{read_edwards_rice}
\title{Read the Edwards Field-Grown Rice Lifecycle Dataset}
\usage{
read_edwards_rice(
  otu_file, metadata_file, taxonomy_file = NULL,
  organelle_file = NULL, min_depth = 5000,
  compartments = c("Rhizosphere", "Endosphere"),
  sites = c("Arbuckle", "Jonesboro"),
  cohort = c("site_year", "site"), rank = NULL,
  taxonomy_feature_col = NULL
)
}
\arguments{
  \item{otu_file}{Path to the gzipped OTU table.}
  \item{metadata_file}{Path to sample metadata.}
  \item{taxonomy_file}{Optional taxonomy RDS.}
  \item{organelle_file}{Optional organelle OTU RDS.}
  \item{min_depth}{Minimum sample depth.}
  \item{compartments}{Compartments to retain.}
  \item{sites}{Field sites to retain.}
  \item{cohort}{Site-year or site cohort definition.}
  \item{rank}{Optional taxonomic aggregation rank.}
  \item{taxonomy_feature_col}{Optional taxonomy feature ID column.}
}
\value{A list with counts, metadata, taxonomy, and source information.}
\description{
Reads locally downloaded files from Dryad DOI 10.5061/dryad.7q7k1 and prepares
the non-Vaccinium generalization dataset.
}
