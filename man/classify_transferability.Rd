\name{classify_transferability}
\alias{classify_transferability}
\title{Assign a Conservative Transferability Evidence Label}
\usage{
classify_transferability(result, min_fraction = 0.7, alpha = 0.05)
}
\arguments{
  \item{result}{A validation result from \pkg{microTransfer}.}
  \item{min_fraction}{Minimum fraction of cohorts required.}
  \item{alpha}{Empirical or nominal threshold.}
}
\value{A one-row evidence summary data frame.}
\description{
Summarizes prespecified fold-level rules. The label does not replace effect
sizes, uncertainty intervals, or cohort-level results.
}
