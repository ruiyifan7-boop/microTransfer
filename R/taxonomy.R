#' Collapse feature counts to a taxonomic rank
#'
#' @param counts Sample-by-feature count matrix.
#' @param taxonomy Data frame with one row per feature.
#' @param rank Taxonomy column to aggregate.
#' @param feature_col Optional taxonomy feature-ID column. Row names are used
#'   when omitted.
#' @param unclassified_label Label used for missing assignments.
#' @return Sample-by-taxon count matrix.
#' @export
collapse_counts_by_taxonomy <- function(
    counts,
    taxonomy,
    rank,
    feature_col = NULL,
    unclassified_label = "Unclassified") {
  counts <- as.matrix(counts)
  ids <- if (is.null(feature_col)) {
    rownames(taxonomy)
  } else {
    as.character(taxonomy[[feature_col]])
  }
  if (is.null(ids)) {
    stop("Taxonomy must have row names or `feature_col`.", call. = FALSE)
  }
  if (!rank %in% names(taxonomy)) {
    stop("Taxonomy rank not found: ", rank, call. = FALSE)
  }
  shared <- intersect(colnames(counts), ids)
  if (!length(shared)) {
    stop("No features are shared by counts and taxonomy.", call. = FALSE)
  }
  tax <- as.character(taxonomy[[rank]][match(shared, ids)])
  tax[is.na(tax) | !nzchar(trimws(tax))] <- unclassified_label
  collapsed <- rowsum(
    t(counts[, shared, drop = FALSE]),
    group = tax,
    reorder = FALSE
  )
  t(collapsed)
}
