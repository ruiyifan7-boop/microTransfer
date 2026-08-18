#' Validate and align a microbiome count matrix and metadata
#'
#' @param counts Numeric sample-by-feature matrix or data frame.
#' @param metadata Sample metadata.
#' @param sample_col Metadata sample identifier column.
#' @param study_col Metadata cohort column.
#' @param outcome_col Optional outcome column.
#' @param allow_zero_depth Retain zero-depth samples when `TRUE`.
#' @return A list containing aligned `counts` and `metadata`.
#' @export
validate_microbiome_input <- function(
    counts,
    metadata,
    sample_col,
    study_col,
    outcome_col = NULL,
    allow_zero_depth = FALSE) {
  counts <- as.matrix(counts)
  storage.mode(counts) <- "double"
  if (is.null(rownames(counts)) || any(!nzchar(rownames(counts)))) {
    stop("`counts` must have unique sample IDs as row names.", call. = FALSE)
  }
  if (anyDuplicated(rownames(counts))) {
    stop("Sample IDs in `counts` are duplicated.", call. = FALSE)
  }
  if (is.null(colnames(counts)) || any(!nzchar(colnames(counts))) ||
      anyDuplicated(colnames(counts))) {
    stop("`counts` must have unique, non-empty feature names.", call. = FALSE)
  }
  if (any(!is.finite(counts)) || any(counts < 0)) {
    stop("`counts` must contain finite, non-negative values.", call. = FALSE)
  }
  needed <- c(sample_col, study_col, outcome_col)
  needed <- needed[!is.na(needed) & nzchar(needed)]
  missing_cols <- setdiff(needed, names(metadata))
  if (length(missing_cols)) {
    stop("Missing metadata column(s): ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  sample_ids <- as.character(metadata[[sample_col]])
  if (anyNA(sample_ids) || any(!nzchar(sample_ids)) ||
      anyDuplicated(sample_ids)) {
    stop("Metadata sample IDs must be unique and non-missing.",
         call. = FALSE)
  }
  shared <- intersect(sample_ids, rownames(counts))
  if (!length(shared)) {
    stop("No sample IDs are shared by counts and metadata.", call. = FALSE)
  }
  metadata <- metadata[match(shared, sample_ids), , drop = FALSE]
  counts <- counts[shared, , drop = FALSE]
  rownames(metadata) <- shared
  keep <- !is.na(metadata[[study_col]]) &
    nzchar(trimws(as.character(metadata[[study_col]])))
  if (!isTRUE(allow_zero_depth)) keep <- keep & rowSums(counts) > 0
  counts <- counts[keep, , drop = FALSE]
  metadata <- metadata[keep, , drop = FALSE]
  if (nrow(counts) < 2L) {
    stop("Fewer than two aligned, non-empty samples remain.", call. = FALSE)
  }
  list(counts = counts, metadata = metadata)
}

#' Centered log-ratio transform
#'
#' @param counts Sample-by-feature non-negative matrix.
#' @param pseudocount Positive pseudocount.
#' @return Numeric CLR matrix.
#' @export
clr_transform <- function(counts, pseudocount = 0.5) {
  if (!is.numeric(pseudocount) || length(pseudocount) != 1L ||
      !is.finite(pseudocount) || pseudocount <= 0) {
    stop("`pseudocount` must be a positive number.", call. = FALSE)
  }
  z <- log(as.matrix(counts) + pseudocount)
  sweep(z, 1L, rowMeans(z), FUN = "-")
}

.transform_counts <- function(counts, transform, pseudocount = 0.5) {
  transform <- match.arg(transform, c("clr", "log1p", "relative", "none"))
  if (transform == "clr") return(clr_transform(counts, pseudocount))
  if (transform == "log1p") return(log1p(counts))
  if (transform == "relative") {
    totals <- rowSums(counts)
    return(sweep(counts, 1L, totals, FUN = "/"))
  }
  as.matrix(counts)
}
