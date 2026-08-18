#' Declare a transferability claim
#'
#' The declaration is deliberately separate from model fitting. It records
#' which property is being claimed and therefore which validation design is
#' appropriate.
#'
#' @param property One of `"membership"`, `"response"`, or `"prediction"`.
#' @param study_col Metadata column identifying independent cohorts or blocks.
#' @param outcome_col Metadata outcome column for response or prediction.
#' @param feature_level Human-readable feature level, for example `"genus"`.
#' @param validation Optional validation design. A property-specific default is
#'   used when omitted.
#' @param provisional Whether the outcome or cohort labels are provisional.
#' @param notes Optional free-text notes.
#' @return An object of class `mt_claim`.
#' @export
declare_transferability <- function(
    property = c("membership", "response", "prediction"),
    study_col,
    outcome_col = NULL,
    feature_level = "feature",
    validation = NULL,
    provisional = FALSE,
    notes = NULL) {
  property <- match.arg(property)
  if (!is.character(study_col) || length(study_col) != 1L ||
      !nzchar(study_col)) {
    stop("`study_col` must be one non-empty column name.", call. = FALSE)
  }
  if (property != "membership" &&
      (is.null(outcome_col) || length(outcome_col) != 1L ||
       !nzchar(outcome_col))) {
    stop("`outcome_col` is required for response and prediction claims.",
         call. = FALSE)
  }
  defaults <- c(
    membership = "leave-one-study-out",
    response = "discovery-projection",
    prediction = "held-out-study"
  )
  if (is.null(validation)) validation <- unname(defaults[property])
  structure(
    list(
      property = property,
      study_col = study_col,
      outcome_col = outcome_col,
      feature_level = feature_level,
      validation = validation,
      provisional = isTRUE(provisional),
      notes = notes
    ),
    class = "mt_claim"
  )
}

#' @export
print.mt_claim <- function(x, ...) {
  cat("<mt_claim>\n")
  cat("  property:   ", x$property, "\n", sep = "")
  cat("  validation: ", x$validation, "\n", sep = "")
  cat("  study:      ", x$study_col, "\n", sep = "")
  if (!is.null(x$outcome_col)) {
    cat("  outcome:    ", x$outcome_col, "\n", sep = "")
  }
  cat("  feature:    ", x$feature_level, "\n", sep = "")
  if (isTRUE(x$provisional)) cat("  labels:      provisional\n")
  invisible(x)
}
