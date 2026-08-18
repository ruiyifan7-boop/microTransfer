#' Assign a conservative transferability evidence label
#'
#' Labels are descriptive summaries of prespecified rules, not replacements for
#' effect estimates, confidence intervals, or cohort-level results.
#'
#' @param result Validation result from this package.
#' @param min_fraction Minimum fraction of cohorts required for a transferable
#'   label.
#' @param alpha Empirical or nominal threshold.
#' @return A one-row data frame with the evidence label and diagnostics.
#' @export
classify_transferability <- function(
    result,
    min_fraction = 0.70,
    alpha = 0.05) {
  property <- result$property
  if (identical(property, "membership")) {
    z <- result$summary
    passed <- is.finite(z$empirical_p) & z$empirical_p <= alpha &
      z$above_null_q95
    fraction <- if (length(passed)) mean(passed) else 0
    label <- if (fraction >= min_fraction &&
                 length(result$stable_features) > 0L) {
      "transferable membership"
    } else if (any(passed)) {
      "context-dependent membership"
    } else {
      "membership not established"
    }
    return(data.frame(
      property = property,
      label = label,
      cohorts_passing = sum(passed),
      cohorts_tested = length(passed),
      passing_fraction = fraction,
      stable_features = length(result$stable_features)
    ))
  }
  if (identical(property, "response")) {
    z <- result$summary
    fraction <- if (nrow(z)) mean(z$replication_fraction > 0) else 0
    overall <- if (nrow(z) && sum(z$candidates_tested) > 0) {
      sum(z$replicated) / sum(z$candidates_tested)
    } else {
      0
    }
    label <- if (fraction >= min_fraction && overall >= min_fraction) {
      "transferable response"
    } else if (overall > 0) {
      "context-dependent response"
    } else {
      "response not established"
    }
    return(data.frame(
      property = property,
      label = label,
      cohorts_with_replication_fraction = fraction,
      candidate_replication_fraction = overall
    ))
  }
  if (identical(property, "prediction")) {
    z <- result$summary
    passed <- z$balanced_accuracy > z$chance
    fraction <- if (length(passed)) mean(passed) else 0
    label <- if (fraction >= min_fraction) {
      "prediction transfers descriptively"
    } else if (any(passed)) {
      "context-dependent prediction"
    } else {
      "prediction not established"
    }
    return(data.frame(
      property = property,
      label = label,
      cohorts_above_chance = sum(passed),
      cohorts_tested = length(passed),
      above_chance_fraction = fraction,
      note = "Permutation or bootstrap uncertainty should be reported separately."
    ))
  }
  stop("Unsupported result property.", call. = FALSE)
}

#' Run the declared transferability workflow
#'
#' @param claim An `mt_claim` object.
#' @param counts Sample-by-feature counts.
#' @param metadata Sample metadata.
#' @param sample_col Sample identifier column.
#' @param ... Property-specific arguments.
#' @return A validation result.
#' @export
run_transferability_workflow <- function(
    claim,
    counts,
    metadata,
    sample_col,
    ...) {
  if (!inherits(claim, "mt_claim")) {
    stop("`claim` must be created by `declare_transferability()`.",
         call. = FALSE)
  }
  if (claim$property == "membership") {
    return(validate_membership(
      counts = counts,
      metadata = metadata,
      sample_col = sample_col,
      study_col = claim$study_col,
      ...
    ))
  }
  if (claim$property == "response") {
    return(validate_response(
      counts = counts,
      metadata = metadata,
      sample_col = sample_col,
      study_col = claim$study_col,
      outcome_col = claim$outcome_col,
      ...
    ))
  }
  validate_prediction(
    counts = counts,
    metadata = metadata,
    sample_col = sample_col,
    study_col = claim$study_col,
    outcome_col = claim$outcome_col,
    ...
  )
}
