.feature_effects <- function(
    counts,
    metadata,
    outcome_col,
    study_col,
    covariates = NULL,
    transform = "clr",
    pseudocount = 0.5,
    min_prevalence = 0.10,
    min_total = 20) {
  if (!is.numeric(metadata[[outcome_col]])) {
    stop("Response validation currently requires a numeric outcome.",
         call. = FALSE)
  }
  z <- .transform_counts(counts, transform, pseudocount)
  studies <- unique(as.character(metadata[[study_col]]))
  rows <- list()
  for (study in studies) {
    ids <- as.character(metadata[[study_col]]) == study &
      is.finite(metadata[[outcome_col]])
    if (sum(ids) < 5L) next
    keep <- colMeans(counts[ids, , drop = FALSE] > 0) >= min_prevalence &
      colSums(counts[ids, , drop = FALSE]) >= min_total
    features <- colnames(counts)[keep]
    if (!length(features)) next
    md <- metadata[ids, , drop = FALSE]
    md$.mt_outcome <- as.numeric(scale(md[[outcome_col]]))
    usable_covariates <- covariates[
      covariates %in% names(md) &
        vapply(covariates, function(v) {
          length(unique(md[[v]][!is.na(md[[v]])])) > 1L
        }, logical(1))
    ]
    rhs <- c(".mt_outcome", usable_covariates)
    for (feature in features) {
      md$.mt_response <- z[ids, feature]
      fit <- try(stats::lm(
        stats::reformulate(rhs, response = ".mt_response"),
        data = md
      ), silent = TRUE)
      if (inherits(fit, "try-error")) next
      co <- summary(fit)$coefficients
      if (!".mt_outcome" %in% rownames(co)) next
      rows[[length(rows) + 1L]] <- data.frame(
        study = study,
        feature = feature,
        n = sum(ids),
        estimate = co[".mt_outcome", "Estimate"],
        standard_error = co[".mt_outcome", "Std. Error"],
        statistic = co[".mt_outcome", "t value"],
        p = co[".mt_outcome", "Pr(>|t|)"],
        adjustment = paste(usable_covariates, collapse = ";")
      )
    }
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

.pool_random_effects <- function(x) {
  x <- x[is.finite(x$estimate) & is.finite(x$standard_error) &
           x$standard_error > 0, , drop = FALSE]
  k <- nrow(x)
  if (k < 2L) return(NULL)
  wi <- 1 / x$standard_error^2
  fixed <- sum(wi * x$estimate) / sum(wi)
  q <- sum(wi * (x$estimate - fixed)^2)
  c_term <- sum(wi) - sum(wi^2) / sum(wi)
  tau2 <- max(0, (q - (k - 1)) / c_term)
  wr <- 1 / (x$standard_error^2 + tau2)
  estimate <- sum(wr * x$estimate) / sum(wr)
  se <- sqrt(1 / sum(wr))
  z <- estimate / se
  data.frame(
    k = k,
    estimate = estimate,
    standard_error = se,
    ci_lower = estimate - stats::qnorm(0.975) * se,
    ci_upper = estimate + stats::qnorm(0.975) * se,
    p = 2 * stats::pnorm(abs(z), lower.tail = FALSE),
    tau2 = tau2,
    I2 = if (q > 0) max(0, (q - (k - 1)) / q) * 100 else 0,
    Q = q,
    Q_p = stats::pchisq(q, df = k - 1, lower.tail = FALSE)
  )
}

.pool_effect_table <- function(effects) {
  if (!nrow(effects)) return(data.frame())
  groups <- split(effects, effects$feature)
  pooled <- lapply(names(groups), function(feature) {
    ans <- .pool_random_effects(groups[[feature]])
    if (is.null(ans)) return(NULL)
    data.frame(feature = feature, ans, stringsAsFactors = FALSE)
  })
  pooled <- pooled[!vapply(pooled, is.null, logical(1))]
  if (!length(pooled)) return(data.frame())
  out <- do.call(rbind, pooled)
  out$q <- stats::p.adjust(out$p, method = "BH")
  rownames(out) <- NULL
  out
}

#' Discover study-aware response candidates
#'
#' @param counts Sample-by-feature counts.
#' @param metadata Sample metadata.
#' @param sample_col Sample identifier column.
#' @param study_col Cohort column.
#' @param outcome_col Numeric outcome column.
#' @param discovery_studies Cohorts used for discovery.
#' @param covariates Optional model covariates.
#' @param transform Feature transform.
#' @param pseudocount CLR pseudocount.
#' @param min_prevalence Minimum within-study prevalence.
#' @param min_total Minimum within-study total count.
#' @param alpha FDR threshold.
#' @return A list with study effects, pooled effects, and candidates.
#' @export
discover_response <- function(
    counts,
    metadata,
    sample_col,
    study_col,
    outcome_col,
    discovery_studies = NULL,
    covariates = NULL,
    transform = "clr",
    pseudocount = 0.5,
    min_prevalence = 0.10,
    min_total = 20,
    alpha = 0.05) {
  dat <- validate_microbiome_input(
    counts, metadata, sample_col, study_col, outcome_col
  )
  counts <- dat$counts
  metadata <- dat$metadata
  if (is.null(discovery_studies)) {
    discovery_studies <- unique(as.character(metadata[[study_col]]))
  }
  ids <- as.character(metadata[[study_col]]) %in% discovery_studies
  effects <- .feature_effects(
    counts[ids, , drop = FALSE],
    metadata[ids, , drop = FALSE],
    outcome_col,
    study_col,
    covariates,
    transform,
    pseudocount,
    min_prevalence,
    min_total
  )
  pooled <- .pool_effect_table(effects)
  candidates <- if (nrow(pooled)) pooled$feature[pooled$q <= alpha] else character()
  structure(
    list(
      study_effects = effects,
      pooled_effects = pooled,
      candidates = candidates,
      discovery_studies = discovery_studies,
      alpha = alpha,
      settings = list(
        transform = transform,
        pseudocount = pseudocount,
        min_prevalence = min_prevalence,
        min_total = min_total,
        covariates = covariates
      )
    ),
    class = "mt_response_discovery"
  )
}

#' Validate response candidates by independent projection
#'
#' @inheritParams discover_response
#' @param validation_studies Cohorts reserved for validation.
#' @param nominal_alpha Per-study validation threshold.
#' @return An `mt_response_validation` object.
#' @export
validate_response <- function(
    counts,
    metadata,
    sample_col,
    study_col,
    outcome_col,
    discovery_studies,
    validation_studies,
    covariates = NULL,
    transform = "clr",
    pseudocount = 0.5,
    min_prevalence = 0.10,
    min_total = 20,
    alpha = 0.05,
    nominal_alpha = 0.05) {
  discovery <- discover_response(
    counts, metadata, sample_col, study_col, outcome_col,
    discovery_studies, covariates, transform, pseudocount,
    min_prevalence, min_total, alpha
  )
  dat <- validate_microbiome_input(
    counts, metadata, sample_col, study_col, outcome_col
  )
  ids <- as.character(dat$metadata[[study_col]]) %in% validation_studies
  validation_effects <- .feature_effects(
    dat$counts[ids, , drop = FALSE],
    dat$metadata[ids, , drop = FALSE],
    outcome_col,
    study_col,
    covariates,
    transform,
    pseudocount,
    min_prevalence,
    min_total
  )
  validation_effects <- validation_effects[
    validation_effects$feature %in% discovery$candidates, , drop = FALSE
  ]
  direction <- if (nrow(discovery$pooled_effects)) {
    stats::setNames(sign(discovery$pooled_effects$estimate),
                    discovery$pooled_effects$feature)
  } else {
    numeric()
  }
  if (nrow(validation_effects)) {
    validation_effects$discovery_direction <- direction[
      validation_effects$feature
    ]
    validation_effects$direction_concordant <-
      sign(validation_effects$estimate) ==
      validation_effects$discovery_direction
    validation_effects$replicated <-
      validation_effects$p <= nominal_alpha &
      validation_effects$direction_concordant
  }
  per_study <- if (nrow(validation_effects)) {
    do.call(rbind, lapply(
      split(validation_effects, validation_effects$study),
      function(z) data.frame(
        study = z$study[1],
        candidates_tested = nrow(z),
        replicated = sum(z$replicated),
        replication_fraction = mean(z$replicated)
      )
    ))
  } else {
    data.frame()
  }
  structure(
    list(
      property = "response",
      discovery = discovery,
      validation_effects = validation_effects,
      summary = per_study,
      validation_studies = validation_studies,
      nominal_alpha = nominal_alpha
    ),
    class = "mt_response_validation"
  )
}
