.prevalence_matrix <- function(counts, studies) {
  lev <- unique(as.character(studies))
  out <- vapply(
    lev,
    function(s) colMeans(counts[as.character(studies) == s, , drop = FALSE] > 0),
    numeric(ncol(counts))
  )
  if (is.null(dim(out))) out <- matrix(out, ncol = 1L)
  rownames(out) <- colnames(counts)
  colnames(out) <- lev
  out
}

#' Discover study-aware membership candidates
#'
#' @param counts Sample-by-feature counts.
#' @param metadata Sample metadata.
#' @param sample_col Sample identifier column.
#' @param study_col Cohort column.
#' @param heldout Optional cohort excluded from discovery.
#' @param prevalence Minimum within-cohort prevalence.
#' @param min_total Minimum total count across discovery samples.
#' @param min_studies Number of discovery cohorts required to pass. The default
#'   requires all discovery cohorts.
#' @return An `mt_membership_discovery` object.
#' @export
discover_membership <- function(
    counts,
    metadata,
    sample_col,
    study_col,
    heldout = NULL,
    prevalence = 0.50,
    min_total = 0,
    min_studies = NULL) {
  dat <- validate_microbiome_input(
    counts, metadata, sample_col, study_col, allow_zero_depth = FALSE
  )
  counts <- dat$counts
  metadata <- dat$metadata
  use <- if (is.null(heldout)) {
    rep(TRUE, nrow(metadata))
  } else {
    as.character(metadata[[study_col]]) != as.character(heldout)
  }
  if (!any(use)) stop("No discovery samples remain.", call. = FALSE)
  training_studies <- unique(as.character(metadata[[study_col]][use]))
  if (length(training_studies) < 1L) {
    stop("No discovery cohort remains.", call. = FALSE)
  }
  prev <- .prevalence_matrix(counts[use, , drop = FALSE],
                             metadata[[study_col]][use])
  if (is.null(min_studies)) min_studies <- ncol(prev)
  min_studies <- as.integer(min_studies)
  if (min_studies < 1L || min_studies > ncol(prev)) {
    stop("`min_studies` must be between 1 and the number of discovery cohorts.",
         call. = FALSE)
  }
  totals <- colSums(counts[use, , drop = FALSE])
  pass <- rowSums(prev >= prevalence, na.rm = TRUE) >= min_studies &
    totals[rownames(prev)] >= min_total
  structure(
    list(
      heldout = heldout,
      candidates = rownames(prev)[pass],
      prevalence = prev,
      min_prevalence = apply(prev, 1L, min),
      mean_prevalence = rowMeans(prev),
      total_counts = totals[rownames(prev)],
      threshold = prevalence,
      min_studies = min_studies,
      training_studies = training_studies
    ),
    class = "mt_membership_discovery"
  )
}

.quantile_bin <- function(x, n = 4L) {
  ok <- is.finite(x)
  if (sum(ok) < 2L || length(unique(x[ok])) <= 1L) {
    return(rep("all", length(x)))
  }
  breaks <- unique(as.numeric(stats::quantile(
    x[ok], probs = seq(0, 1, length.out = n + 1L),
    na.rm = TRUE, names = FALSE, type = 7
  )))
  if (length(breaks) < 3L) return(rep("all", length(x)))
  as.character(cut(x, breaks = breaks, include.lowest = TRUE))
}

.sample_matched_set <- function(core, universe, strata) {
  available <- setdiff(universe, core)
  if (length(available) < length(core)) available <- universe
  available_strata <- as.character(strata[available])
  core_strata <- as.character(strata[core])
  selected <- rep(NA_character_, length(core))

  # Draw only as many features as each stratum needs. Sequential uniform
  # sampling without replacement and a single sample(pool, n) have the same
  # joint distribution, while the latter avoids repeated universe scans.
  for (stratum in unique(core_strata)) {
    positions <- which(core_strata == stratum)
    pool <- available[available_strata == stratum]
    n_take <- min(length(positions), length(pool))
    if (n_take > 0L) {
      selected[positions[seq_len(n_take)]] <-
        sample(pool, n_take, replace = FALSE)
    }
  }

  missing <- which(is.na(selected))
  if (length(missing)) {
    fallback <- setdiff(available, selected[!is.na(selected)])
    if (length(fallback) >= length(missing)) {
      selected[missing] <- sample(fallback, length(missing), replace = FALSE)
    } else {
      selected[missing] <- sample(universe, length(missing), replace = TRUE)
    }
  }
  selected
}

#' Leave-one-study-out validation of membership
#'
#' @param counts Sample-by-feature counts.
#' @param metadata Sample metadata.
#' @param sample_col Sample identifier column.
#' @param study_col Cohort column.
#' @param prevalence Minimum prevalence in discovery and held-out cohorts.
#' @param min_total Minimum total count in discovery samples.
#' @param n_null Number of random-subset draws. Set to zero to skip nulls.
#' @param matched Include prevalence/abundance-matched nulls.
#' @param seed Random seed.
#' @return An `mt_membership_validation` object.
#' @export
validate_membership <- function(
    counts,
    metadata,
    sample_col,
    study_col,
    prevalence = 0.50,
    min_total = 0,
    n_null = 999,
    matched = TRUE,
    seed = 1L) {
  dat <- validate_microbiome_input(
    counts, metadata, sample_col, study_col, allow_zero_depth = FALSE
  )
  counts <- dat$counts
  metadata <- dat$metadata
  studies <- sort(unique(as.character(metadata[[study_col]])))
  if (length(studies) < 2L) {
    stop("Membership validation requires at least two cohorts.",
         call. = FALSE)
  }
  n_null <- as.integer(n_null)
  set.seed(seed)
  summary_rows <- list()
  feature_rows <- list()
  null_rows <- list()

  for (heldout in studies) {
    discovery <- discover_membership(
      counts, metadata, sample_col, study_col,
      heldout = heldout, prevalence = prevalence, min_total = min_total
    )
    core <- discovery$candidates
    test_ids <- as.character(metadata[[study_col]]) == heldout
    test_prev <- colMeans(counts[test_ids, , drop = FALSE] > 0)
    validated <- core[test_prev[core] >= prevalence]
    core_n <- length(core)
    observed_n <- length(validated)
    observed_rate <- if (core_n) observed_n / core_n else NA_real_

    simple_rate <- matched_rate <- numeric()
    if (n_null > 0L && core_n > 0L) {
      universe <- colnames(counts)
      prev_bin <- .quantile_bin(discovery$min_prevalence[universe])
      abundance_bin <- .quantile_bin(
        log10(discovery$total_counts[universe] + 1)
      )
      strata <- paste(prev_bin, abundance_bin, sep = "|")
      names(strata) <- universe
      simple_rate <- numeric(n_null)
      matched_rate <- numeric(n_null)
      for (i in seq_len(n_null)) {
        simple <- sample(universe, core_n, replace = FALSE)
        simple_rate[i] <- mean(test_prev[simple] >= prevalence)
        if (isTRUE(matched)) {
          matched_set <- .sample_matched_set(core, universe, strata)
          matched_rate[i] <- mean(test_prev[matched_set] >= prevalence)
        }
      }
      null_rows[[length(null_rows) + 1L]] <- rbind(
        data.frame(
          heldout_study = heldout,
          permutation = seq_len(n_null),
          null_type = "simple_random",
          validation_rate = simple_rate
        ),
        if (isTRUE(matched)) data.frame(
          heldout_study = heldout,
          permutation = seq_len(n_null),
          null_type = "prevalence_abundance_matched",
          validation_rate = matched_rate
        ) else NULL
      )
    }

    null_for_test <- if (isTRUE(matched) && length(matched_rate)) {
      matched_rate
    } else {
      simple_rate
    }
    null_q95 <- if (length(null_for_test)) {
      unname(stats::quantile(null_for_test, 0.95))
    } else {
      NA_real_
    }
    empirical_p <- if (length(null_for_test) && is.finite(observed_rate)) {
      (sum(null_for_test >= observed_rate) + 1) /
        (length(null_for_test) + 1)
    } else {
      NA_real_
    }
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      heldout_study = heldout,
      n_train = sum(!test_ids),
      n_test = sum(test_ids),
      features = ncol(counts),
      candidates = core_n,
      validated = observed_n,
      validation_rate = observed_rate,
      null_q95 = null_q95,
      empirical_p = empirical_p,
      above_null_q95 = is.finite(observed_rate) &&
        is.finite(null_q95) && observed_rate > null_q95
    )
    if (core_n) {
      feature_rows[[length(feature_rows) + 1L]] <- data.frame(
        heldout_study = heldout,
        feature = core,
        training_min_prevalence = discovery$min_prevalence[core],
        training_mean_prevalence = discovery$mean_prevalence[core],
        heldout_prevalence = test_prev[core],
        validated = test_prev[core] >= prevalence
      )
    }
  }

  all_prev <- .prevalence_matrix(counts, metadata[[study_col]])
  stable <- rownames(all_prev)[apply(all_prev >= prevalence, 1L, all)]
  structure(
    list(
      property = "membership",
      summary = do.call(rbind, summary_rows),
      feature_results = if (length(feature_rows)) {
        do.call(rbind, feature_rows)
      } else {
        data.frame()
      },
      null = if (length(null_rows)) do.call(rbind, null_rows) else data.frame(),
      stable_features = stable,
      prevalence = all_prev,
      threshold = prevalence,
      settings = list(
        min_total = min_total,
        n_null = n_null,
        matched = matched,
        seed = seed
      )
    ),
    class = "mt_membership_validation"
  )
}

#' @export
print.mt_membership_validation <- function(x, ...) {
  cat("<mt_membership_validation>\n")
  cat("  cohorts:          ", nrow(x$summary), "\n", sep = "")
  cat("  stable features:  ", length(x$stable_features), "\n", sep = "")
  cat("  prevalence:       ", x$threshold, "\n", sep = "")
  if (nrow(x$summary)) {
    cat("  validation range: ",
        paste(round(range(x$summary$validation_rate, na.rm = TRUE), 3),
              collapse = " - "), "\n", sep = "")
  }
  invisible(x)
}

#' Compare pooled and study-aware membership definitions
#'
#' @param counts Sample-by-feature counts.
#' @param metadata Sample metadata.
#' @param sample_col Sample identifier column.
#' @param study_col Cohort column.
#' @param thresholds Prevalence thresholds.
#' @return A data frame quantifying pooled-only candidates.
#' @export
pooled_membership_contrast <- function(
    counts,
    metadata,
    sample_col,
    study_col,
    thresholds = c(0.10, 0.30, 0.50)) {
  dat <- validate_microbiome_input(
    counts, metadata, sample_col, study_col, allow_zero_depth = FALSE
  )
  counts <- dat$counts
  metadata <- dat$metadata
  pooled_prev <- colMeans(counts > 0)
  study_prev <- .prevalence_matrix(counts, metadata[[study_col]])
  out <- lapply(thresholds, function(threshold) {
    pooled <- names(pooled_prev)[pooled_prev >= threshold]
    aware <- rownames(study_prev)[apply(
      study_prev >= threshold, 1L, all
    )]
    pooled_only <- setdiff(pooled, aware)
    data.frame(
      threshold = threshold,
      pooled_candidates = length(pooled),
      study_aware_candidates = length(aware),
      pooled_only_candidates = length(pooled_only),
      pooled_only_fraction = if (length(pooled)) {
        length(pooled_only) / length(pooled)
      } else {
        NA_real_
      }
    )
  })
  do.call(rbind, out)
}

#' Membership sensitivity across prevalence and depth choices
#'
#' @param counts Sample-by-feature counts.
#' @param metadata Sample metadata.
#' @param sample_col Sample identifier column.
#' @param study_col Cohort column.
#' @param prevalence Prevalence thresholds.
#' @param min_depth Minimum sample-depth thresholds.
#' @param min_total Minimum discovery total-count thresholds.
#' @param n_null Number of null draws per setting.
#' @param matched Use prevalence/abundance-matched nulls.
#' @param seed Random seed.
#' @return A data frame with one row per setting and held-out cohort.
#' @export
membership_sensitivity <- function(
    counts,
    metadata,
    sample_col,
    study_col,
    prevalence = c(0.30, 0.50, 0.70),
    min_depth = c(0, 5000),
    min_total = c(0, 20),
    n_null = 0,
    matched = TRUE,
    seed = 1L) {
  counts <- as.matrix(counts)
  settings <- expand.grid(
    prevalence = prevalence,
    min_depth = min_depth,
    min_total = min_total,
    stringsAsFactors = FALSE
  )
  rows <- list()
  for (i in seq_len(nrow(settings))) {
    setting <- settings[i, , drop = FALSE]
    keep_samples <- rowSums(counts) >= setting$min_depth
    if (sum(keep_samples) < 2L) next
    fit <- try(validate_membership(
      counts[keep_samples, , drop = FALSE],
      metadata,
      sample_col,
      study_col,
      prevalence = setting$prevalence,
      min_total = setting$min_total,
      n_null = n_null,
      matched = matched,
      seed = seed + i
    ), silent = TRUE)
    if (inherits(fit, "try-error")) next
    z <- fit$summary
    z$prevalence_threshold <- setting$prevalence
    z$minimum_sample_depth <- setting$min_depth
    z$minimum_training_total <- setting$min_total
    z$stable_features <- length(fit$stable_features)
    rows[[length(rows) + 1L]] <- z
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}
