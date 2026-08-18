#' Balanced accuracy
#'
#' @param truth Observed class labels.
#' @param predicted Predicted class labels.
#' @return Mean class-specific recall.
#' @export
balanced_accuracy <- function(truth, predicted) {
  truth <- factor(truth)
  predicted <- factor(predicted, levels = levels(truth))
  mean(vapply(levels(truth), function(level) {
    ids <- truth == level
    if (!any(ids)) return(NA_real_)
    mean(predicted[ids] == level)
  }, numeric(1)), na.rm = TRUE)
}

#' Macro-averaged F1 score
#'
#' @inheritParams balanced_accuracy
#' @return Mean class-specific F1 score.
#' @export
macro_f1 <- function(truth, predicted) {
  truth <- factor(truth)
  predicted <- factor(predicted, levels = levels(truth))
  mean(vapply(levels(truth), function(level) {
    tp <- sum(truth == level & predicted == level)
    fp <- sum(truth != level & predicted == level)
    fn <- sum(truth == level & predicted != level)
    precision <- if ((tp + fp) == 0) 0 else tp / (tp + fp)
    recall <- if ((tp + fn) == 0) 0 else tp / (tp + fn)
    if ((precision + recall) == 0) 0 else {
      2 * precision * recall / (precision + recall)
    }
  }, numeric(1)))
}

#' Fit a nearest-centroid classifier
#'
#' This dependency-free classifier is provided for examples and diagnostics.
#' Publication analyses may supply any `fit_fun` and `predict_fun`, including
#' elastic-net or random-forest models.
#'
#' @param x Numeric training matrix.
#' @param y Training labels.
#' @param ... Unused.
#' @return A nearest-centroid model.
#' @export
fit_nearest_centroid <- function(x, y, ...) {
  y <- factor(y)
  centers <- do.call(rbind, lapply(levels(y), function(level) {
    colMeans(x[y == level, , drop = FALSE])
  }))
  rownames(centers) <- levels(y)
  structure(list(centers = centers, levels = levels(y)),
            class = "mt_nearest_centroid")
}

#' Predict from a nearest-centroid model
#'
#' @param model Model returned by `fit_nearest_centroid()`.
#' @param newx Numeric test matrix.
#' @param ... Unused.
#' @return Predicted class labels.
#' @export
predict_nearest_centroid <- function(model, newx, ...) {
  distance <- vapply(seq_len(nrow(model$centers)), function(i) {
    rowSums((sweep(newx, 2L, model$centers[i, ], FUN = "-"))^2)
  }, numeric(nrow(newx)))
  if (is.null(dim(distance))) distance <- matrix(distance, ncol = 1L)
  model$levels[max.col(-distance, ties.method = "first")]
}

.swap_labels_preserve_counts <- function(y, fraction, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  y <- as.character(y)
  if (fraction <= 0) return(y)
  target <- max(2L, as.integer(round(length(y) * fraction)))
  best <- y
  best_diff <- Inf
  for (attempt in seq_len(200L)) {
    size <- min(length(y), max(2L, target))
    ids <- sample(seq_along(y), size)
    proposal <- y
    proposal[ids] <- sample(y[ids], replace = FALSE)
    changed <- sum(proposal != y)
    diff <- abs(changed - target)
    if (diff < best_diff) {
      best <- proposal
      best_diff <- diff
    }
    if (best_diff == 0L) break
  }
  best
}

.prediction_once <- function(
    counts,
    metadata,
    outcome_col,
    study_col,
    features,
    transform,
    pseudocount,
    min_prevalence,
    min_total,
    fit_fun,
    predict_fun,
    noise_fraction = 0,
    seed = 1L,
    ...) {
  studies <- sort(unique(as.character(metadata[[study_col]])))
  rows <- list()
  predictions <- list()
  for (i in seq_along(studies)) {
    heldout <- studies[i]
    test <- as.character(metadata[[study_col]]) == heldout
    train <- !test & !is.na(metadata[[outcome_col]])
    test <- test & !is.na(metadata[[outcome_col]])
    y_train <- factor(metadata[[outcome_col]][train])
    y_test <- factor(metadata[[outcome_col]][test], levels = levels(y_train))
    if (length(unique(y_train)) < 2L || anyNA(y_test) || !sum(test)) next
    use_features <- if (is.null(features)) {
      keep <- colMeans(counts[train, , drop = FALSE] > 0) >= min_prevalence &
        colSums(counts[train, , drop = FALSE]) >= min_total
      colnames(counts)[keep]
    } else {
      intersect(features, colnames(counts))
    }
    if (!length(use_features)) next
    x <- .transform_counts(
      counts[, use_features, drop = FALSE], transform, pseudocount
    )
    noisy <- .swap_labels_preserve_counts(
      y_train, noise_fraction, seed + i
    )
    model <- fit_fun(x[train, , drop = FALSE], factor(noisy), ...)
    pred <- predict_fun(model, x[test, , drop = FALSE], ...)
    pred <- factor(pred, levels = levels(y_train))
    rows[[length(rows) + 1L]] <- data.frame(
      heldout_study = heldout,
      n_train = sum(train),
      n_test = sum(test),
      features = length(use_features),
      balanced_accuracy = balanced_accuracy(y_test, pred),
      macro_f1 = macro_f1(y_test, pred),
      chance = 1 / nlevels(y_train),
      target_noise = noise_fraction,
      realized_noise = mean(noisy != as.character(y_train))
    )
    predictions[[length(predictions) + 1L]] <- data.frame(
      heldout_study = heldout,
      sample_id = rownames(counts)[test],
      observed = as.character(y_test),
      predicted = as.character(pred)
    )
  }
  list(
    summary = if (length(rows)) do.call(rbind, rows) else data.frame(),
    predictions = if (length(predictions)) {
      do.call(rbind, predictions)
    } else {
      data.frame()
    }
  )
}

#' Held-out-study prediction validation
#'
#' @param counts Sample-by-feature counts.
#' @param metadata Sample metadata.
#' @param sample_col Sample identifier column.
#' @param study_col Cohort column.
#' @param outcome_col Categorical outcome column.
#' @param features Optional prespecified features.
#' @param transform Feature transform.
#' @param pseudocount CLR pseudocount.
#' @param min_prevalence Training prevalence filter.
#' @param min_total Training total-count filter.
#' @param fit_fun Function with arguments `(x, y, ...)`.
#' @param predict_fun Function with arguments `(model, newx, ...)`.
#' @param ... Additional arguments passed to model functions.
#' @return An `mt_prediction_validation` object.
#' @export
validate_prediction <- function(
    counts,
    metadata,
    sample_col,
    study_col,
    outcome_col,
    features = NULL,
    transform = "clr",
    pseudocount = 0.5,
    min_prevalence = 0.10,
    min_total = 20,
    fit_fun = fit_nearest_centroid,
    predict_fun = predict_nearest_centroid,
    ...) {
  dat <- validate_microbiome_input(
    counts, metadata, sample_col, study_col, outcome_col
  )
  ans <- .prediction_once(
    dat$counts, dat$metadata, outcome_col, study_col, features,
    transform, pseudocount, min_prevalence, min_total,
    fit_fun, predict_fun, noise_fraction = 0, seed = 1L, ...
  )
  structure(
    c(
      list(property = "prediction"),
      ans,
      list(settings = list(
        features = features,
        transform = transform,
        pseudocount = pseudocount,
        min_prevalence = min_prevalence,
        min_total = min_total
      ))
    ),
    class = "mt_prediction_validation"
  )
}

#' Within-cohort permutation test for held-out prediction
#'
#' Outcome labels are permuted separately within each cohort. The complete
#' held-out workflow is then repeated against the permuted labels, preserving
#' cohort sizes and cohort-specific class counts.
#'
#' @inheritParams validate_prediction
#' @param permutations Number of label permutations.
#' @param seed Random seed.
#' @return A list with observed performance, null distributions, and empirical
#'   P values.
#' @export
prediction_permutation_test <- function(
    counts,
    metadata,
    sample_col,
    study_col,
    outcome_col,
    features = NULL,
    transform = "clr",
    pseudocount = 0.5,
    min_prevalence = 0.10,
    min_total = 20,
    fit_fun = fit_nearest_centroid,
    predict_fun = predict_nearest_centroid,
    permutations = 999L,
    seed = 1L,
    ...) {
  dat <- validate_microbiome_input(
    counts, metadata, sample_col, study_col, outcome_col
  )
  observed <- validate_prediction(
    dat$counts, dat$metadata, sample_col, study_col, outcome_col,
    features, transform, pseudocount, min_prevalence, min_total,
    fit_fun, predict_fun, ...
  )
  observed_ba <- mean(observed$summary$balanced_accuracy)
  observed_f1 <- mean(observed$summary$macro_f1)
  permutations <- as.integer(permutations)
  null <- vector("list", permutations)
  set.seed(seed)
  study <- as.character(dat$metadata[[study_col]])
  original_outcome <- dat$metadata[[outcome_col]]
  for (i in seq_len(permutations)) {
    permuted <- dat$metadata
    shuffled <- as.character(original_outcome)
    for (cohort in unique(study)) {
      ids <- which(study == cohort)
      shuffled[ids] <- sample(shuffled[ids], replace = FALSE)
    }
    permuted[[outcome_col]] <- if (is.factor(original_outcome)) {
      factor(shuffled, levels = levels(original_outcome))
    } else {
      shuffled
    }
    fit <- validate_prediction(
      dat$counts, permuted, sample_col, study_col, outcome_col,
      features, transform, pseudocount, min_prevalence, min_total,
      fit_fun, predict_fun, ...
    )
    null[[i]] <- data.frame(
      permutation = i,
      balanced_accuracy = mean(fit$summary$balanced_accuracy),
      macro_f1 = mean(fit$summary$macro_f1)
    )
  }
  null <- if (length(null)) do.call(rbind, null) else data.frame()
  list(
    property = "prediction",
    observed = observed,
    null = null,
    observed_balanced_accuracy = observed_ba,
    observed_macro_f1 = observed_f1,
    balanced_accuracy_p = if (nrow(null)) {
      (sum(null$balanced_accuracy >= observed_ba) + 1) / (nrow(null) + 1)
    } else {
      NA_real_
    },
    macro_f1_p = if (nrow(null)) {
      (sum(null$macro_f1 >= observed_f1) + 1) / (nrow(null) + 1)
    } else {
      NA_real_
    }
  )
}

#' @export
print.mt_prediction_validation <- function(x, ...) {
  cat("<mt_prediction_validation>\n")
  cat("  held-out cohorts: ", nrow(x$summary), "\n", sep = "")
  if (nrow(x$summary)) {
    cat("  mean balanced accuracy: ",
        round(mean(x$summary$balanced_accuracy), 3), "\n", sep = "")
    cat("  mean macro-F1: ",
        round(mean(x$summary$macro_f1), 3), "\n", sep = "")
  }
  invisible(x)
}

#' Class-preserving label-noise sensitivity
#'
#' @inheritParams validate_prediction
#' @param noise Fractions of training labels targeted for alteration.
#' @param replicates Replicates per nonzero noise level.
#' @param seed Random seed.
#' @return A list with replicate-level and summarized performance.
#' @export
label_noise_sensitivity <- function(
    counts,
    metadata,
    sample_col,
    study_col,
    outcome_col,
    features = NULL,
    transform = "clr",
    pseudocount = 0.5,
    min_prevalence = 0.10,
    min_total = 20,
    fit_fun = fit_nearest_centroid,
    predict_fun = predict_nearest_centroid,
    noise = c(0, 0.05, 0.10, 0.15, 0.20),
    replicates = 100L,
    seed = 1L,
    ...) {
  dat <- validate_microbiome_input(
    counts, metadata, sample_col, study_col, outcome_col
  )
  rows <- list()
  index <- 0L
  for (fraction in noise) {
    reps <- if (fraction == 0) 1L else as.integer(replicates)
    for (replicate in seq_len(reps)) {
      index <- index + 1L
      ans <- .prediction_once(
        dat$counts, dat$metadata, outcome_col, study_col, features,
        transform, pseudocount, min_prevalence, min_total,
        fit_fun, predict_fun, fraction, seed + index, ...
      )
      if (!nrow(ans$summary)) next
      rows[[length(rows) + 1L]] <- data.frame(
        target_noise = fraction,
        replicate = replicate,
        realized_noise = mean(ans$summary$realized_noise),
        balanced_accuracy = mean(ans$summary$balanced_accuracy),
        macro_f1 = mean(ans$summary$macro_f1),
        chance = mean(ans$summary$chance)
      )
    }
  }
  raw <- if (length(rows)) do.call(rbind, rows) else data.frame()
  summary <- if (nrow(raw)) {
    do.call(rbind, lapply(split(raw, raw$target_noise), function(z) {
      data.frame(
        target_noise = z$target_noise[1],
        mean_realized_noise = mean(z$realized_noise),
        n_replicates = nrow(z),
        mean_balanced_accuracy = mean(z$balanced_accuracy),
        BA_q025 = unname(stats::quantile(z$balanced_accuracy, 0.025)),
        BA_median = stats::median(z$balanced_accuracy),
        BA_q975 = unname(stats::quantile(z$balanced_accuracy, 0.975)),
        mean_macro_f1 = mean(z$macro_f1),
        fraction_BA_above_chance = mean(z$balanced_accuracy > z$chance)
      )
    }))
  } else {
    data.frame()
  }
  list(raw = raw, summary = summary)
}
