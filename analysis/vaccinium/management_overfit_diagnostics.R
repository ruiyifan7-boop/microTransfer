## Management overfit diagnostics for PRJEB110492.
## This script does not modify the existing transferability results.
## It evaluates nested training-only feature selection for k = 3,5,10,20,50,101,
## reports training-vs-held-out gaps, and treats the named three predictors as
## post-hoc diagnostics rather than unbiased performance estimates.

suppressPackageStartupMessages(library(glmnet))
options(stringsAsFactors = FALSE, width = 180)

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[[1]] else Sys.getenv("GLOBAL_HARMONIZED_ROOT", "global_harmonized")
root <- normalizePath(root, mustWork = TRUE)
outdir <- if (length(args) >= 2) args[[2]] else Sys.getenv(
  "MANAGEMENT_OVERFIT_OUTDIR",
  file.path(root, "enhancement", "management_overfit_diagnostics")
)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

ENET_ALPHA <- 0.5
CV_FOLDS <- as.integer(Sys.getenv("MANAGEMENT_CV_FOLDS", "9"))
BOOT_N <- as.integer(Sys.getenv("MANAGEMENT_OVERFIT_BOOT", "2000"))
if (CV_FOLDS < 2) stop("MANAGEMENT_CV_FOLDS must be >= 2")

project_root <- dirname(root)
metadata_path <- Sys.getenv(
  "MANAGEMENT_METADATA",
  file.path(project_root, "PRJEB110492", "metadata", "PRJEB110492_master_metadata_provisional.tsv")
)

pick_counts <- function(kind) {
  candidates <- c(
    file.path(root, "taxa_tables", paste0(kind, "_genus_counts_primary211.rds")),
    file.path(root, "analysis_ready", paste0(kind, "_genus_broad_counts.rds")),
    file.path(root, "analysis_ready", paste0(kind, "_genus_model_counts.rds"))
  )
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) stop("No count matrix found for ", kind)
  as.matrix(readRDS(hit[[1]]))
}

meta_path <- file.path(root, "global_metadata_primary_211.tsv")
if (!file.exists(meta_path)) stop("Missing metadata: ", meta_path)
if (!file.exists(metadata_path)) stop("Missing management metadata: ", metadata_path)

meta <- read.delim(meta_path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA"))
B0 <- pick_counts("bacteria")
F0 <- pick_counts("fungi")
shared <- Reduce(intersect, list(meta$sample_uid, rownames(B0), rownames(F0)))
meta <- meta[match(shared, meta$sample_uid), , drop = FALSE]
B0 <- B0[shared, , drop = FALSE]
F0 <- F0[shared, , drop = FALSE]
stopifnot(identical(meta$sample_uid, rownames(B0)), identical(rownames(B0), rownames(F0)))

coreB <- read.delim(file.path(root, "core_microbiome", "bacteria_core_audit.tsv"), check.names = FALSE)

prevalence_cols <- function(x) grep("^prevalence_", names(x), value = TRUE)
loso_core <- function(x, heldout) {
  cols <- setdiff(prevalence_cols(x), paste0("prevalence_", heldout))
  if (!length(cols)) stop("No discovery prevalence columns for ", heldout)
  keep <- apply(x[, cols, drop = FALSE], 1, function(z) all(!is.na(z) & z >= 0.50))
  unique(x$feature[keep])
}

clr_transform <- function(x) {
  z <- log(as.matrix(x) + 0.5)
  sweep(z, 1, rowMeans(z), "-")
}

class_weights <- function(y) {
  tab <- table(y)
  as.numeric(length(y) / (length(tab) * tab[as.character(y)]))
}

stratified_foldid <- function(y, k = 3, seed = 1) {
  set.seed(seed)
  y <- factor(y)
  k <- min(k, min(table(y)))
  if (k < 2) stop("Too few samples per class for inner CV")
  foldid <- integer(length(y))
  for (lev in levels(y)) {
    idx <- sample(which(y == lev))
    foldid[idx] <- rep(seq_len(k), length.out = length(idx))
  }
  foldid
}

balanced_accuracy <- function(truth, pred) {
  truth <- factor(truth)
  pred <- factor(pred, levels = levels(truth))
  mean(vapply(levels(truth), function(z) mean(pred[truth == z] == z), numeric(1)))
}

macro_f1 <- function(truth, pred) {
  truth <- factor(truth)
  pred <- factor(pred, levels = levels(truth))
  mean(vapply(levels(truth), function(z) {
    tp <- sum(truth == z & pred == z)
    fp <- sum(truth != z & pred == z)
    fn <- sum(truth == z & pred != z)
    precision <- if ((tp + fp) == 0) 0 else tp / (tp + fp)
    recall <- if ((tp + fn) == 0) 0 else tp / (tp + fn)
    if ((precision + recall) == 0) 0 else 2 * precision * recall / (precision + recall)
  }, numeric(1)))
}

## Management target: labels and blocks are used only here.
mm <- read.delim(metadata_path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA"))
mm$sample_uid <- paste("PRJEB110492", mm$biological_sample, sep = "__")
mm <- mm[
  mm$design_family == "transplanted_field_system" &
    mm$treatment_code_raw %in% c("HD", "LD", "wild") &
    mm$sample_uid %in% shared,
  , drop = FALSE
]
mm$treatment <- factor(mm$treatment_code_raw, levels = c("HD", "LD", "wild"))
mm$block <- factor(mm$block_code_raw)
if (nrow(mm) == 0 || nlevels(mm$block) < 2) stop("No usable management samples/blocks")
mm <- mm[order(mm$sample_uid), , drop = FALSE]

## The 101-feature pool is defined from the three non-target cohorts.
external_core <- loso_core(coreB, "PRJEB110492")
external_core <- intersect(external_core, colnames(B0))
rank_tab <- coreB[match(external_core, coreB$feature), , drop = FALSE]
rank_tab$minimum_prevalence <- as.numeric(rank_tab$minimum_prevalence)
rank_tab$geometric_mean_RA <- as.numeric(rank_tab$geometric_mean_RA)
rank_tab <- rank_tab[order(-rank_tab$minimum_prevalence, -rank_tab$geometric_mean_RA, rank_tab$feature), , drop = FALSE]
external_rank <- rank_tab$feature

## Named three predictors were selected after all four outer folds and are therefore post-hoc.
named3 <- c(
  coreB$feature[grep("g__Conexibacter$", coreB$feature)][1],
  coreB$feature[grep("g__Candidatus Solibacter$", coreB$feature)][1],
  coreB$feature[grep("f__67-14;g__Unclassified$", coreB$feature)][1]
)
named3 <- named3[!is.na(named3) & named3 %in% external_core]
if (length(named3) != 3) warning("Named three predictors not all found in external core: ", paste(named3, collapse = " | "))

## Training-only ranking for a genuinely nested sparse model.
rank_training_features <- function(z_train, y_train) {
  scores <- vapply(seq_len(ncol(z_train)), function(j) {
    d <- data.frame(value = z_train[, j], y = factor(y_train))
    fit <- tryCatch(aov(value ~ y, data = d), error = function(e) NULL)
    if (is.null(fit)) return(NA_real_)
    tab <- tryCatch(summary(fit)[[1]], error = function(e) NULL)
    if (is.null(tab) || !"F value" %in% colnames(tab)) return(NA_real_)
    as.numeric(tab[1, "F value"])
  }, numeric(1))
  scores[!is.finite(scores)] <- -Inf
  ord <- order(-scores, colnames(z_train))
  colnames(z_train)[ord]
}

fit_fold <- function(train_ids, test_ids, heldout, k, mode, seed) {
  y_train <- droplevels(mm$treatment[match(train_ids, mm$sample_uid)])
  y_test <- factor(mm$treatment[match(test_ids, mm$sample_uid)], levels = levels(mm$treatment))
  all_ids <- c(train_ids, test_ids)
  z_all <- clr_transform(B0[all_ids, external_rank, drop = FALSE])
  rownames(z_all) <- all_ids

  if (mode == "nested_topk") {
    ordered <- rank_training_features(z_all[train_ids, , drop = FALSE], y_train)
    selected <- ordered[seq_len(min(k, length(ordered)))]
    scope <- "nested_training_only"
  } else if (mode == "external_fixed") {
    selected <- external_rank[seq_len(min(k, length(external_rank)))]
    scope <- "external_prevalence_rank"
  } else if (mode == "named3_posthoc") {
    selected <- named3
    scope <- "posthoc_all_blocks_named3"
  } else stop("Unknown mode: ", mode)

  x_train <- z_all[train_ids, selected, drop = FALSE]
  x_test <- z_all[test_ids, selected, drop = FALSE]
  foldid <- stratified_foldid(y_train, k = CV_FOLDS, seed = seed)
  cvfit <- cv.glmnet(
    x = x_train, y = y_train, family = "multinomial",
    alpha = ENET_ALPHA, foldid = foldid,
    weights = class_weights(y_train), standardize = TRUE,
    type.measure = "class", nlambda = 80,
    control = list(maxit = 1e6)
  )
  lambda <- cvfit$lambda.1se
  pred_train <- as.vector(predict(cvfit, newx = x_train, s = "lambda.1se", type = "class"))
  pred_test <- as.vector(predict(cvfit, newx = x_test, s = "lambda.1se", type = "class"))
  cc <- coef(cvfit$glmnet.fit, s = lambda)
  nz <- sum(vapply(cc, function(m) {
    mmx <- as.matrix(m)
    sum(rownames(mmx) != "(Intercept)" & mmx[, 1] != 0)
  }, numeric(1)))
  list(
    fold = heldout, mode = mode, scope = scope, k = length(selected),
    selected = selected, lambda = lambda, nonzero = nz,
    train_ba = balanced_accuracy(y_train, pred_train),
    train_f1 = macro_f1(y_train, pred_train),
    test_ba = balanced_accuracy(y_test, pred_test),
    test_f1 = macro_f1(y_test, pred_test),
    train_pred = data.frame(sample_uid = train_ids, truth = as.character(y_train), predicted = pred_train),
    test_pred = data.frame(sample_uid = test_ids, truth = as.character(y_test), predicted = pred_test)
  )
}

ks <- c(3L, 5L, 10L, 20L, 50L, 101L)
fold_metrics <- list()
predictions <- list()
model_results <- list()
counter <- 0L

for (mode in c("nested_topk", "external_fixed")) {
  mode_ks <- if (mode == "nested_topk") ks else 101L
  for (k in mode_ks) {
    key <- paste(mode, k, sep = "__")
    fold_out <- list()
    pred_out <- list()
    for (heldout in levels(mm$block)) {
      train_ids <- mm$sample_uid[mm$block != heldout]
      test_ids <- mm$sample_uid[mm$block == heldout]
      fit <- fit_fold(train_ids, test_ids, heldout, k, mode, seed = 7000 + which(levels(mm$block) == heldout) + k)
      counter <- counter + 1L
      fold_out[[length(fold_out) + 1L]] <- data.frame(
        model = key, selection_scope = fit$scope, k = fit$k, heldout_block = heldout,
        n_train = length(train_ids), n_test = length(test_ids),
        train_balanced_accuracy = fit$train_ba, train_macro_f1 = fit$train_f1,
        heldout_balanced_accuracy = fit$test_ba, heldout_macro_f1 = fit$test_f1,
        train_test_BA_gap = fit$train_ba - fit$test_ba,
        lambda = fit$lambda, nonzero_coefficients = fit$nonzero,
        selected_features = paste(fit$selected, collapse = "||")
      )
      pred_out[[length(pred_out) + 1L]] <- rbind(
        data.frame(model = key, selection_scope = fit$scope, k = fit$k, split = "train", fold = heldout, fit$train_pred),
        data.frame(model = key, selection_scope = fit$scope, k = fit$k, split = "heldout", fold = heldout, fit$test_pred)
      )
    }
    fold_df <- do.call(rbind, fold_out)
    pred_df <- do.call(rbind, pred_out)
    held <- pred_df[pred_df$split == "heldout", , drop = FALSE]
    truth <- factor(held$truth, levels = levels(mm$treatment))
    pred <- factor(held$predicted, levels = levels(mm$treatment))
    model_results[[key]] <- list(folds = fold_df, pred = held)
    fold_metrics[[length(fold_metrics) + 1L]] <- fold_df
    predictions[[key]] <- held
    cat(key, "heldout_BA=", sprintf("%.3f", balanced_accuracy(truth, pred)),
        "heldout_F1=", sprintf("%.3f", macro_f1(truth, pred)),
        "mean_gap=", sprintf("%.3f", mean(fold_df$train_test_BA_gap)), "\n")
  }
}

## Named three-feature model: useful as a post-hoc diagnostic only.
key <- "named3_posthoc__3"
fold_out <- list()
pred_out <- list()
for (heldout in levels(mm$block)) {
  train_ids <- mm$sample_uid[mm$block != heldout]
  test_ids <- mm$sample_uid[mm$block == heldout]
  fit <- fit_fold(train_ids, test_ids, heldout, 3L, "named3_posthoc", seed = 9000 + which(levels(mm$block) == heldout))
  fold_out[[length(fold_out) + 1L]] <- data.frame(
    model = key, selection_scope = fit$scope, k = fit$k, heldout_block = heldout,
    n_train = length(train_ids), n_test = length(test_ids),
    train_balanced_accuracy = fit$train_ba, train_macro_f1 = fit$train_f1,
    heldout_balanced_accuracy = fit$test_ba, heldout_macro_f1 = fit$test_f1,
    train_test_BA_gap = fit$train_ba - fit$test_ba,
    lambda = fit$lambda, nonzero_coefficients = fit$nonzero,
    selected_features = paste(fit$selected, collapse = "||")
  )
  pred_out[[length(pred_out) + 1L]] <- data.frame(
    model = key, selection_scope = fit$scope, k = fit$k, split = "heldout", fold = heldout,
    fit$test_pred
  )
}
named_fold <- do.call(rbind, fold_out)
named_pred <- do.call(rbind, pred_out)
fold_metrics[[length(fold_metrics) + 1L]] <- named_fold
predictions[[key]] <- named_pred

fold_metrics_out <- do.call(rbind, fold_metrics)
summary_out <- do.call(rbind, lapply(names(predictions), function(key) {
  d <- predictions[[key]]
  truth <- factor(d$truth, levels = levels(mm$treatment))
  pred <- factor(d$predicted, levels = levels(mm$treatment))
  f <- fold_metrics_out[fold_metrics_out$model == key, , drop = FALSE]
  data.frame(
    model = key, selection_scope = unique(f$selection_scope), k = unique(f$k),
    n = nrow(d), heldout_balanced_accuracy = balanced_accuracy(truth, pred),
    heldout_macro_f1 = macro_f1(truth, pred),
    mean_train_balanced_accuracy = mean(f$train_balanced_accuracy),
    mean_train_macro_f1 = mean(f$train_macro_f1),
    mean_train_test_BA_gap = mean(f$train_test_BA_gap),
    max_train_test_BA_gap = max(f$train_test_BA_gap),
    mean_nonzero_coefficients = mean(f$nonzero_coefficients),
    mean_lambda = mean(f$lambda)
  )
}))

## Paired stratified bootstrap difference versus nested 101-feature model.
bootstrap_diff <- function(a, b, B = 2000, seed = 1) {
  set.seed(seed)
  z <- merge(a[, c("sample_uid", "truth", "predicted")],
             b[, c("sample_uid", "truth", "predicted")],
             by = "sample_uid", suffixes = c("_a", "_b"))
  stopifnot(all(z$truth_a == z$truth_b))
  truth <- factor(z$truth_a, levels = levels(mm$treatment))
  pa <- factor(z$predicted_a, levels = levels(mm$treatment))
  pb <- factor(z$predicted_b, levels = levels(mm$treatment))
  observed <- balanced_accuracy(truth, pa) - balanced_accuracy(truth, pb)
  vals <- numeric(B)
  for (i in seq_len(B)) {
    idx <- unlist(lapply(levels(truth), function(cl) {
      ids <- which(truth == cl)
      sample(ids, length(ids), replace = TRUE)
    }))
    vals[i] <- balanced_accuracy(truth[idx], pa[idx]) - balanced_accuracy(truth[idx], pb[idx])
  }
  data.frame(
    comparison = "model_a_minus_nested_topk_101",
    model_a = NA_character_, model_b = "nested_topk__101",
    observed_delta_BA = observed,
    delta_BA_lo = unname(quantile(vals, 0.025, na.rm = TRUE)),
    delta_BA_hi = unname(quantile(vals, 0.975, na.rm = TRUE)),
    bootstrap_n = B
  )
}

base <- predictions[["nested_topk__101"]]
boot_out <- list()
for (key in setdiff(names(predictions), "nested_topk__101")) {
  if (key == "external_fixed__101") next
  z <- bootstrap_diff(predictions[[key]], base, B = BOOT_N, seed = 12000 + match(key, names(predictions)))
  z$model_a <- key
  boot_out[[length(boot_out) + 1L]] <- z
}
boot_df <- if (length(boot_out)) do.call(rbind, boot_out) else data.frame()

write.table(fold_metrics_out, file.path(outdir, "management_overfit_fold_metrics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary_out, file.path(outdir, "management_overfit_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, predictions), file.path(outdir, "management_overfit_heldout_predictions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(boot_df, file.path(outdir, "management_overfit_paired_bootstrap.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(named3_feature = named3), file.path(outdir, "management_named3_posthoc_features.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nWrote management overfit diagnostics to: ", outdir, "\n", sep = "")
print(summary_out)
