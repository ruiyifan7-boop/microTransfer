#!/usr/bin/env Rscript

## Re-run the PRJEB110492 leave-one-block-out management classifier.
## The observed models and permutation-null models use the same ENET_ALPHA.
## Usage from the project directory:
##   Rscript rerun_management_alpha05.R global_harmonized
## Optional second argument: output directory.

suppressPackageStartupMessages(library(glmnet))

options(stringsAsFactors = FALSE, width = 180)
set.seed(20260706)

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[[1]] else Sys.getenv("GLOBAL_HARMONIZED_ROOT", "global_harmonized")
root <- normalizePath(root, mustWork = TRUE)
outdir <- if (length(args) >= 2) args[[2]] else Sys.getenv("MANAGEMENT_OUTDIR", file.path(root, "enhancement", "transferability"))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

ENET_ALPHA <- 0.5
N_PERM <- as.integer(Sys.getenv("MANAGEMENT_NPERM", "999"))
PERM_CORES <- as.integer(Sys.getenv("MANAGEMENT_PERM_CORES", "4"))
CV_FOLDS <- as.integer(Sys.getenv("MANAGEMENT_CV_FOLDS", "9"))
if (N_PERM < 1) stop("MANAGEMENT_NPERM must be >= 1")
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
  if (!length(hit)) stop("No count matrix found for ", kind, ":\n", paste(candidates, collapse = "\n"))
  message(kind, " counts: ", hit[[1]])
  as.matrix(readRDS(hit[[1]]))
}

meta_path <- file.path(root, "global_metadata_primary_211.tsv")
if (!file.exists(meta_path)) stop("Missing metadata: ", meta_path)
if (!file.exists(metadata_path)) stop("Missing management metadata: ", metadata_path)

meta <- read.delim(meta_path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA"))
B0 <- pick_counts("bacteria")
F0 <- pick_counts("fungi")

shared <- Reduce(intersect, list(meta$sample_uid, rownames(B0), rownames(F0)))
if (!length(shared)) stop("No shared sample IDs across metadata, bacterial counts, and fungal counts")
meta <- meta[match(shared, meta$sample_uid), , drop = FALSE]
B0 <- B0[shared, , drop = FALSE]
F0 <- F0[shared, , drop = FALSE]
stopifnot(identical(meta$sample_uid, rownames(B0)), identical(rownames(B0), rownames(F0)))

coreB_path <- file.path(root, "core_microbiome", "bacteria_core_audit.tsv")
coreF_path <- file.path(root, "core_microbiome", "fungi_core_audit.tsv")
if (!file.exists(coreB_path) || !file.exists(coreF_path)) stop("Missing core audit table(s) under ", file.path(root, "core_microbiome"))
coreB <- read.delim(coreB_path, check.names = FALSE)
coreF <- read.delim(coreF_path, check.names = FALSE)

prevalence_cols <- function(x) grep("^prevalence_", names(x), value = TRUE)
loso_core <- function(x, heldout) {
  cols <- setdiff(prevalence_cols(x), paste0("prevalence_", heldout))
  if (!length(cols)) stop("No discovery prevalence columns for ", heldout)
  keep <- apply(x[, cols, drop = FALSE], 1, function(z) all(!is.na(z) & z >= 0.50))
  unique(x$feature[keep])
}

train_filter <- function(x, ids, min_prev = 0.10, min_total = 20) {
  z <- x[ids, , drop = FALSE]
  colMeans(z > 0) >= min_prev & colSums(z) >= min_total
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

nonzero_coefficients <- function(fit, lambda, kind, variant, fold) {
  cc <- coef(fit, s = lambda)
  do.call(rbind, lapply(names(cc), function(cl) {
    m <- as.matrix(cc[[cl]])
    keep <- rownames(m) != "(Intercept)" & m[, 1] != 0
    if (!any(keep)) return(NULL)
    data.frame(
      task = "management", direction = "leave_one_block_out",
      kind = kind, variant = variant, fold = fold, class = cl,
      feature = rownames(m)[keep], coefficient = m[keep, 1],
      alpha = ENET_ALPHA, stringsAsFactors = FALSE
    )
  }))
}

prepare_variant <- function(kind, variant, train_ids, test_ids, heldout) {
  add_one <- function(x, core, prefix) {
    features <- if (variant == "full") {
      colnames(x)[train_filter(x, train_ids)]
    } else {
      intersect(loso_core(core, heldout), colnames(x))
    }
    if (!length(features)) return(NULL)
    z <- rbind(x[train_ids, features, drop = FALSE], x[test_ids, features, drop = FALSE])
    z <- clr_transform(z)
    colnames(z) <- paste0(prefix, colnames(z))
    z
  }
  blocks <- list()
  if (kind %in% c("bacteria", "combined")) blocks[["B"]] <- add_one(B0, coreB, "B__")
  if (kind %in% c("fungi", "combined")) blocks[["F"]] <- add_one(F0, coreF, "F__")
  blocks <- blocks[!vapply(blocks, is.null, logical(1))]
  if (!length(blocks)) return(NULL)
  z <- do.call(cbind, blocks)
  list(train = z[train_ids, , drop = FALSE], test = z[test_ids, , drop = FALSE])
}

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
if (nrow(mm) == 0 || nlevels(mm$block) < 2) stop("No usable management samples/blocks after filtering")
print(table(mm$block, mm$treatment))

management_results <- list()
management_predictions <- list()
management_coefficients <- list()
management_nulls <- list()
mcounter <- 0

for (kind in c("bacteria", "fungi", "combined")) {
  for (variant in c("full", "loso_core")) {
    fold_predictions <- list()
    fold_coefficients <- list()
    fold_cache <- list()
    lambdas <- numeric()
    feature_counts <- integer()
    fcounter <- 0

    for (held_block in levels(mm$block)) {
      train_ids <- mm$sample_uid[mm$block != held_block]
      test_ids <- mm$sample_uid[mm$block == held_block]
      dat <- prepare_variant(kind, variant, train_ids, test_ids, heldout = "PRJEB110492")
      if (is.null(dat)) next
      y_train <- droplevels(mm$treatment[match(train_ids, mm$sample_uid)])
      y_test <- factor(mm$treatment[match(test_ids, mm$sample_uid)], levels = levels(mm$treatment))
      foldid <- stratified_foldid(y_train, k = CV_FOLDS, seed = 5000 + which(levels(mm$block) == held_block))

      cvfit <- cv.glmnet(
        x = dat$train, y = y_train, family = "multinomial",
        alpha = ENET_ALPHA, foldid = foldid, weights = class_weights(y_train),
        standardize = TRUE, type.measure = "class", nlambda = 80, control = list(maxit = 1e6)
      )
      lambda <- cvfit$lambda.1se
      pred <- as.vector(predict(cvfit, newx = dat$test, s = "lambda.1se", type = "class"))

      fcounter <- fcounter + 1
      lambdas[fcounter] <- lambda
      feature_counts[fcounter] <- ncol(dat$train)
      fold_predictions[[fcounter]] <- data.frame(
        sample_uid = test_ids, truth = as.character(y_test), predicted = pred,
        fold = held_block, stringsAsFactors = FALSE
      )
      fold_cache[[fcounter]] <- list(train_ids = train_ids, test_ids = test_ids,
                                      train = dat$train, test = dat$test, lambda = lambda)
      fold_coefficients[[fcounter]] <- nonzero_coefficients(
        cvfit$glmnet.fit, lambda, kind, variant, held_block
      )
    }

    pp <- do.call(rbind, fold_predictions)
    truth <- factor(pp$truth, levels = levels(mm$treatment))
    pred <- factor(pp$predicted, levels = levels(mm$treatment))
    obs_bal <- balanced_accuracy(truth, pred)

    one_management_permutation <- function(b) {
      set.seed(20260706 + 100000L * mcounter + b)
      permuted <- as.character(mm$treatment)
      for (block_name in levels(mm$block)) {
        idx <- which(mm$block == block_name)
        permuted[idx] <- sample(permuted[idx])
      }
      null_truth <- character()
      null_pred <- character()
      for (fc in fold_cache) {
        y_train_perm <- factor(permuted[match(fc$train_ids, mm$sample_uid)], levels = levels(mm$treatment))
        fit_null <- glmnet(
          x = fc$train, y = y_train_perm, family = "multinomial",
          alpha = ENET_ALPHA, lambda = fc$lambda,
          weights = class_weights(y_train_perm), standardize = TRUE,
          control = list(maxit = 1e6)
        )
        pred_null <- as.vector(predict(fit_null, newx = fc$test, s = fc$lambda, type = "class"))
        null_truth <- c(null_truth, permuted[match(fc$test_ids, mm$sample_uid)])
        null_pred <- c(null_pred, pred_null)
      }
      balanced_accuracy(
        factor(null_truth, levels = levels(mm$treatment)),
        factor(null_pred, levels = levels(mm$treatment))
      )
    }

    if (.Platform$OS.type == "unix" && PERM_CORES > 1) {
      null_bal <- unlist(parallel::mclapply(seq_len(N_PERM), one_management_permutation,
                                            mc.cores = PERM_CORES, mc.set.seed = FALSE))
    } else {
      null_bal <- vapply(seq_len(N_PERM), one_management_permutation, numeric(1))
    }
    perm_p <- (1 + sum(null_bal >= obs_bal, na.rm = TRUE)) / (N_PERM + 1)
    mcounter <- mcounter + 1

    management_results[[mcounter]] <- data.frame(
      task = "management", direction = "leave_one_block_out",
      train_study = "PRJEB110492", test_study = "heldout_block",
      kind = kind, variant = variant,
      n_train = nrow(mm) - nrow(mm) / nlevels(mm$block), n_test = nrow(mm),
      features = round(mean(feature_counts)), balanced_accuracy = obs_bal,
      macro_f1 = macro_f1(truth, pred), AUROC = NA_real_, permutation_p = perm_p,
      lambda = mean(lambdas), alpha = ENET_ALPHA, n_perm = N_PERM,
      stringsAsFactors = FALSE
    )
    management_nulls[[mcounter]] <- data.frame(
      task = "management", direction = "leave_one_block_out",
      kind = kind, variant = variant, permutation = seq_len(N_PERM),
      balanced_accuracy = null_bal, alpha = ENET_ALPHA,
      stringsAsFactors = FALSE
    )
    pp$task <- "management"; pp$direction <- "leave_one_block_out"
    pp$kind <- kind; pp$variant <- variant; pp$alpha <- ENET_ALPHA
    management_predictions[[mcounter]] <- pp
    management_coefficients[[mcounter]] <- do.call(rbind, fold_coefficients)
    cat("management", kind, variant, "BA=", sprintf("%.3f", obs_bal),
        "macroF1=", sprintf("%.3f", macro_f1(truth, pred)),
        "permP=", sprintf("%.4f", perm_p), "alpha=", ENET_ALPHA, "\n")
  }
}

fill_bind <- function(a, b) {
  if (is.null(a) || !nrow(a)) return(b)
  if (is.null(b) || !nrow(b)) return(a)
  cols <- union(names(a), names(b))
  for (nm in setdiff(cols, names(a))) a[[nm]] <- NA
  for (nm in setdiff(cols, names(b))) b[[nm]] <- NA
  rbind(a[, cols, drop = FALSE], b[, cols, drop = FALSE])
}

merge_replace_management <- function(path, new_rows) {
  old <- if (file.exists(path)) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE) else NULL
  if (!is.null(old) && "task" %in% names(old)) old <- old[old$task != "management", , drop = FALSE]
  fill_bind(old, new_rows)
}

summary_new <- do.call(rbind, management_results)
pred_new <- do.call(rbind, management_predictions)
coef_new <- do.call(rbind, management_coefficients)
null_new <- do.call(rbind, management_nulls)

summary_out <- merge_replace_management(file.path(outdir, "transferability_summary.tsv"), summary_new)
pred_out <- merge_replace_management(file.path(outdir, "transferability_predictions.tsv"), pred_new)
coef_out <- merge_replace_management(file.path(outdir, "transferability_nonzero_coefficients.tsv"), coef_new)
null_out <- merge_replace_management(file.path(outdir, "transferability_permutation_null.tsv"), null_new)

write.table(summary_out, file.path(outdir, "transferability_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pred_out, file.path(outdir, "transferability_predictions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(coef_out, file.path(outdir, "transferability_nonzero_coefficients.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(null_out, file.path(outdir, "transferability_permutation_null.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary_new, file.path(outdir, "management_alpha05_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(null_new, file.path(outdir, "management_alpha05_permutation_null.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nWrote alpha=", ENET_ALPHA, " management results to: ", outdir, "\n", sep = "")
cat("Existing non-management transferability rows were preserved.\n")
print(summary_new)
