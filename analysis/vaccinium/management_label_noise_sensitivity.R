## Label-noise sensitivity for the exploratory PRJEB110492 management classifier.
## The held-out labels remain unchanged. Noise is injected only into the training
## labels within each leave-one-block-out split. Labels are exchanged in pairs
## between different classes so that the training class counts remain unchanged.
## The feature set and inner-CV seed are fixed so that the curve isolates label
## noise rather than class imbalance or feature-selection variability. For each
## perturbed training set, inner folds are re-stratified to the perturbed labels.

suppressPackageStartupMessages(library(glmnet))
options(stringsAsFactors = FALSE, width = 180)

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args) >= 1) args[[1]] else "global_harmonized", mustWork = TRUE)
outdir <- if (length(args) >= 2) args[[2]] else file.path(root, "enhancement", "management_label_noise")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

ENET_ALPHA <- 0.5
CV_FOLDS <- as.integer(Sys.getenv("MANAGEMENT_CV_FOLDS", "9"))
REPS <- as.integer(Sys.getenv("MANAGEMENT_NOISE_REPS", "500"))
if (REPS < 50) stop("MANAGEMENT_NOISE_REPS must be >= 50")

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

coreB <- read.delim(file.path(root, "core_microbiome", "bacteria_core_audit.tsv"), check.names = FALSE)
prevalence_cols <- grep("^prevalence_", names(coreB), value = TRUE)
target_prev <- "prevalence_PRJEB110492"
discovery_cols <- setdiff(prevalence_cols, target_prev)
keep <- apply(coreB[, discovery_cols, drop = FALSE], 1, function(z) all(!is.na(z) & z >= 0.50))
external_core <- unique(coreB$feature[keep])
external_core <- intersect(external_core, colnames(B0))
rank_tab <- coreB[match(external_core, coreB$feature), , drop = FALSE]
rank_tab$minimum_prevalence <- as.numeric(rank_tab$minimum_prevalence)
rank_tab$geometric_mean_RA <- as.numeric(rank_tab$geometric_mean_RA)
rank_tab <- rank_tab[order(-rank_tab$minimum_prevalence, -rank_tab$geometric_mean_RA, rank_tab$feature), , drop = FALSE]
external_rank <- rank_tab$feature
if (length(external_rank) < 101) stop("Expected at least 101 external bacterial core features; found ", length(external_rank))
external_rank <- external_rank[seq_len(101)]

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
mm <- mm[order(mm$sample_uid), , drop = FALSE]
if (nrow(mm) != 36 || nlevels(mm$block) != 4) stop("Expected 36 samples in 4 blocks; found ", nrow(mm), " samples and ", nlevels(mm$block), " blocks")

## Fixed inner-CV seeds for each outer block. The folds themselves are
## re-stratified after label perturbation so that every class is represented.
inner_fold_seed <- setNames(
  7000 + 101 + seq_along(levels(mm$block)),
  levels(mm$block)
)

flip_training_labels <- function(y, fraction, seed) {
  y <- factor(y, levels = levels(mm$treatment))
  target <- round(length(y) * fraction)
  if (target == 0) return(list(labels = y, n_flipped = 0L))
  ## Preserve all class counts while approximating the requested number of
  ## flips. Three-way cycles realize 3 flips; pair swaps realize 2 or 4.
  n_triplets <- target %/% 3L
  remainder <- target - 3L * n_triplets
  n_pairs <- 0L
  if (remainder == 1L) {
    if (n_triplets >= 1L) {
      n_triplets <- n_triplets - 1L
      n_pairs <- 2L
    } else {
      n_pairs <- 1L
    }
  } else if (remainder == 2L) {
    n_pairs <- 1L
  }
  set.seed(seed)
  out <- y
  available <- seq_along(y)
  for (triplet_i in seq_len(n_triplets)) {
    chosen <- vapply(levels(y), function(lev) {
      eligible <- which(y == lev & seq_along(y) %in% available)
      if (!length(eligible)) stop("No eligible sample for class-preserving triplet")
      sample(eligible, 1)
    }, integer(1))
    out[chosen[1]] <- y[chosen[2]]
    out[chosen[2]] <- y[chosen[3]]
    out[chosen[3]] <- y[chosen[1]]
    available <- setdiff(available, chosen)
  }
  for (pair_i in seq_len(n_pairs)) {
    found <- FALSE
    for (attempt in seq_len(1000)) {
      if (length(available) < 2) break
      pair <- sample(available, 2)
      if (as.character(y[pair[1]]) != as.character(y[pair[2]])) {
        out[pair[1]] <- y[pair[2]]
        out[pair[2]] <- y[pair[1]]
        available <- setdiff(available, pair)
        found <- TRUE
        break
      }
    }
    if (!found) break
  }
  out <- factor(out, levels = levels(mm$treatment))
  if (!all(as.integer(table(out)) == as.integer(table(y)))) {
    stop("Class-preserving label swap failed")
  }
  list(labels = out, n_flipped = sum(as.character(out) != as.character(y)))
}

fit_outer <- function(train_ids, test_ids, noisy_y, heldout) {
  y_test <- factor(mm$treatment[match(test_ids, mm$sample_uid)], levels = levels(mm$treatment))
  all_ids <- c(train_ids, test_ids)
  z_all <- clr_transform(B0[all_ids, external_rank, drop = FALSE])
  rownames(z_all) <- all_ids
  x_train <- z_all[train_ids, , drop = FALSE]
  x_test <- z_all[test_ids, , drop = FALSE]
  foldid <- stratified_foldid(
    noisy_y,
    k = CV_FOLDS,
    seed = inner_fold_seed[[as.character(heldout)]]
  )
  cvfit <- cv.glmnet(
    x = x_train, y = noisy_y, family = "multinomial", alpha = ENET_ALPHA,
    foldid = foldid, weights = class_weights(noisy_y), standardize = TRUE,
    type.measure = "class", nlambda = 80, control = list(maxit = 1e6)
  )
  pred <- as.vector(predict(cvfit, newx = x_test, s = "lambda.1se", type = "class"))
  data.frame(sample_uid = test_ids, truth = as.character(y_test), predicted = pred, heldout_block = as.character(heldout))
}

noise_levels <- c(0, 0.05, 0.10, 0.15, 0.20)
rows <- vector("list", length(noise_levels) * REPS)
counter <- 0L
for (rate in noise_levels) {
  for (rep in seq_len(REPS)) {
    pred_rows <- list()
    realized <- numeric(length(levels(mm$block)))
    for (j in seq_along(levels(mm$block))) {
      heldout <- levels(mm$block)[j]
      train_ids <- mm$sample_uid[mm$block != heldout]
      test_ids <- mm$sample_uid[mm$block == heldout]
      y_train <- factor(mm$treatment[match(train_ids, mm$sample_uid)], levels = levels(mm$treatment))
      noisy <- flip_training_labels(y_train, rate, seed = 310000 + rep * 100 + j + round(rate * 1000))
      realized[j] <- noisy$n_flipped / length(y_train)
      pred_rows[[j]] <- fit_outer(train_ids, test_ids, noisy$labels, heldout)
    }
    pred_all <- do.call(rbind, pred_rows)
    truth <- factor(pred_all$truth, levels = levels(mm$treatment))
    pred <- factor(pred_all$predicted, levels = levels(mm$treatment))
    counter <- counter + 1L
    rows[[counter]] <- data.frame(
      target_noise = rate, replicate = rep,
      realized_training_noise = mean(realized),
      heldout_balanced_accuracy = balanced_accuracy(truth, pred),
      heldout_macro_f1 = macro_f1(truth, pred)
    )
  }
  cat("completed target noise", rate, "with", REPS, "replicates\n")
}

replicate_out <- do.call(rbind, rows)
baseline <- mean(replicate_out$heldout_balanced_accuracy[replicate_out$target_noise == 0])
summary_out <- do.call(rbind, lapply(split(replicate_out, replicate_out$target_noise), function(d) {
  data.frame(
    target_noise = unique(d$target_noise),
    mean_realized_training_noise = mean(d$realized_training_noise),
    n_replicates = nrow(d),
    mean_heldout_balanced_accuracy = mean(d$heldout_balanced_accuracy),
    BA_q025 = unname(quantile(d$heldout_balanced_accuracy, 0.025)),
    BA_median = median(d$heldout_balanced_accuracy),
    BA_q975 = unname(quantile(d$heldout_balanced_accuracy, 0.975)),
    mean_heldout_macro_F1 = mean(d$heldout_macro_f1),
    F1_q025 = unname(quantile(d$heldout_macro_f1, 0.025)),
    F1_q975 = unname(quantile(d$heldout_macro_f1, 0.975)),
    fraction_BA_above_three_class_chance = mean(d$heldout_balanced_accuracy > (1 / 3)),
    mean_BA_change_from_zero_noise = mean(d$heldout_balanced_accuracy) - baseline
  )
}))

write.table(replicate_out, file.path(outdir, "management_label_noise_replicates.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary_out, file.path(outdir, "management_label_noise_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nWrote label-noise sensitivity results to: ", outdir, "\n", sep = "")
print(summary_out)
