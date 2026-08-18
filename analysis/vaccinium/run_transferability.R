suppressPackageStartupMessages(library(glmnet))

options(stringsAsFactors=FALSE, width=180)
set.seed(20260706)

# Keep the regularization mixing parameter identical for observed and
# permutation-null fits. The manuscript reports this value explicitly.
ENET_ALPHA <- 0.5

root <- "global_harmonized"
outdir <- file.path(root, "enhancement", "transferability")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

pick_counts <- function(kind) {
  candidates <- c(
    file.path(root, "taxa_tables",
              paste0(kind, "_genus_counts_primary211.rds")),
    file.path(root, "analysis_ready",
              paste0(kind, "_genus_broad_counts.rds")),
    file.path(root, "analysis_ready",
              paste0(kind, "_genus_model_counts.rds"))
  )
  hit <- candidates[file.exists(candidates)]
  if (!length(hit))
    stop("No count matrix found for ", kind)
  message(kind, " counts: ", hit[1])
  readRDS(hit[1])
}

meta <- read.delim(
  file.path(root, "global_metadata_primary_211.tsv"),
  check.names=FALSE, stringsAsFactors=FALSE,
  na.strings=c("", "NA")
)
B0 <- pick_counts("bacteria")
F0 <- pick_counts("fungi")

shared <- Reduce(intersect, list(
  meta$sample_uid, rownames(B0), rownames(F0)
))
meta <- meta[match(shared, meta$sample_uid), , drop=FALSE]
B0 <- B0[shared, , drop=FALSE]
F0 <- F0[shared, , drop=FALSE]
stopifnot(
  identical(meta$sample_uid, rownames(B0)),
  identical(rownames(B0), rownames(F0))
)

norm_compartment <- function(x) {
  z <- tolower(trimws(as.character(x)))
  ifelse(
    grepl("rhizo", z), "Rhizosphere",
    ifelse(
      grepl("bulk", z), "Bulk",
      ifelse(grepl("endo|root", z), "Endosphere", NA_character_)
    )
  )
}

meta$compartment_model <- norm_compartment(meta$compartment)
meta$pH_model <- suppressWarnings(as.numeric(meta$pH))
meta$pH_binary <- ifelse(
  is.na(meta$pH_model), NA_character_,
  ifelse(
    meta$pH_model >= 4.0 & meta$pH_model <= 5.5,
    "Favorable", "Outside"
  )
)

coreB <- read.delim(
  file.path(root, "core_microbiome", "bacteria_core_audit.tsv"),
  check.names=FALSE
)
coreF <- read.delim(
  file.path(root, "core_microbiome", "fungi_core_audit.tsv"),
  check.names=FALSE
)
pHB <- read.delim(
  file.path(root, "pH_analysis",
            "bacteria_genus_pH_taxa_validation.tsv"),
  check.names=FALSE
)
pHF <- read.delim(
  file.path(root, "pH_analysis",
            "fungi_genus_pH_taxa_validation.tsv"),
  check.names=FALSE
)

prevalence_cols <- function(x)
  grep("^prevalence_", names(x), value=TRUE)

loso_core <- function(x, heldout) {
  cols <- setdiff(
    prevalence_cols(x), paste0("prevalence_", heldout)
  )
  if (!length(cols))
    stop("No discovery prevalence columns for ", heldout)
  keep <- apply(
    x[, cols, drop=FALSE], 1,
    function(z) all(!is.na(z) & z >= 0.50)
  )
  unique(x$feature[keep])
}

cat("Exact feature mapping:\n")
cat("Bacterial pH candidates:",
    sum(pHB$feature %in% colnames(B0)), "/", nrow(pHB), "\n")
cat("Fungal pH candidates:",
    sum(pHF$feature %in% colnames(F0)), "/", nrow(pHF), "\n")
cat("Bacterial primary core:",
    sum(coreB$strict_core_50 & coreB$feature %in% colnames(B0)),
    "/", sum(coreB$strict_core_50), "\n")
cat("Fungal primary core:",
    sum(coreF$strict_core_50 & coreF$feature %in% colnames(F0)),
    "/", sum(coreF$strict_core_50), "\n")

train_filter <- function(x, ids, min_prev=0.10, min_total=20) {
  z <- x[ids, , drop=FALSE]
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

stratified_foldid <- function(y, k=5, seed=1) {
  set.seed(seed)
  y <- factor(y)
  k <- min(k, min(table(y)))
  if (k < 2) stop("Too few samples per class for inner CV")
  foldid <- integer(length(y))
  for (lev in levels(y)) {
    idx <- which(y == lev)
    idx <- sample(idx)
    foldid[idx] <- rep(seq_len(k), length.out=length(idx))
  }
  foldid
}

balanced_accuracy <- function(truth, pred) {
  truth <- factor(truth)
  pred <- factor(pred, levels=levels(truth))
  mean(vapply(
    levels(truth),
    function(z) mean(pred[truth == z] == z),
    numeric(1)
  ))
}

macro_f1 <- function(truth, pred) {
  truth <- factor(truth)
  pred <- factor(pred, levels=levels(truth))
  mean(vapply(levels(truth), function(z) {
    tp <- sum(truth == z & pred == z)
    fp <- sum(truth != z & pred == z)
    fn <- sum(truth == z & pred != z)
    precision <- if ((tp + fp) == 0) 0 else tp / (tp + fp)
    recall <- if ((tp + fn) == 0) 0 else tp / (tp + fn)
    if ((precision + recall) == 0) 0 else
      2 * precision * recall / (precision + recall)
  }, numeric(1)))
}

binary_auc <- function(truth, score, positive) {
  y <- as.integer(truth == positive)
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  if (!n1 || !n0) return(NA_real_)
  r <- rank(score, ties.method="average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

safe_glmnet <- function(x, y, family, foldid, weights,
                        alpha=ENET_ALPHA) {
  cv.glmnet(
    x=x, y=y, family=family, alpha=alpha,
    foldid=foldid, weights=weights,
    standardize=TRUE, type.measure="class",
    nlambda=80, control=list(maxit=1e6)
  )
}

prepare_variant <- function(kind, variant, train_ids, test_ids,
                            heldout, allow_candidates=FALSE) {
  matrices <- list()

  add_one <- function(prefix, x, core, candidates) {
    if (variant == "full") {
      keep <- train_filter(x, train_ids)
      features <- colnames(x)[keep]
    } else if (variant == "loso_core") {
      features <- intersect(loso_core(core, heldout), colnames(x))
    } else if (variant == "pH_candidates") {
      if (!allow_candidates)
        return(NULL)
      features <- intersect(candidates$feature, colnames(x))
    } else {
      stop("Unknown variant: ", variant)
    }
    if (!length(features)) return(NULL)
    z <- rbind(
      x[train_ids, features, drop=FALSE],
      x[test_ids, features, drop=FALSE]
    )
    z <- clr_transform(z)
    colnames(z) <- paste0(prefix, colnames(z))
    z
  }

  if (kind %in% c("bacteria", "combined")) {
    matrices[["B"]] <- add_one("B__", B0, coreB, pHB)
  }
  if (kind %in% c("fungi", "combined")) {
    matrices[["F"]] <- add_one("F__", F0, coreF, pHF)
  }
  matrices <- matrices[!vapply(matrices, is.null, logical(1))]
  if (!length(matrices)) return(NULL)
  z <- do.call(cbind, matrices)
  list(
    train=z[train_ids, , drop=FALSE],
    test=z[test_ids, , drop=FALSE]
  )
}

nonzero_coefficients <- function(fit, family, lambda,
                                 task, direction, kind, variant,
                                 fold=NA_character_) {
  cc <- coef(fit, s=lambda)
  if (family == "binomial") cc <- list(class=cc)
  out <- do.call(rbind, lapply(names(cc), function(cl) {
    m <- as.matrix(cc[[cl]])
    keep <- rownames(m) != "(Intercept)" & m[, 1] != 0
    if (!any(keep)) return(NULL)
    data.frame(
      task=task, direction=direction, kind=kind,
      variant=variant, fold=fold, class=cl,
      feature=rownames(m)[keep],
      coefficient=m[keep, 1],
      stringsAsFactors=FALSE
    )
  }))
  out
}

fit_external_binary <- function(task, train_study, test_study,
                                train_ids, test_ids, truth,
                                negative, positive,
                                kinds=c("bacteria", "fungi", "combined"),
                                variants=c("full", "loso_core"),
                                candidates=FALSE,
                                nperm=100) {
  direction <- paste0(train_study, "_to_", test_study)
  result <- list()
  predictions <- list()
  coefficients <- list()
  nulls <- list()
  counter <- 0

  if (candidates)
    variants <- unique(c(variants, "pH_candidates"))

  y_all <- factor(truth, levels=c(negative, positive))
  names(y_all) <- c(train_ids, test_ids)
  y_train <- y_all[train_ids]
  y_test <- y_all[test_ids]

  for (kind in kinds) for (variant in variants) {
    allow_candidates <- candidates &&
      train_study == "PRJNA1156347" &&
      test_study == "PRJEB98254"

    dat <- prepare_variant(
      kind, variant, train_ids, test_ids,
      heldout=test_study,
      allow_candidates=allow_candidates
    )
    if (is.null(dat)) next

    foldid <- stratified_foldid(
      y_train, k=5,
      seed=1000 + counter
    )
    weights <- class_weights(y_train)
    cvfit <- safe_glmnet(
      dat$train, y_train, "binomial",
      foldid, weights
    )
    lambda <- cvfit$lambda.1se
    score <- as.numeric(predict(
      cvfit, newx=dat$test,
      s="lambda.1se", type="response"
    ))
    pred <- factor(
      ifelse(score >= 0.5, positive, negative),
      levels=levels(y_test)
    )

    obs_bal <- balanced_accuracy(y_test, pred)
    obs_f1 <- macro_f1(y_test, pred)
    obs_auc <- binary_auc(y_test, score, positive)

    null_bal <- numeric(nperm)
    for (b in seq_len(nperm)) {
      yp <- sample(y_train)
      pf <- glmnet(
        x=dat$train, y=yp, family="binomial",
        alpha=ENET_ALPHA, lambda=lambda,
        weights=class_weights(yp),
        standardize=TRUE, control=list(maxit=1e6)
      )
      ps <- as.numeric(predict(
        pf, newx=dat$test,
        s=lambda, type="response"
      ))
      pp <- factor(
        ifelse(ps >= 0.5, positive, negative),
        levels=levels(y_test)
      )
      null_bal[b] <- balanced_accuracy(y_test, pp)
    }
    perm_p <- (1 + sum(null_bal >= obs_bal)) / (nperm + 1)

    counter <- counter + 1
    result[[counter]] <- data.frame(
      task=task, direction=direction,
      train_study=train_study, test_study=test_study,
      kind=kind, variant=variant,
      n_train=length(train_ids), n_test=length(test_ids),
      features=ncol(dat$train),
      balanced_accuracy=obs_bal,
      macro_f1=obs_f1, AUROC=obs_auc,
      permutation_p=perm_p,
      lambda=lambda
    )
    predictions[[counter]] <- data.frame(
      task=task, direction=direction,
      kind=kind, variant=variant,
      sample_uid=test_ids,
      truth=as.character(y_test),
      predicted=as.character(pred),
      score_positive=score
    )
    coefficients[[counter]] <- nonzero_coefficients(
      cvfit$glmnet.fit, "binomial", lambda,
      task, direction, kind, variant
    )
    nulls[[counter]] <- data.frame(
      task=task, direction=direction,
      kind=kind, variant=variant,
      permutation=seq_len(nperm),
      balanced_accuracy=null_bal
    )
    cat(
      task, direction, kind, variant,
      "features=", ncol(dat$train),
      "BA=", sprintf("%.3f", obs_bal),
      "AUC=", sprintf("%.3f", obs_auc),
      "permP=", sprintf("%.3f", perm_p), "\n"
    )
  }
  list(
    summary=do.call(rbind, result),
    predictions=do.call(rbind, predictions),
    coefficients=do.call(rbind, coefficients),
    null=do.call(rbind, nulls)
  )
}

append_result <- function(acc, x) {
  for (nm in names(acc))
    acc[[nm]] <- rbind(acc[[nm]], x[[nm]])
  acc
}

rbind_fill <- function(...) {
  xs <- list(...)
  xs <- xs[!vapply(xs, is.null, logical(1))]
  if (!length(xs)) return(NULL)
  all_names <- Reduce(union, lapply(xs, names))
  xs <- lapply(xs, function(x) {
    missing <- setdiff(all_names, names(x))
    for (nm in missing) x[[nm]] <- NA
    x[, all_names, drop=FALSE]
  })
  do.call(rbind, xs)
}

all_results <- list(
  summary=NULL, predictions=NULL,
  coefficients=NULL, null=NULL
)

checkpoint <- file.path(outdir, "external_transfer_checkpoint.rds")
if (file.exists(checkpoint)) {
  cat("\n===== RESUME EXTERNAL TRANSFER CHECKPOINT =====\n")
  all_results <- readRDS(checkpoint)
} else {
cat("\n===== EXTERNAL COMPARTMENT TRANSFER =====\n")
task_md <- meta[
  meta$study %in% c("PRJNA1156347", "PRJEB98254") &
    meta$compartment_model %in% c("Bulk", "Rhizosphere"),
  , drop=FALSE
]

for (train_study in c("PRJNA1156347", "PRJEB98254")) {
  test_study <- setdiff(
    c("PRJNA1156347", "PRJEB98254"), train_study
  )
  train_ids <- task_md$sample_uid[task_md$study == train_study]
  test_ids <- task_md$sample_uid[task_md$study == test_study]
  truth <- c(
    task_md$compartment_model[
      match(train_ids, task_md$sample_uid)
    ],
    task_md$compartment_model[
      match(test_ids, task_md$sample_uid)
    ]
  )
  x <- fit_external_binary(
    task="compartment",
    train_study=train_study,
    test_study=test_study,
    train_ids=train_ids, test_ids=test_ids,
    truth=truth,
    negative="Bulk", positive="Rhizosphere",
    candidates=FALSE, nperm=100
  )
  all_results <- append_result(all_results, x)
}

cat("\n===== EXTERNAL pH-CLASS TRANSFER =====\n")
task_md <- meta[
  meta$study %in% c("PRJNA1156347", "PRJEB98254") &
    meta$compartment_model == "Rhizosphere" &
    !is.na(meta$pH_binary),
  , drop=FALSE
]

for (train_study in c("PRJNA1156347", "PRJEB98254")) {
  test_study <- setdiff(
    c("PRJNA1156347", "PRJEB98254"), train_study
  )
  train_ids <- task_md$sample_uid[task_md$study == train_study]
  test_ids <- task_md$sample_uid[task_md$study == test_study]
  truth <- c(
    task_md$pH_binary[
      match(train_ids, task_md$sample_uid)
    ],
    task_md$pH_binary[
      match(test_ids, task_md$sample_uid)
    ]
  )
  x <- fit_external_binary(
    task="pH_4.0_5.5",
    train_study=train_study,
    test_study=test_study,
    train_ids=train_ids, test_ids=test_ids,
    truth=truth,
    negative="Favorable", positive="Outside",
    candidates=TRUE, nperm=100
  )
  all_results <- append_result(all_results, x)
}

saveRDS(
  all_results,
  checkpoint
)
}

cat("\n===== LEAVE-ONE-BLOCK-OUT MANAGEMENT =====\n")
mm <- read.delim(
  "PRJEB110492/metadata/PRJEB110492_master_metadata_provisional.tsv",
  check.names=FALSE, stringsAsFactors=FALSE,
  na.strings=c("", "NA")
)
mm$sample_uid <- paste(
  "PRJEB110492", mm$biological_sample, sep="__"
)
mm <- mm[
  mm$design_family == "transplanted_field_system" &
    mm$treatment_code_raw %in% c("HD", "LD", "wild") &
    mm$sample_uid %in% shared,
  , drop=FALSE
]
mm$treatment <- factor(
  mm$treatment_code_raw, levels=c("HD", "LD", "wild")
)
mm$block <- factor(mm$block_code_raw)
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
      dat <- prepare_variant(
        kind, variant, train_ids, test_ids,
        heldout="PRJEB110492",
        allow_candidates=FALSE
      )
      if (is.null(dat)) next
      y_train <- droplevels(mm$treatment[
        match(train_ids, mm$sample_uid)
      ])
      y_test <- factor(
        mm$treatment[match(test_ids, mm$sample_uid)],
        levels=levels(mm$treatment)
      )
      foldid <- stratified_foldid(
        y_train, k=3,
        seed=5000 + which(levels(mm$block) == held_block)
      )
      cvfit <- safe_glmnet(
        dat$train, y_train, "multinomial",
        foldid, class_weights(y_train)
      )
      lambda <- cvfit$lambda.1se
      pred <- as.vector(predict(
        cvfit, newx=dat$test,
        s="lambda.1se", type="class"
      ))
      fcounter <- fcounter + 1
      lambdas[fcounter] <- lambda
      feature_counts[fcounter] <- ncol(dat$train)
      fold_predictions[[fcounter]] <- data.frame(
        sample_uid=test_ids,
        truth=as.character(y_test),
        predicted=pred,
        fold=held_block
      )
      fold_cache[[fcounter]] <- list(
        train_ids=train_ids,
        test_ids=test_ids,
        train=dat$train,
        test=dat$test,
        lambda=lambda
      )
      fold_coefficients[[fcounter]] <- nonzero_coefficients(
        cvfit$glmnet.fit, "multinomial", lambda,
        "management", "leave_one_block_out",
        kind, variant, fold=held_block
      )
    }
    pp <- do.call(rbind, fold_predictions)
    truth <- factor(pp$truth, levels=levels(mm$treatment))
    pred <- factor(pp$predicted, levels=levels(mm$treatment))
    obs_bal <- balanced_accuracy(truth, pred)
    nperm_management <- as.integer(
      Sys.getenv("MANAGEMENT_NPERM", "999")
    )
    permutation_cores <- as.integer(
      Sys.getenv("MANAGEMENT_PERM_CORES", "4")
    )
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
        y_train_perm <- factor(
          permuted[match(fc$train_ids, mm$sample_uid)],
          levels=levels(mm$treatment)
        )
        fit_null <- glmnet(
          x=fc$train, y=y_train_perm,
          family="multinomial", alpha=ENET_ALPHA,
          lambda=fc$lambda,
          weights=class_weights(y_train_perm),
          standardize=TRUE
        )
        pred_null <- as.vector(predict(
          fit_null, newx=fc$test,
          s=fc$lambda, type="class"
        ))
        null_truth <- c(
          null_truth,
          permuted[
            match(fc$test_ids, mm$sample_uid)
          ]
        )
        null_pred <- c(null_pred, pred_null)
      }
      balanced_accuracy(
        factor(null_truth, levels=levels(mm$treatment)),
        factor(null_pred, levels=levels(mm$treatment))
      )
    }
    if (.Platform$OS.type == "unix" && permutation_cores > 1) {
      null_bal <- unlist(parallel::mclapply(
        seq_len(nperm_management),
        one_management_permutation,
        mc.cores=permutation_cores,
        mc.set.seed=FALSE
      ))
    } else {
      null_bal <- vapply(
        seq_len(nperm_management),
        one_management_permutation,
        numeric(1)
      )
    }
    perm_p <- (
      1 + sum(null_bal >= obs_bal, na.rm=TRUE)
    ) / (nperm_management + 1)
    mcounter <- mcounter + 1
    management_results[[mcounter]] <- data.frame(
      task="management",
      direction="leave_one_block_out",
      train_study="PRJEB110492",
      test_study="heldout_block",
      kind=kind, variant=variant,
      n_train=nrow(mm) - nrow(mm) / nlevels(mm$block),
      n_test=nrow(mm),
      features=round(mean(feature_counts)),
      balanced_accuracy=obs_bal,
      macro_f1=macro_f1(truth, pred),
      AUROC=NA_real_,
      permutation_p=perm_p,
      lambda=mean(lambdas)
    )
    management_nulls[[mcounter]] <- data.frame(
      task="management",
      direction="leave_one_block_out",
      kind=kind,
      variant=variant,
      permutation=seq_len(nperm_management),
      balanced_accuracy=null_bal
    )
    pp$task <- "management"
    pp$direction <- "leave_one_block_out"
    pp$kind <- kind
    pp$variant <- variant
    management_predictions[[mcounter]] <- pp
    management_coefficients[[mcounter]] <- do.call(
      rbind, fold_coefficients
    )
    cat(
      "management", kind, variant,
      "BA=", sprintf(
        "%.3f",
        management_results[[mcounter]]$balanced_accuracy
      ),
      "macroF1=", sprintf(
        "%.3f",
        management_results[[mcounter]]$macro_f1
      ),
      "permP=", sprintf("%.3f", perm_p), "\n"
    )
  }
}

all_results$summary <- rbind(
  all_results$summary, do.call(rbind, management_results)
)
all_results$predictions <- rbind_fill(
  all_results$predictions, do.call(rbind, management_predictions)
)
all_results$coefficients <- rbind_fill(
  all_results$coefficients, do.call(rbind, management_coefficients)
)
all_results$null <- rbind_fill(
  all_results$null, do.call(rbind, management_nulls)
)

write.table(
  all_results$summary,
  file.path(outdir, "transferability_summary.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)
write.table(
  all_results$predictions,
  file.path(outdir, "transferability_predictions.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)
write.table(
  all_results$coefficients,
  file.path(outdir, "transferability_nonzero_coefficients.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)
write.table(
  all_results$null,
  file.path(outdir, "transferability_permutation_null.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)

cat("\n===== SUMMARY =====\n")
print(all_results$summary)
if (!is.null(warnings())) {
  capture.output(
    warnings(),
    file=file.path(outdir, "warnings.txt")
  )
  cat("Warnings written to:", file.path(outdir, "warnings.txt"), "\n")
}
cat("\nTRANSFERABILITY ANALYSIS DONE\n")
