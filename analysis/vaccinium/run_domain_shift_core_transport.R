#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

root <- Sys.getenv("IMETA_EXTENDED_ROOT", "global_harmonized_extended")
outdir <- file.path(root, "domain_shift_core_transport")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

metadata_file <- file.path(root, "global_metadata_primary.tsv")
if (!file.exists(metadata_file))
  stop("Missing metadata: ", metadata_file)

meta <- read.delim(
  metadata_file, check.names=FALSE, stringsAsFactors=FALSE
)

pick_column <- function(x, choices) {
  hit <- choices[choices %in% names(x)]
  if (length(hit)) hit[1] else NA_character_
}

sample_col <- pick_column(meta, c("sample_uid", "sample_id"))
study_col <- pick_column(meta, c("study", "project", "project_accession"))
if (anyNA(c(sample_col, study_col)))
  stop("Metadata must contain sample_uid and study")

meta$sample_uid_model <- as.character(meta[[sample_col]])
meta$study_model <- as.character(meta[[study_col]])

count_files <- c(
  bacteria=file.path(root, "analysis_ready", "bacteria_genus_model_counts.rds"),
  fungi=file.path(root, "analysis_ready", "fungi_genus_model_counts.rds")
)

make_clr <- function(x) {
  z <- log(x + 0.5)
  sweep(z, 1, rowMeans(z), "-")
}

js_divergence <- function(p, q) {
  p <- p / sum(p)
  q <- q / sum(q)
  m <- (p + q) / 2
  kl <- function(a, b) {
    keep <- a > 0 & b > 0
    sum(a[keep] * log(a[keep] / b[keep]))
  }
  (kl(p, m) + kl(q, m)) / 2
}

centroid <- function(x) colMeans(x)

distance_to <- function(x, center) {
  sqrt(rowSums((sweep(x, 2, center, "-"))^2))
}

results <- list()
core_details <- list()

for (marker in names(count_files)) {
  file <- count_files[[marker]]
  if (!file.exists(file)) stop("Missing count table: ", file)
  counts <- readRDS(file)

  shared <- intersect(rownames(counts), meta$sample_uid_model)
  md <- meta[match(shared, meta$sample_uid_model), , drop=FALSE]
  counts <- counts[shared, , drop=FALSE]
  stopifnot(identical(rownames(counts), md$sample_uid_model))

  studies <- sort(unique(md$study_model))
  if (length(studies) < 4)
    stop("Domain-shift audit requires at least four studies")

  for (heldout in studies) {
    test_idx <- which(md$study_model == heldout)
    train_idx <- which(md$study_model != heldout)
    train_studies <- setdiff(studies, heldout)

    keep <- colMeans(counts[train_idx, , drop=FALSE] > 0) >= 0.10 &
      colSums(counts[train_idx, , drop=FALSE]) >= 20
    x <- counts[, keep, drop=FALSE]
    clr <- make_clr(x)

    train_centroids <- do.call(rbind, lapply(
      split(train_idx, md$study_model[train_idx]),
      function(i) centroid(clr[i, , drop=FALSE])
    ))
    heldout_centroid <- centroid(clr[test_idx, , drop=FALSE])
    centroid_distances <- sqrt(rowSums(
      (sweep(train_centroids, 2, heldout_centroid, "-"))^2
    )) / sqrt(ncol(clr))
    nearest_centroid_shift <- min(centroid_distances)

    train_dispersion <- unlist(lapply(
      split(train_idx, md$study_model[train_idx]),
      function(i) {
        c0 <- centroid(clr[i, , drop=FALSE])
        distance_to(clr[i, , drop=FALSE], c0) / sqrt(ncol(clr))
      }
    ))
    heldout_dispersion <- distance_to(
      clr[test_idx, , drop=FALSE], heldout_centroid
    ) / sqrt(ncol(clr))
    shift_ratio <- nearest_centroid_shift / median(train_dispersion)

    train_mean <- colSums(x[train_idx, , drop=FALSE]) + 0.5
    test_mean <- colSums(x[test_idx, , drop=FALSE]) + 0.5
    jsd <- js_divergence(train_mean, test_mean)

    train_prev <- colMeans(x[train_idx, , drop=FALSE] > 0)
    test_prev <- colMeans(x[test_idx, , drop=FALSE] > 0)
    prevalence_drift <- mean(abs(train_prev - test_prev))

    prevalence_by_train_study <- do.call(rbind, lapply(
      train_studies,
      function(s) colMeans(
        x[md$study_model == s, , drop=FALSE] > 0
      )
    ))
    rownames(prevalence_by_train_study) <- train_studies
    core <- colnames(x)[
      apply(prevalence_by_train_study >= 0.50, 2, all)
    ]
    heldout_core_prevalence <- if (length(core))
      colMeans(x[test_idx, core, drop=FALSE] > 0) else numeric()
    validated <- names(heldout_core_prevalence)[
      heldout_core_prevalence >= 0.50
    ]

    results[[length(results) + 1]] <- data.frame(
      marker=marker,
      heldout_study=heldout,
      training_studies=length(train_studies),
      n_train=length(train_idx),
      n_test=length(test_idx),
      features=ncol(x),
      nearest_centroid_shift=nearest_centroid_shift,
      median_training_dispersion=median(train_dispersion),
      median_heldout_dispersion=median(heldout_dispersion),
      centroid_shift_ratio=shift_ratio,
      Jensen_Shannon_divergence=jsd,
      mean_absolute_prevalence_drift=prevalence_drift,
      discovery_core=length(core),
      validated_core=length(validated),
      core_validation_rate=if (length(core))
        length(validated) / length(core) else NA_real_
    )

    if (length(core)) {
      core_details[[length(core_details) + 1]] <- data.frame(
        marker=marker,
        heldout_study=heldout,
        feature=core,
        minimum_training_prevalence=apply(
          prevalence_by_train_study[, core, drop=FALSE], 2, min
        ),
        heldout_prevalence=heldout_core_prevalence[core],
        validated=heldout_core_prevalence[core] >= 0.50
      )
    }
  }
}

summary <- do.call(rbind, results)
details <- do.call(rbind, core_details)

write.table(
  summary, file.path(outdir, "domain_shift_core_transport_summary.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)
write.table(
  details, file.path(outdir, "leave_one_study_out_core_details.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)

for (marker in unique(summary$marker)) {
  z <- summary[summary$marker == marker, , drop=FALSE]
  if (nrow(z) >= 4) {
    test <- cor.test(
      z$centroid_shift_ratio, z$core_validation_rate,
      method="spearman", exact=FALSE
    )
    cat(
      marker,
      "domain-shift/core-retention rho=",
      sprintf("%.3f", unname(test$estimate)),
      "P=", sprintf("%.4f", test$p.value), "\n"
    )
  }
}

cat("Held-out study evaluations:", nrow(summary), "\n")
cat("DOMAIN SHIFT AND CORE TRANSPORT AUDIT DONE\n")
