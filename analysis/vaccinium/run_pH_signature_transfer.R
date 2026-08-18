options(stringsAsFactors=FALSE, width=180)

root <- "global_harmonized"
outdir <- file.path(root, "enhancement", "pH_signature_transfer")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

meta <- read.delim(
  file.path(root, "global_metadata_primary_211.tsv"),
  check.names=FALSE, na.strings=c("", "NA")
)
B <- readRDS(file.path(
  root, "analysis_ready", "bacteria_genus_model_counts.rds"
))
F <- readRDS(file.path(
  root, "analysis_ready", "fungi_genus_model_counts.rds"
))

candidate_files <- c(
  bacteria=file.path(
    root, "pH_analysis", "bacteria_genus_pH_taxa_validation.tsv"
  ),
  fungi=file.path(
    root, "pH_analysis", "fungi_genus_pH_taxa_validation.tsv"
  )
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

shared <- Reduce(intersect, list(
  meta$sample_uid, rownames(B), rownames(F)
))
meta <- meta[match(shared, meta$sample_uid), , drop=FALSE]
B <- B[shared, , drop=FALSE]
F <- F[shared, , drop=FALSE]

discovery <- meta[
  meta$study == "PRJNA1156347" &
    meta$compartment_model == "Rhizosphere" &
    !is.na(meta$pH_model),
  , drop=FALSE
]
validation <- meta[
  meta$study == "PRJEB98254" &
    meta$compartment_model == "Rhizosphere" &
    !is.na(meta$pH_model),
  , drop=FALSE
]

if (nrow(discovery) != 15)
  warning("Expected 15 discovery rhizosphere samples; found ", nrow(discovery))
if (nrow(validation) != 29)
  warning("Expected 29 field rhizosphere samples; found ", nrow(validation))

clr_transform <- function(x) {
  z <- log(as.matrix(x) + 0.5)
  sweep(z, 1, rowMeans(z), "-")
}

find_feature_column <- function(x) {
  hit <- c("feature", "Feature", "taxon", "genus")
  hit <- hit[hit %in% names(x)]
  if (!length(hit))
    stop("No feature column in candidate table")
  hit[1]
}

score_one_marker <- function(marker, counts, candidate_path) {
  candidates <- read.delim(
    candidate_path, check.names=FALSE, stringsAsFactors=FALSE
  )
  feature_col <- find_feature_column(candidates)
  requested <- unique(as.character(candidates[[feature_col]]))
  features <- intersect(requested, colnames(counts))
  if (!length(features))
    stop("No ", marker, " pH candidates mapped to count table")

  clr <- clr_transform(counts)
  d <- clr[discovery$sample_uid, features, drop=FALSE]
  v <- clr[validation$sample_uid, features, drop=FALSE]

  means <- colMeans(d)
  sds <- apply(d, 2, sd)
  usable <- is.finite(sds) & sds > 0
  d <- d[, usable, drop=FALSE]
  v <- v[, usable, drop=FALSE]
  means <- means[usable]
  sds <- sds[usable]

  rho <- vapply(colnames(d), function(feature) {
    suppressWarnings(cor(
      d[, feature], discovery$pH_model,
      method="spearman", use="complete.obs"
    ))
  }, numeric(1))
  direction <- sign(rho)
  usable <- is.finite(direction) & direction != 0
  d <- d[, usable, drop=FALSE]
  v <- v[, usable, drop=FALSE]
  means <- means[usable]
  sds <- sds[usable]
  rho <- rho[usable]
  direction <- direction[usable]

  d_scaled <- sweep(d, 2, means, "-")
  d_scaled <- sweep(d_scaled, 2, sds, "/")
  v_scaled <- sweep(v, 2, means, "-")
  v_scaled <- sweep(v_scaled, 2, sds, "/")

  discovery_score <- rowMeans(
    sweep(d_scaled, 2, direction, "*")
  )
  validation_score <- rowMeans(
    sweep(v_scaled, 2, direction, "*")
  )

  field <- data.frame(
    sample_uid=validation$sample_uid,
    marker=marker,
    pH=validation$pH_model,
    signature_score=as.numeric(validation_score)
  )
  field$pH_z <- as.numeric(scale(field$pH))
  fit <- lm(signature_score ~ pH_z + I(pH_z^2), data=field)
  reduced <- lm(signature_score ~ 1, data=field)
  observed_r2 <- summary(fit)$r.squared
  observed_rho <- suppressWarnings(cor(
    field$signature_score, field$pH, method="spearman"
  ))

  set.seed(if (marker == "bacteria") 202607061 else 202607062)
  nperm <- 999
  null_r2 <- numeric(nperm)
  null_rho <- numeric(nperm)
  for (i in seq_len(nperm)) {
    ph <- sample(field$pH)
    phz <- as.numeric(scale(ph))
    null_fit <- lm(field$signature_score ~ phz + I(phz^2))
    null_r2[i] <- summary(null_fit)$r.squared
    null_rho[i] <- suppressWarnings(cor(
      field$signature_score, ph, method="spearman"
    ))
  }
  p_r2 <- (1 + sum(null_r2 >= observed_r2)) / (nperm + 1)
  p_rho <- (1 + sum(abs(null_rho) >= abs(observed_rho))) /
    (nperm + 1)
  model_cmp <- anova(reduced, fit)

  weights <- data.frame(
    marker=marker,
    feature=colnames(d),
    discovery_spearman_rho=as.numeric(rho),
    direction=as.numeric(direction),
    discovery_mean=as.numeric(means),
    discovery_sd=as.numeric(sds)
  )
  scores <- rbind(
    data.frame(
      sample_uid=discovery$sample_uid,
      study="PRJNA1156347",
      marker=marker,
      pH=discovery$pH_model,
      signature_score=as.numeric(discovery_score)
    ),
    data.frame(
      sample_uid=validation$sample_uid,
      study="PRJEB98254",
      marker=marker,
      pH=validation$pH_model,
      signature_score=as.numeric(validation_score)
    )
  )
  summary_row <- data.frame(
    marker=marker,
    discovery_n=nrow(discovery),
    field_n=nrow(validation),
    requested_candidates=length(requested),
    mapped_candidates=length(features),
    usable_candidates=ncol(d),
    field_spearman_rho=observed_rho,
    field_spearman_permutation_p=p_rho,
    field_quadratic_R2=observed_r2,
    field_quadratic_F=model_cmp$F[2],
    field_quadratic_parametric_p=model_cmp$`Pr(>F)`[2],
    field_quadratic_permutation_p=p_r2
  )
  null <- data.frame(
    marker=marker,
    permutation=seq_len(nperm),
    spearman_rho=null_rho,
    quadratic_R2=null_r2
  )
  list(
    summary=summary_row,
    scores=scores,
    weights=weights,
    null=null
  )
}

results <- list(
  score_one_marker("bacteria", B, candidate_files[["bacteria"]]),
  score_one_marker("fungi", F, candidate_files[["fungi"]])
)

summary_tab <- do.call(rbind, lapply(results, `[[`, "summary"))
scores <- do.call(rbind, lapply(results, `[[`, "scores"))
weights <- do.call(rbind, lapply(results, `[[`, "weights"))
null <- do.call(rbind, lapply(results, `[[`, "null"))

write.table(
  summary_tab, file.path(outdir, "pH_signature_transfer_summary.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)
write.table(
  scores, file.path(outdir, "pH_signature_scores.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)
write.table(
  weights, file.path(outdir, "pH_signature_weights.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)
write.table(
  null, file.path(outdir, "pH_signature_permutation_null.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)

cat("Discovery samples:", nrow(discovery), "\n")
cat("Field samples:", nrow(validation), "\n")
print(summary_tab, row.names=FALSE)
cat("PH SIGNATURE TRANSFER ANALYSIS DONE\n")
