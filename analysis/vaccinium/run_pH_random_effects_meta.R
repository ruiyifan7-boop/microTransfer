#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

root <- Sys.getenv("IMETA_EXTENDED_ROOT", "global_harmonized_extended")
outdir <- file.path(root, "pH_random_effects_meta")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

if (!requireNamespace("metafor", quietly=TRUE))
  stop("Package 'metafor' is required. Install with install.packages('metafor').")
if (!requireNamespace("vegan", quietly=TRUE))
  stop("Package 'vegan' is required.")

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
pH_col <- pick_column(meta, c("pH", "pH_value", "pH_model", "ph"))
compartment_col <- pick_column(
  meta, c("compartment", "compartment_model", "sample_type")
)

needed <- c(sample_col, study_col, pH_col)
if (anyNA(needed))
  stop("Metadata must contain sample_uid, study, and pH columns")

meta$sample_uid_model <- as.character(meta[[sample_col]])
meta$study_model <- as.character(meta[[study_col]])
meta$pH_model <- suppressWarnings(as.numeric(meta[[pH_col]]))
meta$compartment_model <- if (!is.na(compartment_col))
  as.character(meta[[compartment_col]]) else "Rhizosphere"

meta <- meta[
  is.finite(meta$pH_model) &
    tolower(meta$compartment_model) == "rhizosphere",
  , drop=FALSE
]

min_n <- as.integer(Sys.getenv("PH_META_MIN_N", "8"))
min_unique_pH <- as.integer(Sys.getenv("PH_META_MIN_UNIQUE_PH", "3"))
min_pH_range <- as.numeric(Sys.getenv("PH_META_MIN_RANGE", "0.3"))
adjustment_requested <- trimws(strsplit(
  Sys.getenv("PH_META_ADJUST", ""), ",", fixed=TRUE
)[[1]])
adjustment_requested <- adjustment_requested[nzchar(adjustment_requested)]

study_ok <- vapply(
  split(meta, meta$study_model),
  function(z)
    nrow(z) >= min_n &&
    length(unique(z$pH_model)) >= min_unique_pH &&
    diff(range(z$pH_model)) >= min_pH_range,
  logical(1)
)
eligible_studies <- names(study_ok)[study_ok]
meta <- meta[meta$study_model %in% eligible_studies, , drop=FALSE]

if (length(eligible_studies) < 3)
  stop(
    "Random-effects pH meta-analysis needs at least 3 eligible studies; found ",
    length(eligible_studies), ": ", paste(eligible_studies, collapse=", ")
  )

study_audit <- do.call(rbind, lapply(
  split(meta, meta$study_model),
  function(z) data.frame(
    study=z$study_model[1],
    n=nrow(z),
    unique_pH=length(unique(z$pH_model)),
    minimum_pH=min(z$pH_model),
    maximum_pH=max(z$pH_model),
    pH_range=diff(range(z$pH_model))
  )
))
write.table(
  study_audit, file.path(outdir, "eligible_pH_studies.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)

count_files <- c(
  bacteria=file.path(root, "analysis_ready", "bacteria_genus_model_counts.rds"),
  fungi=file.path(root, "analysis_ready", "fungi_genus_model_counts.rds")
)

short_taxon <- function(x) sub(".*;g__", "", x)

make_clr <- function(x) {
  z <- log(x + 0.5)
  sweep(z, 1, rowMeans(z), "-")
}

shannon <- function(x) {
  p <- x / rowSums(x)
  p[p <= 0] <- NA_real_
  -rowSums(p * log(p), na.rm=TRUE)
}

study_formula <- function(z, response) {
  terms <- c("pH_z", "I(pH_z^2)")
  available <- intersect(adjustment_requested, names(z))
  usable <- available[vapply(
    available,
    function(v) length(unique(z[[v]][!is.na(z[[v]])])) > 1,
    logical(1)
  )]
  reformulate(c(terms, usable), response=response)
}

coefficient_rows <- function(y, md, feature, marker) {
  z <- md
  z$response <- as.numeric(y)
  z$pH_z <- as.numeric(scale(z$pH_model))
  fit <- try(lm(study_formula(z, "response"), data=z), silent=TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  co <- summary(fit)$coefficients
  wanted <- intersect(c("pH_z", "I(pH_z^2)"), rownames(co))
  if (!length(wanted)) return(NULL)
  do.call(rbind, lapply(wanted, function(term) data.frame(
    marker=marker,
    feature=feature,
    genus=short_taxon(feature),
    study=z$study_model[1],
    n=nrow(z),
    pH_min=min(z$pH_model),
    pH_max=max(z$pH_model),
    term=ifelse(term == "pH_z", "linear", "quadratic"),
    estimate=co[term, "Estimate"],
    standard_error=co[term, "Std. Error"],
    statistic=co[term, "t value"],
    p=co[term, "Pr(>|t|)"],
    model_adjustment=paste(
      setdiff(attr(terms(fit), "term.labels"), c("pH_z", "I(pH_z^2)")),
      collapse=";"
    )
  )))
}

meta_one <- function(z) {
  z <- z[
    is.finite(z$estimate) & is.finite(z$standard_error) &
      z$standard_error > 0,
    , drop=FALSE
  ]
  if (nrow(z) < 3) return(NULL)
  fit <- try(
    metafor::rma.uni(
      yi=z$estimate, sei=z$standard_error,
      method="REML", test="knha"
    ),
    silent=TRUE
  )
  if (inherits(fit, "try-error")) return(NULL)
  pred <- predict(fit)
  data.frame(
    marker=z$marker[1],
    feature=z$feature[1],
    genus=z$genus[1],
    term=z$term[1],
    k=nrow(z),
    estimate=as.numeric(fit$b),
    standard_error=fit$se,
    ci_lower=fit$ci.lb,
    ci_upper=fit$ci.ub,
    p=fit$pval,
    tau2=fit$tau2,
    I2=fit$I2,
    Q=fit$QE,
    Q_p=fit$QEp,
    prediction_lower=pred$pi.lb,
    prediction_upper=pred$pi.ub
  )
}

all_study_effects <- list()
all_meta_effects <- list()
all_community <- list()
all_alpha <- list()

for (marker in names(count_files)) {
  file <- count_files[[marker]]
  if (!file.exists(file)) stop("Missing count table: ", file)
  counts <- readRDS(file)

  shared <- intersect(rownames(counts), meta$sample_uid_model)
  md <- meta[match(shared, meta$sample_uid_model), , drop=FALSE]
  counts <- counts[shared, , drop=FALSE]
  stopifnot(identical(rownames(counts), md$sample_uid_model))

  prevalence_by_study <- do.call(rbind, lapply(
    split(seq_len(nrow(counts)), md$study_model),
    function(i) colMeans(counts[i, , drop=FALSE] > 0)
  ))
  keep <- colMeans(counts > 0) >= 0.10 &
    colSums(prevalence_by_study >= 0.10) >= 2 &
    colSums(counts) >= 20
  counts <- counts[, keep, drop=FALSE]
  clr <- make_clr(counts)

  marker_effects <- list()
  for (study in eligible_studies) {
    idx <- which(md$study_model == study)
    for (j in seq_len(ncol(clr))) {
      marker_effects[[length(marker_effects) + 1]] <- coefficient_rows(
        clr[idx, j], md[idx, , drop=FALSE], colnames(clr)[j], marker
      )
    }

    z <- md[idx, , drop=FALSE]
    z$pH_z <- as.numeric(scale(z$pH_model))
    d <- dist(clr[idx, , drop=FALSE])
    set.seed(20260706)
    ad <- vegan::adonis2(
      d ~ pH_z + I(pH_z^2), data=z,
      permutations=999, by="margin"
    )
    for (term in intersect(c("pH_z", "I(pH_z^2)"), rownames(ad))) {
      all_community[[length(all_community) + 1]] <- data.frame(
        marker=marker, study=study, n=length(idx),
        term=ifelse(term == "pH_z", "linear", "quadratic"),
        R2=ad[term, "R2"], F=ad[term, "F"], p=ad[term, "Pr(>F)"]
      )
    }

    alpha_rows <- coefficient_rows(
      shannon(counts[idx, , drop=FALSE]),
      md[idx, , drop=FALSE], "Shannon", marker
    )
    all_alpha[[length(all_alpha) + 1]] <- alpha_rows
  }

  marker_effects <- do.call(rbind, marker_effects)
  all_study_effects[[marker]] <- marker_effects
  split_key <- interaction(
    marker_effects$feature, marker_effects$term, drop=TRUE
  )
  marker_meta <- do.call(rbind, lapply(
    split(marker_effects, split_key), meta_one
  ))
  if (!is.null(marker_meta)) {
    marker_meta$q <- p.adjust(marker_meta$p, method="BH")
    all_meta_effects[[marker]] <- marker_meta
  }
}

study_effects <- do.call(rbind, all_study_effects)
meta_effects <- do.call(rbind, all_meta_effects)
community <- do.call(rbind, all_community)
alpha_study <- do.call(rbind, all_alpha)
alpha_meta <- do.call(rbind, lapply(
  split(
    alpha_study,
    interaction(alpha_study$marker, alpha_study$term, drop=TRUE)
  ),
  meta_one
))
if (!is.null(alpha_meta))
  alpha_meta$q <- p.adjust(alpha_meta$p, method="BH")

if (is.null(meta_effects) || !nrow(meta_effects))
  stop("No taxon-level random-effects models could be fitted")

write.table(
  study_effects, file.path(outdir, "taxon_study_effects.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)
write.table(
  meta_effects, file.path(outdir, "taxon_random_effects_meta.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)
write.table(
  community, file.path(outdir, "community_pH_PERMANOVA_by_study.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)
write.table(
  alpha_study, file.path(outdir, "alpha_diversity_study_effects.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)
write.table(
  alpha_meta, file.path(outdir, "alpha_diversity_random_effects_meta.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)

cat("Eligible pH studies:", paste(eligible_studies, collapse=", "), "\n")
cat("Study-level taxon coefficients:", nrow(study_effects), "\n")
cat("Random-effects taxon models:", nrow(meta_effects), "\n")
cat("FDR-significant meta-analytic effects:", sum(meta_effects$q < 0.05), "\n")
cat("PH RANDOM-EFFECTS META-ANALYSIS DONE\n")
