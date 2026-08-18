#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop(
    "Usage: Rscript run_edwards_rice_generalization.R ",
    "lc_study_otu_table.tsv.gz lc_study_mapping_file.tsv output_dir ",
    "[gg_otus_tax.rds] [organelle.rds]"
  )
}

otu_file <- args[[1]]
metadata_file <- args[[2]]
outdir <- args[[3]]
taxonomy_file <- if (length(args) >= 4L) args[[4]] else NULL
organelle_file <- if (length(args) >= 5L) args[[5]] else NULL
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(microTransfer))

rice <- read_edwards_rice(
  otu_file = otu_file,
  metadata_file = metadata_file,
  taxonomy_file = taxonomy_file,
  organelle_file = organelle_file,
  min_depth = 5000,
  compartments = c("Rhizosphere", "Endosphere"),
  sites = c("Arbuckle", "Jonesboro"),
  cohort = "site_year"
)

cohort_size <- as.data.frame(table(rice$metadata$cohort))
names(cohort_size) <- c("cohort", "n")
utils::write.table(
  cohort_size,
  file.path(outdir, "rice_cohort_sizes.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Require enough held-out observations for interpretable fold-level estimates.
eligible <- cohort_size$cohort[cohort_size$n >= 20]
keep <- rice$metadata$cohort %in% eligible
counts <- rice$counts[keep, , drop = FALSE]
metadata <- rice$metadata[keep, , drop = FALSE]

contrast <- pooled_membership_contrast(
  counts, metadata, "SampleID", "cohort",
  thresholds = c(0.10, 0.30, 0.50)
)
utils::write.table(
  contrast,
  file.path(outdir, "rice_pooled_vs_study_aware_membership.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

membership <- validate_membership(
  counts, metadata, "SampleID", "cohort",
  prevalence = 0.50,
  min_total = 20,
  n_null = 999,
  matched = TRUE,
  seed = 20260722
)
utils::write.table(
  membership$summary,
  file.path(outdir, "rice_membership_loso_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
utils::write.table(
  membership$feature_results,
  file.path(outdir, "rice_membership_loso_features.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
utils::write.table(
  data.frame(feature = membership$stable_features),
  file.path(outdir, "rice_stable_members.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
utils::write.table(
  classify_transferability(membership),
  file.path(outdir, "rice_membership_evidence_label.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Prediction tests whether compartment labels transfer across site-year cohorts.
if ("Compartment" %in% names(metadata) &&
    length(unique(metadata$Compartment)) >= 2L &&
    length(membership$stable_features) >= 2L) {
  prediction <- validate_prediction(
    counts, metadata, "SampleID", "cohort", "Compartment",
    features = membership$stable_features,
    transform = "clr"
  )
  utils::write.table(
    prediction$summary,
    file.path(outdir, "rice_compartment_prediction_loso.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  utils::write.table(
    prediction$predictions,
    file.path(outdir, "rice_compartment_predictions.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

# Age response: use the first three eligible cohorts for discovery and reserve
# the final eligible cohort for independent projection.
if ("Age" %in% names(metadata)) {
  metadata$Age <- suppressWarnings(as.numeric(metadata$Age))
  age_ok <- tapply(
    is.finite(metadata$Age), metadata$cohort,
    function(x) sum(x) >= 20L
  )
  age_cohorts <- sort(names(age_ok)[age_ok])
  if (length(age_cohorts) >= 3L) {
    validation_study <- tail(age_cohorts, 1L)
    discovery_studies <- setdiff(age_cohorts, validation_study)
    response <- validate_response(
      counts, metadata, "SampleID", "cohort", "Age",
      discovery_studies = discovery_studies,
      validation_studies = validation_study,
      covariates = "Compartment",
      min_prevalence = 0.05,
      min_total = 20,
      alpha = 0.05
    )
    utils::write.table(
      response$discovery$pooled_effects,
      file.path(outdir, "rice_age_discovery_meta.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    utils::write.table(
      response$validation_effects,
      file.path(outdir, "rice_age_validation_effects.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    utils::write.table(
      response$summary,
      file.path(outdir, "rice_age_validation_summary.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
  }
}

writeLines(
  c(
    paste("Dataset:", rice$source$dataset),
    paste("Dryad DOI:", rice$source$dryad_doi),
    paste("Eligible cohorts:", paste(eligible, collapse = ", ")),
    paste("Samples analyzed:", nrow(counts)),
    paste("Features analyzed:", ncol(counts))
  ),
  file.path(outdir, "rice_generalization_run_info.txt")
)

cat("Rice generalization outputs written to:", normalizePath(outdir), "\n")
