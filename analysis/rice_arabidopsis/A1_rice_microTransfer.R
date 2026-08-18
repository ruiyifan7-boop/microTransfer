#!/usr/bin/env Rscript
# =============================================================================
# A1 — Official rice four-property reproduction via the microTransfer package
# Produces the submission-grade numbers for §2.9 / Figure 7.
# Usage:  Rscript A1_rice_microTransfer.R <rice_data_dir> <out_dir>
#   e.g.  Rscript A1_rice_microTransfer.R data/rice results/rice_official
# Requires: microTransfer installed; data files present in <rice_data_dir>.
# =============================================================================
suppressPackageStartupMessages(library(microTransfer))
args <- commandArgs(trailingOnly = TRUE)
rice_dir <- if (length(args) >= 1) args[1] else "data/rice"
out      <- if (length(args) >= 2) args[2] else "results/rice_official"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
set.seed(1)

## ---- load via the package adapter (endo+rhizo, depth>=5000, site_year cohorts) ----
rice <- read_edwards_rice(
  otu_file       = file.path(rice_dir, "lc_study_otu_table.tsv.gz"),
  metadata_file  = file.path(rice_dir, "lc_study_mapping_file.tsv"),
  taxonomy_file  = file.path(rice_dir, "gg_otus_tax.rds"),
  organelle_file = file.path(rice_dir, "organelle.rds"),
  min_depth      = 5000,
  compartments   = c("Rhizosphere", "Endosphere"),
  cohort         = "site_year"
)
counts <- rice$counts
md     <- rice$metadata
md$SampleID <- rownames(counts)

## >>> THIS IS THE OFFICIAL SAMPLE COUNT to use throughout the manuscript <<<
cat(sprintf("RICE included samples = %d ; OTUs = %d ; cohorts = %s\n",
            nrow(counts), ncol(counts), paste(sort(unique(md$cohort)), collapse = ", ")))
writeLines(sprintf("rice_included_samples\t%d\nrice_features\t%d",
                   nrow(counts), ncol(counts)),
           file.path(out, "SAMPLE_COUNT.txt"))

## ---- (1) MEMBERSHIP: pooled-vs-study-aware + LOSO with matched null ----
contrast <- pooled_membership_contrast(counts, md, "SampleID", "cohort",
                                       thresholds = c(0.10, 0.30, 0.50))
mem <- validate_membership(counts, md, "SampleID", "cohort",
                           prevalence = 0.50, min_total = 20,
                           n_null = 999, matched = TRUE, seed = 1)
write.table(contrast, file.path(out, "membership_contrast.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(mem$summary, file.path(out, "membership_loso.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("\n[MEMBERSHIP] classify:\n"); print(classify_transferability(mem))

## ---- (2) RESPONSE (plant Age): leave-one-cohort-out ----
studies <- sort(unique(md$cohort)); resp <- list()
for (h in studies) {
  vr <- tryCatch(
    validate_response(counts, md, "SampleID", "cohort", outcome_col = "Age",
                      discovery_studies = setdiff(studies, h),
                      validation_studies = h,
                      covariates = c("Compartment", "Cultivar")),
    error = function(e) { message("response fold ", h, ": ", conditionMessage(e)); NULL })
  if (!is.null(vr)) { s <- vr$summary; s$heldout <- h; resp[[h]] <- s }
}
if (length(resp)) {
  resp <- do.call(rbind, resp)
  write.table(resp, file.path(out, "response_loso.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  print(resp)
}

## ---- (3) PREDICTION (Endosphere vs Rhizosphere): held-out + permutation null ----
pred <- prediction_permutation_test(counts, md, "SampleID", "cohort",
                                    outcome_col = "Compartment",
                                    permutations = 999, seed = 1)
write.table(pred$observed$summary, file.path(out, "prediction_folds.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
writeLines(c(sprintf("observed_mean_ba\t%.4f", pred$observed_balanced_accuracy),
             sprintf("observed_mean_macro_f1\t%.4f", pred$observed_macro_f1),
             sprintf("balanced_accuracy_p\t%.4f", pred$balanced_accuracy_p)),
           file.path(out, "prediction_overall.tsv"))
cat(sprintf("\n[PREDICTION] mean BA = %.3f ; permutation P = %.3f\n",
            pred$observed_balanced_accuracy, pred$balanced_accuracy_p))

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("\nDONE. Official rice results in: ", out, "\n")
