#!/usr/bin/env Rscript

# Export the two minimal, portable inputs required for sample-level FAPROTAX
# and seven-fold leave-one-study-out validation.
#
# Run on the analysis server from the project root:
#   IMETA_EXTENDED_ROOT=global_harmonized_extended \
#   FAPROTAX_EXPORT_DIR=faprotax_server_export \
#   Rscript 06_Data_Code_Archive/scripts/export_faprotax_server_inputs.R

options(stringsAsFactors = FALSE, width = 160)

root <- Sys.getenv("IMETA_EXTENDED_ROOT", "global_harmonized_extended")
outdir <- Sys.getenv("FAPROTAX_EXPORT_DIR", "faprotax_server_export")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

pick_existing <- function(paths, label) {
  hit <- paths[file.exists(paths)]
  if (!length(hit)) {
    stop(
      "Missing ", label, ". Checked:\n  ",
      paste(paths, collapse = "\n  ")
    )
  }
  normalizePath(hit[[1]], winslash = "/", mustWork = TRUE)
}

metadata_file <- pick_existing(
  c(
    file.path(root, "global_metadata_primary.tsv"),
    file.path(root, "global_metadata_primary_265.tsv"),
    file.path(root, "global_metadata_extended_265.tsv")
  ),
  "extended seven-cohort metadata"
)

counts_file <- pick_existing(
  c(
    file.path(root, "analysis_ready", "bacteria_genus_model_counts.rds"),
    file.path(root, "taxa_tables", "bacteria_genus_counts_primary265.rds"),
    file.path(root, "analysis_ready", "bacteria_genus_broad_counts.rds")
  ),
  "extended bacterial genus count matrix"
)

meta <- read.delim(
  metadata_file,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)
counts <- as.matrix(readRDS(counts_file))
storage.mode(counts) <- "numeric"

pick_column <- function(dat, choices, label) {
  hit <- choices[choices %in% names(dat)]
  if (!length(hit)) stop("Metadata is missing ", label)
  hit[[1]]
}

sample_col <- pick_column(
  meta,
  c("sample_uid", "sample_uid_model", "sample_id"),
  "sample identifier"
)
study_col <- pick_column(
  meta,
  c("study", "study_model", "project", "project_accession"),
  "study identifier"
)

meta$sample_uid_export <- as.character(meta[[sample_col]])
meta$study_export <- as.character(meta[[study_col]])

if (anyDuplicated(meta$sample_uid_export)) {
  stop("Metadata contains duplicated sample identifiers")
}
if (anyDuplicated(rownames(counts))) {
  stop("Count matrix contains duplicated sample identifiers")
}
if (any(!is.finite(counts)) || any(counts < 0)) {
  stop("Count matrix contains non-finite or negative values")
}

shared <- intersect(meta$sample_uid_export, rownames(counts))
meta <- meta[match(shared, meta$sample_uid_export), , drop = FALSE]
counts <- counts[shared, , drop = FALSE]
stopifnot(identical(meta$sample_uid_export, rownames(counts)))

studies <- sort(unique(meta$study_export))
if (length(studies) != 7L) {
  stop("Expected seven cohorts, found ", length(studies), ": ",
       paste(studies, collapse = ", "))
}
if (nrow(meta) != 265L) {
  warning("Expected 265 paired samples; export contains ", nrow(meta))
}
if (is.null(colnames(counts)) || any(!nzchar(colnames(counts)))) {
  stop("Count matrix feature names/taxonomic paths are missing")
}

metadata_out <- file.path(outdir, "vaccinium_sample_metadata_7cohort.tsv")
counts_out <- file.path(outdir, "vaccinium_bacteria_genus_counts_7cohort.tsv.gz")

metadata_export <- data.frame(
  sample_uid = meta$sample_uid_export,
  study = meta$study_export,
  stringsAsFactors = FALSE
)
extra_cols <- setdiff(
  intersect(
    c(
      "compartment", "compartment_model", "host", "species", "cultivar",
      "pH", "pH_model", "treatment", "block"
    ),
    names(meta)
  ),
  names(metadata_export)
)
for (v in extra_cols) metadata_export[[v]] <- meta[[v]]

write.table(
  metadata_export,
  metadata_out,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

taxon_by_sample <- data.frame(
  taxonomy = colnames(counts),
  t(counts),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
con <- gzfile(counts_out, open = "wt")
on.exit(close(con), add = TRUE)
write.table(
  taxon_by_sample,
  con,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "0"
)
close(con)
on.exit(NULL, add = FALSE)

audit <- data.frame(
  object = c("metadata", "bacterial_genus_counts"),
  source = c(metadata_file, counts_file),
  output = c(metadata_out, counts_out),
  samples = c(nrow(metadata_export), nrow(metadata_export)),
  studies = c(length(studies), length(studies)),
  features = c(NA_integer_, ncol(counts)),
  stringsAsFactors = FALSE
)
write.table(
  audit,
  file.path(outdir, "export_audit.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

cat("Exported", nrow(metadata_export), "samples from", length(studies), "cohorts\n")
print(table(metadata_export$study))
cat("Bacterial genus features:", ncol(counts), "\n")
cat("Metadata:", metadata_out, "\n")
cat("Counts:", counts_out, "\n")

