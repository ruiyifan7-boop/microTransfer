#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 180)
suppressPackageStartupMessages(library(vegan))

root <- normalizePath(Sys.getenv("IMETA_ROOT", "global_harmonized"), mustWork = TRUE)
outdir <- file.path(root, "enhancement", "submission_gap_exports")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

pick_file <- function(candidates) {
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) stop("Missing required file:\n", paste(candidates, collapse = "\n"))
  hit[[1]]
}

pick_col <- function(x, choices, label) {
  hit <- choices[choices %in% names(x)]
  if (!length(hit)) stop("Could not identify ", label, ". Available columns:\n", paste(names(x), collapse = ", "))
  hit[[1]]
}

clr <- function(x) {
  z <- log(as.matrix(x) + 0.5)
  sweep(z, 1, rowMeans(z), "-")
}

norm_compartment <- function(x) {
  z <- tolower(trimws(as.character(x)))
  ifelse(grepl("rhizo", z), "Rhizosphere",
         ifelse(grepl("bulk|soil", z), "Bulk",
                ifelse(grepl("endo|root", z), "Endosphere", NA_character_)))
}

meta_file <- pick_file(c(
  file.path(root, "global_metadata_primary_211.tsv"),
  file.path(root, "global_metadata_primary.tsv")
))
B_file <- pick_file(c(
  file.path(root, "taxa_tables", "bacteria_genus_counts_primary211.rds"),
  file.path(root, "analysis_ready", "bacteria_genus_model_counts.rds"),
  file.path(root, "analysis_ready", "bacteria_genus_broad_counts.rds")
))
F_file <- pick_file(c(
  file.path(root, "taxa_tables", "fungi_genus_counts_primary211.rds"),
  file.path(root, "analysis_ready", "fungi_genus_model_counts.rds"),
  file.path(root, "analysis_ready", "fungi_genus_broad_counts.rds")
))

meta <- read.delim(meta_file, check.names = FALSE, na.strings = c("", "NA"))
B0 <- as.matrix(readRDS(B_file))
F0 <- as.matrix(readRDS(F_file))

sample_col <- pick_col(meta, c("sample_uid", "sample_id"), "sample ID")
study_col <- pick_col(meta, c("study", "project", "project_accession"), "study")
comp_col <- pick_col(meta, c("compartment", "compartment_normalized", "sample_type"), "compartment")
treat_col <- pick_col(meta, c("pH_treatment", "pH_treatment_raw", "treatment", "group", "pH_group"), "pH treatment")

shared <- Reduce(intersect, list(as.character(meta[[sample_col]]), rownames(B0), rownames(F0)))
meta <- meta[match(shared, meta[[sample_col]]), , drop = FALSE]
B0 <- B0[shared, , drop = FALSE]
F0 <- F0[shared, , drop = FALSE]
meta$study_model <- as.character(meta[[study_col]])
meta$compartment_model <- norm_compartment(meta[[comp_col]])
meta$pH_treatment_model <- factor(as.character(meta[[treat_col]]))

controlled_rows <- list()
for (kingdom in c("Bacteria", "Fungi")) {
  counts <- if (kingdom == "Bacteria") B0 else F0
  for (comp in c("Bulk", "Rhizosphere", "Endosphere")) {
    z <- meta[meta$study_model == "PRJNA1156347" & meta$compartment_model == comp, , drop = FALSE]
    z <- z[!is.na(z$pH_treatment_model), , drop = FALSE]
    idx <- match(z[[sample_col]], rownames(counts))
    ok <- !is.na(idx)
    z <- z[ok, , drop = FALSE]
    idx <- idx[ok]
    if (nrow(z) < 4 || nlevels(droplevels(z$pH_treatment_model)) < 2) next
    d <- dist(clr(counts[idx, , drop = FALSE]))
    fit <- vegan::adonis2(d ~ pH_treatment_model, data = z, permutations = 999, by = "margin")
    controlled_rows[[length(controlled_rows) + 1L]] <- data.frame(
      kingdom = kingdom, compartment = comp, n = nrow(z),
      R2 = unname(fit["pH_treatment_model", "R2"]),
      P = unname(fit["pH_treatment_model", "Pr(>F)"]), stringsAsFactors = FALSE
    )
  }
}
controlled <- do.call(rbind, controlled_rows)
write.table(controlled, file.path(outdir, "controlled_pH_exact_R2.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

management_file <- file.path(dirname(root), "PRJEB110492", "metadata", "PRJEB110492_master_metadata_provisional.tsv")
if (!file.exists(management_file)) stop("Missing management metadata: ", management_file)
mm <- read.delim(management_file, check.names = FALSE, na.strings = c("", "NA"))
mm$sample_uid <- paste("PRJEB110492", mm$biological_sample, sep = "__")
mm <- mm[mm$design_family == "transplanted_field_system" &
           mm$treatment_code_raw %in% c("HD", "LD", "wild") & mm$sample_uid %in% shared, , drop = FALSE]
mm$treatment <- factor(mm$treatment_code_raw, levels = c("HD", "LD", "wild"))
mm$block <- factor(mm$block_code_raw)

pair_rows <- list()
for (kingdom in c("Bacteria", "Fungi")) {
  counts <- if (kingdom == "Bacteria") B0 else F0
  for (pair in list(c("HD", "LD"), c("HD", "wild"), c("LD", "wild"))) {
    z <- mm[mm$treatment %in% pair, , drop = FALSE]
    idx <- match(z$sample_uid, rownames(counts))
    ok <- !is.na(idx)
    z <- z[ok, , drop = FALSE]
    idx <- idx[ok]
    z$treatment <- droplevels(z$treatment)
    if (nrow(z) < 6 || nlevels(z$treatment) < 2) next
    d <- dist(clr(counts[idx, , drop = FALSE]))
    fit <- vegan::adonis2(d ~ treatment + block, data = z, permutations = 999, by = "margin")
    pair_rows[[length(pair_rows) + 1L]] <- data.frame(
      kingdom = kingdom, contrast = paste(pair, collapse = " vs "), n = nrow(z),
      R2 = unname(fit["treatment", "R2"]), P = unname(fit["treatment", "Pr(>F)"]),
      q_within_kingdom = NA_real_, stringsAsFactors = FALSE
    )
  }
}
pairwise <- do.call(rbind, pair_rows)
pairwise$q_within_kingdom <- ave(pairwise$P, pairwise$kingdom, FUN = function(x) p.adjust(x, "BH"))
write.table(pairwise, file.path(outdir, "management_pairwise_exact_R2.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## Optional global FDR sensitivity. Supply the two existing p-value tables explicitly when needed.
cross_file <- Sys.getenv("CROSS_KINGDOM_P_FILE", "")
taxon_file <- Sys.getenv("PH_TAXON_P_FILE", "")
read_p <- function(path, family, filter_mantel = FALSE) {
  if (!nzchar(path) || !file.exists(path)) return(NULL)
  if (family == "cross-kingdom Mantel") {
    raw <- readLines(path, warn = FALSE)
    fields <- strsplit(raw, "\t", fixed = TRUE)
    if (length(fields) >= 2L && any(fields[[1]] == "mantel_rmantel_p") &&
        max(lengths(fields)) >= 9L) {
      hdr <- c("group", "n", "metric", "bacterial_features",
               "fungal_features", "mantel_r", "mantel_p",
               "procrustes_correlation", "procrustes_p")
      rows <- lapply(fields[-1L], function(z) {
        length(z) <- length(hdr)
        z
      })
      x <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
      names(x) <- hdr
    } else {
      x <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
                      fill = TRUE, quote = "", comment.char = "")
    }
  } else {
    x <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
                    fill = TRUE, quote = "", comment.char = "")
  }
  pcol <- intersect(c("p", "P", "permutation_p", "p_value", "P_value",
                      "mantel_p", "mantel_P", "mantel.p", "mantel_p_value"), names(x))
  if (!length(pcol)) return(NULL)
  if (filter_mantel) {
    mcol <- intersect(c("metric", "test", "method", "statistic"), names(x))
    if (length(mcol)) x <- x[grepl("mantel", x[[mcol[[1]]]], ignore.case = TRUE), , drop = FALSE]
  }
  data.frame(family = family, p = suppressWarnings(as.numeric(x[[pcol[[1]]]])), source = path)
}

fdr_inputs <- list(
  data.frame(family = "pH PERMANOVA", p = controlled$P, source = "controlled_pH_exact_R2.tsv"),
  data.frame(family = "management PERMANOVA", p = pairwise$P, source = "management_pairwise_exact_R2.tsv"),
  read_p(cross_file, "cross-kingdom Mantel", FALSE),
  read_p(taxon_file, "pH taxon meta-analysis", FALSE)
)
fdr_inputs <- do.call(rbind, fdr_inputs[!vapply(fdr_inputs, is.null, logical(1))])
fdr_inputs <- fdr_inputs[is.finite(fdr_inputs$p), , drop = FALSE]
if (nrow(fdr_inputs)) {
  fdr_inputs$cross_family_q <- p.adjust(fdr_inputs$p, "BH")
  summary <- do.call(rbind, lapply(split(fdr_inputs, fdr_inputs$family), function(x) data.frame(
    family = x$family[1], n_tests = nrow(x), min_raw_P = min(x$p),
    min_cross_family_q = min(x$cross_family_q), n_cross_family_q_lt_0.05 = sum(x$cross_family_q < 0.05)
  )))
  write.table(fdr_inputs, file.path(outdir, "cross_family_fdr_all_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(summary, file.path(outdir, "cross_family_fdr_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  writeLines("Set CROSS_KINGDOM_P_FILE and PH_TAXON_P_FILE, then rerun to compute cross-family FDR.",
             file.path(outdir, "cross_family_fdr_MISSING_INPUTS.txt"))
}

cat("Wrote Table S3 gap exports to: ", outdir, "\n", sep = "")
