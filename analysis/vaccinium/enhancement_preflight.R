options(stringsAsFactors=FALSE, width=160)

root <- "global_harmonized"
out <- file.path(root, "enhancement")
dir.create(out, recursive=TRUE, showWarnings=FALSE)

needed <- c(
  metadata=file.path(root, "global_metadata_primary_211.tsv"),
  bacteria=file.path(root, "analysis_ready", "bacteria_genus_model_counts.rds"),
  fungi=file.path(root, "analysis_ready", "fungi_genus_model_counts.rds"),
  bacteria_core=file.path(root, "core_microbiome", "bacteria_core_audit.tsv"),
  fungi_core=file.path(root, "core_microbiome", "fungi_core_audit.tsv"),
  bacteria_pH=file.path(root, "pH_analysis", "bacteria_genus_pH_taxa_validation.tsv"),
  fungi_pH=file.path(root, "pH_analysis", "fungi_genus_pH_taxa_validation.tsv"),
  management="PRJEB110492/metadata/PRJEB110492_master_metadata_provisional.tsv"
)

cat("===== FILE AUDIT =====\n")
file_audit <- data.frame(
  object=names(needed),
  path=unname(needed),
  exists=file.exists(needed),
  size=ifelse(file.exists(needed), file.info(needed)$size, NA_real_)
)
print(file_audit, row.names=FALSE)
write.table(
  file_audit, file.path(out, "preflight_files.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)

essential <- c("metadata", "bacteria", "fungi")
if (!all(file_audit$exists[match(essential, file_audit$object)])) {
  stop("Essential metadata/count files are missing")
}

packages <- c(
  "vegan", "glmnet", "ranger", "pROC", "mgcv",
  "ANCOMBC", "ALDEx2", "Maaslin2", "LinDA",
  "SpiecEasi", "NetCoMi", "propr"
)
pkg_audit <- data.frame(
  package=packages,
  available=vapply(packages, requireNamespace, logical(1), quietly=TRUE),
  version=vapply(packages, function(p) {
    if (requireNamespace(p, quietly=TRUE))
      as.character(utils::packageVersion(p))
    else NA_character_
  }, character(1))
)
cat("\n===== PACKAGE AUDIT =====\n")
print(pkg_audit, row.names=FALSE)
write.table(
  pkg_audit, file.path(out, "preflight_packages.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)

meta <- read.delim(
  needed[["metadata"]], check.names=FALSE,
  stringsAsFactors=FALSE, na.strings=c("", "NA")
)
B <- readRDS(needed[["bacteria"]])
F <- readRDS(needed[["fungi"]])

cat("\n===== MATRIX AUDIT =====\n")
cat("Metadata:", nrow(meta), "rows x", ncol(meta), "columns\n")
cat("Bacteria:", nrow(B), "samples x", ncol(B), "features\n")
cat("Fungi:", nrow(F), "samples x", ncol(F), "features\n")

if (!"sample_uid" %in% names(meta))
  stop("sample_uid column is missing from global metadata")
if (anyDuplicated(meta$sample_uid))
  stop("Duplicated sample_uid values in metadata")
if (anyDuplicated(rownames(B)) || anyDuplicated(rownames(F)))
  stop("Duplicated sample IDs in count matrices")

cat("Bacteria IDs absent from metadata:",
    sum(!rownames(B) %in% meta$sample_uid), "\n")
cat("Fungi IDs absent from metadata:",
    sum(!rownames(F) %in% meta$sample_uid), "\n")
cat("Primary metadata absent from bacteria:",
    sum(!meta$sample_uid %in% rownames(B)), "\n")
cat("Primary metadata absent from fungi:",
    sum(!meta$sample_uid %in% rownames(F)), "\n")

shared <- Reduce(intersect, list(meta$sample_uid, rownames(B), rownames(F)))
cat("Exactly shared paired samples:", length(shared), "\n")

md <- meta[match(shared, meta$sample_uid), , drop=FALSE]
B <- B[shared, , drop=FALSE]
F <- F[shared, , drop=FALSE]
stopifnot(identical(md$sample_uid, rownames(B)),
          identical(rownames(B), rownames(F)))

cat("\n===== METADATA COLUMNS =====\n")
print(names(md))

candidate_columns <- grep(
  "study|compartment|sample_type|pH|treatment|management|block|group|host|species|cultivar|location|country",
  names(md), ignore.case=TRUE, value=TRUE
)
cat("\n===== CANDIDATE DESIGN FIELDS =====\n")
for (v in candidate_columns) {
  x <- md[[v]]
  nonmissing <- sum(!is.na(x) & trimws(as.character(x)) != "")
  values <- unique(as.character(x[!is.na(x) & trimws(as.character(x)) != ""]))
  cat("\n--", v, "-- nonmissing:", nonmissing,
      "unique:", length(values), "\n")
  if (length(values) <= 30) print(sort(values))
}

pick_column <- function(dat, choices) {
  hit <- choices[choices %in% names(dat)]
  if (length(hit)) hit[1] else NA_character_
}

study_col <- pick_column(md, c("study", "project", "project_accession"))
comp_col <- pick_column(
  md, c("compartment", "compartment_normalized", "sample_type")
)
pH_col <- pick_column(md, c("pH", "pH_value", "ph"))
treat_col <- pick_column(
  md, c("pH_treatment", "treatment", "treatment_normalized")
)
block_col <- pick_column(md, c("block", "block_code_raw"))

cat("\nSelected fields:\n")
print(c(
  study=study_col, compartment=comp_col, pH=pH_col,
  treatment=treat_col, block=block_col
))

if (is.na(study_col))
  stop("Could not identify study column")

norm_compartment <- function(x) {
  z <- tolower(trimws(as.character(x)))
  ifelse(
    grepl("rhizo", z), "Rhizosphere",
    ifelse(
      grepl("bulk|soil", z), "Bulk",
      ifelse(grepl("endo|root", z), "Endosphere", NA_character_)
    )
  )
}

md$study_model <- as.character(md[[study_col]])
md$compartment_model <- if (!is.na(comp_col))
  norm_compartment(md[[comp_col]]) else NA_character_
md$pH_model <- if (!is.na(pH_col))
  suppressWarnings(as.numeric(md[[pH_col]])) else NA_real_
md$pH_binary_4_0_5_5 <- ifelse(
  is.na(md$pH_model), NA_character_,
  ifelse(md$pH_model >= 4.0 & md$pH_model <= 5.5,
         "Favorable_4.0_5.5", "Outside_4.0_5.5")
)

cat("\n===== STUDY COUNTS =====\n")
print(table(md$study_model, useNA="ifany"))

cat("\n===== SHARED COMPARTMENT TRANSFER TASK =====\n")
comp_task <- md[
  md$study_model %in% c("PRJNA1156347", "PRJEB98254") &
    md$compartment_model %in% c("Bulk", "Rhizosphere"),
  , drop=FALSE
]
print(with(
  comp_task,
  table(study_model, compartment_model, useNA="ifany")
))

cat("\n===== SHARED RHIZOSPHERE pH TRANSFER TASK =====\n")
pH_task <- md[
  md$study_model %in% c("PRJNA1156347", "PRJEB98254") &
    md$compartment_model == "Rhizosphere" &
    !is.na(md$pH_binary_4_0_5_5),
  , drop=FALSE
]
print(with(
  pH_task,
  table(study_model, pH_binary_4_0_5_5, useNA="ifany")
))
cat("pH summaries:\n")
print(tapply(pH_task$pH_model, pH_task$study_model, summary))

cat("\n===== CONTROLLED pH TASK =====\n")
controlled <- md[
  md$study_model == "PRJNA1156347" &
    !is.na(md$compartment_model),
  , drop=FALSE
]
if (!is.na(treat_col)) {
  controlled$pH_treatment_model <- as.character(controlled[[treat_col]])
  print(with(
    controlled,
    table(compartment_model, pH_treatment_model, useNA="ifany")
  ))
} else {
  cat("No pH treatment column detected\n")
}

cat("\n===== MANAGEMENT SOURCE AUDIT =====\n")
if (file.exists(needed[["management"]])) {
  mm <- read.delim(
    needed[["management"]], check.names=FALSE,
    stringsAsFactors=FALSE, na.strings=c("", "NA")
  )
  cat("Rows:", nrow(mm), "\n")
  print(names(mm))
  management_fields <- grep(
    "treat|block|design|group|sample",
    names(mm), ignore.case=TRUE, value=TRUE
  )
  for (v in management_fields) {
    values <- unique(na.omit(as.character(mm[[v]])))
    cat("\n--", v, "--\n")
    if (length(values) <= 40) print(sort(values))
  }
}

terminal_genus <- function(x) {
  z <- sub("^.*;g__", "", x)
  sub(";.*$", "", z)
}

feature_audit <- data.frame(
  kingdom=c("Bacteria", "Fungi"),
  features=c(ncol(B), ncol(F)),
  terminal_genera=c(
    length(unique(terminal_genus(colnames(B)))),
    length(unique(terminal_genus(colnames(F))))
  ),
  duplicated_terminal_labels=c(
    sum(duplicated(terminal_genus(colnames(B)))),
    sum(duplicated(terminal_genus(colnames(F))))
  )
)
cat("\n===== FEATURE AUDIT =====\n")
print(feature_audit, row.names=FALSE)
write.table(
  feature_audit, file.path(out, "preflight_features.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)

audit_feature_file <- function(path, matrix_columns, label) {
  if (!file.exists(path)) {
    cat(label, ": missing\n")
    return(NULL)
  }
  x <- read.delim(path, check.names=FALSE, stringsAsFactors=FALSE)
  cat("\n", label, ":", nrow(x), "rows\n", sep="")
  print(names(x))
  genus_col <- pick_column(
    x, c("genus", "Genus", "taxon", "feature", "Feature")
  )
  if (!is.na(genus_col)) {
    target <- unique(sub("^g__", "", as.character(x[[genus_col]])))
    available <- unique(terminal_genus(matrix_columns))
    cat("Mapped terminal genus labels:",
        sum(target %in% available), "/", length(target), "\n")
    cat("Unmapped examples:\n")
    print(head(setdiff(target, available), 20))
  }
  invisible(x)
}

audit_feature_file(
  needed[["bacteria_core"]], colnames(B), "Bacterial core audit"
)
audit_feature_file(
  needed[["fungi_core"]], colnames(F), "Fungal core audit"
)
audit_feature_file(
  needed[["bacteria_pH"]], colnames(B), "Bacterial pH candidates"
)
audit_feature_file(
  needed[["fungi_pH"]], colnames(F), "Fungal pH candidates"
)

write.table(
  md[, c(
    "sample_uid", "study_model", "compartment_model",
    "pH_model", "pH_binary_4_0_5_5"
  )],
  file.path(out, "preflight_model_metadata.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE, na=""
)

cat("\n===== DECISION FLAGS =====\n")
cat("Elastic-net transfer models:", pkg_audit$available[pkg_audit$package=="glmnet"], "\n")
cat("GAM nonlinear sensitivity:", pkg_audit$available[pkg_audit$package=="mgcv"], "\n")
cat("ANCOM-BC2 differential abundance:", pkg_audit$available[pkg_audit$package=="ANCOMBC"], "\n")
cat("ALDEx2 differential abundance:", pkg_audit$available[pkg_audit$package=="ALDEx2"], "\n")
cat("SPIEC-EASI network:", pkg_audit$available[pkg_audit$package=="SpiecEasi"], "\n")
cat("\nPREFLIGHT DONE\n")
