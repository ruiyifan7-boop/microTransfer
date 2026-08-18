#' Read the Edwards et al. field-grown rice lifecycle dataset
#'
#' The public dataset is archived at Dryad DOI
#' `10.5061/dryad.7q7k1`. The function reads the downloaded OTU table and
#' mapping file; it does not download data implicitly.
#'
#' @param otu_file Path to `lc_study_otu_table.tsv.gz`.
#' @param metadata_file Path to `lc_study_mapping_file.tsv`.
#' @param taxonomy_file Optional path to `gg_otus_tax.rds`.
#' @param organelle_file Optional path to `organelle.rds`; matching OTUs are
#'   removed before analysis.
#' @param min_depth Minimum sample read depth.
#' @param compartments Optional compartments to retain.
#' @param sites Field sites to retain. The default reproduces the two-site
#'   lifecycle study rather than auxiliary reference datasets.
#' @param cohort One of `"site_year"` or `"site"`.
#' @param rank Optional taxonomy rank for aggregation.
#' @param taxonomy_feature_col Optional taxonomy feature-ID column.
#' @return A list containing counts, metadata, taxonomy, and source details.
#' @export
read_edwards_rice <- function(
    otu_file,
    metadata_file,
    taxonomy_file = NULL,
    organelle_file = NULL,
    min_depth = 5000,
    compartments = c("Rhizosphere", "Endosphere"),
    sites = c("Arbuckle", "Jonesboro"),
    cohort = c("site_year", "site"),
    rank = NULL,
    taxonomy_feature_col = NULL) {
  cohort <- match.arg(cohort)
  metadata <- utils::read.delim(
    metadata_file,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA")
  )
  required <- c("SampleID", "Site")
  missing <- setdiff(required, names(metadata))
  if (length(missing)) {
    stop("Rice metadata is missing: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  con <- if (grepl("\\.gz$", otu_file, ignore.case = TRUE)) {
    gzfile(otu_file, open = "rt")
  } else {
    otu_file
  }
  otu <- utils::read.delim(
    con,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (inherits(con, "connection")) close(con)
  if (!"OTUID" %in% names(otu)) {
    stop("Rice OTU table must contain an `OTUID` column.", call. = FALSE)
  }
  feature_ids <- as.character(otu$OTUID)
  otu$OTUID <- NULL
  counts <- t(as.matrix(otu))
  storage.mode(counts) <- "double"
  colnames(counts) <- feature_ids
  if (!is.null(organelle_file)) {
    organelle <- as.character(readRDS(organelle_file))
    counts <- counts[, setdiff(colnames(counts), organelle), drop = FALSE]
  }
  shared <- intersect(as.character(metadata$SampleID), rownames(counts))
  metadata <- metadata[match(shared, metadata$SampleID), , drop = FALSE]
  counts <- counts[shared, , drop = FALSE]
  depth <- rowSums(counts)
  keep <- depth >= min_depth
  if (!is.null(compartments) && "Compartment" %in% names(metadata)) {
    keep <- keep & metadata$Compartment %in% compartments
  }
  if (!is.null(sites)) keep <- keep & metadata$Site %in% sites
  metadata <- metadata[keep, , drop = FALSE]
  counts <- counts[keep, , drop = FALSE]
  metadata$Depth_recomputed <- rowSums(counts)
  if (cohort == "site_year") {
    if (!"Season" %in% names(metadata)) {
      stop("`Season` is required for site-year cohorts.", call. = FALSE)
    }
    metadata$cohort <- paste(metadata$Site, metadata$Season, sep = "_")
  } else {
    metadata$cohort <- as.character(metadata$Site)
  }
  taxonomy <- NULL
  if (!is.null(taxonomy_file)) {
    taxonomy <- readRDS(taxonomy_file)
    taxonomy <- as.data.frame(taxonomy, stringsAsFactors = FALSE)
    if (!is.null(rank)) {
      counts <- collapse_counts_by_taxonomy(
        counts, taxonomy, rank, taxonomy_feature_col
      )
    }
  } else if (!is.null(rank)) {
    stop("`taxonomy_file` is required when `rank` is requested.",
         call. = FALSE)
  }
  list(
    counts = counts,
    metadata = metadata,
    taxonomy = taxonomy,
    source = list(
      dataset = "Edwards et al. field-grown rice lifecycle microbiome",
      dryad_doi = "10.5061/dryad.7q7k1",
      article_doi = "10.1371/journal.pbio.2003862",
      feature_definition = if (is.null(rank)) "GreenGenes 13_8 97% OTU" else rank,
      cohort_definition = cohort
    )
  )
}
