options(stringsAsFactors=FALSE, width=180)

root <- "global_harmonized"
outdir <- file.path(root, "enhancement", "management_abundance")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

counts <- readRDS(file.path(
  root, "taxa_tables", "bacteria_genus_counts_primary211.rds"
))
core <- read.delim(
  file.path(root, "core_microbiome", "bacteria_core_audit.tsv"),
  check.names=FALSE
)
meta <- read.delim(
  "PRJEB110492/metadata/PRJEB110492_master_metadata_provisional.tsv",
  check.names=FALSE, na.strings=c("", "NA")
)

meta$sample_uid <- paste(
  "PRJEB110492", meta$biological_sample, sep="__"
)
meta <- meta[
  meta$design_family == "transplanted_field_system" &
    meta$treatment_code_raw %in% c("HD", "LD", "wild") &
    meta$sample_uid %in% rownames(counts),
  , drop=FALSE
]
meta$treatment <- factor(
  meta$treatment_code_raw, levels=c("HD", "LD", "wild")
)
meta$block <- factor(meta$block_code_raw)
counts <- counts[meta$sample_uid, , drop=FALSE]

prevalence_cols <- grep(
  "^prevalence_", names(core), value=TRUE
)
discovery_cols <- setdiff(
  prevalence_cols, "prevalence_PRJEB110492"
)
keep_core <- apply(
  core[, discovery_cols, drop=FALSE], 1,
  function(z) all(!is.na(z) & z >= 0.50)
)
features <- intersect(core$feature[keep_core], colnames(counts))
if (!length(features))
  stop("No externally defined bacterial core features mapped")

x <- counts[, features, drop=FALSE]
clr <- log(x + 0.5)
clr <- sweep(clr, 1, rowMeans(clr), "-")

contrast_test <- function(fit, weights, label, feature) {
  beta <- coef(fit)
  v <- vcov(fit)
  common <- intersect(names(weights), names(beta))
  w <- weights[common]
  estimate <- sum(w * beta[common])
  se <- sqrt(as.numeric(t(w) %*% v[common, common, drop=FALSE] %*% w))
  statistic <- estimate / se
  p <- 2 * pt(abs(statistic), df=df.residual(fit), lower.tail=FALSE)
  data.frame(
    feature=feature, contrast=label,
    clr_difference=estimate, standard_error=se,
    statistic=statistic, df=df.residual(fit), p=p
  )
}

omnibus <- vector("list", length(features))
contrasts <- vector("list", length(features))

for (i in seq_along(features)) {
  feature <- features[i]
  dat <- data.frame(
    y=clr[, feature],
    treatment=meta$treatment,
    block=meta$block
  )
  reduced <- lm(y ~ block, data=dat)
  full <- lm(y ~ treatment + block, data=dat)
  cmp <- anova(reduced, full)
  full_anova <- anova(full)
  treatment_ss <- full_anova["treatment", "Sum Sq"]
  residual_ss <- full_anova["Residuals", "Sum Sq"]
  omnibus[[i]] <- data.frame(
    feature=feature,
    treatment_F=cmp$F[2],
    treatment_p=cmp$`Pr(>F)`[2],
    partial_eta2=treatment_ss / (treatment_ss + residual_ss),
    mean_RA_HD=mean(x[meta$treatment == "HD", feature] /
                      rowSums(counts[meta$treatment == "HD", ])),
    mean_RA_LD=mean(x[meta$treatment == "LD", feature] /
                      rowSums(counts[meta$treatment == "LD", ])),
    mean_RA_wild=mean(x[meta$treatment == "wild", feature] /
                        rowSums(counts[meta$treatment == "wild", ]))
  )
  contrasts[[i]] <- rbind(
    contrast_test(
      full, c(treatmentLD=1), "LD_minus_HD", feature
    ),
    contrast_test(
      full, c(treatmentwild=1), "wild_minus_HD", feature
    ),
    contrast_test(
      full,
      c(treatmentwild=1, treatmentLD=-1),
      "wild_minus_LD", feature
    )
  )
}

omnibus <- do.call(rbind, omnibus)
omnibus$treatment_q <- p.adjust(
  omnibus$treatment_p, method="BH"
)
omnibus <- omnibus[order(
  omnibus$treatment_q, -omnibus$partial_eta2
), ]

contrasts <- do.call(rbind, contrasts)
contrasts$contrast_q <- ave(
  contrasts$p, contrasts$contrast,
  FUN=function(z) p.adjust(z, method="BH")
)
contrasts <- contrasts[order(
  contrasts$contrast, contrasts$contrast_q
), ]

stable_path <- file.path(
  root, "enhancement", "transferability",
  "management_stable_features_with_direction.tsv"
)
stable <- read.delim(stable_path, check.names=FALSE)
stable <- stable[
  stable$kind == "bacteria" & stable$selected_blocks == 4,
  , drop=FALSE
]
stable$feature <- sub("^B__", "", stable$feature)
stable_validation <- merge(
  stable, omnibus, by="feature", all.x=TRUE
)
stable_validation <- merge(
  stable_validation, contrasts, by="feature", all.x=TRUE
)

write.table(
  omnibus,
  file.path(outdir, "core_bacteria_treatment_omnibus.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)
write.table(
  contrasts,
  file.path(outdir, "core_bacteria_treatment_contrasts.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)
write.table(
  stable_validation,
  file.path(outdir, "stable_predictor_abundance_validation.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)

cat("Samples:", nrow(meta), "\n")
print(table(meta$block, meta$treatment))
cat("Externally defined core bacterial genera:", length(features), "\n")
cat("Omnibus q<0.05:", sum(omnibus$treatment_q < 0.05), "\n")
cat("Pairwise q<0.05:", sum(contrasts$contrast_q < 0.05), "\n")
cat("\nStable predictor validation:\n")
print(stable_validation, row.names=FALSE)
cat("\nBLOCKED ABUNDANCE ANALYSIS DONE\n")
