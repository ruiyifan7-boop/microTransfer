#!/usr/bin/env Rscript
# =============================================================================
# A2 (self-contained) — batch-correction specificity control WITHOUT MMUPHin.
# Implements parametric empirical-Bayes ComBat (the engine underlying MMUPHin)
# in base R, so it needs NO Bioconductor packages — only microTransfer.
# Usage:  Rscript A2_mmuphin_comparison.R <rice_data_dir> <out_dir>
# =============================================================================
suppressPackageStartupMessages(library(microTransfer))
args <- commandArgs(trailingOnly = TRUE)
rice_dir <- if (length(args) >= 1) args[1] else "data/rice"
out      <- if (length(args) >= 2) args[2] else "results/combat_official"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
set.seed(3)

rice <- read_edwards_rice(
  otu_file = file.path(rice_dir, "lc_study_otu_table.tsv.gz"),
  metadata_file = file.path(rice_dir, "lc_study_mapping_file.tsv"),
  organelle_file = file.path(rice_dir, "organelle.rds"),
  min_depth = 5000, compartments = c("Rhizosphere", "Endosphere"),
  cohort = "site_year")
counts <- rice$counts; md <- rice$metadata; md$SampleID <- rownames(counts)
batch <- as.character(md$cohort)
comp  <- as.character(md$Compartment)
cult  <- as.character(md$Cultivar)
age   <- suppressWarnings(as.numeric(md$Age)); age[is.na(age)] <- mean(age, na.rm = TRUE)
cult[is.na(cult)] <- "NA"

prev <- colMeans(counts > 0); feat <- which(prev >= 0.10)
clr  <- function(M) { M <- M + 0.5; L <- log(M); sweep(L, 1, rowMeans(L), "-") }
Y    <- clr(counts[, feat, drop = FALSE])                 # n x F
dmat <- function(x) { x <- as.factor(x); if (nlevels(x) < 2) return(matrix(0, length(x), 0)); model.matrix(~x)[, -1, drop = FALSE] }

## ---- parametric empirical-Bayes ComBat (base R) ----
combat <- function(dat, batch, mod) {           # dat: n x F
  n <- nrow(dat); Fn <- ncol(dat)
  b <- as.factor(batch); lev <- levels(b); nb <- length(lev)
  Bat <- model.matrix(~ b - 1)                  # n x nb
  design <- cbind(Bat, mod)
  B <- solve(crossprod(design), crossprod(design, dat))     # (nb+p) x F
  ni <- colSums(Bat); grand <- as.numeric((ni / n) %*% B[1:nb, , drop = FALSE])
  stand <- matrix(grand, n, Fn, byrow = TRUE)
  if (ncol(mod) > 0) stand <- stand + mod %*% B[(nb + 1):nrow(B), , drop = FALSE]
  resid <- dat - design %*% B
  vp <- colMeans(resid^2) + 1e-12
  s <- (dat - stand) / matrix(sqrt(vp), n, Fn, byrow = TRUE)
  gamma <- matrix(0, nb, Fn); delta <- matrix(0, nb, Fn)
  for (i in 1:nb) { m <- b == lev[i]
    gamma[i, ] <- colMeans(s[m, , drop = FALSE])
    delta[i, ] <- apply(s[m, , drop = FALSE], 2, var) + 1e-8 }
  gbar <- rowMeans(gamma); tau2 <- apply(gamma, 1, var)
  dbar <- rowMeans(delta); dvar <- apply(delta, 1, var)
  apr <- (2 * dvar + dbar^2) / dvar; bpr <- (dbar * dvar + dbar^3) / dvar
  gstar <- gamma; dstar <- delta
  for (i in 1:nb) { m <- b == lev[i]; ni_ <- sum(m); sdat <- s[m, , drop = FALSE]
    g <- gamma[i, ]; dl <- delta[i, ]
    for (it in 1:100) {
      gnew <- (ni_ * tau2[i] * gamma[i, ] + dl * gbar[i]) / (ni_ * tau2[i] + dl)
      ss <- colSums(sweep(sdat, 2, gnew, "-")^2)
      dnew <- (0.5 * ss + bpr[i]) / (ni_ / 2 + apr[i] - 1)
      if (max(abs(gnew - g)) < 1e-4 && max(abs(dnew - dl)) < 1e-4) { g <- gnew; dl <- dnew; break }
      g <- gnew; dl <- dnew }
    gstar[i, ] <- g; dstar[i, ] <- dl }
  adj <- s
  for (i in 1:nb) { m <- b == lev[i]
    adj[m, ] <- sweep(sweep(s[m, , drop = FALSE], 2, gstar[i, ], "-"), 2, sqrt(dstar[i, ]), "/") }
  adj * matrix(sqrt(vp), n, Fn, byrow = TRUE) + stand
}
mod <- cbind(dmat(comp), age - mean(age))
CB  <- combat(Y, batch, mod)

## ---- study variance explained (Aitchison PERMANOVA R^2, subsampled) ----
permR2 <- function(X, lab, nsub = 300) {
  idx <- sample(nrow(X), min(nsub, nrow(X))); Xs <- X[idx, ]; g <- as.factor(lab[idx])
  D2 <- as.matrix(dist(Xs))^2; N <- nrow(Xs); sst <- sum(D2) / (2 * N); sw <- 0
  for (l in levels(g)) { ii <- which(g == l); ng <- length(ii); sw <- sw + sum(D2[ii, ii]) / (2 * ng) }
  (sst - sw) / sst
}
r2_raw <- permR2(Y, batch); r2_cb <- permR2(CB, batch)

## ---- pooled age-response significance: raw vs ComBat + out-of-study replication ----
resid_on <- function(Y, D) { D <- cbind(1, D); Y - D %*% solve(crossprod(D), crossprod(D, Y)) }
pool_sig <- function(X) {
  D <- cbind(dmat(comp), dmat(cult)); Yr <- resid_on(X, D)
  ar <- as.numeric(resid_on(matrix(age, ncol = 1), D)); arn <- (ar - mean(ar)) / sd(ar)
  Yn <- scale(Yr); r <- as.numeric(crossprod(Yn, arn)) / (nrow(X) - 1)
  tval <- r * sqrt((nrow(X) - 2) / (1 - r^2 + 1e-9)); p <- 2 * pt(-abs(tval), nrow(X) - 2)
  list(r = r, q = p.adjust(p, "BH"))
}
raw <- pool_sig(Y); cb <- pool_sig(CB)
studies <- sort(unique(batch))
heldrepl <- function(cand, refsign) {
  fr <- c()
  for (h in studies) { m <- batch == h; D <- cbind(dmat(comp[m]), dmat(cult[m]))
    Yr <- resid_on(Y[m, cand, drop = FALSE], D)
    ar <- as.numeric(resid_on(matrix(age[m], ncol = 1), D)); arn <- (ar - mean(ar)) / sd(ar)
    Yn <- scale(Yr); rr <- as.numeric(crossprod(Yn, arn)) / (sum(m) - 1)
    tt <- rr * sqrt((sum(m) - 2) / (1 - rr^2 + 1e-9)); pp <- 2 * pt(-abs(tt), sum(m) - 2)
    fr <- c(fr, mean(pp < 0.05 & sign(rr) == refsign)) }
  mean(fr)
}
cand_cb <- which(cb$q <= 0.05)
repl_cb <- heldrepl(cand_cb, sign(cb$r[cand_cb]))

contrast <- pooled_membership_contrast(counts, md, "SampleID", "cohort", thresholds = 0.50)

res <- c(
  sprintf("study_R2_raw\t%.4f", r2_raw),
  sprintf("study_R2_combat\t%.4f", r2_cb),
  sprintf("age_sig_raw\t%d", sum(raw$q <= 0.05, na.rm = TRUE)),
  sprintf("age_sig_combat\t%d", sum(cb$q <= 0.05, na.rm = TRUE)),
  sprintf("age_combat_heldout_replication\t%.3f", repl_cb),
  sprintf("membership_pooled_core_50\t%d", contrast$pooled_candidates[1]),
  sprintf("membership_studyaware_core_50\t%d", contrast$study_aware_candidates[1]))
writeLines(res, file.path(out, "combat_summary.tsv"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("DONE (self-contained ComBat; no MMUPHin needed)\n"); print(res)
