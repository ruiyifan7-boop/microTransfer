## run_pH_design_aware_robustness.R
## Design-aware robustness of the field bacterial pH signature (reviewer 1-item #5)
## Stress-tests the quadratic signature~pH fit (reported R2=0.242, P=0.025) and the
## community-level pH PERMANOVA. Writes a tidy table for Supplementary Table S3 (3g).
## Run on the analysis server (needs the frozen matrices/metadata). set.seed(20260706).

suppressMessages({library(vegan); library(splines); library(dplyr)})
set.seed(20260706)
outdir <- "global_harmonized/enhancement/review"; dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

## ---------- INPUTS (export these from the frozen analysis) ----------
## 1) Field signature scores in PRJEB98254 rhizosphere, projected WITHOUT refit,
##    one row per sample: columns must include score, pH, and any covariates
##    (cultivar, site, block, SOM, TC, TN, ...). Extra columns auto-detected.
SIG <- "global_harmonized/pH_analysis/pH_signature_field_scores_PRJEB98254.tsv"
## 2) (community part) CLR bacterial genus matrix + metadata for the same cohort
CLR <- "global_harmonized/analysis_ready/PRJEB98254_rhizo_bacteria_clr.rds"   # samples x features
MET <- "global_harmonized/analysis_ready/PRJEB98254_rhizo_metadata.tsv"       # sample_uid, pH, covars

dat <- read.delim(SIG, check.names=FALSE)
stopifnot(all(c("score","pH") %in% names(dat)))
dat$pHz <- as.numeric(scale(dat$pH))
covars <- setdiff(names(dat), c("sample","sample_uid","score","pH","pHz"))
cat("Detected covariates:", if(length(covars)) paste(covars,collapse=", ") else "NONE (only pH available)","\n")

AICc <- function(m){k<-length(coef(m))+1; n<-length(resid(m)); AIC(m)+(2*k*(k+1))/(n-k-1)}
loocv_mse <- function(form,d){
  pr<-numeric(nrow(d)); for(i in seq_len(nrow(d))){fit<-lm(form,d[-i,]); pr[i]<-predict(fit,d[i,])}
  mean((d$score-pr)^2)}

res <- list()

## ---- Check 1: model comparison M0/M1/M2/M3 ----
cov_str <- if(length(covars)) paste("+",paste(sprintf("`%s`",covars),collapse="+")) else ""
f0<-as.formula(paste("score ~ 1", cov_str))
f1<-as.formula(paste("score ~ pHz", cov_str))
f2<-as.formula(paste("score ~ pHz + I(pHz^2)", cov_str))
f3<-as.formula(paste("score ~ ns(pHz,3)", cov_str))
mods<-list(M0_covars=lm(f0,dat),M1_linear=lm(f1,dat),M2_quadratic=lm(f2,dat),M3_spline=lm(f3,dat))
for(nm in names(mods)){m<-mods[[nm]]
  res[[nm]]<-data.frame(model=nm, adjR2=summary(m)$adj.r.squared, AICc=AICc(m),
    LOOCV_MSE=loocv_mse(formula(m),dat))}
cat("\n== Check 1: model comparison ==\n"); print(do.call(rbind,res),row.names=FALSE)

## ---- Check 3: covariate-adjusted quadratic term ----
q_p_adj <- tryCatch(summary(mods$M2_quadratic)$coefficients["I(pHz^2)","Pr(>|t|)"], error=function(e) NA)
cat("\n== Check 3: quadratic-term P after covariate adjustment ==", q_p_adj, "\n")

## ---- Check 2: restricted permutation of pH within block/site ----
grp <- if("block" %in% covars) dat$block else if("site" %in% covars) dat$site else NULL
obsR2 <- summary(mods$M2_quadratic)$r.squared
perm <- replicate(999, {
  d2<-dat
  if(is.null(grp)) d2$pHz<-sample(d2$pHz) else
    d2$pHz<-ave(d2$pHz, grp, FUN=function(x) sample(x))
  summary(lm(f2,d2))$r.squared })
p_restr <- (sum(perm>=obsR2)+1)/1000
cat("== Check 2: restricted-permutation P (within", if(is.null(grp))"NONE (unrestricted)" else "block/site", ") ==", p_restr, "\n")

## ---- Check 4: endpoint exclusion (drop top/bottom 5% pH) ----
lo<-quantile(dat$pH,.05); hi<-quantile(dat$pH,.95)
dd<-dat[dat$pH>lo & dat$pH<hi,]
m_end<-lm(f2,dd); cat("== Check 4: after dropping pH <5%/>95% (n=",nrow(dd),") quad adjR2=",
  summary(m_end)$adj.r.squared," quadP=",summary(m_end)$coefficients["I(pHz^2)","Pr(>|t|)"],"\n")

## ---- Check 5: influence diagnostics ----
ck<-cooks.distance(mods$M2_quadratic); infl<-which(ck > 4/nrow(dat))
cat("== Check 5: high-influence samples (Cook's D > 4/n):", length(infl),
    if(length(infl)) paste0("[rows: ",paste(infl,collapse=","),"]") else "", "\n")
m_noinfl<-lm(f2, dat[setdiff(seq_len(nrow(dat)),infl),])
cat("   quad adjR2 without influential points =", summary(m_noinfl)$adj.r.squared, "\n")

## ---- Check 6: leave-one-group-out ----
if(!is.null(grp)){
  cat("== Check 6: leave-one-group-out quadratic adjR2 ==\n")
  for(g in unique(grp)){mm<-lm(f2, dat[grp!=g,]); cat("  drop",g,":",summary(mm)$adj.r.squared,"\n")}
} else cat("== Check 6: no block/site column -> leave-one-group-out not possible ==\n")

## ---- Community-level design-aware PERMANOVA (marginal SS, restricted perm) ----
if(file.exists(CLR) && file.exists(MET)){
  X<-readRDS(CLR); md<-read.delim(MET,check.names=FALSE)
  md<-md[match(rownames(X),md$sample_uid),]; md$pHz<-as.numeric(scale(md$pH))
  D<-dist(X)  # Aitchison = Euclidean on CLR
  hh<-if("block"%in%names(md)) how(blocks=factor(md$block),nperm=999) else how(nperm=999)
  cov2<-intersect(c("cultivar","site","SOM","TC","TN"),names(md))
  f_comm<-as.formula(paste("D ~", if(length(cov2))paste(cov2,collapse="+")else"1","+ pHz + I(pHz^2)"))
  a<-adonis2(f_comm, data=md, by="margin", permutations=hh)  # marginal (type III)
  cat("\n== Community-level design-aware PERMANOVA (marginal SS, restricted perm) ==\n"); print(a)
  capture.output(print(a), file=file.path(outdir,"pH_community_designaware.txt"))
} else cat("\n[community part skipped: CLR/metadata not found]\n")

## ---- write tidy S3-3g table ----
s3 <- rbind(
  do.call(rbind,res),
  data.frame(model="quad_P_covariate_adjusted", adjR2=NA, AICc=NA, LOOCV_MSE=NA),
  data.frame(model="quad_P_restricted_perm",    adjR2=NA, AICc=NA, LOOCV_MSE=NA))
write.table(do.call(rbind,res), file.path(outdir,"S3_3g_pH_model_comparison.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE)
saveRDS(list(model_cmp=do.call(rbind,res), q_p_adj=q_p_adj, p_restricted=p_restr,
  endpoint=summary(m_end), n_influential=length(infl)),
  file.path(outdir,"pH_designaware_summary.rds"))
cat("\nDONE -> ", outdir, "\n")
