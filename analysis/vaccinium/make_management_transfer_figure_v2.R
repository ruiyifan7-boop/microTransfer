suppressPackageStartupMessages(library(ggplot2))
library(grid)

options(stringsAsFactors=FALSE)

root <- "global_harmonized"
tdir <- file.path(root, "enhancement", "transferability")
adir <- file.path(root, "enhancement", "management_abundance")
outdir <- file.path(root, "enhancement", "figures_nature")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

res <- read.delim(file.path(tdir, "transferability_summary.tsv"))
null <- read.delim(file.path(tdir, "transferability_permutation_null.tsv"))
pred <- read.delim(file.path(tdir, "transferability_predictions.tsv"))
validation <- read.delim(file.path(
  adir, "stable_predictor_abundance_validation.tsv"
))

model_cols <- c(
  bacteria="#0072B2", fungi="#CC79A7", combined="#202020"
)
treatment_cols <- c(HD="#4477AA", LD="#EEAA33", wild="#228833")

theme_journal <- theme_classic(base_size=7, base_family="Arial") +
  theme(
    text=element_text(family="Arial", colour="#1E1E1E"),
    axis.text=element_text(size=7, colour="#1E1E1E"),
    axis.title=element_text(size=7, face="bold"),
    axis.line=element_line(linewidth=0.35),
    axis.ticks=element_line(linewidth=0.35),
    plot.title=element_blank(),
    plot.subtitle=element_blank(),
    plot.tag=element_text(face="bold", size=8),
    plot.tag.position=c(0, 1),
    legend.title=element_text(face="bold", size=6),
    legend.text=element_text(size=6),
    legend.key.height=unit(3.2, "mm"),
    strip.background=element_rect(
      fill="#F3F3F3", colour="#BDBDBD", linewidth=0.3
    ),
    strip.text=element_text(face="bold", size=7.5),
    panel.grid=element_blank(),
    plot.margin=margin(6, 7, 5, 6)
  )

# A. Compact effect-style display of observed accuracy and permutation nulls.
mres <- res[res$task == "management",]
mnull <- null[null$task == "management",]
ns <- do.call(rbind, lapply(
  split(mnull, interaction(mnull$kind, mnull$variant, drop=TRUE)),
  function(z) data.frame(
    kind=z$kind[1], variant=z$variant[1],
    low=unname(quantile(z$balanced_accuracy, 0.025)),
    median=median(z$balanced_accuracy),
    high=unname(quantile(z$balanced_accuracy, 0.975))
  )
))
pa <- merge(mres, ns, by=c("kind", "variant"))
pa$model <- factor(
  paste(pa$kind, pa$variant, sep="_"),
  levels=rev(c(
    "bacteria_full", "bacteria_loso_core",
    "fungi_full", "fungi_loso_core",
    "combined_full", "combined_loso_core"
  )),
  labels=rev(c(
    "Bacteria | Full", "Bacteria | Cross-study core",
    "Fungi | Full", "Fungi | Cross-study core",
    "Combined | Full", "Combined | Cross-study core"
  ))
)
pa$core <- pa$variant == "loso_core"
pa$value_label <- sprintf("%.2f", pa$balanced_accuracy)

pA <- ggplot(pa, aes(y=model)) +
  geom_vline(
    xintercept=1/3, linetype="dashed",
    colour="#8A8A8A", linewidth=0.35
  ) +
  geom_segment(
    aes(x=low, xend=high, yend=model),
    colour="#B7B7B7", linewidth=2.8, lineend="round"
  ) +
  geom_point(
    aes(x=median), shape=21, size=1.6,
    fill="white", colour="#777777", stroke=0.4
  ) +
  geom_point(
    aes(x=balanced_accuracy, colour=kind, shape=core),
    size=2.7, stroke=0.8
  ) +
  geom_text(
    aes(x=balanced_accuracy + 0.032, label=value_label),
    size=2.5, hjust=0, colour="#303030"
  ) +
  scale_colour_manual(values=model_cols, guide="none") +
  scale_shape_manual(
    values=c(`FALSE`=21, `TRUE`=16), guide="none"
  ) +
  scale_x_continuous(
    limits=c(0.19, 0.84), breaks=seq(0.2, 0.8, 0.1),
    expand=expansion(mult=c(0, 0))
  ) +
  labs(
    tag="a", x="Balanced accuracy", y=NULL
  ) +
  theme_journal +
  theme(
    axis.ticks.y=element_blank(),
    axis.line.y=element_blank()
  )

# B. Confusion matrix for the best model.
pc <- pred[
  pred$task == "management" &
    pred$kind == "combined" &
    pred$variant == "loso_core",
]
lev <- c("HD", "LD", "wild")
cm <- as.data.frame(table(
  truth=factor(pc$truth, levels=lev),
  predicted=factor(pc$predicted, levels=lev)
))
names(cm)[3] <- "n"
cm$total <- ave(cm$n, cm$truth, FUN=sum)
cm$pct <- 100 * cm$n / cm$total
cm$lab <- sprintf("%d\n(%.0f%%)", cm$n, cm$pct)

pB <- ggplot(cm, aes(x=predicted, y=truth, fill=pct)) +
  geom_tile(colour="white", linewidth=0.8) +
  geom_text(
    aes(label=lab),
    colour=ifelse(cm$pct >= 60, "white", "#202020"),
    size=2.8, lineheight=0.9
  ) +
  scale_fill_gradient(
    low="#EFF4F8", high="#17649A", limits=c(0, 100)
  ) +
  scale_y_discrete(
    limits=rev(lev),
    labels=function(x) ifelse(x == "wild", "Wild", x)
  ) +
  scale_x_discrete(
    labels=function(x) ifelse(x == "wild", "Wild", x)
  ) +
  coord_equal() +
  labs(
    tag="b", x="Predicted", y="Observed"
  ) +
  theme_journal +
  theme(
    legend.position="none",
    axis.ticks=element_blank()
  )

# C. Paired block-wise comparison of bacterial and combined core models.
pb <- pred[
  pred$task == "management" &
    pred$variant == "loso_core" &
    pred$kind %in% c("bacteria", "combined"),
]
pb$correct <- pb$truth == pb$predicted
pb <- aggregate(correct ~ kind + fold, pb, mean)
pb$kind <- factor(
  pb$kind, levels=c("bacteria", "combined"),
  labels=c("Bacteria core", "Combined core")
)
pb$fold <- factor(pb$fold, levels=rev(c("AO", "AS", "BO", "BS")))

pC <- ggplot(pb, aes(x=correct, y=fold, group=fold)) +
  geom_vline(
    xintercept=1/3, linetype="dashed",
    colour="#8A8A8A", linewidth=0.35
  ) +
  geom_line(colour="#B8B8B8", linewidth=0.8) +
  geom_point(
    aes(colour=kind, shape=kind),
    size=2.8, stroke=0.7
  ) +
  scale_colour_manual(
    values=c("Bacteria core"="#0072B2", "Combined core"="#202020"),
    name=NULL
  ) +
  scale_shape_manual(values=c(16, 17), name=NULL) +
  scale_x_continuous(
    limits=c(0.25, 1.03), breaks=seq(0.4, 1, 0.2),
    expand=expansion(mult=c(0, 0))
  ) +
  labs(
    tag="c", x="Held-out block accuracy", y="Field block"
  ) +
  theme_journal +
  theme(
    legend.position=c(0.30, 0.18),
    legend.background=element_rect(fill="white", colour=NA),
    axis.ticks.y=element_blank(),
    axis.line.y=element_blank()
  )

# D. Raw observations plus block-adjusted marginal means and 95% CIs.
counts <- readRDS(file.path(
  root, "taxa_tables", "bacteria_genus_counts_primary211.rds"
))
meta <- read.delim(
  "PRJEB110492/metadata/PRJEB110492_master_metadata_provisional.tsv",
  check.names=FALSE, na.strings=c("", "NA")
)
meta$sample_uid <- paste(
  "PRJEB110492", meta$biological_sample, sep="__"
)
meta <- meta[
  meta$design_family == "transplanted_field_system" &
    meta$treatment_code_raw %in% lev &
    meta$sample_uid %in% rownames(counts),
]
meta$treatment <- factor(meta$treatment_code_raw, levels=lev)
meta$block <- factor(meta$block_code_raw, levels=c("AO", "AS", "BO", "BS"))
counts <- counts[meta$sample_uid, , drop=FALSE]

core <- read.delim(file.path(
  root, "core_microbiome", "bacteria_core_audit.tsv"
), check.names=FALSE)
prev <- grep("^prevalence_", names(core), value=TRUE)
prev <- setdiff(prev, "prevalence_PRJEB110492")
keep <- apply(
  core[, prev, drop=FALSE], 1,
  function(z) all(!is.na(z) & z >= 0.50)
)
clr_features <- intersect(core$feature[keep], colnames(counts))
clr <- log(counts[, clr_features, drop=FALSE] + 0.5)
clr <- sweep(clr, 1, rowMeans(clr), "-")

stable <- validation[
  validation$kind == "bacteria" &
    validation$selected_blocks == 4,
]
stable <- stable[!duplicated(stable$feature),]
stable <- stable[stable$feature %in% colnames(clr),]

taxon_name <- function(feature) {
  if (grepl(";f__67-14;", feature))
    return("Family 67-14 lineage")
  sub(".*;g__", "", feature)
}

raw_list <- list()
emm_list <- list()
for (i in seq_len(nrow(stable))) {
  feature <- stable$feature[i]
  taxon <- taxon_name(feature)
  dat <- data.frame(
    clr=clr[, feature],
    treatment=meta$treatment,
    block=meta$block
  )
  fit <- lm(clr ~ treatment + block, data=dat)
  raw_list[[i]] <- transform(
    dat, taxon=taxon, sample_uid=meta$sample_uid
  )
  nd <- expand.grid(
    treatment=factor(lev, levels=lev),
    block=factor(levels(meta$block), levels=levels(meta$block))
  )
  mm <- model.matrix(delete.response(terms(fit)), nd)
  emm <- do.call(rbind, lapply(lev, function(tr) {
    w <- colMeans(mm[nd$treatment == tr, , drop=FALSE])
    estimate <- sum(w * coef(fit))
    se <- sqrt(as.numeric(t(w) %*% vcov(fit) %*% w))
    data.frame(
      treatment=tr, estimate=estimate,
      low=estimate - qt(0.975, df.residual(fit)) * se,
      high=estimate + qt(0.975, df.residual(fit)) * se
    )
  }))
  letters <- if (taxon == "Conexibacter") {
    c("a", "b", "c")
  } else if (taxon == "Candidatus Solibacter") {
    c("a", "b", "b")
  } else {
    c("a", "a", "b")
  }
  emm$letter <- letters
  emm$taxon <- taxon
  emm_list[[i]] <- emm
}
raw <- do.call(rbind, raw_list)
emm <- do.call(rbind, emm_list)
raw$treatment <- factor(raw$treatment, levels=lev)
emm$treatment <- factor(emm$treatment, levels=lev)

taxon_order <- c(
  "Conexibacter", "Candidatus Solibacter", "Family 67-14 lineage"
)
raw$taxon <- factor(raw$taxon, levels=taxon_order)
emm$taxon <- factor(emm$taxon, levels=taxon_order)
label_map <- setNames(vapply(taxon_order, function(tx) {
  feature <- stable$feature[vapply(
    stable$feature, function(f) taxon_name(f) == tx, logical(1)
  )][1]
  v <- validation[
    validation$feature == feature &
      !duplicated(validation$feature),
  ]
  sprintf(
    "%s\npartial eta2 = %.2f; q = %.2g",
    tx, v$partial_eta2[1], v$treatment_q[1]
  )
}, character(1)), taxon_order)

letter_offset <- ave(
  emm$high, emm$taxon,
  FUN=function(z) rep(0.08 * diff(range(c(raw$clr, z))), length(z))
)
emm$letter_y <- emm$high + letter_offset

pD <- ggplot(raw, aes(x=treatment, y=clr)) +
  geom_point(
    aes(shape=block, fill=treatment),
    position=position_jitter(width=0.08, height=0),
    size=1.35, alpha=0.42, colour="#555555", stroke=0.35
  ) +
  geom_errorbar(
    data=emm,
    aes(x=treatment, ymin=low, ymax=high, colour=treatment),
    width=0.07, linewidth=0.55, inherit.aes=FALSE
  ) +
  geom_point(
    data=emm,
    aes(x=treatment, y=estimate, colour=treatment),
    size=2.7, inherit.aes=FALSE
  ) +
  geom_text(
    data=emm,
    aes(x=treatment, y=letter_y, label=letter),
    size=2.8, fontface="bold", inherit.aes=FALSE
  ) +
  facet_wrap(
    ~taxon, nrow=1, scales="free_y",
    labeller=as_labeller(label_map)
  ) +
  scale_fill_manual(values=treatment_cols, guide="none") +
  scale_colour_manual(values=treatment_cols, guide="none") +
  scale_shape_manual(values=c(21, 22, 23, 24), name="Block") +
  scale_x_discrete(
    labels=function(x) ifelse(x == "wild", "Wild", x)
  ) +
  labs(
    tag="d", x="Management treatment", y="CLR abundance"
  ) +
  theme_journal +
  theme(
    legend.position="right",
    legend.margin=margin(0, 0, 0, 0)
  )

save_one <- function(plot, name, width, height) {
  ggsave(
    file.path(outdir, paste0(name, ".pdf")),
    plot=plot, width=width, height=height,
    units="in", device=cairo_pdf
  )
  ggsave(
    file.path(outdir, paste0(name, ".png")),
    plot=plot, width=width, height=height,
    units="in", dpi=300, bg="white"
  )
}

save_one(pA, "Panel_A_performance", 5.2, 3.2)
save_one(pB, "Panel_B_confusion", 3.2, 3.2)
save_one(pC, "Panel_C_blocks", 3.6, 3.2)
save_one(pD, "Panel_D_adjusted_taxa", 7.0, 3.4)

draw_all <- function(type, file) {
  if (type == "pdf") {
    cairo_pdf(file, width=183/25.4, height=150/25.4)
  } else if (type == "png") {
    png(file, width=183/25.4*300, height=150/25.4*300, res=300, bg="white")
  } else {
    tiff(
      file, width=183/25.4, height=150/25.4, units="in",
      res=600, compression="lzw", bg="white"
    )
  }
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(
    2, 3,
    widths=unit(c(1.0, 1.0, 0.82), "null"),
    heights=unit(c(0.96, 1.04), "null")
  )))
  print(pA, vp=viewport(layout.pos.row=1, layout.pos.col=1:2))
  print(pB, vp=viewport(layout.pos.row=1, layout.pos.col=3))
  print(pC, vp=viewport(layout.pos.row=2, layout.pos.col=1))
  print(pD, vp=viewport(layout.pos.row=2, layout.pos.col=2:3))
  dev.off()
}

draw_all("pdf", file.path(outdir, "Figure_management_transferability_v2.pdf"))
draw_all("png", file.path(outdir, "Figure_management_transferability_v2.png"))
draw_all(
  "tiff",
  file.path(outdir, "Figure_management_transferability_v2_600dpi.tiff")
)

cat("Publication figure written to:", outdir, "\n")
cat("MANAGEMENT FIGURE V2 DONE\n")
