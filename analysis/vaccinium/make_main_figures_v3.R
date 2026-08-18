library(ggplot2)
library(grid)

root <- "global_harmonized"
out <- file.path(root, "figures_nature")
dir.create(out, recursive=TRUE, showWarnings=FALSE)

study_cols <- c(
  PRJEB35843="#0072B2",
  PRJNA1156347="#D55E00",
  PRJEB98254="#009E73",
  PRJEB110492="#CC79A7"
)
kingdom_cols <- c(Bacteria="#0072B2", Fungi="#D55E00")
ph_cols <- c(Low="#3B82F6", Optimum="#009E73", High="#D55E00")

theme_pub <- function() {
  theme_classic(base_size=7, base_family="Arial") +
    theme(
      text=element_text(family="Arial", color="black"),
      axis.title=element_text(size=7, color="black"),
      axis.text=element_text(size=7, color="black"),
      legend.title=element_text(size=6),
      legend.text=element_text(size=6),
      plot.title=element_blank(),
      plot.subtitle=element_blank(),
      plot.tag=element_text(size=8, face="bold"),
      plot.margin=margin(5, 7, 5, 5)
    )
}

save_plot <- function(p, name, width_mm=90, height_mm=75) {
  ggsave(
    file.path(out, paste0(name, ".pdf")),
    plot=p, width=width_mm, height=height_mm, units="mm", device=cairo_pdf
  )
  ggsave(
    file.path(out, paste0(name, ".png")),
    plot=p, width=width_mm, height=height_mm, units="mm", dpi=300
  )
}

save_grid <- function(plots, name, nrow, ncol, width_mm=183,
                      height_mm=120, widths=NULL, heights=NULL) {
  draw <- function() {
    if (is.null(widths))
      widths <- unit(rep(1, ncol), "null")
    if (is.null(heights))
      heights <- unit(rep(1, nrow), "null")
    grid.newpage()
    lay <- grid.layout(nrow, ncol, widths=widths, heights=heights)
    pushViewport(viewport(layout=lay))
    for (i in seq_along(plots)) {
      r <- ceiling(i/ncol)
      c <- i - (r-1)*ncol
      print(plots[[i]], vp=viewport(layout.pos.row=r, layout.pos.col=c))
    }
  }
  cairo_pdf(
    file.path(out, paste0(name, ".pdf")),
    width=width_mm/25.4, height=height_mm/25.4
  )
  draw()
  dev.off()
  png(
    file.path(out, paste0(name, ".png")),
    width=width_mm/25.4*300, height=height_mm/25.4*300, res=300
  )
  draw()
  dev.off()
}

fmt_p <- function(p) {
  ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
}

# Figure 1: cohort design and sequencing depth
meta_all <- read.delim(
  file.path(root, "global_metadata_UNITE_phylum_QC.tsv"),
  check.names=FALSE, stringsAsFactors=FALSE
)

study_order <- c("PRJEB110492", "PRJEB35843", "PRJEB98254", "PRJNA1156347")
qc <- do.call(rbind, lapply(study_order, function(s) {
  z <- meta_all[meta_all$study == s, , drop=FALSE]
  data.frame(
    study=s,
    tier=c("Available pairs", "Primary: UNITE phylum", "Strict: EUKARYOME"),
    n=c(
      nrow(z),
      sum(z$primary_UNITE_phylum_QC %in% c(TRUE, "TRUE", 1, "1")),
      sum(z$strict_pair_QC %in% c(TRUE, "TRUE", 1, "1"))
    )
  )
}))
qc$study <- factor(qc$study, levels=rev(study_order))
qc$tier <- factor(
  qc$tier,
  levels=c("Available pairs", "Primary: UNITE phylum", "Strict: EUKARYOME")
)

p1a <- ggplot(qc, aes(n, study, fill=tier)) +
  geom_col(position=position_dodge(width=0.72), width=0.65) +
  geom_text(
    aes(label=n),
    position=position_dodge(width=0.72), hjust=-0.25, size=2.4
  ) +
  scale_fill_manual(
    values=c("#9CA3AF", "#2563EB", "#E69F00"),
    labels=c("Available", "Primary", "Strict")
  ) +
  scale_x_continuous(expand=expansion(mult=c(0, 0.12))) +
  labs(x="Paired biological samples", y=NULL, fill=NULL, tag="b") +
  theme_pub() +
  theme(legend.position="bottom")

primary <- meta_all[
  meta_all$primary_UNITE_phylum_QC %in% c(TRUE, "TRUE", 1, "1"), ,
  drop=FALSE
]
depth <- rbind(
  data.frame(
    study=primary$study, kingdom="Bacteria",
    reads=primary$bacterial_reads
  ),
  data.frame(
    study=primary$study, kingdom="Fungi",
    reads=primary$fungal_reads_UNITE_phylum
  )
)
depth$study <- factor(depth$study, levels=study_order)

p1b <- ggplot(depth, aes(study, reads, fill=kingdom)) +
  geom_boxplot(
    width=0.58, outlier.shape=NA, linewidth=0.35,
    position=position_dodge(width=0.65)
  ) +
  geom_point(
    aes(color=kingdom),
    alpha=0.22, size=0.65,
    position=position_jitterdodge(
      jitter.width=0.14, dodge.width=0.65
    )
  ) +
  scale_y_log10(
    labels=scales::label_number(scale_cut=scales::cut_si(""))
  ) +
  scale_fill_manual(values=kingdom_cols) +
  scale_color_manual(values=kingdom_cols) +
  labs(x=NULL, y="High-confidence reads per sample", fill=NULL, color=NULL, tag="c") +
  theme_pub() +
  theme(
    axis.text.x=element_text(angle=25, hjust=1),
    legend.position="bottom"
  )

# A compact study-aware entry point for the paper.  This is deliberately a
# workflow panel rather than a decorative title: it makes the discovery /
# external-validation split explicit before the QC plots are read.
p1workflow <- ggplot() +
  annotate(
    "text", x=0.48, y=1.40, label="a", family="Arial",
    fontface="bold", size=3.0, hjust=0
  ) +
  annotate(
    "label", x=1.0, y=0.92,
    label="7 paired cohorts\n16S + ITS\n265 primary-QC pairs",
    family="Arial", fontface="bold", size=2.6, lineheight=0.92,
    label.size=0.35, color="#111827", fill="#F3F4F6"
  ) +
  annotate(
    "segment", x=1.62, xend=2.28, y=0.92, yend=0.92,
    linewidth=0.5, color="#6B7280",
    arrow=arrow(length=unit(2.0, "mm"), type="closed")
  ) +
  annotate(
    "label", x=2.95, y=0.92,
    label="Discovery\n4 cohorts\ncore definition",
    family="Arial", fontface="bold", size=2.6, lineheight=0.92,
    label.size=0.35, color="#111827", fill="#E8F1FA"
  ) +
  annotate(
    "segment", x=3.57, xend=4.23, y=0.92, yend=0.92,
    linewidth=0.5, color="#6B7280",
    arrow=arrow(length=unit(2.0, "mm"), type="closed")
  ) +
  annotate(
    "label", x=4.90, y=0.92,
    label="Study-held-out\nbenchmark\nLOSO transfer",
    family="Arial", fontface="bold", size=2.6, lineheight=0.92,
    label.size=0.35, color="#111827", fill="#EAF5F1"
  ) +
  annotate(
    "segment", x=5.52, xend=6.18, y=0.92, yend=0.92,
    linewidth=0.5, color="#6B7280",
    arrow=arrow(length=unit(2.0, "mm"), type="closed")
  ) +
  annotate(
    "label", x=6.85, y=0.92,
    label="External validation\n3 cohorts\nindependent test",
    family="Arial", fontface="bold", size=2.6, lineheight=0.92,
    label.size=0.35, color="#111827", fill="#FFF4E5"
  ) +
  scale_x_continuous(limits=c(0.40, 7.60), expand=c(0, 0)) +
  scale_y_continuous(limits=c(0.48, 1.46), expand=c(0, 0)) +
  theme_void(base_family="Arial") +
  theme(plot.margin=margin(2, 5, 1, 5))

save_figure1 <- function(p_top, p_left, p_right, name,
                         width_mm=183, height_mm=103) {
  draw <- function() {
    grid.newpage()
    lay <- grid.layout(
      nrow=2, ncol=2,
      heights=unit(c(0.58, 1.42), "null"),
      widths=unit(c(1, 1), "null")
    )
    pushViewport(viewport(layout=lay))
    print(p_top,   vp=viewport(layout.pos.row=1, layout.pos.col=1:2))
    print(p_left,  vp=viewport(layout.pos.row=2, layout.pos.col=1))
    print(p_right, vp=viewport(layout.pos.row=2, layout.pos.col=2))
    upViewport()
  }
  cairo_pdf(
    file.path(out, paste0(name, ".pdf")),
    width=width_mm/25.4, height=height_mm/25.4
  )
  draw()
  dev.off()
  png(
    file.path(out, paste0(name, ".png")),
    width=width_mm/25.4*300, height=height_mm/25.4*300, res=300
  )
  draw()
  dev.off()
}

save_figure1(
  p1workflow, p1a, p1b, "Figure1_study_design_QC_nature",
  width_mm=183, height_mm=103
)

# Figure 2: controlled pH response shown within each compartment
meta_primary <- read.delim(
  file.path(root, "global_metadata_primary_211.tsv"),
  check.names=FALSE, stringsAsFactors=FALSE
)

make_ph_pca <- function(marker, compartment, tag) {
  clr <- readRDS(file.path(
    root, "analysis_ready", paste0(marker, "_clr.rds")
  ))
  md <- meta_primary[match(rownames(clr), meta_primary$sample_uid), ]
  keep <- md$study == "PRJNA1156347" & md$compartment == compartment
  z <- clr[keep, , drop=FALSE]
  m <- md[keep, , drop=FALSE]
  fit <- prcomp(z, center=TRUE, scale.=FALSE)
  vv <- 100 * fit$sdev^2 / sum(fit$sdev^2)
  d <- data.frame(
    PC1=fit$x[, 1], PC2=fit$x[, 2],
    pH_treatment=factor(
      m$pH_treatment, levels=c("Low", "Optimum", "High")
    )
  )
  cen <- aggregate(cbind(PC1, PC2) ~ pH_treatment, d, mean)
  ggplot(d, aes(PC1, PC2, color=pH_treatment)) +
    geom_point(size=1.9, alpha=0.78) +
    geom_path(
      data=cen, aes(group=1), color="#374151", linewidth=0.45,
      arrow=arrow(length=unit(1.4, "mm"), type="closed")
    ) +
    geom_point(
      data=cen, size=3.2, shape=21, fill="white", stroke=0.9
    ) +
    geom_text(
      data=cen, aes(label=substr(pH_treatment, 1, 1)),
      size=2.1, fontface="bold",
      show.legend=FALSE
    ) +
    scale_color_manual(values=ph_cols) +
    labs(
      x=sprintf("PC1 (%.1f%%)", vv[1]),
      y=sprintf("PC2 (%.1f%%)", vv[2]),
      color="pH treatment", tag=tag
    ) +
    theme_pub() +
    theme(legend.position="none")
}

ph_plots <- list()
tags <- letters[1:6]
k <- 0
for (marker in c("bacteria_genus", "fungi_genus")) {
  for (comp in c("Bulk", "Rhizosphere", "Endosphere")) {
    k <- k + 1
    ph_plots[[k]] <- make_ph_pca(marker, comp, tags[k])
  }
}
save_grid(
  ph_plots, "Figure2_controlled_pH_ordination_nature",
  nrow=2, ncol=3, width_mm=183, height_mm=120
)

# Figure 3: pH effect sizes and independent validation
ph <- read.delim(
  file.path(root, "pH_analysis", "stratified_pH_summary.tsv"),
  stringsAsFactors=FALSE
)
controlled <- ph[
  ph$study == "PRJNA1156347" &
    ph$metric == "Aitchison" &
    ph$model == "categorical_pH", ,
  drop=FALSE
]
controlled$context <- controlled$subset
controlled$evidence <- "Controlled experiment"

field <- ph[
  ph$study == "PRJEB98254" &
    ph$metric == "Aitchison" &
    ph$model == "quadratic_pH" &
    ph$term == "pH_z2", ,
  drop=FALSE
]
field$context <- "Field rhizosphere"
field$evidence <- "Independent field cohort"

ph_eff <- rbind(
  controlled[, c("marker", "context", "evidence", "R2", "p")],
  field[, c("marker", "context", "evidence", "R2", "p")]
)
ph_eff$kingdom <- ifelse(
  ph_eff$marker == "bacteria_genus", "Bacteria", "Fungi"
)
ph_eff$context <- factor(
  ph_eff$context,
  levels=rev(c("Bulk", "Rhizosphere", "Endosphere", "Field rhizosphere"))
)

p3a <- ggplot(ph_eff, aes(R2, context, color=kingdom, shape=evidence)) +
  geom_point(size=2.8, position=position_dodge(width=0.45)) +
  geom_text(
    aes(label=paste0("P ", fmt_p(p))),
    position=position_dodge(width=0.45),
    hjust=-0.15, size=2.2, show.legend=FALSE
  ) +
  scale_color_manual(values=kingdom_cols) +
  scale_shape_manual(values=c(
    "Controlled experiment"=16,
    "Independent field cohort"=17
  )) +
  scale_x_continuous(expand=expansion(mult=c(0.02, 0.28))) +
  labs(
    x=expression(Aitchison~R^2), y=NULL,
    color=NULL, shape=NULL, tag="a"
  ) +
  theme_pub() +
  guides(shape="none") +
  theme(legend.position="bottom")

validation <- do.call(rbind, lapply(
  c(Bacteria="bacteria_genus", Fungi="fungi_genus"),
  function(marker) {
    x <- read.delim(file.path(
      root, "pH_analysis", paste0(marker, "_pH_taxa_validation.tsv")
    ))
    data.frame(
      kingdom=ifelse(marker == "bacteria_genus", "Bacteria", "Fungi"),
      tested=nrow(x),
      validated=sum(x$overall_p < 0.05)
    )
  }
))
validation$rate <- validation$validated / validation$tested
validation$expected <- 0.05
validation$label <- paste0(
  validation$validated, "/", validation$tested
)
ci <- t(vapply(seq_len(nrow(validation)), function(i) {
  binom.test(validation$validated[i], validation$tested[i])$conf.int
}, numeric(2)))
validation$lower <- ci[, 1]
validation$upper <- ci[, 2]

p3b <- ggplot(validation, aes(kingdom, rate, color=kingdom)) +
  geom_hline(
    yintercept=0.05, linetype=2, linewidth=0.45, color="#6B7280"
  ) +
  geom_errorbar(
    aes(ymin=lower, ymax=upper), width=0.12, linewidth=0.6
  ) +
  geom_point(size=3.2) +
  geom_text(aes(label=label), vjust=-1.0, size=2.5, color="black") +
  scale_color_manual(values=kingdom_cols) +
  scale_y_continuous(
    labels=scales::percent_format(accuracy=1),
    limits=c(0, max(validation$upper) * 1.08)
  ) +
  labs(
    x=NULL, y="Nominal validation rate (95% CI)",
    tag="b"
  ) +
  annotate(
    "text", x=1.5, y=0.052, label="5% null expectation",
    vjust=-0.4, size=2.3, color="#6B7280"
  ) +
  theme_pub() +
  theme(legend.position="none")

save_grid(
  list(p3a, p3b), "Figure3_pH_effect_validation_nature",
  nrow=1, ncol=2, width_mm=183, height_mm=82
)

# Figure 4: management ordination and pairwise effects
management_meta <- read.delim(
  "PRJEB110492/metadata/PRJEB110492_master_metadata_provisional.tsv",
  stringsAsFactors=FALSE
)
management_meta$sample_uid <- paste(
  "PRJEB110492", management_meta$biological_sample, sep="__"
)
management_meta <- management_meta[
  management_meta$design_family == "transplanted_field_system" &
    management_meta$treatment_code_raw != "Bulk", ,
  drop=FALSE
]
management_meta$treatment <- factor(
  management_meta$treatment_code_raw, levels=c("HD", "LD", "wild")
)
management_meta$block <- factor(management_meta$block_code_raw)

make_management_pca <- function(marker, tag) {
  clr <- readRDS(file.path(
    root, "analysis_ready", paste0(marker, "_clr.rds")
  ))
  z <- clr[management_meta$sample_uid, , drop=FALSE]
  fit <- prcomp(z, center=TRUE, scale.=FALSE)
  vv <- 100 * fit$sdev^2 / sum(fit$sdev^2)
  d <- data.frame(
    PC1=fit$x[, 1], PC2=fit$x[, 2],
    treatment=management_meta$treatment,
    block=management_meta$block
  )
  cen <- aggregate(cbind(PC1, PC2) ~ treatment, d, mean)
  ggplot(d, aes(PC1, PC2, color=treatment, shape=block)) +
    geom_point(size=2.0, alpha=0.8) +
    geom_point(
      data=cen, aes(shape=NULL), size=4, shape=21,
      fill="white", stroke=1
    ) +
    geom_text(
      data=cen, aes(label=treatment, shape=NULL),
      size=2.3, fontface="bold", nudge_y=0.35,
      show.legend=FALSE
    ) +
    scale_color_manual(values=c(HD="#0072B2", LD="#E69F00", wild="#009E73")) +
    labs(
      x=sprintf("PC1 (%.1f%%)", vv[1]),
      y=sprintf("PC2 (%.1f%%)", vv[2]),
      color="Treatment code", shape="Block", tag=tolower(tag)
    ) +
    guides(color="none", shape="none") +
    theme_pub() +
    theme(legend.position="none")
}

pairwise <- read.delim(
  file.path(root, "management_analysis", "pairwise_treatments.tsv")
)
pairwise <- pairwise[pairwise$metric == "Aitchison", , drop=FALSE]
pairwise$kingdom <- ifelse(
  pairwise$marker == "bacteria", "Bacteria", "Fungi"
)
pairwise$contrast <- factor(
  pairwise$contrast,
  levels=rev(c("HD_vs_LD", "HD_vs_wild", "LD_vs_wild"))
)

p4c <- ggplot(pairwise, aes(R2, contrast, color=kingdom)) +
  geom_point(size=2.8, position=position_dodge(width=0.42)) +
  scale_color_manual(values=kingdom_cols) +
  scale_x_continuous(expand=expansion(mult=c(0.08, 0.08))) +
  labs(
    x=expression(Aitchison~R^2), y=NULL,
    color=NULL, tag="c"
  ) +
  theme_pub() +
  theme(legend.position="bottom")

save_grid(
  list(
    make_management_pca("bacteria_genus", "A"),
    make_management_pca("fungi_genus", "B"),
    p4c
  ),
  "Figure4_management_response_nature",
  nrow=1, ncol=3, width_mm=183, height_mm=72
)

# Figure 5: core microbiome, split by kingdom and annotated by sensitivity
core_long <- function(marker, top_n=NULL) {
  x <- read.delim(file.path(
    root, "core_microbiome", paste0(marker, "_core_audit.tsv")
  ))
  x <- x[x$strict_core_50, , drop=FALSE]
  x <- x[order(-x$minimum_prevalence, -x$geometric_mean_RA), ]
  if (!is.null(top_n)) x <- head(x, top_n)
  pc <- grep("^prevalence_", names(x), value=TRUE)
  d <- do.call(rbind, lapply(pc, function(v) {
    data.frame(
      genus=x$genus,
      study=sub("^prevalence_", "", v),
      prevalence=x[[v]]
    )
  }))
  d$genus <- factor(d$genus, levels=rev(x$genus))
  d$study <- factor(d$study, levels=study_order)
  d$text_col <- ifelse(d$prevalence >= 0.75, "white", "black")
  d
}

make_core_heatmap <- function(marker, top_n, tag) {
  d <- core_long(marker, top_n)
  if (marker == "fungi") {
    lev <- levels(d$genus)
    lev[lev == "Penicillium"] <- "Penicillium*"
    levels(d$genus) <- lev
  }
  ggplot(d, aes(study, genus, fill=prevalence)) +
    geom_tile(color="white", linewidth=0.4) +
    geom_text(
      aes(label=sprintf("%.2f", prevalence), color=text_col),
      size=2.0
    ) +
    scale_color_identity() +
    scale_fill_gradientn(
      colors=c("#F7FBFF", "#6BAED6", "#08306B"),
      limits=c(0.5, 1), oob=scales::squish
    ) +
    labs(
      x=NULL, y=NULL, fill="Prevalence", tag=tolower(tag)
    ) +
    theme_pub() +
    theme(
      axis.text.x=element_text(angle=25, hjust=1),
      legend.position="none"
    )
}

save_grid(
  list(
    make_core_heatmap("bacteria", 15, "a"),
    make_core_heatmap("fungi", NULL, "b")
  ),
  "Figure5_core_microbiome_nature",
  nrow=1, ncol=2, width_mm=183, height_mm=105,
  widths=unit(c(1.45, 1), "null")
)

# Figure 6: cross-kingdom concordance without connecting independent studies
main_c <- read.delim(
  file.path(root, "residual_cross_kingdom_coupling.tsv")
)
strict_c <- read.delim(
  file.path(root, "strict175_coupling_sensitivity.tsv")
)

main_long <- rbind(
  data.frame(
    study=main_c$study, adjustment="Raw",
    r=main_c$raw_mantel_r, p=main_c$raw_mantel_p
  ),
  data.frame(
    study=main_c$study, adjustment="Adjusted",
    r=main_c$residual_mantel_r, p=main_c$residual_mantel_p
  )
)
strict_long <- rbind(
  data.frame(
    study=strict_c$study, adjustment="Raw",
    r=strict_c$raw_Aitchison_mantel_r,
    p=strict_c$raw_Aitchison_mantel_p
  ),
  data.frame(
    study=strict_c$study, adjustment="Adjusted",
    r=strict_c$residual_mantel_r,
    p=strict_c$residual_mantel_p
  )
)

make_coupling_plot <- function(d, title, tag) {
  d$study <- factor(d$study, levels=rev(study_order))
  d$significant <- factor(
    ifelse(d$p < 0.05, "P < 0.05", "Not significant"),
    levels=c("P < 0.05", "Not significant")
  )
  ggplot(d, aes(r, study, color=adjustment, shape=significant)) +
    geom_vline(xintercept=0, color="#9CA3AF", linewidth=0.4) +
    geom_point(size=3.0, position=position_dodge(width=0.48)) +
    scale_color_manual(values=c(Raw="#2563EB", Adjusted="#D55E00")) +
    scale_shape_manual(values=c("P < 0.05"=16, "Not significant"=1)) +
    coord_cartesian(xlim=c(-0.08, 0.75)) +
    labs(
      x="Cross-kingdom Mantel correlation",
      y=NULL, color=NULL, shape=NULL, tag=tolower(tag)
    ) +
    theme_pub() +
    theme(legend.position="none")
}

save_grid(
  list(
    make_coupling_plot(main_long, "Primary high-confidence analysis", "a"),
    make_coupling_plot(strict_long, "Strict lower-bound sensitivity", "b")
  ),
  "Figure6_cross_kingdom_concordance_nature",
  nrow=1, ncol=2, width_mm=183, height_mm=78
)

cat("Publication-style v3 figures written to:", out, "\n")
