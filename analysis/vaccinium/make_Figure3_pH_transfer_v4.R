suppressPackageStartupMessages(library(ggplot2))
library(grid)

root <- "global_harmonized"
outdir <- file.path(root, "figures_nature")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

kingdom_cols <- c(Bacteria="#0072B2", Fungi="#D55E00")

theme_pub <- function() {
  theme_classic(base_size=7, base_family="Arial") +
    theme(
      text=element_text(family="Arial", color="black"),
      axis.title=element_text(size=7, color="black"),
      axis.text=element_text(size=7, color="black"),
      legend.title=element_text(size=6),
      legend.text=element_text(size=6),
      plot.tag=element_text(size=8, face="bold"),
      plot.margin=margin(5, 7, 5, 5),
      strip.background=element_rect(
        fill="#F3F3F3", colour="#B8B8B8", linewidth=0.3
      ),
      strip.text=element_text(size=7.5, face="bold")
    )
}

fmt_p <- function(p) {
  ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
}

# A. Community-level pH effect sizes.
ph <- read.delim(
  file.path(root, "pH_analysis", "stratified_pH_summary.tsv"),
  stringsAsFactors=FALSE
)
controlled <- ph[
  ph$study == "PRJNA1156347" &
    ph$metric == "Aitchison" &
    ph$model == "categorical_pH",
  , drop=FALSE
]
controlled$context <- controlled$subset
controlled$evidence <- "Controlled experiment"

field <- ph[
  ph$study == "PRJEB98254" &
    ph$metric == "Aitchison" &
    ph$model == "quadratic_pH" &
    ph$term == "pH_z2",
  , drop=FALSE
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

pA <- ggplot(
  ph_eff,
  aes(R2, context, color=kingdom, shape=evidence)
) +
  geom_point(size=2.7, position=position_dodge(width=0.45)) +
  geom_text(
    aes(label=paste0("P ", fmt_p(p))),
    position=position_dodge(width=0.45),
    hjust=-0.12, size=2.2, show.legend=FALSE
  ) +
  scale_color_manual(values=kingdom_cols) +
  scale_shape_manual(values=c(16, 17)) +
  scale_x_continuous(
    limits=c(0, max(ph_eff$R2) * 1.25),
    expand=expansion(mult=c(0, 0))
  ) +
  labs(
    x=expression(Aitchison~R^2),
    y=NULL, color=NULL, shape=NULL, tag="a"
  ) +
  theme_pub() +
  guides(shape="none") +
  theme(legend.position="bottom")

# B. Candidate-level nominal validation.
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
validation$label <- paste0(
  validation$validated, "/", validation$tested
)
ci <- t(vapply(seq_len(nrow(validation)), function(i) {
  binom.test(validation$validated[i], validation$tested[i])$conf.int
}, numeric(2)))
validation$lower <- ci[, 1]
validation$upper <- ci[, 2]

pB <- ggplot(validation, aes(kingdom, rate, color=kingdom)) +
  geom_hline(
    yintercept=0.05, linetype=2,
    linewidth=0.4, color="#6B7280"
  ) +
  geom_errorbar(
    aes(ymin=lower, ymax=upper),
    width=0.12, linewidth=0.55
  ) +
  geom_point(size=3) +
  geom_text(
    aes(label=label), vjust=-1,
    size=2.4, color="black"
  ) +
  scale_color_manual(values=kingdom_cols) +
  scale_y_continuous(
    labels=scales::percent_format(accuracy=1),
    limits=c(0, max(validation$upper) * 1.08)
  ) +
  labs(
    x=NULL, y="Nominal validation rate (95% CI)", tag="b"
  ) +
  annotate(
    "text", x=1.5, y=0.052,
    label="5% null expectation",
    vjust=-0.4, size=2.2, color="#6B7280"
  ) +
  theme_pub() +
  theme(legend.position="none")

# C. Locked controlled-to-field pH signature transfer.
sig_dir <- file.path(root, "enhancement", "pH_signature_transfer")
scores <- read.delim(file.path(sig_dir, "pH_signature_scores.tsv"))
scores <- scores[scores$study == "PRJEB98254", , drop=FALSE]
scores$kingdom <- ifelse(
  scores$marker == "bacteria", "Bacteria", "Fungi"
)

curves <- list()
for (marker in c("bacteria", "fungi")) {
  z <- scores[scores$marker == marker, , drop=FALSE]
  fit <- lm(signature_score ~ pH + I(pH^2), data=z)
  nd <- data.frame(pH=seq(min(z$pH), max(z$pH), length.out=200))
  pr <- predict(fit, newdata=nd, interval="confidence")
  curves[[marker]] <- data.frame(
    marker=marker,
    kingdom=ifelse(marker == "bacteria", "Bacteria", "Fungi"),
    pH=nd$pH,
    fit=pr[, "fit"],
    low=pr[, "lwr"],
    high=pr[, "upr"]
  )
}
curves <- do.call(rbind, curves)

pC <- ggplot(
  scores,
  aes(pH, signature_score, color=kingdom)
) +
  geom_ribbon(
    data=curves,
    aes(x=pH, ymin=low, ymax=high, fill=kingdom),
    alpha=0.16, colour=NA, inherit.aes=FALSE
  ) +
  geom_line(
    data=curves,
    aes(x=pH, y=fit, color=kingdom),
    linewidth=0.7, inherit.aes=FALSE
  ) +
  geom_point(size=1.8, alpha=0.82) +
  facet_wrap(
    ~kingdom, nrow=1,
    scales="free_y",
    labeller=label_value
  ) +
  scale_color_manual(values=kingdom_cols, guide="none") +
  scale_fill_manual(values=kingdom_cols, guide="none") +
  labs(
    x="Field rhizosphere pH",
    y="Transferred pH signature score",
    tag="c"
  ) +
  theme_pub()

draw <- function() {
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(
    2, 2,
    widths=unit(c(1.25, 0.75), "null"),
    heights=unit(c(0.95, 1.05), "null")
  )))
  print(pA, vp=viewport(layout.pos.row=1, layout.pos.col=1))
  print(pB, vp=viewport(layout.pos.row=1, layout.pos.col=2))
  print(pC, vp=viewport(layout.pos.row=2, layout.pos.col=1:2))
}

cairo_pdf(
  file.path(outdir, "Figure3_pH_external_signature_transfer.pdf"),
  width=183/25.4, height=138/25.4
)
draw()
dev.off()

png(
  file.path(outdir, "Figure3_pH_external_signature_transfer.png"),
  width=183/25.4*300, height=138/25.4*300, res=300
)
draw()
dev.off()

tiff(
  file.path(outdir, "Figure3_pH_external_signature_transfer_600dpi.tiff"),
  width=183/25.4, height=138/25.4,
  units="in", res=600, compression="lzw"
)
draw()
dev.off()

cat("FIGURE 3 PH SIGNATURE TRANSFER DONE\n")
