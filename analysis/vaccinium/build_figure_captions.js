const fs = require("fs");
const path = require("path");
const {
  Document,
  Packer,
  Paragraph,
  TextRun,
  HeadingLevel,
  AlignmentType,
  Footer,
  PageNumber,
} = require("docx");

const outputDir =
  "C:\\Users\\admin\\Desktop\\figures_v2\\publication_vector_title_free";

const captions = [
  {
    number: "Figure 1.",
    text:
      "Cohort construction and sequencing depth. (A) Available paired biological samples, samples retained in the primary analysis using UNITE phylum-confident fungal reads, and samples retained in the strict EUKARYOME fungal sensitivity set. (B) High-confidence bacterial and fungal read depths among primary-analysis samples. Points represent individual samples; boxes show medians and interquartile ranges; the y-axis is log10-scaled.",
  },
  {
    number: "Figure 2.",
    text:
      "Controlled pH responses across compartments. Principal component ordinations of centered log-ratio genus profiles for bacterial (A-C) and fungal (D-F) communities in bulk soil, rhizosphere, and endosphere samples from PRJNA1156347. Filled points are biological samples. Open labeled points are treatment centroids (L, low; O, optimum; H, high). Lines connect ordered treatment centroids to summarize the pH sequence and do not represent longitudinal sample trajectories. Axis labels report explained variance.",
  },
  {
    number: "Figure 3.",
    text:
      "Community-, taxon-, and signature-level pH transferability. (A) Aitchison PERMANOVA effect sizes for categorical pH treatment within compartments of the controlled study (circles) and for the quadratic pH term in the independent field rhizosphere cohort (triangles). Controlled and field R2 values arise from different model structures and should be compared as evidence strength, not as identical estimands. (B) Nominal validation fractions for prespecified pH-associated bacterial and fungal genera. Error bars are exact binomial 95% confidence intervals; the dashed line is the 5% null expectation. (C) Locked bacterial and fungal pH signature scores projected from the controlled rhizosphere discovery set into the independent field rhizosphere cohort. Curves and ribbons show quadratic fits and 95% confidence intervals; facet labels report field R2 and within-field pH-label permutation P values.",
  },
  {
    number: "Figure 4.",
    text:
      "Cross-block transferability of management-associated bacterial signatures. (A) Leave-one-block-out balanced accuracy for full and cross-study core feature sets. Colored points are observed values; gray intervals show the central 95% of 999 within-block label-permutation null distributions; the dashed line marks three-class chance performance. (B) Confusion matrix for the combined bacterial-fungal core classifier. Percentages are calculated within observed treatment. (C) Held-out-block accuracy for bacterial and combined core models. (D) Raw centered log-ratio observations and block-adjusted marginal means with 95% confidence intervals for three bacterial predictors selected with consistent signs in all four blocks. Letters summarize false-discovery-rate-adjusted treatment contrasts. HD, high density; LD, low density.",
  },
  {
    number: "Figure 5.",
    text:
      "Cross-study core genera. Study-specific prevalence of (A) the 15 most prevalent members of the 37-genus bacterial core and (B) all four fungal genera meeting the primary core criterion. A core genus was present in at least 50% of samples in every primary cohort. The asterisk marks Penicillium, the only fungal genus that also met the criterion in the strict EUKARYOME sensitivity analysis. Cell values are prevalence proportions.",
  },
  {
    number: "Figure 6.",
    text:
      "Bacterial-fungal community concordance and sensitivity to fungal filtering. Spearman Mantel correlations between bacterial and fungal Aitchison distance matrices in (A) the 211-sample primary analysis and (B) the 175-sample strict fungal sensitivity set. Raw correlations are blue; correlations after residualization for available study-specific design variables are orange. Filled symbols indicate P < 0.05 and open symbols indicate P >= 0.05. Studies were tested independently; points are not connected.",
  },
  {
    number: "Supplementary Figure S1.",
    text:
      "Management-associated community shifts in PRJEB110492. Aitchison ordinations for (A) bacteria and (B) fungi. Colors denote high-density (HD), low-density (LD), and wild treatments; point symbols denote four experimental blocks. Open labeled points are treatment centroids. (C) Pairwise treatment effect sizes from block-aware Aitchison PERMANOVA. All displayed contrasts had Benjamini-Hochberg-adjusted q <= 0.005.",
  },
];

const children = [
  new Paragraph({
    heading: HeadingLevel.TITLE,
    alignment: AlignmentType.CENTER,
    spacing: { after: 360 },
    children: [new TextRun({ text: "Figure Captions", bold: true })],
  }),
];

for (const caption of captions) {
  children.push(
    new Paragraph({
      keepLines: true,
      spacing: { after: 240, line: 300 },
      children: [
        new TextRun({ text: `${caption.number} `, bold: true }),
        new TextRun({ text: caption.text }),
      ],
    })
  );
}

const doc = new Document({
  styles: {
    default: {
      document: {
        run: { font: "Times New Roman", size: 22 },
        paragraph: { spacing: { line: 300 } },
      },
    },
  },
  sections: [
    {
      properties: {
        page: {
          size: { width: 12240, height: 15840 },
          margin: {
            top: 1080,
            right: 1080,
            bottom: 1080,
            left: 1080,
          },
        },
      },
      footers: {
        default: new Footer({
          children: [
            new Paragraph({
              alignment: AlignmentType.CENTER,
              children: [
                new TextRun("Page "),
                new TextRun({ children: [PageNumber.CURRENT] }),
              ],
            }),
          ],
        }),
      },
      children,
    },
  ],
});

Packer.toBuffer(doc).then((buffer) => {
  const target = path.join(outputDir, "Figure_Captions.docx");
  fs.writeFileSync(target, buffer);
  console.log(`Wrote ${target}`);
});
