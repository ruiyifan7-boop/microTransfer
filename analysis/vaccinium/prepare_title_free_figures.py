from pathlib import Path
import shutil
import xml.etree.ElementTree as ET

import fitz


ROOT = Path(r"C:\Users\admin\Desktop")
SOURCE_V3 = ROOT / "Figures" / "figures_v3"
SOURCE_V4 = ROOT / "Figures" / "figures_v4"
SOURCE_NEW = ROOT / "figures_v2" / "manuscript_ready"
OUT = ROOT / "figures_v2" / "publication_vector_title_free"
OUT.mkdir(parents=True, exist_ok=True)

FIGURES = [
    (
        SOURCE_V3 / "Figure1_study_design_QC_v3.pdf",
        "Figure1_cohort_design_and_QC",
        [],
    ),
    (
        SOURCE_V3 / "Figure2_controlled_pH_ordination_v3.pdf",
        "Figure2_controlled_pH_responses",
        [],
    ),
    (
        SOURCE_V4 / "Figure3_pH_external_signature_transfer.pdf",
        "Figure3_pH_transferability",
        [],
    ),
    (
        SOURCE_V4 / "Figure4_management_transferability_999.pdf",
        "Figure4_management_transferability",
        [],
    ),
    (
        SOURCE_V3 / "Figure5_core_microbiome_v3.pdf",
        "Figure5_cross_study_core_microbiome",
        [],
    ),
    (
        SOURCE_V3 / "Figure6_cross_kingdom_concordance_v3.pdf",
        "Figure6_cross_kingdom_concordance",
        [],
    ),
    (
        SOURCE_NEW / "Supplementary_Figure_S1_management_ordination_effect_sizes.pdf",
        "FigureS1_management_ordination_effect_sizes",
        [],
    ),
]


def prepare_pdf(source: Path, target: Path, remove_text: list[str]) -> None:
    if not source.exists():
        raise FileNotFoundError(source)
    doc = fitz.open(source)
    if len(doc) != 1:
        raise ValueError(f"{source} contains {len(doc)} pages")
    page = doc[0]
    for text in remove_text:
        matches = page.search_for(text)
        if not matches:
            raise ValueError(f"Title not found in {source.name}: {text}")
        for rect in matches:
            rect.x0 -= 1
            rect.y0 -= 1
            rect.x1 += 1
            rect.y1 += 1
            page.add_redact_annot(rect, fill=(1, 1, 1))
    if remove_text:
        page.apply_redactions(images=0, graphics=0, text=0)
    doc.save(target, garbage=4, deflate=True)
    doc.close()


def pdf_to_svg(pdf_path: Path, svg_path: Path) -> None:
    doc = fitz.open(pdf_path)
    svg = doc[0].get_svg_image(text_as_path=True)
    svg_path.write_text(svg, encoding="utf-8")
    doc.close()
    ET.parse(svg_path)


for source, stem, remove_text in FIGURES:
    pdf_out = OUT / f"{stem}.pdf"
    svg_out = OUT / f"{stem}.svg"
    prepare_pdf(source, pdf_out, remove_text)
    pdf_to_svg(pdf_out, svg_out)

readme = """Publication artwork

Figures 1-3, 5-6, and S1 are single-page native vector PDF/SVG exports.
Figure 4 is generated from the locally available 4725 x 3240 PNG after the
999-permutation analysis. Its PDF and SVG preserve the high-resolution image
but are raster-embedded compatibility files, not element-editable vectors.
No overall figure titles are embedded in the artwork.
Necessary panel identifiers, axis labels, legends, and facet labels are retained.
Figure explanations and statistical details are provided separately in Figure_Captions.
"""
(OUT / "README.txt").write_text(readme, encoding="utf-8")

print(f"Wrote {len(FIGURES)} PDF/SVG figure pairs to {OUT}")
