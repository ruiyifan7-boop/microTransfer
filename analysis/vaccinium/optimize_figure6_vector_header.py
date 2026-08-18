from io import BytesIO
from pathlib import Path
import re

from pypdf import PageObject, PdfReader, PdfWriter, Transformation
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas
from pypdfium2 import PdfDocument


ROOT = Path("figures_submission_final_20260719")
SOURCE = ROOT / "Figure6_extended_validation.pdf"
OUT_PDF = ROOT / "Figure6_extended_validation_optimized.pdf"
OUT_PNG = ROOT / "Figure6_extended_validation_optimized.png"
OUT_SVG = ROOT / "Figure6_extended_validation_optimized.svg"
WIDTH, HEIGHT, HEADER = 662.0, 504.0, 48.0
CONTENT_SCALE = 0.90


def make_header_pdf():
    buf = BytesIO()
    font_path = Path(r"C:\Windows\Fonts\arial.ttf")
    bold_path = Path(r"C:\Windows\Fonts\arialbd.ttf")
    if font_path.exists() and bold_path.exists():
        pdfmetrics.registerFont(TTFont("ArialLocal", str(font_path)))
        pdfmetrics.registerFont(TTFont("ArialLocal-Bold", str(bold_path)))
        regular, bold = "ArialLocal", "ArialLocal-Bold"
    else:
        regular, bold = "Helvetica", "Helvetica-Bold"

    c = canvas.Canvas(buf, pagesize=(WIDTH, HEIGHT + HEADER))
    c.setFillColorRGB(1, 1, 1)
    c.rect(0, HEIGHT, WIDTH, HEADER, fill=1, stroke=0)
    c.setStrokeColorRGB(0.82, 0.84, 0.87)
    c.setLineWidth(0.7)
    c.line(20, HEIGHT + 5, WIDTH - 20, HEIGHT + 5)
    c.setFillColorRGB(0.08, 0.10, 0.14)
    c.setFont(bold, 10.0)
    c.drawString(22, HEIGHT + 30, "Study-held-out external validation")
    c.setFillColorRGB(0.25, 0.28, 0.34)
    c.setFont(regular, 6.4)
    c.drawString(
        22, HEIGHT + 16,
        "Bacterial core: broadly transferable   |   Fungal core: narrower and variable   |   pH effects: context dependent",
    )
    c.save()
    buf.seek(0)
    return buf


def build_pdf():
    source_page = PdfReader(str(SOURCE)).pages[0]
    header_page = PdfReader(make_header_pdf()).pages[0]
    page = PageObject.create_blank_page(width=WIDTH, height=HEIGHT + HEADER)
    tx = (WIDTH - WIDTH * CONTENT_SCALE) / 2
    page.merge_transformed_page(
        source_page,
        Transformation().scale(CONTENT_SCALE).translate(tx=tx, ty=HEADER),
    )
    page.merge_page(header_page)
    writer = PdfWriter()
    writer.add_page(page)
    with OUT_PDF.open("wb") as fh:
        writer.write(fh)


def render_png():
    doc = PdfDocument(str(OUT_PDF))
    page = doc[0]
    image = page.render(scale=300 / 72).to_pil()
    image.save(OUT_PNG, dpi=(300, 300))


def build_svg():
    source_svg = (ROOT / "Figure6_extended_validation.svg").read_text(encoding="utf-8")
    source_svg = source_svg.replace(
        'height="504pt" viewBox="0 0 662 504"',
        'height="552pt" viewBox="0 0 662 552"',
        1,
    )
    header = """
<rect x="0" y="0" width="662" height="48" fill="#FFFFFF"/>
<line x1="20" y1="43" x2="642" y2="43" stroke="#D1D5DB" stroke-width="0.7"/>
<text x="22" y="18" font-family="Arial, Liberation Sans, sans-serif" font-size="10" font-weight="700" fill="#141820">Study-held-out external validation</text>
<text x="22" y="33" font-family="Arial, Liberation Sans, sans-serif" font-size="6.4" fill="#404650">Bacterial core: broadly transferable   |   Fungal core: narrower and variable   |   pH effects: context dependent</text>
<g transform="translate(0,48)">
"""
    source_svg = source_svg.replace("</defs>", "</defs>" + header, 1)
    source_svg = source_svg.replace("</svg>", "</g></svg>", 1)
    OUT_SVG.write_text(source_svg, encoding="utf-8", newline="\n")


def main():
    build_pdf()
    render_png()
    build_svg()
    print(f"Wrote {OUT_PDF}")
    print(f"Wrote {OUT_PNG}")
    print(f"Wrote {OUT_SVG}")


if __name__ == "__main__":
    main()
