from io import BytesIO
from pathlib import Path

from pypdf import PageObject, PdfReader, PdfWriter
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas
from pypdfium2 import PdfDocument


ROOT = Path("figures_submission_final_20260719")
SOURCE = ROOT / "Figure4_management_transferability.pdf"
OUT_PDF = ROOT / "Figure4_management_transferability_optimized.pdf"
OUT_PNG = ROOT / "Figure4_management_transferability_optimized.png"
OUT_TIFF = ROOT / "Figure4_management_transferability_optimized_600dpi.tiff"
W, H = 518.0, 425.0


def overlay_pdf():
    buf = BytesIO()
    arial = Path(r"C:\Windows\Fonts\arial.ttf")
    arial_bold = Path(r"C:\Windows\Fonts\arialbd.ttf")
    if arial.exists() and arial_bold.exists():
        pdfmetrics.registerFont(TTFont("ArialLocal", str(arial)))
        pdfmetrics.registerFont(TTFont("ArialLocal-Bold", str(arial_bold)))
        font = "ArialLocal-Bold"
    else:
        font = "Helvetica-Bold"
    c = canvas.Canvas(buf, pagesize=(W, H))
    # Replace only the three facet strips in panel d. Coordinates are in the
    # original 518 x 425 pt PDF coordinate system.
    labels = [
        (205, 284, "Conexibacter"),
        (301, 381, "Candidatus Solibacter"),
        (395, 474, "Family 67-14 lineage"),
    ]
    for x0, x1, label in labels:
        c.setFillColorRGB(1, 1, 1)
        c.rect(x0 - 1, 201, (x1 - x0) + 2, 22, fill=1, stroke=0)
        c.setFillColorRGB(0.08, 0.10, 0.14)
        c.setFont(font, 6.1)
        c.drawCentredString((x0 + x1) / 2, 211.5, label)
    c.save()
    buf.seek(0)
    return buf


def main():
    source = PdfReader(str(SOURCE)).pages[0]
    overlay = PdfReader(overlay_pdf()).pages[0]
    page = PageObject.create_blank_page(width=W, height=H)
    page.merge_page(source)
    page.merge_page(overlay)
    writer = PdfWriter()
    writer.add_page(page)
    with OUT_PDF.open("wb") as fh:
        writer.write(fh)
    doc = PdfDocument(str(OUT_PDF))
    image = doc[0].render(scale=300 / 72).to_pil()
    image.save(OUT_PNG, dpi=(300, 300))
    image_600 = doc[0].render(scale=600 / 72).to_pil()
    image_600.save(OUT_TIFF, dpi=(600, 600), compression="tiff_lzw")
    print(f"Wrote {OUT_PDF}")
    print(f"Wrote {OUT_PNG}")
    print(f"Wrote {OUT_TIFF}")


if __name__ == "__main__":
    main()
