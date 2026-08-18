from io import BytesIO
from pathlib import Path

from pypdf import PageObject, PdfReader, PdfWriter, Transformation
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas
from pypdfium2 import PdfDocument


ROOT = Path("figures_submission_final_20260719")
SOURCE = ROOT / "Figure2_controlled_pH_responses.pdf"
OUT_PDF = ROOT / "Figure2_controlled_pH_responses_optimized.pdf"
OUT_PNG = ROOT / "Figure2_controlled_pH_responses_optimized.png"
W, H, HEADER = 518.0, 340.0, 46.0
SCALE = 0.88


def header_pdf():
    buf = BytesIO()
    regular_path = Path(r"C:\Windows\Fonts\arial.ttf")
    bold_path = Path(r"C:\Windows\Fonts\arialbd.ttf")
    if regular_path.exists() and bold_path.exists():
        pdfmetrics.registerFont(TTFont("ArialLocal", str(regular_path)))
        pdfmetrics.registerFont(TTFont("ArialLocal-Bold", str(bold_path)))
        regular, bold = "ArialLocal", "ArialLocal-Bold"
    else:
        regular, bold = "Helvetica", "Helvetica-Bold"
    c = canvas.Canvas(buf, pagesize=(W, H + HEADER))
    c.setFillColorRGB(1, 1, 1)
    c.rect(0, H, W, HEADER, fill=1, stroke=0)
    c.setStrokeColorRGB(0.82, 0.84, 0.87)
    c.setLineWidth(0.7)
    c.line(18, H + 4, W - 18, H + 4)
    c.setFillColorRGB(0.08, 0.10, 0.14)
    c.setFont(bold, 9.2)
    c.drawString(20, H + 29, "Controlled pH response across compartments")
    c.setFillColorRGB(0.28, 0.30, 0.35)
    c.setFont(regular, 5.9)
    c.drawString(
        20, H + 16,
        "Ordinations show treatment-centroid trajectories from Low to Optimum to High pH; rows separate bacteria and fungi.",
    )
    legend = [("Low", (0.231, 0.510, 0.965)), ("Optimum", (0.0, 0.620, 0.451)), ("High", (0.835, 0.369, 0.0))]
    x = 340
    for label, col in legend:
        c.setFillColorRGB(*col)
        c.circle(x, H + 16, 3.2, fill=1, stroke=0)
        c.setFillColorRGB(0.18, 0.20, 0.24)
        c.setFont(regular, 5.7)
        c.drawString(x + 7, H + 14, label)
        x += 53 if label != "Optimum" else 65
    c.save()
    buf.seek(0)
    return buf


def main():
    source = PdfReader(str(SOURCE)).pages[0]
    overlay = PdfReader(header_pdf()).pages[0]
    page = PageObject.create_blank_page(width=W, height=H + HEADER)
    tx = (W - W * SCALE) / 2
    page.merge_transformed_page(
        source, Transformation().scale(SCALE).translate(tx=tx, ty=HEADER)
    )
    page.merge_page(overlay)
    writer = PdfWriter()
    writer.add_page(page)
    with OUT_PDF.open("wb") as fh:
        writer.write(fh)
    image = PdfDocument(str(OUT_PDF))[0].render(scale=300 / 72).to_pil()
    image.save(OUT_PNG, dpi=(300, 300))
    print(f"Wrote {OUT_PDF}")
    print(f"Wrote {OUT_PNG}")


if __name__ == "__main__":
    main()
