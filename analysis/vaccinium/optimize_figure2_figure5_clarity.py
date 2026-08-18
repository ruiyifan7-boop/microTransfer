from io import BytesIO
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont
from pypdf import PdfReader, PdfWriter
from reportlab.pdfgen import canvas
from pypdfium2 import PdfDocument


ROOT = Path("submission_package_FINAL_20260719_v5")
FIG = ROOT / "02_Main_Figures"


def font(size, bold=False):
    p = Path(r"C:\Windows\Fonts\arialbd.ttf" if bold else r"C:\Windows\Fonts\arial.ttf")
    if p.exists():
        return ImageFont.truetype(str(p), size)
    return ImageFont.load_default()


def optimize_figure2():
    source = FIG / "Figure2_controlled_pH_responses.pdf"
    target_pdf = FIG / "Figure2_controlled_pH_responses.pdf"
    overlay = BytesIO()
    c = canvas.Canvas(overlay, pagesize=(518, 386))
    c.setFillColorRGB(0.08, 0.12, 0.18)
    c.setFont("Helvetica-Bold", 7.2)
    for label, y in (("Bacteria", 258), ("Fungi", 116)):
        c.saveState()
        c.translate(16, y)
        c.rotate(90)
        c.drawCentredString(0, 0, label)
        c.restoreState()
    c.save()
    overlay.seek(0)

    base = PdfReader(str(source))
    over = PdfReader(overlay)
    base.pages[0].merge_page(over.pages[0])
    tmp_pdf = target_pdf.with_suffix(".tmp.pdf")
    with tmp_pdf.open("wb") as fh:
        PdfWriter().write(fh) if False else None
        writer = PdfWriter()
        writer.add_page(base.pages[0])
        writer.write(fh)
    tmp_pdf.replace(target_pdf)

    png = PdfDocument(str(target_pdf))[0].render(scale=300 / 72).to_pil()
    png.save(FIG / "Figure2_controlled_pH_responses.png", dpi=(300, 300))


def optimize_figure5():
    path = FIG / "Figure5_core_and_crosskingdom.png"
    base = Image.open(path).convert("RGB")
    im = Image.new("RGB", (base.width, base.height + 80), "white")
    im.paste(base, (0, 0))
    draw = ImageDraw.Draw(im)
    navy = (20, 40, 77)
    grey = (75, 82, 93)
    blue = (36, 104, 222)
    orange = (209, 95, 0)
    legend_font = font(20, bold=False)

    # Missing legend for the two correlation summaries in panels C and D.
    # Use a dedicated bottom margin so the legend cannot obscure any data panel.
    y = base.height + 40
    draw.ellipse((735, y - 11, 757, y + 11), fill=blue)
    draw.text((770, y - 15), "Raw", fill=grey, font=legend_font)
    draw.ellipse((900, y - 11, 922, y + 11), fill=orange)
    draw.text((935, y - 15), "Design-adjusted", fill=grey, font=legend_font)
    draw.text((1165, y - 15), "filled: P < 0.05; open: P >= 0.05", fill=grey, font=legend_font)

    out = FIG / "Figure5_core_and_crosskingdom.png"
    im.save(out, dpi=(450, 450))


if __name__ == "__main__":
    optimize_figure2()
    optimize_figure5()
    print("Optimized Figure 2 row labels and Figure 5 correlation legend.")
