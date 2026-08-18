from io import BytesIO
from pathlib import Path

from pypdf import PdfReader, PdfWriter, PageObject
from pypdfium2 import PdfDocument
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas


SOURCE = Path(r"C:\Users\admin\Desktop\iMeta\figures_v2\publication_vector_title_free")
OUT = Path("submission_package_FINAL_20260719_v5/02_Main_Figures")
TOP = SOURCE / "Figure5_cross_study_core_microbiome.pdf"
BOTTOM = SOURCE / "Figure6_cross_kingdom_concordance.pdf"
PDF = OUT / "Figure5_core_and_crosskingdom.pdf"
PNG = OUT / "Figure5_core_and_crosskingdom.png"
SVG = OUT / "Figure5_core_and_crosskingdom.svg"
TOP_W, TOP_H = 510.0, 297.0
BOT_W, BOT_H = 510.0, 221.0
GAP = 12.0
TOTAL_H = TOP_H + GAP + BOT_H


def register_arial():
    reg = Path(r"C:\Windows\Fonts\arial.ttf")
    bold = Path(r"C:\Windows\Fonts\arialbd.ttf")
    if reg.exists() and bold.exists():
        pdfmetrics.registerFont(TTFont("ArialLocal", str(reg)))
        pdfmetrics.registerFont(TTFont("ArialLocal-Bold", str(bold)))
        return "ArialLocal", "ArialLocal-Bold"
    return "Helvetica", "Helvetica-Bold"


def label_overlay():
    regular, bold = register_arial()
    buf = BytesIO()
    c = canvas.Canvas(buf, pagesize=(TOP_W, TOTAL_H))
    c.setFillColorRGB(1, 1, 1)
    # Mask the source uppercase panel letters only; titles and plot content remain vector-native.
    masks = [(0, TOTAL_H - 30), (300, TOTAL_H - 30), (0, BOT_H - 30), (250, BOT_H - 30)]
    for x, y in masks:
        c.rect(x, y, 24, 30, fill=1, stroke=0)
    c.setFillColorRGB(0.02, 0.02, 0.02)
    c.setFont(bold, 8.5)
    c.drawString(7, TOTAL_H - 20, "a")
    c.drawString(307, TOTAL_H - 20, "b")
    c.drawString(7, BOT_H - 20, "c")
    c.drawString(257, BOT_H - 20, "d")
    c.save()
    buf.seek(0)
    return PdfReader(buf).pages[0]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    top = PdfReader(str(TOP)).pages[0]
    bottom = PdfReader(str(BOTTOM)).pages[0]
    page = PageObject.create_blank_page(width=TOP_W, height=TOTAL_H)
    page.merge_translated_page(bottom, 0, 0)
    page.merge_translated_page(top, 0, BOT_H + GAP)
    page.merge_page(label_overlay())

    writer = PdfWriter()
    writer.add_page(page)
    with PDF.open("wb") as fh:
        writer.write(fh)

    rendered = PdfDocument(str(PDF))[0].render(scale=450 / 72).to_pil()
    rendered.save(PNG, dpi=(450, 450))
    print(PDF)
    print(PNG)


if __name__ == "__main__":
    main()
