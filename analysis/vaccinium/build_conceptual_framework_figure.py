from io import BytesIO
from pathlib import Path

from pypdf import PdfReader, PdfWriter
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas
from pypdfium2 import PdfDocument


OUT = Path("submission_package_FINAL_20260719_v5/03_Supplementary_Figures")
OUT.mkdir(parents=True, exist_ok=True)
PDF = OUT / "FigureS2_conceptual_framework.pdf"
PNG = OUT / "FigureS2_conceptual_framework.png"
SVG = OUT / "FigureS2_conceptual_framework.svg"
W, H = 518.0, 270.0


def fonts():
    reg = Path(r"C:\Windows\Fonts\arial.ttf")
    bold = Path(r"C:\Windows\Fonts\arialbd.ttf")
    if reg.exists() and bold.exists():
        pdfmetrics.registerFont(TTFont("ArialLocal", str(reg)))
        pdfmetrics.registerFont(TTFont("ArialLocal-Bold", str(bold)))
        return "ArialLocal", "ArialLocal-Bold"
    return "Helvetica", "Helvetica-Bold"


def draw_wrapped(c, text, x, y, width, font, size, leading, color):
    c.setFont(font, size)
    c.setFillColorRGB(*color)
    words = text.split()
    lines, line = [], ""
    for word in words:
        trial = (line + " " + word).strip()
        if c.stringWidth(trial, font, size) <= width or not line:
            line = trial
        else:
            lines.append(line)
            line = word
    if line:
        lines.append(line)
    for i, item in enumerate(lines):
        c.drawCentredString(x + width / 2, y - i * leading, item)


def build_pdf():
    regular, bold = fonts()
    buf = BytesIO()
    c = canvas.Canvas(buf, pagesize=(W, H))
    c.setFillColorRGB(1, 1, 1)
    c.rect(0, 0, W, H, fill=1, stroke=0)

    navy = (0.08, 0.16, 0.30)
    blue = (0.00, 0.45, 0.70)
    orange = (0.83, 0.37, 0.00)
    green = (0.00, 0.55, 0.42)
    grey = (0.27, 0.30, 0.35)

    c.setFillColorRGB(*navy)
    c.roundRect(145, 231, 228, 25, 7, fill=1, stroke=0)
    c.setFillColorRGB(1, 1, 1)
    c.setFont(bold, 9.2)
    c.drawCentredString(W / 2, 240, "Study-aware dual-kingdom benchmark")

    c.setStrokeColorRGB(*grey)
    c.setLineWidth(0.8)
    c.line(W / 2, 231, W / 2, 217)
    c.line(W / 2, 217, 92, 217)
    c.line(W / 2, 217, 259, 217)
    c.line(W / 2, 217, 426, 217)
    for x in (92, 259, 426):
        c.line(x, 217, x, 207)
        c.setFillColorRGB(*grey)
        c.setStrokeColorRGB(*grey)
        c.setLineWidth(0.8)
        c.line(x, 207, x - 3, 212)
        c.line(x, 207, x + 3, 212)

    boxes = [
        (20, blue, "CORE MEMBERSHIP", "Cross-study conserved bacterial backbone", "held-out prevalence validation"),
        (187, orange, "ENVIRONMENTAL RESPONSE", "pH filters communities", "context-dependent; no universal pH biomarker"),
        (354, green, "PREDICTIVE TRANSFERABILITY", "Within-study management stress test", "portable prediction not established"),
    ]
    for x, col, heading, main, foot in boxes:
        c.setFillColorRGB(0.97, 0.98, 0.99)
        c.setStrokeColorRGB(*col)
        c.setLineWidth(1.2)
        c.roundRect(x, 118, 144, 88, 7, fill=1, stroke=1)
        c.setFillColorRGB(*col)
        c.setFont(bold, 6.2)
        c.drawCentredString(x + 72, 190, heading)
        draw_wrapped(c, main, x + 8, 164, 128, bold, 8.0, 10.0, navy)
        draw_wrapped(c, foot, x + 10, 137, 124, regular, 6.5, 8.2, grey)

    c.setFillColorRGB(0.95, 0.96, 0.97)
    c.roundRect(42, 43, 434, 42, 8, fill=1, stroke=0)
    c.setFillColorRGB(*navy)
    c.setFont(bold, 8.0)
    c.drawCentredString(W / 2, 69, "Interpretation")
    c.setFont(regular, 7.0)
    c.drawCentredString(W / 2, 55, "Reproducible membership does not imply reproducible response or portable prediction")
    c.setFillColorRGB(*grey)
    c.setFont(regular, 5.9)
    c.drawCentredString(W / 2, 26, "Study-held-out designs separate what recurs across cohorts from what remains context-specific.")
    c.save()
    buf.seek(0)
    with PDF.open("wb") as fh:
        fh.write(buf.read())


def build_svg():
    svg = f'''<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="{W}pt" height="{H}pt" viewBox="0 0 {W} {H}">
<rect width="100%" height="100%" fill="white"/>
<rect x="145" y="14" width="228" height="25" rx="7" fill="#14284D"/>
<text x="259" y="30" text-anchor="middle" font-family="Arial, sans-serif" font-size="9.2" font-weight="700" fill="white">Study-aware dual-kingdom benchmark</text>
<path d="M259 39 V53 H92 M259 53 H259 M259 53 H426 M92 53 V63 M259 53 V63 M426 53 V63" fill="none" stroke="#454B57" stroke-width="0.8"/>
<g font-family="Arial, sans-serif">
<rect x="20" y="64" width="144" height="88" rx="7" fill="#F8FAFB" stroke="#0073B3" stroke-width="1.2"/>
<rect x="187" y="64" width="144" height="88" rx="7" fill="#F8FAFB" stroke="#D45F00" stroke-width="1.2"/>
<rect x="354" y="64" width="144" height="88" rx="7" fill="#F8FAFB" stroke="#008C6B" stroke-width="1.2"/>
<text x="92" y="80" text-anchor="middle" font-size="6.2" font-weight="700" fill="#0073B3">CORE MEMBERSHIP</text>
<text x="259" y="80" text-anchor="middle" font-size="6.2" font-weight="700" fill="#D45F00">ENVIRONMENTAL RESPONSE</text>
<text x="426" y="80" text-anchor="middle" font-size="6.2" font-weight="700" fill="#008C6B">PREDICTIVE TRANSFERABILITY</text>
<text x="92" y="108" text-anchor="middle" font-size="8" font-weight="700" fill="#14284D">Cross-study conserved</text>
<text x="92" y="119" text-anchor="middle" font-size="8" font-weight="700" fill="#14284D">bacterial backbone</text>
<text x="92" y="139" text-anchor="middle" font-size="6.5" fill="#454B57">held-out prevalence validation</text>
<text x="259" y="108" text-anchor="middle" font-size="8" font-weight="700" fill="#14284D">pH filters communities</text>
<text x="259" y="130" text-anchor="middle" font-size="6.5" fill="#454B57">context-dependent; no universal</text>
<text x="259" y="139" text-anchor="middle" font-size="6.5" fill="#454B57">pH biomarker</text>
<text x="426" y="108" text-anchor="middle" font-size="8" font-weight="700" fill="#14284D">Within-study management</text>
<text x="426" y="119" text-anchor="middle" font-size="8" font-weight="700" fill="#14284D">stress test</text>
<text x="426" y="139" text-anchor="middle" font-size="6.5" fill="#454B57">portable prediction not established</text>
<rect x="42" y="185" width="434" height="42" rx="8" fill="#F2F4F6"/>
<text x="259" y="201" text-anchor="middle" font-size="8" font-weight="700" fill="#14284D">Interpretation</text>
<text x="259" y="215" text-anchor="middle" font-size="7" fill="#14284D">Reproducible membership does not imply reproducible response or portable prediction</text>
<text x="259" y="242" text-anchor="middle" font-size="5.9" fill="#454B57">Study-held-out designs separate what recurs across cohorts from what remains context-specific.</text>
</g></svg>'''
    SVG.write_text(svg, encoding="utf-8")


def main():
    build_pdf()
    build_svg()
    image = PdfDocument(str(PDF))[0].render(scale=300 / 72).to_pil()
    image.save(PNG, dpi=(300, 300))
    print(PDF)
    print(PNG)
    print(SVG)


if __name__ == "__main__":
    main()
