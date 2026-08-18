from io import BytesIO
from pathlib import Path

from pypdf import PageObject, PdfReader, PdfWriter
from pypdf.generic import FloatObject, NameObject
from pypdfium2 import PdfDocument
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas


ROOT = Path("submission_package_FINAL_20260719_v5")
FIG = ROOT / "02_Main_Figures"
FONT = Path(r"C:\Windows\Fonts\arial.ttf")
FONT_BOLD = Path(r"C:\Windows\Fonts\arialbd.ttf")


def fonts():
    if FONT.exists() and FONT_BOLD.exists():
        pdfmetrics.registerFont(TTFont("ArialLocal", str(FONT)))
        pdfmetrics.registerFont(TTFont("ArialLocal-Bold", str(FONT_BOLD)))
        return "ArialLocal", "ArialLocal-Bold"
    return "Helvetica", "Helvetica-Bold"


def overlay_for(width, height, letters, titles):
    regular, bold = fonts()
    buf = BytesIO()
    c = canvas.Canvas(buf, pagesize=(width, height))
    navy = (0.08, 0.10, 0.14)

    for x, top, label in letters:
        c.setFillColorRGB(1, 1, 1)
        c.rect(x - 2, height - top - 18, 27, 19, fill=1, stroke=0)
        c.setFillColorRGB(0.02, 0.02, 0.02)
        c.setFont(bold, 8.2)
        c.drawString(x, height - top - 9, label)

    c.setFillColorRGB(*navy)
    c.setFont(bold, 5.7)
    for x, top, title in titles:
        c.drawString(x, height - top - 6.1, title)
    c.save()
    buf.seek(0)
    return PdfReader(buf).pages[0]


def apply_overlay(name, letters, titles):
    path = FIG / name
    base = PdfReader(str(path)).pages[0]
    width = float(base.mediabox.width)
    height = float(base.mediabox.height)
    page = PageObject.create_blank_page(width=width, height=height)
    page.merge_page(base)
    page.merge_page(overlay_for(width, height, letters, titles))
    tmp = path.with_suffix(".tmp.pdf")
    writer = PdfWriter()
    writer.add_page(page)
    with tmp.open("wb") as fh:
        writer.write(fh)
    tmp.replace(path)
    png = path.with_suffix(".png")
    image = PdfDocument(str(path))[0].render(scale=300 / 72).to_pil()
    image.save(png, dpi=(300, 300))


def figure1_to4():
    apply_overlay(
        "Figure1_cohort_design_and_QC.pdf",
        [(5, 4, "A"), (264, 4, "B")],
        [(22, 16, "Cohort retention"), (282, 16, "Sequencing depth")],
    )
    apply_overlay(
        "Figure2_controlled_pH_responses.pdf",
        [(35, 43, "A"), (187, 43, "B"), (339, 43, "C"),
         (35, 192, "D"), (187, 192, "E"), (339, 192, "F")],
        [],
    )
    apply_overlay(
        "Figure3_pH_transferability.pdf",
        [(5, 4, "A"), (329, 4, "B"), (5, 190, "C")],
        [(20, 16, "Community-level pH effects"),
         (344, 16, "Taxon-level validation")],
    )
    apply_overlay(
        "Figure4_management_transferability.pdf",
        [(3, 2, "A"), (371, 33, "B"), (3, 205, "C"), (187, 205, "D")],
        [],
    )


def figure6_compact_header():
    source_pdf_original = Path(r"C:\Users\admin\Desktop\iMeta\figures_extended_validation\Figure_extended_validation.pdf")
    source_svg = Path(r"C:\Users\admin\Desktop\iMeta\figures_extended_validation\Figure_extended_validation.svg")
    target_pdf = FIG / "Figure6_extended_validation.pdf"
    target_png = FIG / "Figure6_extended_validation.png"
    target_svg = FIG / "Figure6_extended_validation.svg"
    width, content_h, header, scale = 662.0, 504.0, 34.0, 0.93
    regular, bold = fonts()

    # Move the four source panel labels in the source coordinate system,
    # instead of masking them after scaling. This preserves all plotted marks
    # and count labels while placing A-D in the whitespace outside the panels.
    repositioned_pdf = Path("work") / "Figure_extended_validation_repositioned.pdf"
    source_reader = PdfReader(str(source_pdf_original))
    source_page = source_reader.pages[0]
    from pypdf.generic import ContentStream
    content = ContentStream(source_page.get_contents(), source_reader)
    moved_labels = {
        241: (7.3, 11.8),    # A
        585: (297.7, 11.8),  # B
        788: (25.0, 245.0),  # C
        1011: (395.0, 245.0),# D
    }
    for idx, (x, y) in moved_labels.items():
        tm_idx = idx - 2
        if tm_idx < len(content.operations) and content.operations[tm_idx][1] == b"Tm":
            operands = list(content.operations[tm_idx][0])
            operands[-2] = FloatObject(x)
            operands[-1] = FloatObject(y)
            content.operations[tm_idx] = (operands, b"Tm")
    # Remove the original D text object so the replacement can be drawn after
    # the heatmap layer and remain visible in both PDF and SVG outputs.
    content.operations = [
        operation for i, operation in enumerate(content.operations)
        if not 1008 <= i <= 1012
    ]
    source_page.replace_contents(content)
    repositioned_pdf.parent.mkdir(parents=True, exist_ok=True)
    repositioned_writer = PdfWriter()
    repositioned_writer.add_page(source_page)
    with repositioned_pdf.open("wb") as fh:
        repositioned_writer.write(fh)

    buf = BytesIO()
    c = canvas.Canvas(buf, pagesize=(width, content_h + header))
    c.setFillColorRGB(1, 1, 1)
    c.rect(0, content_h, width, header, fill=1, stroke=0)
    c.setStrokeColorRGB(0.82, 0.84, 0.87)
    c.setLineWidth(0.7)
    c.line(20, content_h + 5, width - 20, content_h + 5)
    c.setFillColorRGB(0.08, 0.10, 0.14)
    c.setFont(bold, 9.2)
    c.drawString(22, content_h + 20, "Study-held-out validation")
    c.setFillColorRGB(0.25, 0.28, 0.34)
    c.setFont(regular, 6.1)
    c.drawString(22, content_h + 9, "Core membership, pH response and management discrimination")
    c.save()
    buf.seek(0)

    source_page = PdfReader(str(repositioned_pdf)).pages[0]
    header_page = PdfReader(buf).pages[0]
    page = PageObject.create_blank_page(width=width, height=content_h + header)
    tx = (width - width * scale) / 2
    from pypdf import Transformation
    page.merge_transformed_page(source_page, Transformation().scale(scale).translate(tx=tx, ty=header))
    page.merge_page(header_page)
    d_buf = BytesIO()
    d_canvas = canvas.Canvas(d_buf, pagesize=(width, content_h + header))
    d_canvas.setFillColorRGB(0.02, 0.02, 0.02)
    d_canvas.setFont(bold, 13.2)
    d_canvas.drawString(tx + scale * 395.0, header + scale * (content_h - 245.0), "D")
    d_canvas.save()
    d_buf.seek(0)
    page.merge_page(PdfReader(d_buf).pages[0])

    writer = PdfWriter()
    writer.add_page(page)
    with target_pdf.open("wb") as fh:
        writer.write(fh)
    image = PdfDocument(str(target_pdf))[0].render(scale=300 / 72).to_pil()
    image.save(target_png, dpi=(300, 300))

    text = source_svg.read_text(encoding="utf-8")
    text = text.replace('x="44.75" y="23.081055"', 'x="7.3" y="11.8"', 1)
    text = text.replace('x="375.75" y="23.081055"', 'x="297.7" y="11.8"', 1)
    text = text.replace('x="40.75" y="272.588867"', 'x="25.0" y="245.0"', 1)
    text = text.replace('x="406.632812" y="272.588867"', 'x="395.0" y="245.0"', 1)
    text = text.replace(
        '<g style="fill:rgb(0%,0%,0%);fill-opacity:1;">\n'
        '  <use xlink:href="#glyph2-4" x="395.0" y="245.0"/>\n'
        '</g>\n',
        '',
        1,
    )
    text = text.replace('height="504pt" viewBox="0 0 662 504"', 'height="538pt" viewBox="0 0 662 538"', 1)
    header_svg = (
        '<rect x="0" y="0" width="662" height="34" fill="#FFFFFF"/>'
        '<line x1="20" y1="31" x2="642" y2="31" stroke="#D1D5DB" stroke-width="0.7"/>'
        '<text x="22" y="15" font-family="Arial, Liberation Sans, sans-serif" font-size="9.2" font-weight="700" fill="#141820">Study-held-out validation</text>'
        '<text x="22" y="27" font-family="Arial, Liberation Sans, sans-serif" font-size="6.1" fill="#404650">Core membership, pH response and management discrimination</text>'
        '<g transform="translate(0,34)">'
    )
    d_svg = '<text x="395" y="245" font-family="Arial, Liberation Sans, sans-serif" font-size="14.2" font-weight="700" fill="#000000">D</text>'
    text = text.replace("</defs>", "</defs>" + header_svg, 1).replace("</svg>", d_svg + "</g></svg>", 1)
    target_svg.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    figure1_to4()
    figure6_compact_header()
    print("Applied iMeta-style panel labels and titles to Figures 1-4 and compacted Figure 6 header.")
