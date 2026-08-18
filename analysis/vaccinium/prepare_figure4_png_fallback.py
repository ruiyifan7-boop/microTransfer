from pathlib import Path

from PIL import Image
import fitz


source = Path(r"C:\Users\admin\Desktop\Figures\Figure_management_transferability_v2.png")
outdir = Path(r"C:\Users\admin\Desktop\Figures\figures_v4")
outdir.mkdir(parents=True, exist_ok=True)

png_out = outdir / "Figure4_management_transferability_999.png"
pdf_out = outdir / "Figure4_management_transferability_999.pdf"

image = Image.open(source).convert("RGB")
image.save(png_out, format="PNG", optimize=True)

# PDF compatibility version. The source available locally is a high-resolution
# raster, so this PDF preserves resolution but is not a native editable vector.
page_width = 756
page_height = page_width * image.height / image.width
doc = fitz.open()
page = doc.new_page(width=page_width, height=page_height)
page.insert_image(page.rect, filename=str(png_out))
doc.save(pdf_out, deflate=True)
doc.close()

print(png_out)
print(pdf_out)
