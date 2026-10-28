from PIL import Image, ImageDraw, ImageFont
from pathlib import Path

PROJECT = Path(r"C:\Users\akidi\Desktop\ISO_AK_Docking")

folder = PROJECT / "results" / "AK1_panelH_Zn"

h1_file = folder / "AK1_panelH1_ZN601_clean.png"
h2_file = folder / "AK1_panelH2_ZN641_clean.png"

out_file = folder / "AK1_panelH_Zn_environment_FINAL.png"

# ------------------------------------------------------------
# Load images
# ------------------------------------------------------------

h1 = Image.open(h1_file).convert("RGB")
h2 = Image.open(h2_file).convert("RGB")

# Make both exactly the same size
target_w = max(h1.width, h2.width)
target_h = max(h1.height, h2.height)

def fit_white(img, w, h):
    canvas = Image.new("RGB", (w, h), "white")
    x = (w - img.width) // 2
    y = (h - img.height) // 2
    canvas.paste(img, (x, y))
    return canvas

h1 = fit_white(h1, target_w, target_h)
h2 = fit_white(h2, target_w, target_h)

# ------------------------------------------------------------
# Layout
# ------------------------------------------------------------

top_space = 260
bottom_space = 100
gap = 60
side_margin = 60

canvas_w = side_margin * 2 + target_w * 2 + gap
canvas_h = top_space + target_h + bottom_space

canvas = Image.new(
    "RGB",
    (canvas_w, canvas_h),
    "white"
)

canvas.paste(
    h1,
    (side_margin, top_space)
)

canvas.paste(
    h2,
    (side_margin + target_w + gap, top_space)
)

draw = ImageDraw.Draw(canvas)

# ------------------------------------------------------------
# Fonts
# ------------------------------------------------------------

arial_bold = r"C:\Windows\Fonts\arialbd.ttf"
arial = r"C:\Windows\Fonts\arial.ttf"

panel_font = ImageFont.truetype(
    arial_bold,
    72
)

title_font = ImageFont.truetype(
    arial_bold,
    56
)

subtitle_font = ImageFont.truetype(
    arial,
    42
)

# ------------------------------------------------------------
# Left panel H1
# ------------------------------------------------------------

left_x = side_margin

draw.text(
    (left_x, 35),
    "H1",
    font=panel_font,
    fill="black"
)

draw.text(
    (left_x + 110, 48),
    "AP5-associated ZN601 site",
    font=title_font,
    fill="black"
)

draw.text(
    (left_x + 110, 125),
    "ZN601 coordinates AP5 phosphate oxygens at ~2.05–2.07 Å",
    font=subtitle_font,
    fill="black"
)

# ------------------------------------------------------------
# Right panel H2
# ------------------------------------------------------------

right_x = side_margin + target_w + gap

draw.text(
    (right_x, 35),
    "H2",
    font=panel_font,
    fill="black"
)

draw.text(
    (right_x + 110, 48),
    "Protein-associated ZN641 site",
    font=title_font,
    fill="black"
)

draw.text(
    (right_x + 110, 125),
    "ZN641 coordinates GLU26 and HIS36 at ~1.97 and ~2.32 Å",
    font=subtitle_font,
    fill="black"
)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

canvas.save(
    out_file,
    dpi=(600, 600)
)

print("Saved:")
print(out_file)