from PIL import Image, ImageDraw, ImageFont
from pathlib import Path

PROJECT = Path(r"C:\Users\akidi\Desktop\ISO_AK_Docking")

inp = (
    PROJECT
    / "results"
    / "AK1_panelC"
    / "AK1_panelC_AP5_redocking_validation.png"
)

out = (
    PROJECT
    / "results"
    / "AK1_panelC"
    / "AK1_panelC_AP5_redocking_validation_LEGEND.png"
)

img = Image.open(inp).convert("RGB")
draw = ImageDraw.Draw(img)

# ------------------------------------------------------------
# Fonts
# ------------------------------------------------------------

arial = r"C:\Windows\Fonts\arial.ttf"
arial_bold = r"C:\Windows\Fonts\arialbd.ttf"

title_font = ImageFont.truetype(arial_bold, 58)
text_font = ImageFont.truetype(arial, 48)
rmsd_font = ImageFont.truetype(arial_bold, 54)

# ------------------------------------------------------------
# Legend box
# Put it in the lower-left white area
# ------------------------------------------------------------

x0 = 100
y0 = img.height - 600
x1 = 1250
y1 = img.height - 100

draw.rounded_rectangle(
    [x0, y0, x1, y1],
    radius=28,
    fill=(255, 255, 255),
    outline=(170, 170, 170),
    width=3
)

# ------------------------------------------------------------
# Title
# ------------------------------------------------------------

draw.text(
    (x0 + 50, y0 + 35),
    "AK1 native-ligand redocking",
    font=title_font,
    fill=(0, 0, 0)
)

# ------------------------------------------------------------
# Legend entries
# ------------------------------------------------------------

entries = [
    ("Crystallographic AP5", (0, 220, 220)),
    ("Redocked AP5", (220, 0, 220)),
    ("AK1 protein", (150, 150, 150)),
    ("AP5-site residues", (145, 210, 145)),
]

start_y = y0 + 125
step = 72

for i, (label, color) in enumerate(entries):

    yy = start_y + i * step

    draw.line(
        [(x0 + 60, yy + 22), (x0 + 160, yy + 22)],
        fill=color,
        width=18
    )

    draw.text(
        (x0 + 195, yy),
        label,
        font=text_font,
        fill=(0, 0, 0)
    )

# ------------------------------------------------------------
# RMSD annotation
# ------------------------------------------------------------

draw.text(
    (x0 + 50, y1 - 90),
    "Redocking RMSD = 1.776 Å",
    font=rmsd_font,
    fill=(0, 0, 0)
)

# ------------------------------------------------------------
# Save as NEW file
# ------------------------------------------------------------

img.save(out, dpi=(600, 600))

print("Saved:")
print(out)