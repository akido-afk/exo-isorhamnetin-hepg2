from pathlib import Path
from PIL import Image

frame_dir = Path(r"..\results\AK1_rotate")
frames = sorted(frame_dir.glob("frame_*.png"))

if not frames:
    raise ValueError("No PNG frames found.")

images = [Image.open(f).convert("RGB") for f in frames]

outfile = r"..\results\AK1_ISO_rotation_clean.gif"

images[0].save(
    outfile,
    save_all=True,
    append_images=images[1:],
    duration=80,      # milliseconds per frame
    loop=0,
    disposal=2,       # important: clear previous frame
    optimize=False
)

print("Saved:", outfile)