import imageio.v2 as imageio
from pathlib import Path

frame_dir = Path(r"..\results\AK1_rotate")
frames = sorted(frame_dir.glob("frame_*.png"))

if not frames:
    raise ValueError("No PNG frames found.")

images = [imageio.imread(f) for f in frames]

outfile = r"..\results\AK1_ISO_rotation.gif"
imageio.mimsave(outfile, images, duration=0.08, loop=0)

print("Saved:", outfile)