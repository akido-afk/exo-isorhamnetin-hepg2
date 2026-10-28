from pathlib import Path
import imageio.v2 as imageio

# folder containing the PNG frames
frame_dir = Path(r"C:\Users\akidi\Desktop\ISO_AK_Docking\results\AK1_rotate")

# output GIF path
gif_path = Path(r"C:\Users\akidi\Desktop\ISO_AK_Docking\results\AK1_rotate_v2.gif")

# collect and sort frames
frames = sorted(frame_dir.glob("frame_*.png"))

if not frames:
    raise FileNotFoundError("No PNG frames found in the AK1_rotate folder.")

images = []
for f in frames:
    images.append(imageio.imread(f))

# duration = seconds per frame
# 0.08 = fairly smooth
imageio.mimsave(gif_path, images, duration=0.08, loop=0)

print("GIF saved to:")
print(gif_path)
print(f"Frames used: {len(frames)}")