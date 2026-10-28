pdb = r"..\ligand\AK1_native_AP5_chainA.pdb"

coords = []

with open(pdb) as f:
    for line in f:
        if line.startswith("HETATM"):
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coords.append((x, y, z))

xs = [x[0] for x in coords]
ys = [x[1] for x in coords]
zs = [x[2] for x in coords]

center_x = (min(xs) + max(xs)) / 2
center_y = (min(ys) + max(ys)) / 2
center_z = (min(zs) + max(zs)) / 2

# Add 5 Å of space on every side of the native ligand
size_x = (max(xs) - min(xs)) + 10
size_y = (max(ys) - min(ys)) + 10
size_z = (max(zs) - min(zs)) + 10

print()
print("AK1 docking box")
print("----------------")
print("center_x =", round(center_x, 3))
print("center_y =", round(center_y, 3))
print("center_z =", round(center_z, 3))
print()
print("size_x =", round(size_x, 3))
print("size_y =", round(size_y, 3))
print("size_z =", round(size_z, 3))