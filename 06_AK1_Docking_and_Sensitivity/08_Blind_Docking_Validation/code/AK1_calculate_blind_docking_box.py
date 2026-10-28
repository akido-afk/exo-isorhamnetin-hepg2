from pathlib import Path
import numpy as np

PROJECT = Path(r"C:\Users\akidi\Desktop\ISO_AK_Docking")

RECEPTOR = (
    PROJECT
    / "receptors"
    / "AK1_chainA_FIXED.pdb"
)

OUT_DIR = (
    PROJECT
    / "results"
    / "AK1_blind_docking"
)

OUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

PADDING = 6.0


coords = []


with open(
    RECEPTOR,
    "r",
    encoding="utf-8",
    errors="ignore"
) as handle:

    for line in handle:

        if not line.startswith("ATOM"):
            continue

        element = (
            line[76:78]
            .strip()
            .upper()
            if len(line) >= 78
            else ""
        )

        if element == "H":
            continue

        try:

            x = float(
                line[30:38]
            )

            y = float(
                line[38:46]
            )

            z = float(
                line[46:54]
            )

        except Exception:
            continue

        coords.append(
            [x, y, z]
        )


coords = np.array(
    coords,
    dtype=float
)


if len(coords) == 0:

    raise RuntimeError(
        "No receptor heavy-atom coordinates found."
    )


minimum = np.min(
    coords,
    axis=0
)

maximum = np.max(
    coords,
    axis=0
)

center = (
    minimum
    +
    maximum
) / 2.0


protein_size = (
    maximum
    -
    minimum
)


blind_size = (
    protein_size
    +
    2.0 * PADDING
)


volume = float(
    np.prod(
        blind_size
    )
)


print()
print("=" * 80)
print("AK1 WHOLE-PROTEIN BLIND-DOCKING BOX")
print("=" * 80)

print(
    f"Heavy atoms: {len(coords)}"
)

print()

print(
    f"Protein min XYZ : "
    f"{minimum[0]:.3f}, "
    f"{minimum[1]:.3f}, "
    f"{minimum[2]:.3f}"
)

print(
    f"Protein max XYZ : "
    f"{maximum[0]:.3f}, "
    f"{maximum[1]:.3f}, "
    f"{maximum[2]:.3f}"
)

print()

print(
    f"Protein span XYZ: "
    f"{protein_size[0]:.3f}, "
    f"{protein_size[1]:.3f}, "
    f"{protein_size[2]:.3f} A"
)

print()

print(
    f"Padding per side: {PADDING:.1f} A"
)

print()

print(
    f"center_x = {center[0]:.3f}"
)

print(
    f"center_y = {center[1]:.3f}"
)

print(
    f"center_z = {center[2]:.3f}"
)

print()

print(
    f"size_x = {blind_size[0]:.3f}"
)

print(
    f"size_y = {blind_size[1]:.3f}"
)

print(
    f"size_z = {blind_size[2]:.3f}"
)

print()

print(
    f"Search volume = {volume:.1f} A^3"
)


config_file = (
    OUT_DIR
    / "AK1_blind_box.txt"
)


with open(
    config_file,
    "w"
) as handle:

    handle.write(
        f"center_x = {center[0]:.3f}\n"
    )

    handle.write(
        f"center_y = {center[1]:.3f}\n"
    )

    handle.write(
        f"center_z = {center[2]:.3f}\n"
    )

    handle.write(
        f"size_x = {blind_size[0]:.3f}\n"
    )

    handle.write(
        f"size_y = {blind_size[1]:.3f}\n"
    )

    handle.write(
        f"size_z = {blind_size[2]:.3f}\n"
    )


print()
print(
    "Saved:"
)

print(
    config_file
)

print()
print(
    "Blind-box calculation complete."
)