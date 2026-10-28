from pathlib import Path
import re
import math
import csv

PROJECT = Path(r"C:\Users\akidi\Desktop\ISO_AK_Docking")
DOCKING = PROJECT / "docking"
RESULTS = PROJECT / "results"

FILES = sorted(
    DOCKING.glob("AK1_ISO_seed*.pdbqt"),
    key=lambda p: int(re.search(r"seed(\d+)", p.name).group(1))
)


def read_best_pose_heavy_atoms(path):
    coords = []
    in_model_1 = False
    model_seen = False

    with open(path, "r") as f:
        for line in f:

            if line.startswith("MODEL"):
                model_number = int(line.split()[1])

                if model_number == 1:
                    in_model_1 = True
                    model_seen = True
                    continue

                elif in_model_1:
                    break

            if line.startswith("ENDMDL") and in_model_1:
                break

            # Some Vina files may not explicitly use MODEL
            if not model_seen and (
                line.startswith("ATOM") or line.startswith("HETATM")
            ):
                in_model_1 = True

            if in_model_1 and (
                line.startswith("ATOM") or line.startswith("HETATM")
            ):

                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                # AutoDock atom type is normally the last field
                atom_type = line.split()[-1].upper()

                # exclude hydrogens
                if atom_type.startswith("H"):
                    continue

                coords.append((x, y, z))

    return coords


def rmsd(coords1, coords2):

    if len(coords1) != len(coords2):
        raise ValueError(
            f"Different atom counts: {len(coords1)} vs {len(coords2)}"
        )

    squared = 0.0

    for a, b in zip(coords1, coords2):

        dx = a[0] - b[0]
        dy = a[1] - b[1]
        dz = a[2] - b[2]

        squared += dx*dx + dy*dy + dz*dz

    return math.sqrt(squared / len(coords1))


poses = {}

for file in FILES:

    seed = int(re.search(r"seed(\d+)", file.name).group(1))

    poses[seed] = read_best_pose_heavy_atoms(file)

    print(
        f"Seed {seed}: "
        f"{len(poses[seed])} heavy atoms"
    )


seeds = sorted(poses)

matrix = []

all_pairwise = []

for seed1 in seeds:

    row = []

    for seed2 in seeds:

        value = rmsd(
            poses[seed1],
            poses[seed2]
        )

        row.append(value)

        if seed2 > seed1:
            all_pairwise.append(value)

    matrix.append(row)


# ------------------------------------------------------------
# Average RMSD of each pose to all other poses
# ------------------------------------------------------------

average_to_others = {}

for i, seed in enumerate(seeds):

    values = [
        matrix[i][j]
        for j in range(len(seeds))
        if i != j
    ]

    average_to_others[seed] = sum(values) / len(values)


representative_seed = min(
    average_to_others,
    key=average_to_others.get
)


# ------------------------------------------------------------
# Output CSV matrix
# ------------------------------------------------------------

matrix_file = RESULTS / "AK1_ISO_pose_RMSD_matrix.csv"

with open(matrix_file, "w", newline="") as f:

    writer = csv.writer(f)

    writer.writerow(
        ["Seed"] + seeds
    )

    for seed, row in zip(seeds, matrix):

        writer.writerow(
            [seed] +
            [round(x, 3) for x in row]
        )


# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

mean_pairwise = sum(all_pairwise) / len(all_pairwise)

minimum_pairwise = min(all_pairwise)
maximum_pairwise = max(all_pairwise)


print()
print("=" * 60)
print("AK1 ISO POSE CONVERGENCE")
print("=" * 60)

print(
    "Number of best poses:",
    len(seeds)
)

print(
    "Pairwise comparisons:",
    len(all_pairwise)
)

print(
    f"Mean pairwise RMSD: {mean_pairwise:.3f} A"
)

print(
    f"Minimum pairwise RMSD: {minimum_pairwise:.3f} A"
)

print(
    f"Maximum pairwise RMSD: {maximum_pairwise:.3f} A"
)

print()

print(
    "Most representative pose:"
)

print(
    f"Seed {representative_seed}"
)

print(
    f"Mean RMSD to other poses: "
    f"{average_to_others[representative_seed]:.3f} A"
)

print()

print(
    "RMSD matrix saved to:"
)

print(matrix_file)