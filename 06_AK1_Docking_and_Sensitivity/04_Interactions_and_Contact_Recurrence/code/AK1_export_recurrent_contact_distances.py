from pathlib import Path
from collections import defaultdict
import math
import csv

from rdkit import Chem
import xlsxwriter


# ============================================================
# PATHS
# ============================================================

PROJECT = Path(r"C:\Users\akidi\Desktop\ISO_AK_Docking")

PROTEIN_FILE = PROJECT / "receptors" / "AK1_chainA_FIXED.pdb"

ANALYSIS_DIR = (
    PROJECT
    / "results"
    / "AK1_dominant8_contact_frequency"
)

POSE_DIR = ANALYSIS_DIR / "poses"

FREQUENCY_FILE = (
    ANALYSIS_DIR
    / "AK1_dominant8_contact_frequency.csv"
)

CSV_OUTPUT = (
    ANALYSIS_DIR
    / "AK1_dominant8_recurrent_polar_distances_long.csv"
)

XLSX_OUTPUT = (
    ANALYSIS_DIR
    / "AK1_dominant8_recurrent_polar_distances.xlsx"
)


# ============================================================
# DOMINANT CLUSTER
# ============================================================

SEEDS = [
    1001,
    2002,
    3003,
    4004,
    5005,
    7007,
    9009,
    10101
]


# ============================================================
# DISTANCE FUNCTION
# ============================================================

def distance(a, b):

    return math.sqrt(
        (a[0] - b[0]) ** 2
        + (a[1] - b[1]) ** 2
        + (a[2] - b[2]) ** 2
    )


# ============================================================
# IDENTIFY THE 8/8 POLAR-CONTACT RESIDUES
# FROM OUR PREVIOUS REAL ANALYSIS
# ============================================================

recurrent_residues = []

with open(
    FREQUENCY_FILE,
    newline="",
    encoding="utf-8"
) as f:

    reader = csv.DictReader(f)

    for row in reader:

        if int(row["Polar contact n"]) == 8:

            recurrent_residues.append(
                row["Residue"]
            )


print()
print("Recurrent 8/8 polar-contact residues:")
print(recurrent_residues)


# ============================================================
# READ PROTEIN POLAR ATOMS
# ============================================================

protein_atoms = defaultdict(list)

with open(
    PROTEIN_FILE,
    "r"
) as f:

    for line in f:

        if not line.startswith("ATOM"):
            continue

        atom_name = line[12:16].strip()
        resname = line[17:20].strip()
        resid = line[22:26].strip()

        element = line[76:78].strip()

        if not element:
            element = atom_name[0]

        element = element.upper()

        if element not in {"O", "N", "S"}:
            continue

        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])

        residue = f"{resname}{resid}"

        protein_atoms[residue].append({
            "atom": atom_name,
            "element": element,
            "coord": (x, y, z)
        })


# ============================================================
# CALCULATE MINIMUM POLAR DISTANCE
# FOR EACH RESIDUE IN EACH SEED
# ============================================================

rows = []

for seed in SEEDS:

    sdf_file = (
        POSE_DIR
        / f"AK1_ISO_seed{seed}_best.sdf"
    )

    if not sdf_file.exists():

        raise FileNotFoundError(
            f"Missing pose file: {sdf_file}"
        )

    mol = Chem.MolFromMolFile(
        str(sdf_file),
        removeHs=False
    )

    if mol is None:

        raise RuntimeError(
            f"Could not read: {sdf_file}"
        )

    conf = mol.GetConformer()

    iso_polar_atoms = []

    for atom in mol.GetAtoms():

        element = atom.GetSymbol().upper()

        if element not in {"O", "N", "S"}:
            continue

        idx = atom.GetIdx()

        pos = conf.GetAtomPosition(idx)

        iso_polar_atoms.append({
            "index": idx + 1,
            "element": element,
            "coord": (
                float(pos.x),
                float(pos.y),
                float(pos.z)
            )
        })

    for residue in recurrent_residues:

        best_distance = None
        best_protein_atom = None
        best_iso_atom = None

        for p in protein_atoms[residue]:

            for iso_atom in iso_polar_atoms:

                d = distance(
                    p["coord"],
                    iso_atom["coord"]
                )

                if (
                    best_distance is None
                    or d < best_distance
                ):

                    best_distance = d
                    best_protein_atom = p["atom"]
                    best_iso_atom = (
                        f"{iso_atom['element']}"
                        f"{iso_atom['index']}"
                    )

        rows.append({
            "Seed": seed,
            "Residue": residue,
            "Minimum polar distance (A)": best_distance,
            "Protein atom": best_protein_atom,
            "ISO atom": best_iso_atom,
            "Within 3.6 A cutoff": (
                "YES"
                if best_distance <= 3.6
                else "NO"
            )
        })


# ============================================================
# PRINT TABLE
# ============================================================

print()
print("=" * 95)
print("AK1 DOMINANT-CLUSTER RECURRENT POLAR CONTACT DISTANCES")
print("=" * 95)

print(
    f"{'Seed':<10}"
    f"{'Residue':<12}"
    f"{'Distance':<14}"
    f"{'Protein atom':<15}"
    f"{'ISO atom':<12}"
)

print("-" * 95)

for row in rows:

    print(
        f"{row['Seed']:<10}"
        f"{row['Residue']:<12}"
        f"{row['Minimum polar distance (A)']:<14.3f}"
        f"{row['Protein atom']:<15}"
        f"{row['ISO atom']:<12}"
    )


# ============================================================
# SAVE CSV
# ============================================================

columns = [
    "Seed",
    "Residue",
    "Minimum polar distance (A)",
    "Protein atom",
    "ISO atom",
    "Within 3.6 A cutoff"
]

with open(
    CSV_OUTPUT,
    "w",
    newline="",
    encoding="utf-8"
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=columns
    )

    writer.writeheader()
    writer.writerows(rows)


# ============================================================
# SAVE EXCEL
# ============================================================

workbook = xlsxwriter.Workbook(
    XLSX_OUTPUT
)

ws = workbook.add_worksheet(
    "Per_Seed_Distances"
)

title_fmt = workbook.add_format({
    "bold": True,
    "font_size": 15,
    "font_color": "white",
    "bg_color": "#1F4E78"
})

header_fmt = workbook.add_format({
    "bold": True,
    "font_color": "white",
    "bg_color": "#4472C4",
    "border": 1,
    "align": "center"
})

normal_fmt = workbook.add_format({
    "border": 1
})

distance_fmt = workbook.add_format({
    "border": 1,
    "num_format": "0.000"
})


ws.merge_range(
    "A1:F1",
    "AK1 Dominant-Cluster Recurrent Polar Contact Distances",
    title_fmt
)

for col, header in enumerate(columns):

    ws.write(
        2,
        col,
        header,
        header_fmt
    )


for row_index, row in enumerate(
    rows,
    start=3
):

    ws.write_number(
        row_index,
        0,
        row["Seed"],
        normal_fmt
    )

    ws.write(
        row_index,
        1,
        row["Residue"],
        normal_fmt
    )

    ws.write_number(
        row_index,
        2,
        row["Minimum polar distance (A)"],
        distance_fmt
    )

    ws.write(
        row_index,
        3,
        row["Protein atom"],
        normal_fmt
    )

    ws.write(
        row_index,
        4,
        row["ISO atom"],
        normal_fmt
    )

    ws.write(
        row_index,
        5,
        row["Within 3.6 A cutoff"],
        normal_fmt
    )


ws.set_column("A:A", 12)
ws.set_column("B:B", 14)
ws.set_column("C:C", 25)
ws.set_column("D:F", 18)

workbook.close()


print()
print("=" * 95)
print("FILES SAVED")
print("=" * 95)

print(CSV_OUTPUT)
print(XLSX_OUTPUT)