from pathlib import Path
from collections import defaultdict
import subprocess
import sys
import math
import csv

from rdkit import Chem
import xlsxwriter


# ============================================================
# PROJECT PATHS
# ============================================================

PROJECT = Path(r"C:\Users\akidi\Desktop\ISO_AK_Docking")

DOCKING = PROJECT / "docking"
RECEPTORS = PROJECT / "receptors"

OUTPUT_DIR = (
    PROJECT
    / "results"
    / "AK1_dominant8_contact_frequency"
)

POSE_DIR = OUTPUT_DIR / "poses"

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

POSE_DIR.mkdir(
    parents=True,
    exist_ok=True
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
# DISTANCE CUTOFFS
# ============================================================

ANY_CONTACT_CUTOFF = 4.0
POLAR_CONTACT_CUTOFF = 3.6


# ============================================================
# INPUT PROTEIN
# ============================================================

PROTEIN_FILE = (
    RECEPTORS
    / "AK1_chainA_FIXED.pdb"
)


# ============================================================
# MEEKO EXPORT SCRIPT
# ============================================================

MEEKO_EXPORT = (
    Path(sys.prefix)
    / "Lib"
    / "site-packages"
    / "meeko"
    / "cli"
    / "mk_export.py"
)


# ============================================================
# BASIC DISTANCE FUNCTION
# ============================================================

def distance(a, b):

    return math.sqrt(
        (a[0] - b[0]) ** 2
        + (a[1] - b[1]) ** 2
        + (a[2] - b[2]) ** 2
    )


# ============================================================
# READ AK1 PROTEIN ATOMS
# ============================================================

protein_atoms = []

with open(
    PROTEIN_FILE,
    "r"
) as f:

    for line in f:

        if not line.startswith("ATOM"):
            continue

        atom_name = line[12:16].strip()
        resname = line[17:20].strip()
        chain = line[21].strip()
        resid = line[22:26].strip()

        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])

        element = line[76:78].strip()

        if not element:
            element = atom_name[0]

        # ignore hydrogen
        if element.upper() == "H":
            continue

        protein_atoms.append({
            "atom": atom_name,
            "resname": resname,
            "resid": resid,
            "chain": chain,
            "element": element.upper(),
            "coord": (x, y, z)
        })


# ============================================================
# EXPORT / LOAD THE BEST POSE FROM EACH SEED
# ============================================================

def get_best_pose(seed):

    pdbqt = (
        DOCKING
        / f"AK1_ISO_seed{seed}.pdbqt"
    )

    all_sdf = (
        POSE_DIR
        / f"AK1_ISO_seed{seed}_all.sdf"
    )

    best_sdf = (
        POSE_DIR
        / f"AK1_ISO_seed{seed}_best.sdf"
    )

    if not pdbqt.exists():

        raise FileNotFoundError(
            f"Missing docking file: {pdbqt}"
        )

    # Export all Vina poses if needed
    if not all_sdf.exists():

        print(
            f"Exporting seed {seed}..."
        )

        subprocess.run(
            [
                sys.executable,
                str(MEEKO_EXPORT),
                str(pdbqt),
                "-s",
                str(all_sdf)
            ],
            check=True
        )

    supplier = Chem.SDMolSupplier(
        str(all_sdf),
        removeHs=False
    )

    mol = None

    for candidate in supplier:

        if candidate is not None:
            mol = candidate
            break

    if mol is None:

        raise RuntimeError(
            f"Could not read pose for seed {seed}"
        )

    # Save first pose separately for provenance
    if not best_sdf.exists():

        writer = Chem.SDWriter(
            str(best_sdf)
        )

        writer.write(mol)
        writer.close()

    return mol


# ============================================================
# ANALYZE ONE ISO POSE
# ============================================================

def analyze_pose(mol):

    conf = mol.GetConformer()

    # residue -> shortest heavy-atom distance
    any_contacts = {}

    # residue -> shortest polar-heavy-atom distance
    polar_contacts = {}

    for atom in mol.GetAtoms():

        element = atom.GetSymbol().upper()

        if element == "H":
            continue

        idx = atom.GetIdx()

        pos = conf.GetAtomPosition(idx)

        iso_coord = (
            float(pos.x),
            float(pos.y),
            float(pos.z)
        )

        for p in protein_atoms:

            d = distance(
                iso_coord,
                p["coord"]
            )

            residue = (
                f"{p['resname']}{p['resid']}"
            )

            # --------------------------------------------
            # ANY HEAVY-ATOM CONTACT
            # --------------------------------------------

            if d <= ANY_CONTACT_CUTOFF:

                if (
                    residue not in any_contacts
                    or d < any_contacts[residue]
                ):

                    any_contacts[residue] = d

            # --------------------------------------------
            # POLAR HEAVY-ATOM CONTACT
            # --------------------------------------------

            if (
                element in {"O", "N", "S"}
                and
                p["element"] in {"O", "N", "S"}
                and
                d <= POLAR_CONTACT_CUTOFF
            ):

                if (
                    residue not in polar_contacts
                    or d < polar_contacts[residue]
                ):

                    polar_contacts[residue] = d

    return any_contacts, polar_contacts


# ============================================================
# RUN ALL 8 DOMINANT POSES
# ============================================================

all_results = {}

any_frequency = defaultdict(int)
polar_frequency = defaultdict(int)

any_distances = defaultdict(list)
polar_distances = defaultdict(list)


for seed in SEEDS:

    print()
    print("=" * 70)
    print(
        f"ANALYZING DOMINANT POSE — SEED {seed}"
    )
    print("=" * 70)

    mol = get_best_pose(seed)

    any_contacts, polar_contacts = (
        analyze_pose(mol)
    )

    all_results[seed] = {
        "any": any_contacts,
        "polar": polar_contacts
    }

    for residue, d in any_contacts.items():

        any_frequency[residue] += 1
        any_distances[residue].append(d)

    for residue, d in polar_contacts.items():

        polar_frequency[residue] += 1
        polar_distances[residue].append(d)

    print(
        "Heavy-atom contacting residues:",
        len(any_contacts)
    )

    print(
        "Polar-contact residues:",
        len(polar_contacts)
    )


# ============================================================
# ALL RESIDUES OBSERVED
# ============================================================

all_residues = sorted(
    set(any_frequency.keys())
    | set(polar_frequency.keys()),
    key=lambda r: (
        int(
            ''.join(
                c for c in r
                if c.isdigit()
            )
        )
    )
)


# ============================================================
# BUILD SUMMARY
# ============================================================

summary = []

for residue in all_residues:

    any_n = any_frequency.get(
        residue,
        0
    )

    polar_n = polar_frequency.get(
        residue,
        0
    )

    any_mean = (
        sum(any_distances[residue])
        / len(any_distances[residue])
        if any_distances.get(residue)
        else None
    )

    any_min = (
        min(any_distances[residue])
        if any_distances.get(residue)
        else None
    )

    polar_mean = (
        sum(polar_distances[residue])
        / len(polar_distances[residue])
        if polar_distances.get(residue)
        else None
    )

    polar_min = (
        min(polar_distances[residue])
        if polar_distances.get(residue)
        else None
    )

    summary.append({
        "Residue": residue,

        "Any contact n": any_n,
        "Any contact %": (
            100 * any_n / len(SEEDS)
        ),

        "Any mean min-distance": any_mean,
        "Any shortest distance": any_min,

        "Polar contact n": polar_n,
        "Polar contact %": (
            100 * polar_n / len(SEEDS)
        ),

        "Polar mean min-distance": polar_mean,
        "Polar shortest distance": polar_min
    })


# Sort by:
# 1) polar recurrence
# 2) any-contact recurrence
# 3) shortest distance
summary = sorted(
    summary,
    key=lambda x: (
        -x["Polar contact n"],
        -x["Any contact n"],
        x["Any shortest distance"]
        if x["Any shortest distance"] is not None
        else 999
    )
)


# ============================================================
# PRINT IMPORTANT RESULTS
# ============================================================

print()
print("=" * 90)
print(
    "AK1 — DOMINANT 8-POSE CONTACT RECURRENCE"
)
print("=" * 90)

print()
print(
    f"{'Residue':<12}"
    f"{'Any':<12}"
    f"{'Polar':<12}"
    f"{'Polar %':<12}"
    f"{'Min polar dist.'}"
)

print("-" * 90)


for row in summary:

    if row["Any contact n"] == 0:
        continue

    polar_dist_text = (
        f"{row['Polar shortest distance']:.3f} A"
        if row["Polar shortest distance"]
        is not None
        else "-"
    )

    print(
        f"{row['Residue']:<12}"
        f"{row['Any contact n']}/8{'':<8}"
        f"{row['Polar contact n']}/8{'':<6}"
        f"{row['Polar contact %']:<12.1f}"
        f"{polar_dist_text}"
    )


# ============================================================
# SAVE CSV
# ============================================================

CSV_FILE = (
    OUTPUT_DIR
    / "AK1_dominant8_contact_frequency.csv"
)

with open(
    CSV_FILE,
    "w",
    newline="",
    encoding="utf-8"
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=list(
            summary[0].keys()
        )
    )

    writer.writeheader()
    writer.writerows(summary)


# ============================================================
# CREATE EXCEL
# ============================================================

XLSX_FILE = (
    OUTPUT_DIR
    / "AK1_dominant8_contact_frequency.xlsx"
)

workbook = xlsxwriter.Workbook(
    XLSX_FILE
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
    "align": "center",
    "text_wrap": True
})

normal_fmt = workbook.add_format({
    "border": 1
})

num_fmt = workbook.add_format({
    "border": 1,
    "num_format": "0.000"
})

percent_fmt = workbook.add_format({
    "border": 1,
    "num_format": "0.0"
})

strong_fmt = workbook.add_format({
    "border": 1,
    "bg_color": "#E2F0D9"
})


# ============================================================
# SHEET 1 — FREQUENCY SUMMARY
# ============================================================

ws = workbook.add_worksheet(
    "Frequency_Summary"
)

headers = list(
    summary[0].keys()
)

ws.merge_range(
    0,
    0,
    0,
    len(headers) - 1,
    "AK1–Isorhamnetin Dominant 8-Pose Contact Recurrence",
    title_fmt
)

for col, header in enumerate(headers):

    ws.write(
        2,
        col,
        header,
        header_fmt
    )


for row_index, row in enumerate(
    summary,
    start=3
):

    strongly_recurrent = (
        row["Polar contact n"] >= 6
        or
        row["Any contact n"] == 8
    )

    for col, header in enumerate(headers):

        value = row[header]

        fmt = (
            strong_fmt
            if strongly_recurrent
            else normal_fmt
        )

        if value is None:

            ws.write(
                row_index,
                col,
                "",
                fmt
            )

        elif isinstance(value, float):

            if "%" in header:

                ws.write_number(
                    row_index,
                    col,
                    value,
                    percent_fmt
                )

            else:

                ws.write_number(
                    row_index,
                    col,
                    value,
                    num_fmt
                )

        else:

            ws.write(
                row_index,
                col,
                value,
                fmt
            )


ws.set_column(
    0,
    0,
    15
)

ws.set_column(
    1,
    len(headers) - 1,
    22
)

ws.freeze_panes(
    3,
    0
)


# ============================================================
# SHEET 2 — SEED-BY-SEED PRESENCE MATRIX
# ============================================================

ws2 = workbook.add_worksheet(
    "Seed_Matrix"
)

matrix_headers = [
    "Residue"
]

for seed in SEEDS:

    matrix_headers.append(
        f"{seed} any"
    )

for seed in SEEDS:

    matrix_headers.append(
        f"{seed} polar"
    )


ws2.merge_range(
    0,
    0,
    0,
    len(matrix_headers) - 1,
    "Contact Presence Across the Eight Dominant AK1–ISO Poses",
    title_fmt
)

for col, header in enumerate(
    matrix_headers
):

    ws2.write(
        2,
        col,
        header,
        header_fmt
    )


for row_index, residue in enumerate(
    all_residues,
    start=3
):

    ws2.write(
        row_index,
        0,
        residue,
        normal_fmt
    )

    col = 1

    for seed in SEEDS:

        present = int(
            residue in all_results[
                seed
            ]["any"]
        )

        ws2.write_number(
            row_index,
            col,
            present,
            normal_fmt
        )

        col += 1

    for seed in SEEDS:

        present = int(
            residue in all_results[
                seed
            ]["polar"]
        )

        ws2.write_number(
            row_index,
            col,
            present,
            normal_fmt
        )

        col += 1


ws2.set_column(
    0,
    0,
    14
)

ws2.set_column(
    1,
    len(matrix_headers) - 1,
    12
)


# ============================================================
# SHEET 3 — METHOD
# ============================================================

ws3 = workbook.add_worksheet(
    "Method"
)

method_rows = [
    [
        "Dominant cluster seeds",
        ", ".join(
            str(x)
            for x in SEEDS
        )
    ],

    [
        "Number of poses",
        8
    ],

    [
        "Any-contact definition",
        "At least one AK1 heavy atom within 4.0 Å of an isorhamnetin heavy atom."
    ],

    [
        "Polar-contact definition",
        "At least one AK1 O/N/S atom within 3.6 Å of an isorhamnetin O/N/S atom."
    ],

    [
        "Interpretation",
        "Recurrence indicates repeated geometric engagement across independently generated members of the dominant docking pose family."
    ],

    [
        "Important caution",
        "Polar-contact recurrence does not by itself establish a hydrogen bond; donor/acceptor identity and geometry should be evaluated separately."
    ]
]


ws3.write(
    0,
    0,
    "Parameter",
    header_fmt
)

ws3.write(
    0,
    1,
    "Value",
    header_fmt
)


for row_index, values in enumerate(
    method_rows,
    start=1
):

    ws3.write(
        row_index,
        0,
        values[0],
        normal_fmt
    )

    ws3.write(
        row_index,
        1,
        values[1],
        normal_fmt
    )


ws3.set_column(
    0,
    0,
    28
)

ws3.set_column(
    1,
    1,
    100
)


workbook.close()


# ============================================================
# FINISHED
# ============================================================

print()
print("=" * 90)
print("FILES SAVED")
print("=" * 90)

print()
print(CSV_FILE)
print(XLSX_FILE)