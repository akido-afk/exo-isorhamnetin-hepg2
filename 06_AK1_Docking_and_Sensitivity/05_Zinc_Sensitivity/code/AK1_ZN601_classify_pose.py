from pathlib import Path
import subprocess
import sys
import math

from rdkit import Chem
import xlsxwriter


# ============================================================
# PATHS
# ============================================================

PROJECT = Path(r"C:\Users\akidi\Desktop\ISO_AK_Docking")

DOCKING = PROJECT / "docking"

OUTPUT_DIR = (
    PROJECT
    / "results"
    / "AK1_ZN_sensitivity"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

POSE_DIR = OUTPUT_DIR / "original_10_reference_poses"

POSE_DIR.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# ORIGINAL SEEDS / CLUSTER MEMBERSHIP
# ============================================================

SEEDS = [
    1001,
    2002,
    3003,
    4004,
    5005,
    6006,
    7007,
    8008,
    9009,
    10101
]

DOMINANT = {
    1001,
    2002,
    3003,
    4004,
    5005,
    7007,
    9009,
    10101
}

ALTERNATIVE = {
    6006,
    8008
}


# ============================================================
# Zn601 PILOT POSE
# ============================================================

ZN_POSE_FILE = (
    OUTPUT_DIR
    / "AK1_ISO_seed9009_ZN601_best.sdf"
)


# ============================================================
# MEEKO EXPORT
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
# EXPORT FIRST/BEST ORIGINAL POSE
# ============================================================

def get_original_pose(seed):

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

    if not all_sdf.exists():

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
            f"Could not read seed {seed}"
        )

    if not best_sdf.exists():

        writer = Chem.SDWriter(
            str(best_sdf)
        )

        writer.write(mol)
        writer.close()

    return mol


# ============================================================
# READ Zn601 POSE
# ============================================================

zn_pose = Chem.MolFromMolFile(
    str(ZN_POSE_FILE),
    removeHs=False
)

if zn_pose is None:

    raise RuntimeError(
        "Could not read ZN601 pilot pose."
    )


# ============================================================
# POSITIONAL RMSD
# No superposition: all poses are already in same AK1 frame
# ============================================================

def positional_rmsd(mol1, mol2):

    if (
        mol1.GetNumAtoms()
        !=
        mol2.GetNumAtoms()
    ):

        raise RuntimeError(
            "Different atom counts."
        )

    conf1 = mol1.GetConformer()
    conf2 = mol2.GetConformer()

    values = []

    for i in range(
        mol1.GetNumAtoms()
    ):

        a1 = mol1.GetAtomWithIdx(i)
        a2 = mol2.GetAtomWithIdx(i)

        if (
            a1.GetAtomicNum()
            !=
            a2.GetAtomicNum()
        ):

            raise RuntimeError(
                f"Atom mismatch at {i}"
            )

        if a1.GetAtomicNum() == 1:
            continue

        p1 = conf1.GetAtomPosition(i)
        p2 = conf2.GetAtomPosition(i)

        values.append(
            (p1.x - p2.x) ** 2
            +
            (p1.y - p2.y) ** 2
            +
            (p1.z - p2.z) ** 2
        )

    return math.sqrt(
        sum(values)
        /
        len(values)
    )


# ============================================================
# COMPARE AGAINST ALL ORIGINAL POSES
# ============================================================

rows = []

for seed in SEEDS:

    original = get_original_pose(seed)

    rmsd = positional_rmsd(
        zn_pose,
        original
    )

    if seed in DOMINANT:

        cluster = "Dominant"

    else:

        cluster = "Alternative"

    rows.append({
        "Seed": seed,
        "Original cluster": cluster,
        "RMSD to ZN601 pose (A)": rmsd
    })


rows.sort(
    key=lambda x:
        x["RMSD to ZN601 pose (A)"]
)


# ============================================================
# CLUSTER SUMMARIES
# ============================================================

dominant_values = [
    r["RMSD to ZN601 pose (A)"]
    for r in rows
    if r["Original cluster"] == "Dominant"
]

alternative_values = [
    r["RMSD to ZN601 pose (A)"]
    for r in rows
    if r["Original cluster"] == "Alternative"
]


dominant_mean = (
    sum(dominant_values)
    /
    len(dominant_values)
)

alternative_mean = (
    sum(alternative_values)
    /
    len(alternative_values)
)


nearest = rows[0]


if (
    alternative_mean
    <
    dominant_mean
):

    closest_family = (
        "Original alternative cluster"
    )

else:

    closest_family = (
        "Original dominant cluster"
    )


# ============================================================
# PRINT
# ============================================================

print()
print("=" * 80)
print("ZN601 PILOT POSE vs ORIGINAL 10 Zn-FREE POSES")
print("=" * 80)

print()

print(
    f"{'Seed':<10}"
    f"{'Cluster':<15}"
    f"{'RMSD (A)':<12}"
)

print("-" * 50)

for row in rows:

    print(
        f"{row['Seed']:<10}"
        f"{row['Original cluster']:<15}"
        f"{row['RMSD to ZN601 pose (A)']:<12.3f}"
    )


print()
print(
    "Mean RMSD to dominant cluster:",
    f"{dominant_mean:.3f} A"
)

print(
    "Mean RMSD to alternative cluster:",
    f"{alternative_mean:.3f} A"
)

print()

print(
    "Nearest original pose:",
    f"seed {nearest['Seed']} "
    f"({nearest['Original cluster']}) "
    f"RMSD = "
    f"{nearest['RMSD to ZN601 pose (A)']:.3f} A"
)

print()

print(
    "Closest original pose family:",
    closest_family
)


# ============================================================
# EXCEL
# ============================================================

XLSX = (
    OUTPUT_DIR
    / "AK1_ZN601_pilot_pose_classification.xlsx"
)

workbook = xlsxwriter.Workbook(
    XLSX
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

number_fmt = workbook.add_format({
    "border": 1,
    "num_format": "0.000"
})

dominant_fmt = workbook.add_format({
    "border": 1,
    "bg_color": "#FFF2CC"
})

alternative_fmt = workbook.add_format({
    "border": 1,
    "bg_color": "#DDEBF7"
})


ws = workbook.add_worksheet(
    "Pose_Comparison"
)

ws.merge_range(
    "A1:C1",
    "AK1 ZN601 Pilot Pose vs Original Zn-Free Pose Ensemble",
    title_fmt
)

headers = [
    "Seed",
    "Original cluster",
    "RMSD to ZN601 pose (A)"
]

for col, h in enumerate(headers):

    ws.write(
        2,
        col,
        h,
        header_fmt
    )


for row_i, row in enumerate(
    rows,
    start=3
):

    fmt = (
        dominant_fmt
        if row["Original cluster"]
        == "Dominant"
        else alternative_fmt
    )

    ws.write_number(
        row_i,
        0,
        row["Seed"],
        fmt
    )

    ws.write(
        row_i,
        1,
        row["Original cluster"],
        fmt
    )

    ws.write_number(
        row_i,
        2,
        row["RMSD to ZN601 pose (A)"],
        number_fmt
    )


ws.set_column(
    "A:A",
    12
)

ws.set_column(
    "B:B",
    20
)

ws.set_column(
    "C:C",
    24
)


# ============================================================
# SUMMARY SHEET
# ============================================================

ws2 = workbook.add_worksheet(
    "Summary"
)

summary_rows = [
    [
        "Mean RMSD to original dominant cluster",
        dominant_mean
    ],

    [
        "Mean RMSD to original alternative cluster",
        alternative_mean
    ],

    [
        "Nearest original seed",
        nearest["Seed"]
    ],

    [
        "Nearest original cluster",
        nearest["Original cluster"]
    ],

    [
        "Nearest-pose RMSD (A)",
        nearest[
            "RMSD to ZN601 pose (A)"
        ]
    ],

    [
        "Closest original pose family",
        closest_family
    ]
]


for row_i, values in enumerate(
    summary_rows
):

    ws2.write(
        row_i,
        0,
        values[0],
        header_fmt
    )

    if isinstance(
        values[1],
        float
    ):

        ws2.write_number(
            row_i,
            1,
            values[1],
            number_fmt
        )

    else:

        ws2.write(
            row_i,
            1,
            values[1],
            normal_fmt
        )


ws2.set_column(
    "A:A",
    45
)

ws2.set_column(
    "B:B",
    40
)


workbook.close()


print()
print("=" * 80)
print("EXCEL SAVED")
print("=" * 80)

print(XLSX)