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
RECEPTORS = PROJECT / "receptors"

OUTPUT_DIR = (
    PROJECT
    / "results"
    / "AK1_ZN_sensitivity"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)


ZNFREE_PDBQT = (
    DOCKING
    / "AK1_ISO_seed9009.pdbqt"
)

ZN601_PDBQT = (
    DOCKING
    / "AK1_ZN601_ISO_seed9009.pdbqt"
)

ZN_RECEPTOR_PDB = (
    RECEPTORS
    / "AK1_chainA_FIXED_ZN601.pdb"
)


# ============================================================
# OUTPUT FILES
# ============================================================

ZNFREE_ALL_SDF = (
    OUTPUT_DIR
    / "AK1_ISO_seed9009_ZNfree_allposes.sdf"
)

ZN601_ALL_SDF = (
    OUTPUT_DIR
    / "AK1_ISO_seed9009_ZN601_allposes.sdf"
)

ZNFREE_BEST_SDF = (
    OUTPUT_DIR
    / "AK1_ISO_seed9009_ZNfree_best.sdf"
)

ZN601_BEST_SDF = (
    OUTPUT_DIR
    / "AK1_ISO_seed9009_ZN601_best.sdf"
)

XLSX_OUTPUT = (
    OUTPUT_DIR
    / "AK1_ZN601_seed9009_pilot_comparison.xlsx"
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
# READ BEST VINA SCORE
# ============================================================

def get_best_score(path):

    with open(
        path,
        "r",
        encoding="utf-8"
    ) as f:

        for line in f:

            if line.startswith(
                "REMARK VINA RESULT:"
            ):

                return float(
                    line.split()[3]
                )

    raise RuntimeError(
        f"No Vina score found in {path}"
    )


# ============================================================
# EXPORT PDBQT TO SDF
# ============================================================

def export_to_sdf(
    pdbqt,
    output_sdf
):

    subprocess.run(
        [
            sys.executable,
            str(MEEKO_EXPORT),
            str(pdbqt),
            "-s",
            str(output_sdf)
        ],
        check=True
    )


# ============================================================
# READ FIRST / BEST POSE
# ============================================================

def read_best_pose(
    all_sdf,
    best_sdf
):

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
            f"No valid molecule in {all_sdf}"
        )


    writer = Chem.SDWriter(
        str(best_sdf)
    )

    writer.write(mol)
    writer.close()

    return mol


# ============================================================
# EXPORT BOTH POSES
# ============================================================

print()
print("Exporting Zn-free pose...")

export_to_sdf(
    ZNFREE_PDBQT,
    ZNFREE_ALL_SDF
)

print("Exporting ZN601-retained pose...")

export_to_sdf(
    ZN601_PDBQT,
    ZN601_ALL_SDF
)


znfree = read_best_pose(
    ZNFREE_ALL_SDF,
    ZNFREE_BEST_SDF
)

zn601 = read_best_pose(
    ZN601_ALL_SDF,
    ZN601_BEST_SDF
)


# ============================================================
# VERIFY ATOM ORDER
# ============================================================

if znfree.GetNumAtoms() != zn601.GetNumAtoms():

    raise RuntimeError(
        "The two poses have different atom counts."
    )


for i in range(
    znfree.GetNumAtoms()
):

    a = znfree.GetAtomWithIdx(i)
    b = zn601.GetAtomWithIdx(i)

    if (
        a.GetAtomicNum()
        !=
        b.GetAtomicNum()
    ):

        raise RuntimeError(
            f"Atom mismatch at index {i}"
        )


# ============================================================
# POSITIONAL RMSD
#
# IMPORTANT:
# No ligand superposition is performed.
# Both poses are already in the same AK1 coordinate frame.
# ============================================================

conf1 = znfree.GetConformer()
conf2 = zn601.GetConformer()

squared = []
heavy_atom_count = 0


for i in range(
    znfree.GetNumAtoms()
):

    atom = znfree.GetAtomWithIdx(i)

    if atom.GetAtomicNum() == 1:
        continue

    p1 = conf1.GetAtomPosition(i)
    p2 = conf2.GetAtomPosition(i)

    d2 = (
        (p1.x - p2.x) ** 2
        +
        (p1.y - p2.y) ** 2
        +
        (p1.z - p2.z) ** 2
    )

    squared.append(d2)

    heavy_atom_count += 1


pose_rmsd = math.sqrt(
    sum(squared)
    /
    len(squared)
)


# ============================================================
# READ ZN601 COORDINATE
# ============================================================

zn_coord = None

with open(
    ZN_RECEPTOR_PDB,
    "r"
) as f:

    for line in f:

        if (
            line.startswith("HETATM")
            and
            line[17:20].strip() == "ZN"
            and
            line[22:26].strip() == "601"
        ):

            zn_coord = (
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54])
            )

            break


if zn_coord is None:

    raise RuntimeError(
        "Could not locate ZN601."
    )


# ============================================================
# ZN601 -> ISO MINIMUM DISTANCES
# ============================================================

conf = zn601.GetConformer()

min_all = None
min_all_atom = None

min_polar = None
min_polar_atom = None


for atom in zn601.GetAtoms():

    if atom.GetAtomicNum() == 1:
        continue

    i = atom.GetIdx()

    pos = conf.GetAtomPosition(i)

    d = math.sqrt(
        (pos.x - zn_coord[0]) ** 2
        +
        (pos.y - zn_coord[1]) ** 2
        +
        (pos.z - zn_coord[2]) ** 2
    )


    if (
        min_all is None
        or d < min_all
    ):

        min_all = d
        min_all_atom = (
            f"{atom.GetSymbol()}{i+1}"
        )


    if (
        atom.GetSymbol().upper()
        in {"O", "N", "S"}
    ):

        if (
            min_polar is None
            or d < min_polar
        ):

            min_polar = d
            min_polar_atom = (
                f"{atom.GetSymbol()}{i+1}"
            )


# ============================================================
# SCORES
# ============================================================

score_znfree = get_best_score(
    ZNFREE_PDBQT
)

score_zn601 = get_best_score(
    ZN601_PDBQT
)

score_difference = (
    score_zn601
    -
    score_znfree
)


# ============================================================
# PRINT RESULTS
# ============================================================

print()
print("=" * 75)
print(
    "AK1 ZN601 SENSITIVITY — SEED 9009 PILOT"
)
print("=" * 75)

print()
print(
    f"Zn-free score:        "
    f"{score_znfree:.3f} kcal/mol"
)

print(
    f"ZN601-retained score: "
    f"{score_zn601:.3f} kcal/mol"
)

print(
    f"Score difference:      "
    f"{score_difference:+.3f} kcal/mol"
)

print()

print(
    f"Heavy-atom positional RMSD: "
    f"{pose_rmsd:.3f} A"
)

print()

print(
    f"Closest ISO atom to ZN601: "
    f"{min_all_atom} at {min_all:.3f} A"
)

print(
    f"Closest ISO polar atom to ZN601: "
    f"{min_polar_atom} at {min_polar:.3f} A"
)


# ============================================================
# INTERPRETATION
# ============================================================

if pose_rmsd < 1.0:

    pose_interpretation = (
        "Very similar pose"
    )

elif pose_rmsd < 2.0:

    pose_interpretation = (
        "Similar pose"
    )

elif pose_rmsd < 3.0:

    pose_interpretation = (
        "Moderately changed pose"
    )

else:

    pose_interpretation = (
        "Substantially different pose"
    )


# ============================================================
# EXCEL
# ============================================================

workbook = xlsxwriter.Workbook(
    XLSX_OUTPUT
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
    "border": 1
})

normal_fmt = workbook.add_format({
    "border": 1
})

number_fmt = workbook.add_format({
    "border": 1,
    "num_format": "0.000"
})


ws = workbook.add_worksheet(
    "Pilot_Comparison"
)


ws.merge_range(
    "A1:C1",
    "AK1 ZN601 Sensitivity Pilot — Seed 9009",
    title_fmt
)


headers = [
    "Metric",
    "Value",
    "Interpretation"
]


for col, h in enumerate(headers):

    ws.write(
        2,
        col,
        h,
        header_fmt
    )


rows = [

    [
        "Zn-free best Vina score",
        score_znfree,
        "Original AK1 receptor"
    ],

    [
        "ZN601-retained best Vina score",
        score_zn601,
        "AK1 receptor with crystallographic ZN601"
    ],

    [
        "Score difference",
        score_difference,
        "ZN601 minus Zn-free"
    ],

    [
        "Pose heavy-atom positional RMSD",
        pose_rmsd,
        pose_interpretation
    ],

    [
        "Closest ISO heavy atom to ZN601",
        min_all,
        min_all_atom
    ],

    [
        "Closest ISO polar atom to ZN601",
        min_polar,
        min_polar_atom
    ],

    [
        "Heavy atoms used for RMSD",
        heavy_atom_count,
        ""
    ]
]


for r, values in enumerate(
    rows,
    start=3
):

    ws.write(
        r,
        0,
        values[0],
        normal_fmt
    )

    if isinstance(
        values[1],
        float
    ):

        ws.write_number(
            r,
            1,
            values[1],
            number_fmt
        )

    else:

        ws.write(
            r,
            1,
            values[1],
            normal_fmt
        )

    ws.write(
        r,
        2,
        values[2],
        normal_fmt
    )


ws.set_column(
    "A:A",
    38
)

ws.set_column(
    "B:B",
    20
)

ws.set_column(
    "C:C",
    55
)


# ============================================================
# METHOD NOTES
# ============================================================

ws2 = workbook.add_worksheet(
    "Method_Notes"
)

notes = [

    [
        "Comparison",
        "Zn-free AK1 versus AK1 with crystallographic ZN601 retained."
    ],

    [
        "Seed",
        "9009"
    ],

    [
        "Docking box",
        "Identical to the original AK1 docking experiment."
    ],

    [
        "Exhaustiveness",
        "16"
    ],

    [
        "Number of modes",
        "10"
    ],

    [
        "RMSD definition",
        "Direct positional RMSD of ligand heavy atoms in the common AK1 coordinate frame; no ligand superposition was performed."
    ],

    [
        "Interpretation caution",
        "The retained-Zn run is a receptor-state sensitivity analysis using standard Vina scoring, not a dedicated metallodocking free-energy calculation."
    ]
]


for r, values in enumerate(
    notes
):

    ws2.write(
        r,
        0,
        values[0],
        header_fmt
    )

    ws2.write(
        r,
        1,
        values[1],
        normal_fmt
    )


ws2.set_column(
    "A:A",
    30
)

ws2.set_column(
    "B:B",
    110
)


workbook.close()


print()
print("=" * 75)
print("FILES SAVED")
print("=" * 75)

print()
print(XLSX_OUTPUT)

print()
print(ZNFREE_BEST_SDF)
print(ZN601_BEST_SDF)