from pathlib import Path
import math
import csv

import numpy as np
import xlsxwriter
from scipy.optimize import linear_sum_assignment


# ============================================================
# PATHS
# ============================================================

PROJECT = Path(
    r"C:\Users\akidi\Desktop\ISO_AK_Docking"
)

DOCKING = PROJECT / "docking"
LIGAND = PROJECT / "ligand"

STATE_DIR = (
    LIGAND /
    "ISO_protonation_states"
)

OUT_DIR = (
    PROJECT /
    "results" /
    "AK1_ligand_state_sensitivity" /
    "analysis_v1"
)

OUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# OUTPUT FILES
# ============================================================

XLSX_FILE = (
    OUT_DIR /
    "AK1_DIANION_VS_NEUTRAL_POSE_FAMILIES.xlsx"
)

CSV_PER_SEED = (
    OUT_DIR /
    "AK1_dianion_vs_neutral_pose_families_per_seed.csv"
)

CSV_SUMMARY = (
    OUT_DIR /
    "AK1_dianion_vs_neutral_pose_families_summary.csv"
)

CSV_ALL_PAIRWISE = (
    OUT_DIR /
    "AK1_dianion_vs_all_neutral_poses_RMSD.csv"
)


# ============================================================
# SEEDS
# ============================================================

DIANION_SEEDS = [
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

NEUTRAL_DOMINANT_SEEDS = [
    1001,
    2002,
    3003,
    4004,
    5005,
    7007,
    9009,
    10101
]

NEUTRAL_ALTERNATIVE_SEEDS = [
    6006,
    8008
]

ALL_NEUTRAL_SEEDS = (
    NEUTRAL_DOMINANT_SEEDS +
    NEUTRAL_ALTERNATIVE_SEEDS
)


# ============================================================
# READ PDBQT HEAVY ATOMS
# ============================================================

def read_heavy_atoms(
    path,
    model_number=None
):

    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(path)

    lines = path.read_text(
        errors="ignore"
    ).splitlines()

    has_models = any(
        line.startswith("MODEL")
        for line in lines
    )

    current_model = None
    atoms = []

    for line in lines:

        if line.startswith("MODEL"):

            parts = line.split()

            if len(parts) >= 2:

                try:
                    current_model = int(
                        parts[1]
                    )
                except Exception:
                    current_model = None

            continue

        if line.startswith("ENDMDL"):

            current_model = None
            continue

        if not (
            line.startswith("ATOM")
            or line.startswith("HETATM")
        ):
            continue

        if has_models:

            if current_model != model_number:
                continue

        fields = line.split()

        atom_type = fields[-1]

        if atom_type.upper().startswith("H"):
            continue

        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])

        atoms.append({
            "atom_type": atom_type,
            "coords": np.array(
                [x, y, z],
                dtype=float
            )
        })

    if len(atoms) == 0:

        raise RuntimeError(
            f"No heavy atoms found in {path}"
        )

    return atoms


# ============================================================
# READ BEST VINA SCORE
# ============================================================

def read_best_score(path):

    lines = Path(path).read_text(
        errors="ignore"
    ).splitlines()

    for line in lines:

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
# PREPARED LIGANDS
# ============================================================

neutral_input = read_heavy_atoms(
    LIGAND /
    "isorhamnetin.pdbqt"
)

dianion_input = read_heavy_atoms(
    STATE_DIR /
    "isorhamnetin_state2_7O_4primeO_dianion.pdbqt"
)


neutral_start = np.array(
    [
        atom["coords"]
        for atom in neutral_input
    ]
)

dianion_start = np.array(
    [
        atom["coords"]
        for atom in dianion_input
    ]
)


if len(neutral_start) != len(dianion_start):

    raise RuntimeError(
        "Neutral and dianion heavy-atom counts differ."
    )


# ============================================================
# MAP DIANION ATOMS TO NEUTRAL ORDER
# ============================================================

distance_matrix = np.linalg.norm(
    neutral_start[:, None, :]
    -
    dianion_start[None, :, :],
    axis=2
)

neutral_idx, dianion_idx = (
    linear_sum_assignment(
        distance_matrix
    )
)

mapping = np.zeros(
    len(neutral_start),
    dtype=int
)

mapping[
    neutral_idx
] = dianion_idx


matched_distances = (
    distance_matrix[
        neutral_idx,
        dianion_idx
    ]
)

mapping_mean = float(
    np.mean(
        matched_distances
    )
)

mapping_max = float(
    np.max(
        matched_distances
    )
)


print()
print("=" * 80)
print("ATOM MAPPING CHECK")
print("=" * 80)

print(
    f"Mean starting-coordinate difference: "
    f"{mapping_mean:.4f} A"
)

print(
    f"Maximum starting-coordinate difference: "
    f"{mapping_max:.4f} A"
)


if mapping_max > 0.20:

    raise RuntimeError(
        "Dianion atom mapping appears unreliable."
    )


# ============================================================
# READ DOCKED POSES
# ============================================================

def get_pose(
    path,
    reorder=None
):

    atoms = read_heavy_atoms(
        path,
        model_number=1
    )

    coords = np.array(
        [
            atom["coords"]
            for atom in atoms
        ]
    )

    if reorder is not None:

        coords = coords[
            reorder
        ]

    return coords


neutral_poses = {}
neutral_scores = {}

for seed in ALL_NEUTRAL_SEEDS:

    path = (
        DOCKING /
        f"AK1_ISO_seed{seed}.pdbqt"
    )

    neutral_poses[seed] = (
        get_pose(path)
    )

    neutral_scores[seed] = (
        read_best_score(path)
    )


dianion_poses = {}
dianion_scores = {}

for seed in DIANION_SEEDS:

    path = (
        DOCKING /
        f"AK1_ISO_dianion_seed{seed}.pdbqt"
    )

    dianion_poses[seed] = (
        get_pose(
            path,
            reorder=mapping
        )
    )

    dianion_scores[seed] = (
        read_best_score(path)
    )


# ============================================================
# POSITIONAL RMSD
#
# No superposition.
# All poses occupy the same AK1 receptor coordinate system.
# ============================================================

def rmsd(
    a,
    b
):

    return math.sqrt(
        np.mean(
            np.sum(
                (a - b) ** 2,
                axis=1
            )
        )
    )


# ============================================================
# ALL 10 x 10 COMPARISONS
# ============================================================

pairwise_rows = []


for dianion_seed in DIANION_SEEDS:

    dianion_pose = (
        dianion_poses[
            dianion_seed
        ]
    )

    for neutral_seed in ALL_NEUTRAL_SEEDS:

        family = (
            "Dominant"
            if neutral_seed
            in NEUTRAL_DOMINANT_SEEDS
            else
            "Alternative"
        )

        value = rmsd(
            dianion_pose,
            neutral_poses[
                neutral_seed
            ]
        )

        pairwise_rows.append({
            "Dianion_seed":
                dianion_seed,

            "Neutral_seed":
                neutral_seed,

            "Neutral_family":
                family,

            "RMSD_A":
                value
        })


# ============================================================
# PER-SEED FAMILY ASSIGNMENT
# ============================================================

seed_rows = []


for dianion_seed in DIANION_SEEDS:

    relevant = [
        row
        for row in pairwise_rows
        if row[
            "Dianion_seed"
        ] == dianion_seed
    ]

    dominant_rows = [
        row
        for row in relevant
        if row[
            "Neutral_family"
        ] == "Dominant"
    ]

    alternative_rows = [
        row
        for row in relevant
        if row[
            "Neutral_family"
        ] == "Alternative"
    ]

    nearest_dom = min(
        dominant_rows,
        key=lambda row:
            row["RMSD_A"]
    )

    nearest_alt = min(
        alternative_rows,
        key=lambda row:
            row["RMSD_A"]
    )

    rmsd_6006 = next(
        row["RMSD_A"]
        for row in alternative_rows
        if row["Neutral_seed"] == 6006
    )

    rmsd_8008 = next(
        row["RMSD_A"]
        for row in alternative_rows
        if row["Neutral_seed"] == 8008
    )


    if (
        nearest_dom["RMSD_A"]
        <
        nearest_alt["RMSD_A"]
    ):

        assigned_family = (
            "Neutral dominant family"
        )

    elif (
        nearest_alt["RMSD_A"]
        <
        nearest_dom["RMSD_A"]
    ):

        assigned_family = (
            "Neutral alternative family"
        )

    else:

        assigned_family = "Equal"


    seed_rows.append({

        "Dianion_seed":
            dianion_seed,

        "Dianion_Vina_score":
            dianion_scores[
                dianion_seed
            ],

        "Nearest_dominant_seed":
            nearest_dom[
                "Neutral_seed"
            ],

        "Nearest_dominant_RMSD_A":
            nearest_dom[
                "RMSD_A"
            ],

        "Nearest_alternative_seed":
            nearest_alt[
                "Neutral_seed"
            ],

        "Nearest_alternative_RMSD_A":
            nearest_alt[
                "RMSD_A"
            ],

        "RMSD_to_alt6006_A":
            rmsd_6006,

        "RMSD_to_alt8008_A":
            rmsd_8008,

        "Assigned_pose_family":
            assigned_family,

        "Within_2A_of_dominant":
            nearest_dom[
                "RMSD_A"
            ] <= 2.0,

        "Within_3A_of_dominant":
            nearest_dom[
                "RMSD_A"
            ] <= 3.0,

        "Within_2A_of_alternative":
            nearest_alt[
                "RMSD_A"
            ] <= 2.0,

        "Within_3A_of_alternative":
            nearest_alt[
                "RMSD_A"
            ] <= 3.0
    })


# ============================================================
# SUMMARY
# ============================================================

dom_values = np.array(
    [
        row[
            "Nearest_dominant_RMSD_A"
        ]
        for row in seed_rows
    ]
)

alt_values = np.array(
    [
        row[
            "Nearest_alternative_RMSD_A"
        ]
        for row in seed_rows
    ]
)


dom_n_2A = sum(
    value <= 2.0
    for value in dom_values
)

dom_n_3A = sum(
    value <= 3.0
    for value in dom_values
)

alt_n_2A = sum(
    value <= 2.0
    for value in alt_values
)

alt_n_3A = sum(
    value <= 3.0
    for value in alt_values
)


assigned_dom = sum(
    row[
        "Assigned_pose_family"
    ] ==
    "Neutral dominant family"
    for row in seed_rows
)

assigned_alt = sum(
    row[
        "Assigned_pose_family"
    ] ==
    "Neutral alternative family"
    for row in seed_rows
)


summary_rows = [

    {
        "Comparison":
            "Dianion vs neutral dominant family",

        "Mean_nearest_RMSD_A":
            float(np.mean(dom_values)),

        "SD_nearest_RMSD_A":
            float(
                np.std(
                    dom_values,
                    ddof=1
                )
            ),

        "Minimum_RMSD_A":
            float(np.min(dom_values)),

        "Maximum_RMSD_A":
            float(np.max(dom_values)),

        "n_within_2A":
            int(dom_n_2A),

        "Percent_within_2A":
            100.0 *
            dom_n_2A /
            len(dom_values),

        "n_within_3A":
            int(dom_n_3A),

        "Percent_within_3A":
            100.0 *
            dom_n_3A /
            len(dom_values)
    },

    {
        "Comparison":
            "Dianion vs neutral alternative family",

        "Mean_nearest_RMSD_A":
            float(np.mean(alt_values)),

        "SD_nearest_RMSD_A":
            float(
                np.std(
                    alt_values,
                    ddof=1
                )
            ),

        "Minimum_RMSD_A":
            float(np.min(alt_values)),

        "Maximum_RMSD_A":
            float(np.max(alt_values)),

        "n_within_2A":
            int(alt_n_2A),

        "Percent_within_2A":
            100.0 *
            alt_n_2A /
            len(alt_values),

        "n_within_3A":
            int(alt_n_3A),

        "Percent_within_3A":
            100.0 *
            alt_n_3A /
            len(alt_values)
    }
]


# ============================================================
# CSV WRITER
# ============================================================

def write_csv(
    path,
    rows
):

    headers = list(
        rows[0].keys()
    )

    with open(
        path,
        "w",
        newline="",
        encoding="utf-8"
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=headers
        )

        writer.writeheader()

        writer.writerows(
            rows
        )


write_csv(
    CSV_PER_SEED,
    seed_rows
)

write_csv(
    CSV_SUMMARY,
    summary_rows
)

write_csv(
    CSV_ALL_PAIRWISE,
    pairwise_rows
)


# ============================================================
# EXCEL WORKBOOK
# ============================================================

workbook = xlsxwriter.Workbook(
    XLSX_FILE
)


header_fmt = workbook.add_format({
    "bold": True,
    "font_color": "white",
    "bg_color": "#4472C4",
    "border": 1
})

cell_fmt = workbook.add_format({
    "border": 1
})

number_fmt = workbook.add_format({
    "border": 1,
    "num_format": "0.000"
})


def write_sheet(
    name,
    rows,
    width=25
):

    ws = workbook.add_worksheet(
        name
    )

    headers = list(
        rows[0].keys()
    )

    for col, header in enumerate(
        headers
    ):

        ws.write(
            0,
            col,
            header,
            header_fmt
        )

    for row_i, row in enumerate(
        rows,
        start=1
    ):

        for col, header in enumerate(
            headers
        ):

            value = row[header]

            if isinstance(
                value,
                (
                    float,
                    np.floating
                )
            ):

                ws.write_number(
                    row_i,
                    col,
                    float(value),
                    number_fmt
                )

            elif isinstance(
                value,
                (
                    int,
                    np.integer
                )
            ):

                ws.write_number(
                    row_i,
                    col,
                    int(value),
                    cell_fmt
                )

            elif isinstance(
                value,
                (
                    bool,
                    np.bool_
                )
            ):

                ws.write(
                    row_i,
                    col,
                    str(value),
                    cell_fmt
                )

            else:

                ws.write(
                    row_i,
                    col,
                    value,
                    cell_fmt
                )

    ws.freeze_panes(
        1,
        0
    )

    ws.set_column(
        0,
        len(headers) - 1,
        width
    )


write_sheet(
    "Per_Seed",
    seed_rows,
    width=27
)

write_sheet(
    "Summary",
    summary_rows,
    width=30
)

write_sheet(
    "All_Pairwise_RMSD",
    pairwise_rows,
    width=24
)


# ============================================================
# FAMILY ASSIGNMENT SHEET
# ============================================================

family_rows = [

    {
        "Assigned_family":
            "Neutral dominant family",

        "n":
            assigned_dom,

        "Percent":
            assigned_dom /
            len(seed_rows) *
            100.0
    },

    {
        "Assigned_family":
            "Neutral alternative family",

        "n":
            assigned_alt,

        "Percent":
            assigned_alt /
            len(seed_rows) *
            100.0
    }
]


write_sheet(
    "Family_Assignment",
    family_rows,
    width=30
)


# ============================================================
# ATOM MAPPING
# ============================================================

mapping_rows = []

for neutral_i in range(
    len(mapping)
):

    dianion_i = int(
        mapping[
            neutral_i
        ]
    )

    mapping_rows.append({

        "Neutral_atom_index":
            neutral_i + 1,

        "Dianion_atom_index":
            dianion_i + 1,

        "Starting_coordinate_difference_A":
            float(
                matched_distances[
                    neutral_i
                ]
            )
    })


write_sheet(
    "Atom_Mapping",
    mapping_rows,
    width=30
)


# ============================================================
# METHOD NOTES
# ============================================================

notes_sheet = (
    workbook.add_worksheet(
        "Method_Notes"
    )
)


notes = [

    [
        "Purpose",
        (
            "Determine whether the 7-O-/4'-O- dianion "
            "corresponds to the neutral dominant or alternative "
            "AK1 pose family."
        )
    ],

    [
        "Neutral dominant family",
        (
            "Seeds 1001, 2002, 3003, 4004, 5005, "
            "7007, 9009 and 10101."
        )
    ],

    [
        "Neutral alternative family",
        "Seeds 6006 and 8008."
    ],

    [
        "Dianion",
        (
            "Isorhamnetin with 7-OH and 4'-OH "
            "deprotonated; formal charge -2."
        )
    ],

    [
        "RMSD",
        (
            "Direct positional heavy-atom RMSD in the common "
            "AK1 receptor coordinate frame; no structural "
            "superposition was performed."
        )
    ],

    [
        "Atom mapping",
        (
            f"Mean starting-coordinate difference = "
            f"{mapping_mean:.4f} A; maximum = "
            f"{mapping_max:.4f} A."
        )
    ],

    [
        "Interpretation caution",
        (
            "Pose-family assignment describes computational "
            "docking convergence and should not be interpreted "
            "as an experimental ligand population or "
            "thermodynamic state preference."
        )
    ]
]


for row_i, (
    label,
    description
) in enumerate(notes):

    notes_sheet.write(
        row_i,
        0,
        label,
        header_fmt
    )

    notes_sheet.write(
        row_i,
        1,
        description,
        cell_fmt
    )


notes_sheet.set_column(
    "A:A",
    30
)

notes_sheet.set_column(
    "B:B",
    110
)


workbook.close()


# ============================================================
# CONSOLE OUTPUT
# ============================================================

print()
print("=" * 80)
print("DIANION VS NEUTRAL POSE FAMILIES")
print("=" * 80)


for row in seed_rows:

    print()

    print(
        f"Dianion seed "
        f"{row['Dianion_seed']}"
    )

    print(
        f"  nearest dominant    : "
        f"seed {row['Nearest_dominant_seed']}   "
        f"{row['Nearest_dominant_RMSD_A']:.3f} A"
    )

    print(
        f"  nearest alternative : "
        f"seed {row['Nearest_alternative_seed']}   "
        f"{row['Nearest_alternative_RMSD_A']:.3f} A"
    )

    print(
        f"  assigned family     : "
        f"{row['Assigned_pose_family']}"
    )


print()
print("=" * 80)
print("SUMMARY")
print("=" * 80)

print()
print(
    "Dianion -> neutral dominant family"
)

print(
    f"Mean nearest RMSD : "
    f"{np.mean(dom_values):.3f} A"
)

print(
    f"Range             : "
    f"{np.min(dom_values):.3f} - "
    f"{np.max(dom_values):.3f} A"
)

print(
    f"Within 2 A        : "
    f"{dom_n_2A}/10"
)

print(
    f"Within 3 A        : "
    f"{dom_n_3A}/10"
)


print()
print(
    "Dianion -> neutral alternative family"
)

print(
    f"Mean nearest RMSD : "
    f"{np.mean(alt_values):.3f} A"
)

print(
    f"Range             : "
    f"{np.min(alt_values):.3f} - "
    f"{np.max(alt_values):.3f} A"
)

print(
    f"Within 2 A        : "
    f"{alt_n_2A}/10"
)

print(
    f"Within 3 A        : "
    f"{alt_n_3A}/10"
)


print()
print("Family assignment:")

print(
    f"  Dominant family    : "
    f"{assigned_dom}/10"
)

print(
    f"  Alternative family : "
    f"{assigned_alt}/10"
)


print()
print("=" * 80)
print("SAVED")
print("=" * 80)

print(
    XLSX_FILE
)

print(
    CSV_PER_SEED
)

print(
    CSV_SUMMARY
)

print(
    CSV_ALL_PAIRWISE
)

print()
print("Analysis complete.")