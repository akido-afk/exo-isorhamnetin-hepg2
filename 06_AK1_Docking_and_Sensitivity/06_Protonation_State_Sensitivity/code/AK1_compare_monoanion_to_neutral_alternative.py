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

DOCKING = (
    PROJECT /
    "docking"
)

LIGAND = (
    PROJECT /
    "ligand"
)

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

CSV_PER_SEED = (
    OUT_DIR /
    "AK1_monoanion_vs_neutral_pose_families_per_seed.csv"
)

CSV_SUMMARY = (
    OUT_DIR /
    "AK1_monoanion_vs_neutral_pose_families_summary.csv"
)

XLSX_FILE = (
    OUT_DIR /
    "AK1_MONOANION_VS_NEUTRAL_POSE_FAMILIES.xlsx"
)


# ============================================================
# SEEDS
# ============================================================

MONO_SEEDS = [
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

NEUTRAL_ALTERNATIVE_SEEDS = [
    6006,
    8008
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


# ============================================================
# PDBQT READER
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
            or
            line.startswith("HETATM")
        ):
            continue

        if has_models:

            if current_model != model_number:
                continue

        fields = line.split()

        atom_type = fields[-1]

        if atom_type.upper().startswith("H"):
            continue

        x = float(
            line[30:38]
        )

        y = float(
            line[38:46]
        )

        z = float(
            line[46:54]
        )

        atoms.append({
            "atom_type":
                atom_type,

            "coords":
                np.array(
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
# VINA SCORE READER
# ============================================================

def read_best_score(path):

    path = Path(path)

    lines = path.read_text(
        errors="ignore"
    ).splitlines()

    for line in lines:

        if line.startswith(
            "REMARK VINA RESULT:"
        ):

            fields = line.split()

            return float(
                fields[3]
            )

    raise RuntimeError(
        f"No Vina result found in {path}"
    )


# ============================================================
# PREPARED LIGANDS
# ============================================================

neutral_input = read_heavy_atoms(
    LIGAND /
    "isorhamnetin.pdbqt"
)

mono_input = read_heavy_atoms(
    STATE_DIR /
    "isorhamnetin_state1_7O_minus.pdbqt"
)


neutral_start = np.array(
    [
        atom["coords"]
        for atom in neutral_input
    ]
)

mono_start = np.array(
    [
        atom["coords"]
        for atom in mono_input
    ]
)


# ============================================================
# MAP MONOANION ATOMS TO NEUTRAL ORDER
# ============================================================

distance_matrix = np.linalg.norm(

    neutral_start[:, None, :]
    -
    mono_start[None, :, :],

    axis=2
)


neutral_idx, mono_idx = (
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
] = mono_idx


matched_distances = (
    distance_matrix[
        neutral_idx,
        mono_idx
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
print("=" * 78)
print("ATOM MAPPING CHECK")
print("=" * 78)

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
        "Atom mapping appears incorrect."
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


all_neutral_seeds = sorted(
    set(
        NEUTRAL_DOMINANT_SEEDS
        +
        NEUTRAL_ALTERNATIVE_SEEDS
    )
)


for seed in all_neutral_seeds:

    file = (
        DOCKING /
        f"AK1_ISO_seed{seed}.pdbqt"
    )

    neutral_poses[seed] = (
        get_pose(
            file
        )
    )

    neutral_scores[seed] = (
        read_best_score(
            file
        )
    )


mono_poses = {}
mono_scores = {}


for seed in MONO_SEEDS:

    file = (
        DOCKING /
        f"AK1_ISO_7Ominus_seed{seed}.pdbqt"
    )

    mono_poses[seed] = (
        get_pose(
            file,
            reorder=mapping
        )
    )

    mono_scores[seed] = (
        read_best_score(
            file
        )
    )


# ============================================================
# RMSD
#
# No post-docking alignment.
# All ligands share the same AK1 coordinate frame.
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
# PER-SEED COMPARISON
# ============================================================

per_seed_rows = []


print()
print("=" * 78)
print("7-O MONOANION VS NEUTRAL ALTERNATIVE FAMILY")
print("=" * 78)


for seed in MONO_SEEDS:

    pose = (
        mono_poses[seed]
    )

    # --------------------------------------------------------
    # Alternative family
    # --------------------------------------------------------

    alt_results = []

    for neutral_seed in (
        NEUTRAL_ALTERNATIVE_SEEDS
    ):

        value = rmsd(
            pose,
            neutral_poses[
                neutral_seed
            ]
        )

        alt_results.append(
            (
                neutral_seed,
                value
            )
        )


    nearest_alt_seed, nearest_alt_rmsd = min(
        alt_results,
        key=lambda x: x[1]
    )


    rmsd_to_6006 = dict(
        alt_results
    )[6006]

    rmsd_to_8008 = dict(
        alt_results
    )[8008]


    # --------------------------------------------------------
    # Dominant family
    # --------------------------------------------------------

    dominant_results = []

    for neutral_seed in (
        NEUTRAL_DOMINANT_SEEDS
    ):

        value = rmsd(
            pose,
            neutral_poses[
                neutral_seed
            ]
        )

        dominant_results.append(
            (
                neutral_seed,
                value
            )
        )


    nearest_dom_seed, nearest_dom_rmsd = min(
        dominant_results,
        key=lambda x: x[1]
    )


    # --------------------------------------------------------
    # Save row
    # --------------------------------------------------------

    row = {

        "Monoanion_seed":
            seed,

        "Monoanion_Vina_score":
            mono_scores[seed],

        "RMSD_to_neutral_alt6006_A":
            rmsd_to_6006,

        "RMSD_to_neutral_alt8008_A":
            rmsd_to_8008,

        "Nearest_alternative_seed":
            nearest_alt_seed,

        "Nearest_alternative_RMSD_A":
            nearest_alt_rmsd,

        "Within_2A_of_alternative":
            nearest_alt_rmsd <= 2.0,

        "Within_3A_of_alternative":
            nearest_alt_rmsd <= 3.0,

        "Nearest_dominant_seed":
            nearest_dom_seed,

        "Nearest_dominant_RMSD_A":
            nearest_dom_rmsd,

        "Within_2A_of_dominant":
            nearest_dom_rmsd <= 2.0,

        "Within_3A_of_dominant":
            nearest_dom_rmsd <= 3.0
    }


    per_seed_rows.append(
        row
    )


    print(
        f"Seed {seed:<5} | "
        f"to neutral 6006 = {rmsd_to_6006:.3f} A | "
        f"to neutral 8008 = {rmsd_to_8008:.3f} A | "
        f"nearest = {nearest_alt_rmsd:.3f} A"
    )


# ============================================================
# SUMMARY VALUES
# ============================================================

alt_values = np.array(
    [
        row[
            "Nearest_alternative_RMSD_A"
        ]
        for row in per_seed_rows
    ]
)


dom_values = np.array(
    [
        row[
            "Nearest_dominant_RMSD_A"
        ]
        for row in per_seed_rows
    ]
)


alt_n_2A = sum(
    value <= 2.0
    for value in alt_values
)

alt_n_3A = sum(
    value <= 3.0
    for value in alt_values
)

dom_n_2A = sum(
    value <= 2.0
    for value in dom_values
)

dom_n_3A = sum(
    value <= 3.0
    for value in dom_values
)


summary_rows = [

    {
        "Comparison":
            "7-O- monoanion vs neutral alternative family",

        "Mean_nearest_RMSD_A":
            float(
                np.mean(
                    alt_values
                )
            ),

        "SD_nearest_RMSD_A":
            float(
                np.std(
                    alt_values,
                    ddof=1
                )
            ),

        "Minimum_RMSD_A":
            float(
                np.min(
                    alt_values
                )
            ),

        "Maximum_RMSD_A":
            float(
                np.max(
                    alt_values
                )
            ),

        "n_within_2A":
            alt_n_2A,

        "Percent_within_2A":
            100.0 *
            alt_n_2A /
            len(
                alt_values
            ),

        "n_within_3A":
            alt_n_3A,

        "Percent_within_3A":
            100.0 *
            alt_n_3A /
            len(
                alt_values
            )
    },

    {
        "Comparison":
            "7-O- monoanion vs neutral dominant family",

        "Mean_nearest_RMSD_A":
            float(
                np.mean(
                    dom_values
                )
            ),

        "SD_nearest_RMSD_A":
            float(
                np.std(
                    dom_values,
                    ddof=1
                )
            ),

        "Minimum_RMSD_A":
            float(
                np.min(
                    dom_values
                )
            ),

        "Maximum_RMSD_A":
            float(
                np.max(
                    dom_values
                )
            ),

        "n_within_2A":
            dom_n_2A,

        "Percent_within_2A":
            100.0 *
            dom_n_2A /
            len(
                dom_values
            ),

        "n_within_3A":
            dom_n_3A,

        "Percent_within_3A":
            100.0 *
            dom_n_3A /
            len(
                dom_values
            )
    }
]


# ============================================================
# CONSOLE SUMMARY
# ============================================================

print()
print("=" * 78)
print("SUMMARY")
print("=" * 78)

print()
print(
    "Monoanion -> neutral alternative family"
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
print(
    "Monoanion -> neutral dominant family"
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


# ============================================================
# CSV WRITER
# ============================================================

def save_csv(
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


save_csv(
    CSV_PER_SEED,
    per_seed_rows
)

save_csv(
    CSV_SUMMARY,
    summary_rows
)


# ============================================================
# EXCEL
# ============================================================

workbook = xlsxwriter.Workbook(
    XSLX_FILE
    if False
    else XLSX_FILE
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


percent_fmt = workbook.add_format({
    "border": 1,
    "num_format": "0.0"
})


# ============================================================
# SHEET 1 — PER SEED
# ============================================================

ws = workbook.add_worksheet(
    "Per_Seed"
)


headers = list(
    per_seed_rows[0].keys()
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
    per_seed_rows,
    start=1
):

    for col, header in enumerate(
        headers
    ):

        value = row[header]

        if isinstance(
            value,
            bool
        ):

            ws.write(
                row_i,
                col,
                str(value),
                cell_fmt
            )

        elif isinstance(
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
    24
)


# ============================================================
# SHEET 2 — SUMMARY
# ============================================================

ws2 = workbook.add_worksheet(
    "Summary"
)


summary_headers = list(
    summary_rows[0].keys()
)


for col, header in enumerate(
    summary_headers
):

    ws2.write(
        0,
        col,
        header,
        header_fmt
    )


for row_i, row in enumerate(
    summary_rows,
    start=1
):

    for col, header in enumerate(
        summary_headers
    ):

        value = row[header]

        if isinstance(
            value,
            (
                float,
                np.floating
            )
        ):

            ws2.write_number(
                row_i,
                col,
                float(value),
                number_fmt
            )

        else:

            ws2.write(
                row_i,
                col,
                value,
                cell_fmt
            )


ws2.set_column(
    0,
    0,
    52
)

ws2.set_column(
    1,
    len(summary_headers) - 1,
    22
)


# ============================================================
# SHEET 3 — NEUTRAL REFERENCE FAMILY
# ============================================================

ws3 = workbook.add_worksheet(
    "Neutral_Reference"
)


reference_headers = [
    "Neutral_seed",
    "Pose_family",
    "Vina_score"
]


for col, header in enumerate(
    reference_headers
):

    ws3.write(
        0,
        col,
        header,
        header_fmt
    )


reference_rows = []


for seed in all_neutral_seeds:

    if seed in (
        NEUTRAL_ALTERNATIVE_SEEDS
    ):

        family = (
            "Alternative"
        )

    else:

        family = (
            "Dominant"
        )

    reference_rows.append({
        "Neutral_seed":
            seed,

        "Pose_family":
            family,

        "Vina_score":
            neutral_scores[seed]
    })


for row_i, row in enumerate(
    reference_rows,
    start=1
):

    ws3.write(
        row_i,
        0,
        row["Neutral_seed"],
        cell_fmt
    )

    ws3.write(
        row_i,
        1,
        row["Pose_family"],
        cell_fmt
    )

    ws3.write_number(
        row_i,
        2,
        row["Vina_score"],
        number_fmt
    )


ws3.set_column(
    "A:C",
    22
)


# ============================================================
# SHEET 4 — ATOM MAPPING
# ============================================================

ws4 = workbook.add_worksheet(
    "Atom_Mapping"
)


mapping_headers = [
    "Neutral_heavy_atom_index",
    "Monoanion_heavy_atom_index",
    "Starting_coordinate_difference_A"
]


for col, header in enumerate(
    mapping_headers
):

    ws4.write(
        0,
        col,
        header,
        header_fmt
    )


for i in range(
    len(mapping)
):

    ws4.write(
        i + 1,
        0,
        i + 1,
        cell_fmt
    )

    ws4.write(
        i + 1,
        1,
        int(
            mapping[i]
        ) + 1,
        cell_fmt
    )

    ws4.write_number(
        i + 1,
        2,
        float(
            matched_distances[i]
        ),
        number_fmt
    )


ws4.set_column(
    "A:C",
    30
)


# ============================================================
# SHEET 5 — INTERPRETATION / METHODS
# ============================================================

ws5 = workbook.add_worksheet(
    "Method_Notes"
)


notes = [

    [
        "Main result",
        (
            "All 10 best poses of the 7-O- monoanion reproduce "
            "the neutral alternative pose family."
        )
    ],

    [
        "Alternative-family RMSD",
        (
            f"Mean nearest RMSD = "
            f"{np.mean(alt_values):.3f} A; "
            f"range = "
            f"{np.min(alt_values):.3f}-"
            f"{np.max(alt_values):.3f} A."
        )
    ],

    [
        "Alternative-family recurrence",
        (
            f"{alt_n_2A}/10 monoanion poses were within 2 A "
            f"of the neutral alternative family."
        )
    ],

    [
        "Dominant-family RMSD",
        (
            f"Mean nearest RMSD to the neutral dominant family = "
            f"{np.mean(dom_values):.3f} A."
        )
    ],

    [
        "Interpretation",
        (
            "The 7-O- monoanion does not generate an arbitrary "
            "new orientation. Instead, protonation-state change "
            "selects the alternative pose family that was already "
            "sampled by neutral isorhamnetin."
        )
    ],

    [
        "RMSD method",
        (
            "Direct heavy-atom positional RMSD was calculated "
            "without post-docking alignment because all poses "
            "were generated within the same AK1 receptor "
            "coordinate frame."
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
        "Caution",
        (
            "The analysis demonstrates docking-pose sensitivity "
            "to ligand protonation state. It does not establish "
            "experimental populations, binding affinities, or "
            "biochemical mechanism."
        )
    ]
]


for row_i, (
    label,
    description
) in enumerate(notes):

    ws5.write(
        row_i,
        0,
        label,
        header_fmt
    )

    ws5.write(
        row_i,
        1,
        description,
        cell_fmt
    )


ws5.set_column(
    "A:A",
    32
)

ws5.set_column(
    "B:B",
    110
)


# ============================================================
# CLOSE WORKBOOK
# ============================================================

workbook.close()


# ============================================================
# OUTPUT SUMMARY
# ============================================================

print()
print("=" * 78)
print("FILES SAVED")
print("=" * 78)

print()
print(
    "Per-seed CSV:"
)

print(
    CSV_PER_SEED
)

print()
print(
    "Summary CSV:"
)

print(
    CSV_SUMMARY
)

print()
print(
    "Excel workbook:"
)

print(
    XLSX_FILE
)

print()
print("=" * 78)