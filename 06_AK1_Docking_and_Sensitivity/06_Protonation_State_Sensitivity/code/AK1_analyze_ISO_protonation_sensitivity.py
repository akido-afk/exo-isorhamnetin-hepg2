from pathlib import Path
from collections import Counter
import csv
import math

import numpy as np
import xlsxwriter

from scipy.optimize import linear_sum_assignment
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform


# ============================================================
# PROJECT PATHS
# ============================================================

PROJECT = Path(
    r"C:\Users\akidi\Desktop\ISO_AK_Docking"
)

DOCKING = PROJECT / "docking"
LIGAND = PROJECT / "ligand"

STATE_LIGAND_DIR = (
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
# SEEDS
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


NEUTRAL_REPRESENTATIVE_SEED = 9009


# ============================================================
# STATE DEFINITIONS
# ============================================================

STATES = {

    "Neutral": {
        "input_pdbqt":
            LIGAND /
            "isorhamnetin.pdbqt",

        "dock_template":
            "AK1_ISO_seed{seed}.pdbqt",

        "charge":
            0
    },

    "7-O- monoanion": {
        "input_pdbqt":
            STATE_LIGAND_DIR /
            "isorhamnetin_state1_7O_minus.pdbqt",

        "dock_template":
            "AK1_ISO_7Ominus_seed{seed}.pdbqt",

        "charge":
            -1
    },

    "7-O-/4'-O- dianion": {
        "input_pdbqt":
            STATE_LIGAND_DIR /
            "isorhamnetin_state2_7O_4primeO_dianion.pdbqt",

        "dock_template":
            "AK1_ISO_dianion_seed{seed}.pdbqt",

        "charge":
            -2
    }
}


# ============================================================
# PDBQT PARSING
# ============================================================

def is_hydrogen_atom_type(atom_type):

    atom_type = atom_type.strip().upper()

    return atom_type.startswith("H")


def read_pdbqt_atoms(
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

            if model_number is None:
                continue

            if current_model != model_number:
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

        except Exception as exc:

            raise RuntimeError(
                f"Could not parse coordinates in:\n"
                f"{path}\n"
                f"Line:\n{line}"
            ) from exc

        fields = line.split()

        atom_type = fields[-1]

        if is_hydrogen_atom_type(
            atom_type
        ):
            continue

        atom_name = (
            line[12:16].strip()
        )

        atoms.append({
            "atom_name":
                atom_name,

            "atom_type":
                atom_type,

            "x":
                x,

            "y":
                y,

            "z":
                z
        })

    if len(atoms) == 0:

        raise RuntimeError(
            f"No heavy atoms found in {path}"
        )

    return atoms


def atom_coordinates(atoms):

    return np.array(
        [
            [
                atom["x"],
                atom["y"],
                atom["z"]
            ]
            for atom in atoms
        ],
        dtype=float
    )


def read_best_vina_score(path):

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
        f"No Vina score found in {path}"
    )


# ============================================================
# READ PREPARED INPUT LIGANDS
#
# These share the same original heavy-atom coordinates.
# Therefore we can establish the atom correspondence
# objectively before looking at docked poses.
# ============================================================

input_atoms = {}
input_coords = {}


for state_name, info in STATES.items():

    atoms = read_pdbqt_atoms(
        info["input_pdbqt"]
    )

    coords = atom_coordinates(
        atoms
    )

    input_atoms[state_name] = atoms
    input_coords[state_name] = coords


# ============================================================
# VALIDATE HEAVY-ATOM COUNTS
# ============================================================

heavy_atom_counts = {
    state:
        len(atoms)
    for state, atoms
    in input_atoms.items()
}


print()
print("=" * 80)
print("PREPARED LIGAND HEAVY-ATOM COUNTS")
print("=" * 80)

for state, count in heavy_atom_counts.items():

    print(
        f"{state:<23} : {count}"
    )


if len(
    set(
        heavy_atom_counts.values()
    )
) != 1:

    raise RuntimeError(
        "The protonation states do not contain the same "
        "number of heavy atoms."
    )


N_HEAVY = next(
    iter(
        heavy_atom_counts.values()
    )
)


# ============================================================
# BUILD ATOM CORRESPONDENCE TO NEUTRAL ISO
#
# Heavy atoms originated from the same parent structure.
# Coordinate matching therefore avoids relying on potentially
# different PDBQT atom ordering between protonation states.
# ============================================================

neutral_input_coords = (
    input_coords["Neutral"]
)


def build_coordinate_mapping(
    reference_coords,
    target_coords
):

    distance_matrix = np.linalg.norm(
        reference_coords[:, None, :]
        -
        target_coords[None, :, :],
        axis=2
    )

    ref_indices, target_indices = (
        linear_sum_assignment(
            distance_matrix
        )
    )

    mapping = np.zeros(
        len(reference_coords),
        dtype=int
    )

    mapping[
        ref_indices
    ] = target_indices

    matched_distances = (
        distance_matrix[
            ref_indices,
            target_indices
        ]
    )

    return (
        mapping,
        matched_distances
    )


state_mappings = {}
mapping_distances = {}


for state_name in STATES:

    if state_name == "Neutral":

        state_mappings[state_name] = (
            np.arange(
                N_HEAVY,
                dtype=int
            )
        )

        mapping_distances[state_name] = (
            np.zeros(
                N_HEAVY
            )
        )

    else:

        mapping, distances = (
            build_coordinate_mapping(
                neutral_input_coords,
                input_coords[state_name]
            )
        )

        state_mappings[state_name] = mapping
        mapping_distances[state_name] = distances


# ============================================================
# VALIDATE COORDINATE-BASED MAPPING
# ============================================================

print()
print("=" * 80)
print("HEAVY-ATOM MAPPING VALIDATION")
print("=" * 80)


for state_name in STATES:

    distances = (
        mapping_distances[state_name]
    )

    print(
        f"{state_name:<23} "
        f"mean starting-coordinate difference = "
        f"{np.mean(distances):.4f} A   "
        f"max = {np.max(distances):.4f} A"
    )

    if np.max(distances) > 0.20:

        raise RuntimeError(
            f"Atom mapping for {state_name} is suspicious: "
            f"maximum starting-coordinate mismatch "
            f"{np.max(distances):.3f} A"
        )


# ============================================================
# READ ALL DOCKED BEST POSES
# ============================================================

scores = {}
poses = {}


for state_name, info in STATES.items():

    scores[state_name] = {}
    poses[state_name] = {}

    input_types = [
        atom["atom_type"]
        for atom in input_atoms[state_name]
    ]

    mapping = (
        state_mappings[state_name]
    )

    for seed in SEEDS:

        dock_file = (
            DOCKING /
            info["dock_template"].format(
                seed=seed
            )
        )

        if not dock_file.exists():

            raise FileNotFoundError(
                dock_file
            )

        score = read_best_vina_score(
            dock_file
        )

        best_atoms = read_pdbqt_atoms(
            dock_file,
            model_number=1
        )

        if len(best_atoms) != N_HEAVY:

            raise RuntimeError(
                f"{dock_file.name}: "
                f"expected {N_HEAVY} heavy atoms, "
                f"found {len(best_atoms)}."
            )

        output_types = [
            atom["atom_type"]
            for atom in best_atoms
        ]

        if output_types != input_types:

            print()
            print(
                "WARNING: PDBQT atom-type order differs "
                "between prepared ligand and docking output:"
            )

            print(
                dock_file
            )

        coords_state_order = (
            atom_coordinates(
                best_atoms
            )
        )

        # Reorder every protonation state into the
        # neutral heavy-atom correspondence.
        coords_neutral_order = (
            coords_state_order[
                mapping
            ]
        )

        scores[state_name][seed] = score

        poses[state_name][seed] = (
            coords_neutral_order
        )


# ============================================================
# RMSD FUNCTION
#
# IMPORTANT:
# No structural superposition.
# All ligands were docked into the same AK1 coordinate frame.
# ============================================================

def positional_rmsd(
    coords_a,
    coords_b
):

    difference = (
        coords_a -
        coords_b
    )

    squared_distance = np.sum(
        difference ** 2,
        axis=1
    )

    return math.sqrt(
        np.mean(
            squared_distance
        )
    )


# ============================================================
# PAIRWISE RMSD MATRICES
# ============================================================

pairwise_matrices = {}


for state_name in STATES:

    matrix = np.zeros(
        (
            len(SEEDS),
            len(SEEDS)
        ),
        dtype=float
    )

    for i, seed_i in enumerate(SEEDS):

        for j, seed_j in enumerate(SEEDS):

            matrix[i, j] = (
                positional_rmsd(
                    poses[state_name][seed_i],
                    poses[state_name][seed_j]
                )
            )

    pairwise_matrices[state_name] = (
        matrix
    )


# ============================================================
# COMPLETE-LINKAGE CLUSTERING
# ============================================================

def cluster_pose_matrix(
    matrix,
    threshold
):

    n = matrix.shape[0]

    if n == 1:

        return np.array(
            [1],
            dtype=int
        )

    condensed = squareform(
        matrix,
        checks=False
    )

    hierarchy = linkage(
        condensed,
        method="complete"
    )

    labels = fcluster(
        hierarchy,
        t=threshold,
        criterion="distance"
    )

    return labels


cluster_results = {}


for state_name in STATES:

    matrix = (
        pairwise_matrices[state_name]
    )

    labels_2A = cluster_pose_matrix(
        matrix,
        2.0
    )

    labels_3A = cluster_pose_matrix(
        matrix,
        3.0
    )

    cluster_results[state_name] = {
        "2A":
            labels_2A,

        "3A":
            labels_3A
    }


# ============================================================
# REPRESENTATIVE / MEDOID OF LARGEST 2 A CLUSTER
# ============================================================

representatives = {}


for state_name in STATES:

    labels = (
        cluster_results[state_name]["2A"]
    )

    matrix = (
        pairwise_matrices[state_name]
    )

    counts = Counter(
        labels
    )

    largest_cluster_label = (
        counts.most_common(1)[0][0]
    )

    member_indices = [
        i
        for i, label
        in enumerate(labels)
        if label == largest_cluster_label
    ]

    best_index = None
    best_mean_rmsd = None

    for i in member_indices:

        distances = [
            matrix[i, j]
            for j in member_indices
            if j != i
        ]

        if len(distances) == 0:

            mean_distance = 0.0

        else:

            mean_distance = float(
                np.mean(
                    distances
                )
            )

        if (
            best_mean_rmsd is None
            or
            mean_distance < best_mean_rmsd
        ):

            best_mean_rmsd = (
                mean_distance
            )

            best_index = i

    representative_seed = (
        SEEDS[best_index]
    )

    representatives[state_name] = {
        "seed":
            representative_seed,

        "largest_cluster_n":
            len(member_indices),

        "mean_rmsd_within_cluster":
            best_mean_rmsd,

        "cluster_label":
            largest_cluster_label
    }


# ============================================================
# SCORE / CLUSTER SUMMARY
# ============================================================

state_summary_rows = []


for state_name, info in STATES.items():

    state_scores = np.array(
        [
            scores[state_name][seed]
            for seed in SEEDS
        ],
        dtype=float
    )

    labels_2A = (
        cluster_results[state_name]["2A"]
    )

    labels_3A = (
        cluster_results[state_name]["3A"]
    )

    rep = (
        representatives[state_name]
    )

    state_summary_rows.append({
        "State":
            state_name,

        "Formal_charge":
            info["charge"],

        "n":
            len(state_scores),

        "Mean_score":
            float(
                np.mean(
                    state_scores
                )
            ),

        "Median_score":
            float(
                np.median(
                    state_scores
                )
            ),

        "SD_score":
            float(
                np.std(
                    state_scores,
                    ddof=1
                )
            ),

        "Best_score":
            float(
                np.min(
                    state_scores
                )
            ),

        "Worst_score":
            float(
                np.max(
                    state_scores
                )
            ),

        "Score_range":
            float(
                np.max(
                    state_scores
                )
                -
                np.min(
                    state_scores
                )
            ),

        "Clusters_2A":
            len(
                set(
                    labels_2A
                )
            ),

        "Largest_2A_cluster_n":
            rep[
                "largest_cluster_n"
            ],

        "Clusters_3A":
            len(
                set(
                    labels_3A
                )
            ),

        "Representative_seed":
            rep[
                "seed"
            ],

        "Representative_mean_RMSD":
            rep[
                "mean_rmsd_within_cluster"
            ]
    })


# ============================================================
# SEED-LEVEL SCORE TABLE
# ============================================================

score_rows = []


for state_name, info in STATES.items():

    labels_2A = (
        cluster_results[state_name]["2A"]
    )

    labels_3A = (
        cluster_results[state_name]["3A"]
    )

    for i, seed in enumerate(SEEDS):

        score_rows.append({
            "State":
                state_name,

            "Formal_charge":
                info["charge"],

            "Seed":
                seed,

            "Best_Vina_score":
                scores[state_name][seed],

            "Cluster_2A":
                int(
                    labels_2A[i]
                ),

            "Cluster_3A":
                int(
                    labels_3A[i]
                )
        })


# ============================================================
# CROSS-STATE RMSD
#
# Reference 1:
#   neutral representative seed 9009
#
# Reference 2:
#   closest member of the established
#   eight-pose neutral dominant family
# ============================================================

neutral_rep_pose = (
    poses["Neutral"][
        NEUTRAL_REPRESENTATIVE_SEED
    ]
)


cross_state_rows = []


for state_name in STATES:

    for seed in SEEDS:

        query_pose = (
            poses[state_name][seed]
        )

        rmsd_to_rep = (
            positional_rmsd(
                query_pose,
                neutral_rep_pose
            )
        )

        dominant_rmsds = []

        for neutral_seed in (
            NEUTRAL_DOMINANT_SEEDS
        ):

            value = positional_rmsd(
                query_pose,
                poses["Neutral"][
                    neutral_seed
                ]
            )

            dominant_rmsds.append(
                (
                    neutral_seed,
                    value
                )
            )

        closest_seed, closest_rmsd = min(
            dominant_rmsds,
            key=lambda x: x[1]
        )

        cross_state_rows.append({
            "State":
                state_name,

            "Seed":
                seed,

            "RMSD_to_neutral_rep9009_A":
                rmsd_to_rep,

            "Nearest_neutral_dominant_seed":
                closest_seed,

            "Minimum_RMSD_to_neutral_dominant_A":
                closest_rmsd,

            "Within_2A_of_neutral_dominant":
                closest_rmsd <= 2.0,

            "Within_3A_of_neutral_dominant":
                closest_rmsd <= 3.0
        })


# ============================================================
# CROSS-STATE SUMMARY
# ============================================================

cross_summary_rows = []


for state_name in STATES:

    rows = [
        row
        for row in cross_state_rows
        if row["State"] == state_name
    ]

    rep_values = np.array(
        [
            row[
                "RMSD_to_neutral_rep9009_A"
            ]
            for row in rows
        ]
    )

    closest_values = np.array(
        [
            row[
                "Minimum_RMSD_to_neutral_dominant_A"
            ]
            for row in rows
        ]
    )

    n_2A = sum(
        row[
            "Within_2A_of_neutral_dominant"
        ]
        for row in rows
    )

    n_3A = sum(
        row[
            "Within_3A_of_neutral_dominant"
        ]
        for row in rows
    )

    cross_summary_rows.append({
        "State":
            state_name,

        "Mean_RMSD_to_rep9009_A":
            float(
                np.mean(
                    rep_values
                )
            ),

        "Min_RMSD_to_rep9009_A":
            float(
                np.min(
                    rep_values
                )
            ),

        "Max_RMSD_to_rep9009_A":
            float(
                np.max(
                    rep_values
                )
            ),

        "Mean_min_RMSD_to_neutral_dominant_A":
            float(
                np.mean(
                    closest_values
                )
            ),

        "Min_of_min_RMSD_A":
            float(
                np.min(
                    closest_values
                )
            ),

        "Max_of_min_RMSD_A":
            float(
                np.max(
                    closest_values
                )
            ),

        "n_within_2A_of_neutral_dominant":
            int(
                n_2A
            ),

        "Percent_within_2A":
            100.0 *
            n_2A /
            len(rows),

        "n_within_3A_of_neutral_dominant":
            int(
                n_3A
            ),

        "Percent_within_3A":
            100.0 *
            n_3A /
            len(rows)
    })


# ============================================================
# CSV WRITER
# ============================================================

def write_dict_csv(
    path,
    rows
):

    if not rows:
        return

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


# ============================================================
# SAVE CSV TABLES
# ============================================================

write_dict_csv(
    OUT_DIR /
    "AK1_ISO_protonation_state_scores.csv",
    score_rows
)

write_dict_csv(
    OUT_DIR /
    "AK1_ISO_protonation_state_summary.csv",
    state_summary_rows
)

write_dict_csv(
    OUT_DIR /
    "AK1_ISO_cross_state_RMSD.csv",
    cross_state_rows
)

write_dict_csv(
    OUT_DIR /
    "AK1_ISO_cross_state_RMSD_summary.csv",
    cross_summary_rows
)


# ============================================================
# SAVE RMSD MATRICES AS CSV
# ============================================================

for state_name, matrix in (
    pairwise_matrices.items()
):

    safe_name = (
        state_name
        .replace(" ", "_")
        .replace("/", "_")
        .replace("'", "")
    )

    path = (
        OUT_DIR /
        f"AK1_ISO_{safe_name}_pairwise_RMSD.csv"
    )

    with open(
        path,
        "w",
        newline="",
        encoding="utf-8"
    ) as handle:

        writer = csv.writer(
            handle
        )

        writer.writerow(
            ["Seed"] + SEEDS
        )

        for i, seed in enumerate(SEEDS):

            writer.writerow(
                [seed]
                +
                [
                    float(x)
                    for x in matrix[i]
                ]
            )


# ============================================================
# EXCEL WORKBOOK
# ============================================================

xlsx_file = (
    OUT_DIR /
    "AK1_ISO_PROTONATION_STATE_SENSITIVITY.xlsx"
)


workbook = xlsxwriter.Workbook(
    xlsx_file
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


def write_table_sheet(
    workbook,
    name,
    rows
):

    ws = workbook.add_worksheet(
        name
    )

    if not rows:
        return

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
        22
    )


write_table_sheet(
    workbook,
    "Scores",
    score_rows
)

write_table_sheet(
    workbook,
    "State_Summary",
    state_summary_rows
)

write_table_sheet(
    workbook,
    "Cross_State_RMSD",
    cross_state_rows
)

write_table_sheet(
    workbook,
    "Cross_RMSD_Summary",
    cross_summary_rows
)


# ============================================================
# MAPPING VALIDATION SHEET
# ============================================================

mapping_sheet = (
    workbook.add_worksheet(
        "Atom_Mapping"
    )
)

mapping_headers = [
    "State",
    "Neutral_heavy_atom_index",
    "State_heavy_atom_index",
    "Neutral_atom_type",
    "State_atom_type",
    "Starting_coordinate_difference_A"
]


for col, header in enumerate(
    mapping_headers
):

    mapping_sheet.write(
        0,
        col,
        header,
        header_fmt
    )


excel_row = 1


for state_name in STATES:

    mapping = (
        state_mappings[state_name]
    )

    distances = (
        mapping_distances[state_name]
    )

    for neutral_idx in range(
        N_HEAVY
    ):

        state_idx = int(
            mapping[neutral_idx]
        )

        values = [
            state_name,
            neutral_idx + 1,
            state_idx + 1,
            input_atoms[
                "Neutral"
            ][neutral_idx][
                "atom_type"
            ],
            input_atoms[
                state_name
            ][state_idx][
                "atom_type"
            ],
            float(
                distances[
                    neutral_idx
                ]
            )
        ]

        for col, value in enumerate(
            values
        ):

            if isinstance(
                value,
                float
            ):

                mapping_sheet.write_number(
                    excel_row,
                    col,
                    value,
                    number_fmt
                )

            else:

                mapping_sheet.write(
                    excel_row,
                    col,
                    value,
                    cell_fmt
                )

        excel_row += 1


mapping_sheet.set_column(
    "A:F",
    28
)


# ============================================================
# RMSD MATRIX SHEETS
# ============================================================

for state_name, matrix in (
    pairwise_matrices.items()
):

    sheet_name = {
        "Neutral":
            "RMSD_Neutral",

        "7-O- monoanion":
            "RMSD_7Ominus",

        "7-O-/4'-O- dianion":
            "RMSD_Dianion"
    }[state_name]

    ws = workbook.add_worksheet(
        sheet_name
    )

    ws.write(
        0,
        0,
        "Seed",
        header_fmt
    )

    for col, seed in enumerate(
        SEEDS,
        start=1
    ):

        ws.write(
            0,
            col,
            seed,
            header_fmt
        )

    for row, seed in enumerate(
        SEEDS,
        start=1
    ):

        ws.write(
            row,
            0,
            seed,
            header_fmt
        )

        for col in range(
            len(SEEDS)
        ):

            ws.write_number(
                row,
                col + 1,
                float(
                    matrix[
                        row - 1,
                        col
                    ]
                ),
                number_fmt
            )

    ws.set_column(
        0,
        len(SEEDS),
        13
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
        "Comparison",
        "Neutral ISO, 7-O monoanion, and 7-O/4'-O dianion."
    ],
    [
        "Docking protocol",
        "Same AK1 receptor, search box, exhaustiveness 16, "
        "10 modes and identical 10-seed series."
    ],
    [
        "RMSD atoms",
        f"{N_HEAVY} ligand heavy atoms."
    ],
    [
        "Atom mapping",
        "Charged-state PDBQT heavy atoms were mapped to the "
        "neutral ligand using their identical starting coordinates."
    ],
    [
        "RMSD convention",
        "Direct positional heavy-atom RMSD in the common AK1 "
        "receptor coordinate frame; no post-docking superposition."
    ],
    [
        "Neutral reference",
        "Seed 9009, previously identified as the medoid/representative "
        "of the dominant neutral ISO pose family."
    ],
    [
        "Neutral dominant family",
        "Seeds 1001, 2002, 3003, 4004, 5005, 7007, 9009, 10101."
    ],
    [
        "Clustering",
        "Complete-linkage hierarchical clustering evaluated at "
        "2.0 A and 3.0 A RMSD thresholds."
    ],
    [
        "Interpretation",
        "The principal endpoint is preservation or alteration of "
        "the neutral AK1 pose family across ligand protonation states. "
        "Vina scores are computational docking scores and are not "
        "experimental binding free energies."
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
# CONSOLE SUMMARY
# ============================================================

print()
print("=" * 80)
print("STATE SUMMARY")
print("=" * 80)


for row in state_summary_rows:

    print()
    print(
        row["State"]
    )

    print(
        f"  score mean +/- SD : "
        f"{row['Mean_score']:.3f} +/- "
        f"{row['SD_score']:.3f} kcal/mol"
    )

    print(
        f"  score range       : "
        f"{row['Best_score']:.3f} to "
        f"{row['Worst_score']:.3f}"
    )

    print(
        f"  clusters at 2 A   : "
        f"{row['Clusters_2A']}"
    )

    print(
        f"  largest 2 A group : "
        f"{row['Largest_2A_cluster_n']}/10"
    )

    print(
        f"  clusters at 3 A   : "
        f"{row['Clusters_3A']}"
    )

    print(
        f"  representative    : "
        f"seed {row['Representative_seed']}"
    )

    print(
        f"  representative mean within-cluster RMSD : "
        f"{row['Representative_mean_RMSD']:.3f} A"
    )


print()
print("=" * 80)
print("CROSS-STATE COMPARISON TO NEUTRAL DOMINANT FAMILY")
print("=" * 80)


for row in cross_summary_rows:

    print()
    print(
        row["State"]
    )

    print(
        f"  mean RMSD to neutral seed9009 : "
        f"{row['Mean_RMSD_to_rep9009_A']:.3f} A"
    )

    print(
        f"  RMSD range to seed9009        : "
        f"{row['Min_RMSD_to_rep9009_A']:.3f} - "
        f"{row['Max_RMSD_to_rep9009_A']:.3f} A"
    )

    print(
        f"  mean nearest RMSD to dominant : "
        f"{row['Mean_min_RMSD_to_neutral_dominant_A']:.3f} A"
    )

    print(
        f"  <= 2 A from dominant family   : "
        f"{row['n_within_2A_of_neutral_dominant']}/10 "
        f"({row['Percent_within_2A']:.1f}%)"
    )

    print(
        f"  <= 3 A from dominant family   : "
        f"{row['n_within_3A_of_neutral_dominant']}/10 "
        f"({row['Percent_within_3A']:.1f}%)"
    )


print()
print("=" * 80)
print("OUTPUTS")
print("=" * 80)

print(
    OUT_DIR
)

print()
print(
    xlsx_file
)

print()
print(
    "Analysis complete."
)