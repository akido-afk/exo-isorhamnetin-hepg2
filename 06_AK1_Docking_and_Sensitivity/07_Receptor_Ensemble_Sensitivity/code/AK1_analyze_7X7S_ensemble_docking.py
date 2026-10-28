from pathlib import Path
from collections import Counter
import csv
import math

import numpy as np
import xlsxwriter

from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform


# ============================================================
# PROJECT PATHS
# ============================================================

PROJECT = Path(
    r"C:\Users\akidi\Desktop\ISO_AK_Docking"
)

DOCKING = PROJECT / "docking"
LIGAND_DIR = PROJECT / "ligand"

OUT_DIR = (
    PROJECT
    / "results"
    / "AK1_ensemble_sensitivity"
    / "docking_analysis_v1"
)

OUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# OUTPUT FILES
# ============================================================

XLSX_FILE = (
    OUT_DIR
    / "AK1_7X7S_ENSEMBLE_DOCKING_ANALYSIS.xlsx"
)

SEED_CSV = (
    OUT_DIR
    / "AK1_7X7S_ensemble_seed_results.csv"
)

SUMMARY_CSV = (
    OUT_DIR
    / "AK1_7X7S_ensemble_summary.csv"
)

CROSS_RMSD_CSV = (
    OUT_DIR
    / "AK1_7X7S_cross_family_RMSD.csv"
)

AP5_REGION_CSV = (
    OUT_DIR
    / "AK1_7X7S_AP5_region_metrics.csv"
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


# ============================================================
# ORIGINAL 1Z83 POSE FAMILIES
# ============================================================

CLOSED_DOMINANT_SEEDS = [
    1001,
    2002,
    3003,
    4004,
    5005,
    7007,
    9009,
    10101
]

CLOSED_ALTERNATIVE_SEEDS = [
    6006,
    8008
]


# ============================================================
# RECEPTOR STATES
#
# Opening values from previous ensemble inspection.
# ============================================================

RECEPTORS = {

    "1Z83 closed": {

        "template":
            "AK1_ISO_seed{seed}.pdbqt",

        "opening_A":
            18.77,

        "state":
            "Closed AP5-bound crystal"
    },

    "7X7S least-open": {

        "template":
            "AK1_7X7S_model20_least_open_ISO_seed{seed}.pdbqt",

        "opening_A":
            26.36,

        "state":
            "Apo NMR model 20"
    },

    "7X7S medoid": {

        "template":
            "AK1_7X7S_model8_ensemble_medoid_ISO_seed{seed}.pdbqt",

        "opening_A":
            27.38,

        "state":
            "Apo NMR model 8"
    },

    "7X7S most-open": {

        "template":
            "AK1_7X7S_model1_most_open_ISO_seed{seed}.pdbqt",

        "opening_A":
            29.11,

        "state":
            "Apo NMR model 1"
    }
}


# ============================================================
# REFERENCE FILES
# ============================================================

NEUTRAL_INPUT = (
    LIGAND_DIR
    / "isorhamnetin.pdbqt"
)

AP5_NATIVE = (
    LIGAND_DIR
    / "AK1_native_AP5_chainA.pdb"
)


# ============================================================
# PDBQT HEAVY-ATOM READER
# ============================================================

def pdbqt_element(
    atom_type
):

    atom_type = (
        atom_type
        .strip()
        .upper()
    )

    if atom_type.startswith("H"):
        return "H"

    if atom_type.startswith("O"):
        return "O"

    if atom_type.startswith("N"):
        return "N"

    if atom_type.startswith("S"):
        return "S"

    if atom_type.startswith("P"):
        return "P"

    if atom_type in (
        "C",
        "A"
    ):
        return "C"

    return atom_type[:1]


def read_pdbqt_atoms(
    path,
    model_number=None
):

    path = Path(path)

    if not path.exists():

        raise FileNotFoundError(
            path
        )

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

            fields = line.split()

            try:

                current_model = int(
                    fields[1]
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


        atom_name = (
            line[12:16]
            .strip()
        )

        fields = line.split()

        atom_type = (
            fields[-1]
        )

        element = (
            pdbqt_element(
                atom_type
            )
        )


        if element == "H":
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

            "atom_name":
                atom_name,

            "atom_type":
                atom_type,

            "element":
                element,

            "coord":
                np.array(
                    [x, y, z],
                    dtype=float
                )
        })


    if not atoms:

        raise RuntimeError(
            f"No heavy atoms found in {path}"
        )


    return atoms


# ============================================================
# PDB HEAVY-ATOM READER FOR AP5
# ============================================================

def infer_pdb_element(
    line
):

    element = (
        line[76:78]
        .strip()
        .upper()
        if len(line) >= 78
        else ""
    )

    if element:

        return element


    atom_name = (
        line[12:16]
        .strip()
        .upper()
    )

    letters = "".join(
        c
        for c in atom_name
        if c.isalpha()
    )

    return (
        letters[0]
        if letters
        else ""
    )


def read_pdb_heavy_atoms(
    path
):

    path = Path(path)

    if not path.exists():

        raise FileNotFoundError(
            path
        )

    atoms = []


    with open(
        path,
        "r",
        encoding="utf-8",
        errors="ignore"
    ) as handle:

        for line in handle:

            if not (
                line.startswith("ATOM")
                or
                line.startswith("HETATM")
            ):
                continue


            element = (
                infer_pdb_element(
                    line
                )
            )


            if element == "H":
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


            atoms.append(
                np.array(
                    [x, y, z],
                    dtype=float
                )
            )


    if not atoms:

        raise RuntimeError(
            f"No heavy atoms found in {path}"
        )


    return np.array(
        atoms,
        dtype=float
    )


# ============================================================
# READ VINA SCORE
# ============================================================

def read_best_score(
    path
):

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
# POSITIONAL RMSD
#
# No post-docking structural superposition.
# All receptors were placed in the same 1Z83 coordinate frame.
# ============================================================

def positional_rmsd(
    coords_a,
    coords_b
):

    if coords_a.shape != coords_b.shape:

        raise ValueError(
            "Ligand coordinate arrays have different shapes."
        )


    diff = (
        coords_a -
        coords_b
    )


    return math.sqrt(

        np.mean(

            np.sum(
                diff ** 2,
                axis=1
            )
        )
    )


# ============================================================
# READ PREPARED NEUTRAL LIGAND
# ============================================================

input_atoms = read_pdbqt_atoms(
    NEUTRAL_INPUT
)


INPUT_TYPES = [
    atom[
        "atom_type"
    ]
    for atom in input_atoms
]


N_HEAVY = len(
    input_atoms
)


print()
print("=" * 80)
print("NEUTRAL ISO INPUT")
print("=" * 80)

print(
    "Heavy atoms:",
    N_HEAVY
)


# ============================================================
# READ NATIVE AP5
# ============================================================

ap5_coords = (
    read_pdb_heavy_atoms(
        AP5_NATIVE
    )
)


ap5_centroid = (
    np.mean(
        ap5_coords,
        axis=0
    )
)


print(
    "Native AP5 heavy atoms:",
    len(
        ap5_coords
    )
)


# ============================================================
# READ ALL BEST POSES
# ============================================================

scores = {}
poses = {}


for receptor_name, info in (
    RECEPTORS.items()
):

    scores[
        receptor_name
    ] = {}

    poses[
        receptor_name
    ] = {}


    for seed in SEEDS:

        docking_file = (
            DOCKING
            /
            info[
                "template"
            ].format(
                seed=seed
            )
        )


        if not docking_file.exists():

            raise FileNotFoundError(
                docking_file
            )


        score = read_best_score(
            docking_file
        )


        atoms = read_pdbqt_atoms(
            docking_file,
            model_number=1
        )


        if len(atoms) != N_HEAVY:

            raise RuntimeError(
                f"{docking_file.name}: "
                f"expected {N_HEAVY} heavy atoms, "
                f"found {len(atoms)}."
            )


        output_types = [
            atom[
                "atom_type"
            ]
            for atom in atoms
        ]


        if output_types != INPUT_TYPES:

            raise RuntimeError(
                f"Atom order/type mismatch in:\n"
                f"{docking_file}"
            )


        coords = np.array(
            [
                atom[
                    "coord"
                ]
                for atom in atoms
            ],
            dtype=float
        )


        scores[
            receptor_name
        ][
            seed
        ] = score


        poses[
            receptor_name
        ][
            seed
        ] = coords


# ============================================================
# BUILD ORIGINAL CLOSED REFERENCE FAMILIES
# ============================================================

closed_dominant_poses = [

    poses[
        "1Z83 closed"
    ][
        seed
    ]

    for seed in
    CLOSED_DOMINANT_SEEDS
]


closed_alternative_poses = [

    poses[
        "1Z83 closed"
    ][
        seed
    ]

    for seed in
    CLOSED_ALTERNATIVE_SEEDS
]


# ============================================================
# PAIRWISE RMSD MATRICES
# ============================================================

pairwise_matrices = {}


for receptor_name in RECEPTORS:

    matrix = np.zeros(
        (
            len(SEEDS),
            len(SEEDS)
        ),
        dtype=float
    )


    for i, seed_i in enumerate(
        SEEDS
    ):

        for j, seed_j in enumerate(
            SEEDS
        ):

            matrix[
                i,
                j
            ] = positional_rmsd(

                poses[
                    receptor_name
                ][
                    seed_i
                ],

                poses[
                    receptor_name
                ][
                    seed_j
                ]
            )


    pairwise_matrices[
        receptor_name
    ] = matrix


# ============================================================
# CLUSTERING
# ============================================================

def cluster_matrix(
    matrix,
    threshold
):

    if matrix.shape[0] == 1:

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


    return fcluster(
        hierarchy,
        t=threshold,
        criterion="distance"
    )


cluster_results = {}


for receptor_name in RECEPTORS:

    matrix = (
        pairwise_matrices[
            receptor_name
        ]
    )


    cluster_results[
        receptor_name
    ] = {

        "2A":
            cluster_matrix(
                matrix,
                2.0
            ),

        "3A":
            cluster_matrix(
                matrix,
                3.0
            )
    }


# ============================================================
# REPRESENTATIVE POSE
#
# Medoid of the largest 2 A cluster.
# ============================================================

representatives = {}


for receptor_name in RECEPTORS:

    matrix = (
        pairwise_matrices[
            receptor_name
        ]
    )

    labels = (
        cluster_results[
            receptor_name
        ][
            "2A"
        ]
    )


    counts = Counter(
        labels
    )


    largest_label = (
        counts.most_common(
            1
        )[0][0]
    )


    member_indices = [

        i

        for i, label
        in enumerate(
            labels
        )

        if label ==
        largest_label
    ]


    best_index = None
    best_mean = None


    for i in member_indices:

        values = [

            matrix[
                i,
                j
            ]

            for j in member_indices

            if j != i
        ]


        mean_value = (
            float(
                np.mean(
                    values
                )
            )
            if values
            else 0.0
        )


        if (
            best_mean is None
            or
            mean_value < best_mean
        ):

            best_mean = (
                mean_value
            )

            best_index = i


    representatives[
        receptor_name
    ] = {

        "seed":
            SEEDS[
                best_index
            ],

        "cluster_n":
            len(
                member_indices
            ),

        "mean_within_cluster_RMSD_A":
            best_mean
    }


# ============================================================
# AP5 REGION METRICS
# ============================================================

def ap5_region_metrics(
    iso_coords
):

    iso_centroid = (
        np.mean(
            iso_coords,
            axis=0
        )
    )


    centroid_distance = float(
        np.linalg.norm(
            iso_centroid -
            ap5_centroid
        )
    )


    distance_matrix = np.linalg.norm(

        iso_coords[
            :,
            None,
            :
        ]
        -
        ap5_coords[
            None,
            :,
            :
        ],

        axis=2
    )


    minimum_distance = float(
        np.min(
            distance_matrix
        )
    )


    min_to_ap5_per_iso_atom = (
        np.min(
            distance_matrix,
            axis=1
        )
    )


    n_within_4 = int(
        np.sum(
            min_to_ap5_per_iso_atom <= 4.0
        )
    )


    n_within_6 = int(
        np.sum(
            min_to_ap5_per_iso_atom <= 6.0
        )
    )


    percent_within_4 = (
        100.0 *
        n_within_4 /
        len(
            iso_coords
        )
    )


    percent_within_6 = (
        100.0 *
        n_within_6 /
        len(
            iso_coords
        )
    )


    return {

        "ISO_AP5_centroid_distance_A":
            centroid_distance,

        "Minimum_ISO_AP5_heavy_atom_distance_A":
            minimum_distance,

        "ISO_atoms_within_4A_of_AP5_n":
            n_within_4,

        "ISO_atoms_within_4A_of_AP5_percent":
            percent_within_4,

        "ISO_atoms_within_6A_of_AP5_n":
            n_within_6,

        "ISO_atoms_within_6A_of_AP5_percent":
            percent_within_6
    }


# ============================================================
# PER-SEED RESULTS
# ============================================================

seed_rows = []
cross_rows = []
ap5_rows = []


for receptor_name, info in (
    RECEPTORS.items()
):

    labels_2A = (
        cluster_results[
            receptor_name
        ][
            "2A"
        ]
    )

    labels_3A = (
        cluster_results[
            receptor_name
        ][
            "3A"
        ]
    )


    for i, seed in enumerate(
        SEEDS
    ):

        pose = (
            poses[
                receptor_name
            ][
                seed
            ]
        )


        # ----------------------------------------------------
        # RMSD to closed dominant family
        # ----------------------------------------------------

        dom_values = [

            positional_rmsd(
                pose,
                closed_pose
            )

            for closed_pose
            in closed_dominant_poses
        ]


        alt_values = [

            positional_rmsd(
                pose,
                closed_pose
            )

            for closed_pose
            in closed_alternative_poses
        ]


        min_dom = float(
            np.min(
                dom_values
            )
        )


        min_alt = float(
            np.min(
                alt_values
            )
        )


        if min_dom < min_alt:

            nearest_family = (
                "Closed dominant family"
            )

            nearest_family_rmsd = (
                min_dom
            )

        elif min_alt < min_dom:

            nearest_family = (
                "Closed alternative family"
            )

            nearest_family_rmsd = (
                min_alt
            )

        else:

            nearest_family = (
                "Equal"
            )

            nearest_family_rmsd = (
                min_dom
            )


        # ----------------------------------------------------
        # AP5 spatial region
        # ----------------------------------------------------

        ap5_result = (
            ap5_region_metrics(
                pose
            )
        )


        seed_row = {

            "Receptor":
                receptor_name,

            "Receptor_state":
                info[
                    "state"
                ],

            "Opening_metric_A":
                info[
                    "opening_A"
                ],

            "Seed":
                seed,

            "Best_Vina_score":
                scores[
                    receptor_name
                ][
                    seed
                ],

            "Cluster_2A":
                int(
                    labels_2A[
                        i
                    ]
                ),

            "Cluster_3A":
                int(
                    labels_3A[
                        i
                    ]
                ),

            "Min_RMSD_to_closed_dominant_A":
                min_dom,

            "Min_RMSD_to_closed_alternative_A":
                min_alt,

            "Nearest_closed_pose_family":
                nearest_family,

            "Nearest_closed_family_RMSD_A":
                nearest_family_rmsd,

            "Within_2A_of_any_closed_family":
                nearest_family_rmsd <= 2.0,

            "Within_3A_of_any_closed_family":
                nearest_family_rmsd <= 3.0,

            **ap5_result
        }


        seed_rows.append(
            seed_row
        )


        cross_rows.append({

            "Receptor":
                receptor_name,

            "Seed":
                seed,

            "Min_RMSD_to_closed_dominant_A":
                min_dom,

            "Min_RMSD_to_closed_alternative_A":
                min_alt,

            "Nearest_closed_pose_family":
                nearest_family,

            "Nearest_closed_family_RMSD_A":
                nearest_family_rmsd,

            "Within_2A":
                nearest_family_rmsd <= 2.0,

            "Within_3A":
                nearest_family_rmsd <= 3.0
        })


        ap5_rows.append({

            "Receptor":
                receptor_name,

            "Seed":
                seed,

            **ap5_result
        })


# ============================================================
# RECEPTOR SUMMARY
# ============================================================

summary_rows = []


for receptor_name, info in (
    RECEPTORS.items()
):

    group_rows = [

        row

        for row in seed_rows

        if row[
            "Receptor"
        ] == receptor_name
    ]


    score_values = np.array(
        [
            row[
                "Best_Vina_score"
            ]
            for row in group_rows
        ],
        dtype=float
    )


    nearest_rmsd_values = np.array(
        [
            row[
                "Nearest_closed_family_RMSD_A"
            ]
            for row in group_rows
        ],
        dtype=float
    )


    centroid_values = np.array(
        [
            row[
                "ISO_AP5_centroid_distance_A"
            ]
            for row in group_rows
        ],
        dtype=float
    )


    percent6_values = np.array(
        [
            row[
                "ISO_atoms_within_6A_of_AP5_percent"
            ]
            for row in group_rows
        ],
        dtype=float
    )


    n_dom = sum(

        row[
            "Nearest_closed_pose_family"
        ] ==
        "Closed dominant family"

        for row in group_rows
    )


    n_alt = sum(

        row[
            "Nearest_closed_pose_family"
        ] ==
        "Closed alternative family"

        for row in group_rows
    )


    n_2A = sum(
        row[
            "Within_2A_of_any_closed_family"
        ]
        for row in group_rows
    )


    n_3A = sum(
        row[
            "Within_3A_of_any_closed_family"
        ]
        for row in group_rows
    )


    labels2 = (
        cluster_results[
            receptor_name
        ][
            "2A"
        ]
    )


    labels3 = (
        cluster_results[
            receptor_name
        ][
            "3A"
        ]
    )


    rep = (
        representatives[
            receptor_name
        ]
    )


    summary_rows.append({

        "Receptor":
            receptor_name,

        "Receptor_state":
            info[
                "state"
            ],

        "Opening_metric_A":
            info[
                "opening_A"
            ],

        "n":
            len(
                group_rows
            ),

        "Mean_Vina_score":
            float(
                np.mean(
                    score_values
                )
            ),

        "Median_Vina_score":
            float(
                np.median(
                    score_values
                )
            ),

        "SD_Vina_score":
            float(
                np.std(
                    score_values,
                    ddof=1
                )
            ),

        "Best_Vina_score":
            float(
                np.min(
                    score_values
                )
            ),

        "Worst_Vina_score":
            float(
                np.max(
                    score_values
                )
            ),

        "Vina_score_range":
            float(
                np.max(
                    score_values
                )
                -
                np.min(
                    score_values
                )
            ),

        "Clusters_2A":
            len(
                set(
                    labels2
                )
            ),

        "Largest_2A_cluster_n":
            rep[
                "cluster_n"
            ],

        "Clusters_3A":
            len(
                set(
                    labels3
                )
            ),

        "Representative_seed":
            rep[
                "seed"
            ],

        "Representative_mean_within_cluster_RMSD_A":
            rep[
                "mean_within_cluster_RMSD_A"
            ],

        "Nearest_closed_dominant_n":
            n_dom,

        "Nearest_closed_alternative_n":
            n_alt,

        "Within_2A_of_any_closed_family_n":
            n_2A,

        "Within_2A_of_any_closed_family_percent":
            100.0 *
            n_2A /
            len(
                group_rows
            ),

        "Within_3A_of_any_closed_family_n":
            n_3A,

        "Within_3A_of_any_closed_family_percent":
            100.0 *
            n_3A /
            len(
                group_rows
            ),

        "Mean_nearest_closed_family_RMSD_A":
            float(
                np.mean(
                    nearest_rmsd_values
                )
            ),

        "Mean_ISO_AP5_centroid_distance_A":
            float(
                np.mean(
                    centroid_values
                )
            ),

        "Mean_ISO_atoms_within_6A_of_AP5_percent":
            float(
                np.mean(
                    percent6_values
                )
            )
    })


# ============================================================
# CSV HELPER
# ============================================================

def write_csv(
    path,
    rows
):

    if not rows:
        return


    headers = list(
        rows[
            0
        ].keys()
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
# SAVE CSV
# ============================================================

write_csv(
    SEED_CSV,
    seed_rows
)

write_csv(
    SUMMARY_CSV,
    summary_rows
)

write_csv(
    CROSS_RMSD_CSV,
    cross_rows
)

write_csv(
    AP5_REGION_CSV,
    ap5_rows
)


# ============================================================
# SAVE RMSD MATRICES
# ============================================================

for receptor_name, matrix in (
    pairwise_matrices.items()
):

    safe_name = (
        receptor_name
        .replace(" ", "_")
        .replace("-", "_")
    )


    path = (
        OUT_DIR
        /
        f"AK1_{safe_name}_pairwise_RMSD.csv"
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
            ["Seed"]
            +
            SEEDS
        )


        for i, seed in enumerate(
            SEEDS
        ):

            writer.writerow(

                [seed]
                +
                [
                    float(
                        value
                    )
                    for value
                    in matrix[
                        i,
                        :
                    ]
                ]
            )


# ============================================================
# EXCEL
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


# ============================================================
# GENERIC TABLE WRITER
# ============================================================

def write_sheet(
    name,
    rows,
    width=25
):

    ws = workbook.add_worksheet(
        name
    )


    if not rows:

        return


    headers = list(
        rows[
            0
        ].keys()
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

            value = row.get(
                header,
                ""
            )


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
                    float(
                        value
                    ),
                    number_fmt
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
                    str(
                        bool(
                            value
                        )
                    ),
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

    ws.autofilter(
        0,
        0,
        len(
            rows
        ),
        len(
            headers
        ) - 1
    )

    ws.set_column(
        0,
        len(
            headers
        ) - 1,
        width
    )


# ============================================================
# EXCEL TABLES
# ============================================================

write_sheet(
    "Seed_Results",
    seed_rows,
    width=26
)

write_sheet(
    "Receptor_Summary",
    summary_rows,
    width=28
)

write_sheet(
    "Cross_Family_RMSD",
    cross_rows,
    width=27
)

write_sheet(
    "AP5_Region_Metrics",
    ap5_rows,
    width=28
)


# ============================================================
# RMSD MATRIX SHEETS
# ============================================================

sheet_names = {

    "1Z83 closed":
        "RMSD_1Z83",

    "7X7S least-open":
        "RMSD_LeastOpen",

    "7X7S medoid":
        "RMSD_Medoid",

    "7X7S most-open":
        "RMSD_MostOpen"
}


for receptor_name, matrix in (
    pairwise_matrices.items()
):

    ws = workbook.add_worksheet(
        sheet_names[
            receptor_name
        ]
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
            len(
                SEEDS
            )
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
        len(
            SEEDS
        ),
        13
    )


# ============================================================
# METHOD NOTES
# ============================================================

ws_notes = workbook.add_worksheet(
    "Method_Notes"
)


notes = [

    [
        "Experiment",
        (
            "Neutral isorhamnetin docked against closed 1Z83 "
            "and three selected apo 7X7S conformers."
        )
    ],

    [
        "Variable changed",
        (
            "Receptor conformation only. Ligand preparation, "
            "Vina search box, exhaustiveness, number of modes "
            "and 10-seed series were kept constant."
        )
    ],

    [
        "7X7S models",
        (
            "Model 20 = least-open selected conformer; "
            "model 8 = ensemble medoid; "
            "model 1 = most-open selected conformer."
        )
    ],

    [
        "Opening metric",
        (
            "Mean of V47-R132 and V47-T145 C-alpha distances "
            "from the previous ensemble inspection."
        )
    ],

    [
        "Coordinate frame",
        (
            "7X7S conformers were sequence-mapped and rigid-body "
            "aligned onto the stable 1Z83 CORE before receptor "
            "preparation, allowing direct positional comparison."
        )
    ],

    [
        "Pose RMSD",
        (
            "Direct ligand heavy-atom positional RMSD in the "
            "shared 1Z83 coordinate frame; no post-docking "
            "ligand superposition."
        )
    ],

    [
        "Closed pose families",
        (
            "Original 1Z83 neutral ISO dominant family = "
            "seeds 1001, 2002, 3003, 4004, 5005, 7007, "
            "9009 and 10101; alternative family = seeds "
            "6006 and 8008."
        )
    ],

    [
        "AP5-region metrics",
        (
            "Spatial metrics compare docked ISO positions with "
            "the crystallographic AP5 coordinates. They are used "
            "only to assess localization within the experimentally "
            "defined nucleotide-region neighborhood and do not "
            "imply that ISO chemically mimics AP5."
        )
    ],

    [
        "Vina score caution",
        (
            "Vina scores are docking scores and should not be "
            "interpreted as experimental binding free energies "
            "or direct affinities across receptor conformations."
        )
    ]
]


for row_i, (
    label,
    description
) in enumerate(
    notes
):

    ws_notes.write(
        row_i,
        0,
        label,
        header_fmt
    )

    ws_notes.write(
        row_i,
        1,
        description,
        cell_fmt
    )


ws_notes.set_column(
    "A:A",
    30
)

ws_notes.set_column(
    "B:B",
    115
)


workbook.close()


# ============================================================
# CONSOLE SUMMARY
# ============================================================

print()
print("=" * 80)
print("AK1 RECEPTOR-ENSEMBLE DOCKING SUMMARY")
print("=" * 80)


for row in summary_rows:

    print()
    print(
        row[
            "Receptor"
        ]
    )

    print(
        f"  opening metric       : "
        f"{row['Opening_metric_A']:.2f} A"
    )

    print(
        f"  Vina mean +/- SD     : "
        f"{row['Mean_Vina_score']:.3f} +/- "
        f"{row['SD_Vina_score']:.3f} kcal/mol"
    )

    print(
        f"  score range          : "
        f"{row['Best_Vina_score']:.3f} to "
        f"{row['Worst_Vina_score']:.3f}"
    )

    print(
        f"  clusters at 2 A      : "
        f"{row['Clusters_2A']}"
    )

    print(
        f"  largest 2 A cluster  : "
        f"{row['Largest_2A_cluster_n']}/10"
    )

    print(
        f"  clusters at 3 A      : "
        f"{row['Clusters_3A']}"
    )

    print(
        f"  representative seed  : "
        f"{row['Representative_seed']}"
    )

    print(
        f"  nearest closed dominant : "
        f"{row['Nearest_closed_dominant_n']}/10"
    )

    print(
        f"  nearest closed alternative : "
        f"{row['Nearest_closed_alternative_n']}/10"
    )

    print(
        f"  within 2 A of either closed family : "
        f"{row['Within_2A_of_any_closed_family_n']}/10 "
        f"({row['Within_2A_of_any_closed_family_percent']:.1f}%)"
    )

    print(
        f"  within 3 A of either closed family : "
        f"{row['Within_3A_of_any_closed_family_n']}/10 "
        f"({row['Within_3A_of_any_closed_family_percent']:.1f}%)"
    )

    print(
        f"  mean nearest closed-family RMSD : "
        f"{row['Mean_nearest_closed_family_RMSD_A']:.3f} A"
    )

    print(
        f"  mean ISO-AP5 centroid distance : "
        f"{row['Mean_ISO_AP5_centroid_distance_A']:.3f} A"
    )

    print(
        f"  mean ISO atoms within 6 A of AP5 : "
        f"{row['Mean_ISO_atoms_within_6A_of_AP5_percent']:.1f}%"
    )


print()
print("=" * 80)
print("FILES SAVED")
print("=" * 80)

print(
    XLSX_FILE
)

print(
    SUMMARY_CSV
)

print(
    SEED_CSV
)

print(
    CROSS_RMSD_CSV
)

print(
    AP5_REGION_CSV
)

print()
print(
    "Ensemble docking analysis complete."
)