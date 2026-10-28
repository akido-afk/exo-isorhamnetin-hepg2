from pathlib import Path
from collections import Counter
import csv
import math

import numpy as np
import xlsxwriter


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
    / "AK1_blind_docking"
    / "full_10seed_analysis_v1"
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


# ============================================================
# TARGETED CLOSED-STATE REFERENCE FAMILIES
# ============================================================

DOMINANT_SEEDS = [
    1001,
    2002,
    3003,
    4004,
    5005,
    7007,
    9009,
    10101
]

ALTERNATIVE_SEEDS = [
    6006,
    8008
]


# ============================================================
# REFERENCE FILES
# ============================================================

AP5_FILE = (
    LIGAND_DIR
    / "AK1_native_AP5_chainA.pdb"
)


# ============================================================
# OUTPUT FILES
# ============================================================

ALL_POSES_CSV = (
    OUT_DIR
    / "AK1_BLIND_all_poses.csv"
)

SEED_SUMMARY_CSV = (
    OUT_DIR
    / "AK1_BLIND_seed_summary.csv"
)

OVERALL_SUMMARY_CSV = (
    OUT_DIR
    / "AK1_BLIND_overall_summary.csv"
)

XLSX_FILE = (
    OUT_DIR
    / "AK1_BLIND_10SEED_COMPREHENSIVE_ANALYSIS.xlsx"
)


# ============================================================
# ELEMENT HANDLING
# ============================================================

def pdbqt_element(atom_type):

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


# ============================================================
# READ MULTI-MODE PDBQT
# ============================================================

def read_pdbqt_models(path):

    path = Path(path)

    if not path.exists():

        raise FileNotFoundError(
            path
        )


    lines = path.read_text(
        errors="ignore"
    ).splitlines()


    models = {}
    scores = {}

    current_model = None
    atoms = []


    for line in lines:

        if line.startswith("MODEL"):

            current_model = int(
                line.split()[1]
            )

            atoms = []

            continue


        if line.startswith(
            "REMARK VINA RESULT:"
        ):

            if current_model is not None:

                scores[
                    current_model
                ] = float(
                    line.split()[3]
                )

            continue


        if line.startswith("ENDMDL"):

            if current_model is not None:

                models[
                    current_model
                ] = atoms

            current_model = None
            atoms = []

            continue


        if current_model is None:
            continue


        if not (
            line.startswith("ATOM")
            or
            line.startswith("HETATM")
        ):
            continue


        atom_name = (
            line[12:16]
            .strip()
        )

        atom_type = (
            line.split()[-1]
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


    if not models:

        raise RuntimeError(
            f"No docking models found in {path}"
        )


    return (
        models,
        scores
    )


# ============================================================
# READ TARGETED MODE 1
# ============================================================

def read_targeted_pose(path):

    models, scores = (
        read_pdbqt_models(
            path
        )
    )


    if 1 not in models:

        raise RuntimeError(
            f"No mode 1 found in {path}"
        )


    atoms = (
        models[
            1
        ]
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


    atom_types = [
        atom[
            "atom_type"
        ]
        for atom in atoms
    ]


    return (
        coords,
        atom_types
    )


# ============================================================
# READ AP5 HEAVY ATOMS
# ============================================================

def read_ap5(path):

    coords = []


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


    if not coords:

        raise RuntimeError(
            "No AP5 heavy atoms found."
        )


    return np.array(
        coords,
        dtype=float
    )


# ============================================================
# POSITIONAL RMSD
#
# Same receptor coordinate frame.
# No ligand post-superposition.
# ============================================================

def positional_rmsd(
    coords_a,
    coords_b
):

    if coords_a.shape != coords_b.shape:

        raise RuntimeError(
            "Coordinate arrays have different shapes."
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
# AP5 REGION METRICS
# ============================================================

def ap5_metrics(
    iso_coords,
    ap5_coords
):

    iso_centroid = np.mean(
        iso_coords,
        axis=0
    )

    ap5_centroid = np.mean(
        ap5_coords,
        axis=0
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


    min_to_ap5 = np.min(
        distance_matrix,
        axis=1
    )


    percent_within_4 = float(
        100.0
        *
        np.mean(
            min_to_ap5 <= 4.0
        )
    )


    percent_within_6 = float(
        100.0
        *
        np.mean(
            min_to_ap5 <= 6.0
        )
    )


    return {

        "ISO_AP5_centroid_distance_A":
            centroid_distance,

        "Minimum_ISO_AP5_distance_A":
            minimum_distance,

        "ISO_atoms_within_4A_percent":
            percent_within_4,

        "ISO_atoms_within_6A_percent":
            percent_within_6
    }


# ============================================================
# LOAD TARGETED REFERENCE FAMILIES
# ============================================================

dominant_poses = []
alternative_poses = []

reference_atom_types = None


for seed in DOMINANT_SEEDS:

    path = (
        DOCKING
        /
        f"AK1_ISO_seed{seed}.pdbqt"
    )


    coords, atom_types = (
        read_targeted_pose(
            path
        )
    )


    if reference_atom_types is None:

        reference_atom_types = (
            atom_types
        )


    if atom_types != reference_atom_types:

        raise RuntimeError(
            f"Atom order mismatch in targeted seed {seed}"
        )


    dominant_poses.append(
        coords
    )


for seed in ALTERNATIVE_SEEDS:

    path = (
        DOCKING
        /
        f"AK1_ISO_seed{seed}.pdbqt"
    )


    coords, atom_types = (
        read_targeted_pose(
            path
        )
    )


    if atom_types != reference_atom_types:

        raise RuntimeError(
            f"Atom order mismatch in targeted seed {seed}"
        )


    alternative_poses.append(
        coords
    )


# ============================================================
# LOAD AP5
# ============================================================

ap5_coords = (
    read_ap5(
        AP5_FILE
    )
)


print()
print("=" * 80)
print("AK1 10-SEED WHOLE-PROTEIN BLIND-DOCKING ANALYSIS")
print("=" * 80)

print(
    f"Targeted dominant reference poses : "
    f"{len(dominant_poses)}"
)

print(
    f"Targeted alternative poses        : "
    f"{len(alternative_poses)}"
)

print(
    f"Native AP5 heavy atoms            : "
    f"{len(ap5_coords)}"
)


# ============================================================
# ANALYZE ALL BLIND POSES
# ============================================================

all_pose_rows = []


for seed in SEEDS:

    blind_file = (
        DOCKING
        /
        f"AK1_ISO_BLIND_seed{seed}.pdbqt"
    )


    models, scores = (
        read_pdbqt_models(
            blind_file
        )
    )


    print()
    print(
        f"Seed {seed}: "
        f"{len(models)} modes detected"
    )


    for mode in sorted(
        models
    ):

        atoms = (
            models[
                mode
            ]
        )


        atom_types = [
            atom[
                "atom_type"
            ]
            for atom in atoms
        ]


        if atom_types != reference_atom_types:

            raise RuntimeError(
                f"Atom-order mismatch: "
                f"seed {seed}, mode {mode}"
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


        # ----------------------------------------------------
        # RMSD to original targeted families
        # ----------------------------------------------------

        dominant_rmsds = [

            positional_rmsd(
                coords,
                ref_pose
            )

            for ref_pose
            in dominant_poses
        ]


        alternative_rmsds = [

            positional_rmsd(
                coords,
                ref_pose
            )

            for ref_pose
            in alternative_poses
        ]


        min_dom = float(
            np.min(
                dominant_rmsds
            )
        )


        min_alt = float(
            np.min(
                alternative_rmsds
            )
        )


        if min_dom < min_alt:

            nearest_family = (
                "Dominant"
            )

            nearest_family_rmsd = (
                min_dom
            )

        else:

            nearest_family = (
                "Alternative"
            )

            nearest_family_rmsd = (
                min_alt
            )


        # ----------------------------------------------------
        # AP5 region
        # ----------------------------------------------------

        site_metrics = (
            ap5_metrics(
                coords,
                ap5_coords
            )
        )


        all_pose_rows.append({

            "Seed":
                seed,

            "Mode":
                mode,

            "Vina_score":
                scores[
                    mode
                ],

            "Min_RMSD_to_targeted_dominant_A":
                min_dom,

            "Min_RMSD_to_targeted_alternative_A":
                min_alt,

            "Nearest_targeted_family":
                nearest_family,

            "Nearest_targeted_family_RMSD_A":
                nearest_family_rmsd,

            "Targeted_like_2A":
                nearest_family_rmsd <= 2.0,

            "Targeted_like_3A":
                nearest_family_rmsd <= 3.0,

            "Dominant_like_2A":
                min_dom <= 2.0,

            "Alternative_like_2A":
                min_alt <= 2.0,

            "Dominant_like_3A":
                min_dom <= 3.0,

            "Alternative_like_3A":
                min_alt <= 3.0,

            "AP5_centroid_within_4A":
                site_metrics[
                    "ISO_AP5_centroid_distance_A"
                ] <= 4.0,

            "AP5_centroid_within_6A":
                site_metrics[
                    "ISO_AP5_centroid_distance_A"
                ] <= 6.0,

            **site_metrics
        })


# ============================================================
# ADD SCORE DELTA FROM MODE 1
# ============================================================

for seed in SEEDS:

    seed_rows = [

        row

        for row in all_pose_rows

        if row[
            "Seed"
        ] == seed
    ]


    global_best_score = min(
        row[
            "Vina_score"
        ]
        for row in seed_rows
    )


    for row in seed_rows:

        row[
            "Delta_score_from_seed_best"
        ] = (
            row[
                "Vina_score"
            ]
            -
            global_best_score
        )


# ============================================================
# PER-SEED SUMMARY
# ============================================================

seed_summary_rows = []


for seed in SEEDS:

    seed_rows = sorted(

        [
            row

            for row in all_pose_rows

            if row[
                "Seed"
            ] == seed
        ],

        key=lambda row:
            row[
                "Mode"
            ]
    )


    mode1 = next(

        row

        for row in seed_rows

        if row[
            "Mode"
        ] == 1
    )


    targeted2 = [

        row

        for row in seed_rows

        if row[
            "Targeted_like_2A"
        ]
    ]


    targeted3 = [

        row

        for row in seed_rows

        if row[
            "Targeted_like_3A"
        ]
    ]


    dominant2 = [

        row

        for row in seed_rows

        if row[
            "Dominant_like_2A"
        ]
    ]


    alternative2 = [

        row

        for row in seed_rows

        if row[
            "Alternative_like_2A"
        ]
    ]


    ap5_centroid6 = [

        row

        for row in seed_rows

        if row[
            "AP5_centroid_within_6A"
        ]
    ]


    def first_rank(rows):

        if not rows:
            return None

        return min(
            row[
                "Mode"
            ]
            for row in rows
        )


    def best_score(rows):

        if not rows:
            return None

        return min(
            row[
                "Vina_score"
            ]
            for row in rows
        )


    best2 = (
        best_score(
            targeted2
        )
    )


    best3 = (
        best_score(
            targeted3
        )
    )


    seed_summary_rows.append({

        "Seed":
            seed,

        "Modes_returned":
            len(
                seed_rows
            ),

        "Global_best_score":
            mode1[
                "Vina_score"
            ],

        "Mode1_nearest_family":
            mode1[
                "Nearest_targeted_family"
            ],

        "Mode1_nearest_family_RMSD_A":
            mode1[
                "Nearest_targeted_family_RMSD_A"
            ],

        "Mode1_targeted_like_2A":
            mode1[
                "Targeted_like_2A"
            ],

        "Mode1_targeted_like_3A":
            mode1[
                "Targeted_like_3A"
            ],

        "Mode1_AP5_centroid_distance_A":
            mode1[
                "ISO_AP5_centroid_distance_A"
            ],

        "Any_targeted_like_2A":
            len(
                targeted2
            ) > 0,

        "Any_targeted_like_3A":
            len(
                targeted3
            ) > 0,

        "Targeted_like_2A_n":
            len(
                targeted2
            ),

        "Targeted_like_3A_n":
            len(
                targeted3
            ),

        "First_targeted_like_2A_rank":
            first_rank(
                targeted2
            ),

        "First_targeted_like_3A_rank":
            first_rank(
                targeted3
            ),

        "Best_targeted_like_2A_score":
            best2,

        "Best_targeted_like_3A_score":
            best3,

        "Score_penalty_to_best_targeted_like_2A":
            (
                best2
                -
                mode1[
                    "Vina_score"
                ]
                if best2 is not None
                else None
            ),

        "Score_penalty_to_best_targeted_like_3A":
            (
                best3
                -
                mode1[
                    "Vina_score"
                ]
                if best3 is not None
                else None
            ),

        "Dominant_like_2A_n":
            len(
                dominant2
            ),

        "Alternative_like_2A_n":
            len(
                alternative2
            ),

        "AP5_centroid_within_6A_n":
            len(
                ap5_centroid6
            )
    })


# ============================================================
# OVERALL SUMMARY
# ============================================================

total_poses = len(
    all_pose_rows
)


mode1_rows = [

    row

    for row in all_pose_rows

    if row[
        "Mode"
    ] == 1
]


n_mode1_2A = sum(
    row[
        "Targeted_like_2A"
    ]
    for row in mode1_rows
)


n_mode1_3A = sum(
    row[
        "Targeted_like_3A"
    ]
    for row in mode1_rows
)


n_seed_any_2A = sum(
    row[
        "Any_targeted_like_2A"
    ]
    for row in seed_summary_rows
)


n_seed_any_3A = sum(
    row[
        "Any_targeted_like_3A"
    ]
    for row in seed_summary_rows
)


all_targeted2 = [

    row

    for row in all_pose_rows

    if row[
        "Targeted_like_2A"
    ]
]


all_targeted3 = [

    row

    for row in all_pose_rows

    if row[
        "Targeted_like_3A"
    ]
]


mode1_family_counts = Counter(

    row[
        "Nearest_targeted_family"
    ]

    for row in mode1_rows
)


first_rank_2A_values = [

    row[
        "First_targeted_like_2A_rank"
    ]

    for row in seed_summary_rows

    if row[
        "First_targeted_like_2A_rank"
    ] is not None
]


first_rank_3A_values = [

    row[
        "First_targeted_like_3A_rank"
    ]

    for row in seed_summary_rows

    if row[
        "First_targeted_like_3A_rank"
    ] is not None
]


global_scores = np.array(

    [
        row[
            "Global_best_score"
        ]

        for row in seed_summary_rows
    ],

    dtype=float
)


mode1_rmsds = np.array(

    [
        row[
            "Mode1_nearest_family_RMSD_A"
        ]

        for row in seed_summary_rows
    ],

    dtype=float
)


overall_summary_rows = [

    {
        "Metric":
            "Blind seeds analyzed",

        "Value":
            len(
                SEEDS
            )
    },

    {
        "Metric":
            "Total blind poses analyzed",

        "Value":
            total_poses
    },

    {
        "Metric":
            "Mean global-best Vina score",

        "Value":
            float(
                np.mean(
                    global_scores
                )
            )
    },

    {
        "Metric":
            "SD global-best Vina score",

        "Value":
            float(
                np.std(
                    global_scores,
                    ddof=1
                )
            )
    },

    {
        "Metric":
            "Mode 1 targeted-like <=2A seeds",

        "Value":
            n_mode1_2A
    },

    {
        "Metric":
            "Mode 1 targeted-like <=2A percent",

        "Value":
            100.0
            *
            n_mode1_2A
            /
            len(
                SEEDS
            )
    },

    {
        "Metric":
            "Mode 1 targeted-like <=3A seeds",

        "Value":
            n_mode1_3A
    },

    {
        "Metric":
            "Mode 1 targeted-like <=3A percent",

        "Value":
            100.0
            *
            n_mode1_3A
            /
            len(
                SEEDS
            )
    },

    {
        "Metric":
            "Seeds recovering any targeted-like pose <=2A",

        "Value":
            n_seed_any_2A
    },

    {
        "Metric":
            "Seeds recovering any targeted-like pose <=2A percent",

        "Value":
            100.0
            *
            n_seed_any_2A
            /
            len(
                SEEDS
            )
    },

    {
        "Metric":
            "Seeds recovering any targeted-like pose <=3A",

        "Value":
            n_seed_any_3A
    },

    {
        "Metric":
            "Seeds recovering any targeted-like pose <=3A percent",

        "Value":
            100.0
            *
            n_seed_any_3A
            /
            len(
                SEEDS
            )
    },

    {
        "Metric":
            "All blind poses targeted-like <=2A",

        "Value":
            len(
                all_targeted2
            )
    },

    {
        "Metric":
            "All blind poses targeted-like <=2A percent",

        "Value":
            100.0
            *
            len(
                all_targeted2
            )
            /
            total_poses
    },

    {
        "Metric":
            "All blind poses targeted-like <=3A",

        "Value":
            len(
                all_targeted3
            )
    },

    {
        "Metric":
            "All blind poses targeted-like <=3A percent",

        "Value":
            100.0
            *
            len(
                all_targeted3
            )
            /
            total_poses
    },

    {
        "Metric":
            "Mean mode-1 nearest targeted-family RMSD",

        "Value":
            float(
                np.mean(
                    mode1_rmsds
                )
            )
    },

    {
        "Metric":
            "Median first <=2A targeted-like rank",

        "Value":
            (
                float(
                    np.median(
                        first_rank_2A_values
                    )
                )
                if first_rank_2A_values
                else None
            )
    },

    {
        "Metric":
            "Median first <=3A targeted-like rank",

        "Value":
            (
                float(
                    np.median(
                        first_rank_3A_values
                    )
                )
                if first_rank_3A_values
                else None
            )
    },

    {
        "Metric":
            "Mode 1 nearest dominant-family seeds",

        "Value":
            mode1_family_counts[
                "Dominant"
            ]
    },

    {
        "Metric":
            "Mode 1 nearest alternative-family seeds",

        "Value":
            mode1_family_counts[
                "Alternative"
            ]
    }
]


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
# SAVE CSVs
# ============================================================

write_csv(
    ALL_POSES_CSV,
    all_pose_rows
)

write_csv(
    SEED_SUMMARY_CSV,
    seed_summary_rows
)

write_csv(
    OVERALL_SUMMARY_CSV,
    overall_summary_rows
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


def write_sheet(
    name,
    rows,
    width=25
):

    ws = workbook.add_worksheet(
        name
    )


    if not rows:

        ws.write(
            0,
            0,
            "No data"
        )

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


            if value is None:

                ws.write_blank(
                    row_i,
                    col,
                    None,
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


write_sheet(
    "All_Blind_Poses",
    all_pose_rows,
    width=27
)

write_sheet(
    "Seed_Summary",
    seed_summary_rows,
    width=28
)

write_sheet(
    "Overall_Summary",
    overall_summary_rows,
    width=42
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
            "Whole-protein blind docking of neutral isorhamnetin "
            "against closed human AK1."
        )
    ],

    [
        "Search space",
        (
            "Whole AK1 receptor bounding box with 6 A padding "
            "on each side; box selection did not use AP5."
        )
    ],

    [
        "Docking",
        (
            "10 independent seeds, AutoDock Vina 1.2.7, "
            "exhaustiveness 128, up to 20 modes per seed, "
            "energy range 10 kcal/mol."
        )
    ],

    [
        "Targeted reference",
        (
            "Original AP5-site targeted neutral-ISO docking: "
            "dominant family = 8 seeds; alternative family = 2 seeds."
        )
    ],

    [
        "Pose similarity",
        (
            "Direct ligand heavy-atom positional RMSD in the common "
            "1Z83 receptor coordinate frame; no post-docking "
            "ligand superposition."
        )
    ],

    [
        "Targeted-like thresholds",
        (
            "<=2 A used as stringent recovery and <=3 A as "
            "a broader pose-family recovery threshold."
        )
    ],

    [
        "AP5 metrics",
        (
            "AP5-relative localization was evaluated only after "
            "blind docking and therefore did not influence the "
            "whole-protein search."
        )
    ],

    [
        "Interpretation",
        (
            "Recovery of targeted pose families by unconstrained "
            "whole-protein docking supports that the AP5/nucleotide "
            "region was not selected solely because the original "
            "targeted docking box was centered there."
        )
    ],

    [
        "Caution",
        (
            "Vina scores are computational docking scores, not "
            "experimental binding affinities or free energies."
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
# CONSOLE: PER-SEED RESULTS
# ============================================================

print()
print("=" * 80)
print("BLIND DOCKING PER-SEED RECOVERY")
print("=" * 80)


for row in seed_summary_rows:

    print()

    print(
        f"Seed {row['Seed']}"
    )

    print(
        f"  global best score       : "
        f"{row['Global_best_score']:.3f}"
    )

    print(
        f"  mode 1 nearest family   : "
        f"{row['Mode1_nearest_family']}"
    )

    print(
        f"  mode 1 family RMSD      : "
        f"{row['Mode1_nearest_family_RMSD_A']:.3f} A"
    )

    print(
        f"  mode 1 <=2 A recovery   : "
        f"{row['Mode1_targeted_like_2A']}"
    )

    print(
        f"  mode 1 <=3 A recovery   : "
        f"{row['Mode1_targeted_like_3A']}"
    )

    print(
        f"  first <=2 A rank        : "
        f"{row['First_targeted_like_2A_rank']}"
    )

    print(
        f"  first <=3 A rank        : "
        f"{row['First_targeted_like_3A_rank']}"
    )

    print(
        f"  <=2 A poses             : "
        f"{row['Targeted_like_2A_n']}"
    )

    print(
        f"  dominant-like <=2 A     : "
        f"{row['Dominant_like_2A_n']}"
    )

    print(
        f"  alternative-like <=2 A  : "
        f"{row['Alternative_like_2A_n']}"
    )


# ============================================================
# CONSOLE: OVERALL RESULTS
# ============================================================

print()
print("=" * 80)
print("BLIND DOCKING OVERALL RECOVERY SUMMARY")
print("=" * 80)


for row in overall_summary_rows:

    value = row[
        "Value"
    ]


    if isinstance(
        value,
        float
    ):

        print(
            f"{row['Metric']}: "
            f"{value:.3f}"
        )

    else:

        print(
            f"{row['Metric']}: "
            f"{value}"
        )


print()
print("=" * 80)
print("FILES SAVED")
print("=" * 80)

print(
    XLSX_FILE
)

print(
    ALL_POSES_CSV
)

print(
    SEED_SUMMARY_CSV
)

print(
    OVERALL_SUMMARY_CSV
)

print()
print(
    "10-seed blind-docking analysis complete."
)