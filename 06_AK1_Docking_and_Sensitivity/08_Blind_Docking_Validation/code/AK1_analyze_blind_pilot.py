from pathlib import Path
import csv
import math

import numpy as np
import xlsxwriter


# ============================================================
# PATHS
# ============================================================

PROJECT = Path(
    r"C:\Users\akidi\Desktop\ISO_AK_Docking"
)

BLIND_FILE = (
    PROJECT
    / "docking"
    / "AK1_ISO_BLIND_seed1001.pdbqt"
)

AP5_FILE = (
    PROJECT
    / "ligand"
    / "AK1_native_AP5_chainA.pdb"
)

DOCKING_DIR = (
    PROJECT
    / "docking"
)

OUT_DIR = (
    PROJECT
    / "results"
    / "AK1_blind_docking"
    / "pilot_seed1001"
)

OUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)


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
# PDBQT ELEMENT
# ============================================================

def pdbqt_element(atom_type):

    atom_type = atom_type.strip().upper()

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

    if atom_type in ("C", "A"):
        return "C"

    return atom_type[:1]


# ============================================================
# READ ALL MODELS FROM PDBQT
# ============================================================

def read_pdbqt_models(path):

    lines = Path(path).read_text(
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


        if line.startswith("REMARK VINA RESULT:"):

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


    return (
        models,
        scores
    )


# ============================================================
# READ BEST MODE FROM TARGETED PDBQT
# ============================================================

def read_best_targeted_pose(path):

    models, scores = read_pdbqt_models(
        path
    )

    if 1 not in models:

        raise RuntimeError(
            f"No mode 1 in {path}"
        )

    coords = np.array(
        [
            atom[
                "coord"
            ]
            for atom in models[1]
        ],
        dtype=float
    )

    atom_types = [
        atom[
            "atom_type"
        ]
        for atom in models[1]
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


            x = float(
                line[30:38]
            )

            y = float(
                line[38:46]
            )

            z = float(
                line[46:54]
            )

            coords.append(
                [x, y, z]
            )


    return np.array(
        coords,
        dtype=float
    )


# ============================================================
# RMSD
# ============================================================

def positional_rmsd(
    coords_a,
    coords_b
):

    if coords_a.shape != coords_b.shape:

        raise RuntimeError(
            "Coordinate arrays differ in shape."
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


    minimum_distance = float(
        np.min(
            distance_matrix
        )
    )


    return (
        centroid_distance,
        minimum_distance,
        percent_within_4,
        percent_within_6
    )


# ============================================================
# LOAD BLIND MODES
# ============================================================

blind_models, blind_scores = (
    read_pdbqt_models(
        BLIND_FILE
    )
)


print()
print("=" * 80)
print("AK1 BLIND-DOCKING PILOT ANALYSIS")
print("=" * 80)

print(
    "Blind modes detected:",
    len(
        blind_models
    )
)


# ============================================================
# LOAD AP5
# ============================================================

ap5_coords = (
    read_ap5(
        AP5_FILE
    )
)

print(
    "Native AP5 heavy atoms:",
    len(
        ap5_coords
    )
)


# ============================================================
# LOAD CLOSED TARGETED REFERENCE FAMILIES
# ============================================================

dominant_poses = []
alternative_poses = []

reference_atom_types = None


for seed in DOMINANT_SEEDS:

    path = (
        DOCKING_DIR
        /
        f"AK1_ISO_seed{seed}.pdbqt"
    )

    coords, atom_types = (
        read_best_targeted_pose(
            path
        )
    )

    if reference_atom_types is None:

        reference_atom_types = (
            atom_types
        )

    dominant_poses.append(
        coords
    )


for seed in ALTERNATIVE_SEEDS:

    path = (
        DOCKING_DIR
        /
        f"AK1_ISO_seed{seed}.pdbqt"
    )

    coords, atom_types = (
        read_best_targeted_pose(
            path
        )
    )

    alternative_poses.append(
        coords
    )


# ============================================================
# ANALYZE ALL 20 BLIND MODES
# ============================================================

rows = []


for mode in sorted(
    blind_models
):

    atoms = (
        blind_models[
            mode
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


    if (
        atom_types !=
        reference_atom_types
    ):

        raise RuntimeError(
            f"Atom order mismatch in blind mode {mode}"
        )


    (
        centroid_distance,
        minimum_distance,
        percent_within_4,
        percent_within_6
    ) = ap5_metrics(
        coords,
        ap5_coords
    )


    dom_rmsds = [
        positional_rmsd(
            coords,
            ref
        )
        for ref in dominant_poses
    ]


    alt_rmsds = [
        positional_rmsd(
            coords,
            ref
        )
        for ref in alternative_poses
    ]


    min_dom = float(
        np.min(
            dom_rmsds
        )
    )


    min_alt = float(
        np.min(
            alt_rmsds
        )
    )


    nearest_family = (
        "Dominant"
        if min_dom < min_alt
        else "Alternative"
    )


    nearest_family_rmsd = min(
        min_dom,
        min_alt
    )


    rows.append({

        "Mode":
            mode,

        "Vina_score":
            blind_scores[
                mode
            ],

        "ISO_AP5_centroid_distance_A":
            centroid_distance,

        "Minimum_ISO_AP5_distance_A":
            minimum_distance,

        "ISO_atoms_within_4A_percent":
            percent_within_4,

        "ISO_atoms_within_6A_percent":
            percent_within_6,

        "Min_RMSD_to_targeted_dominant_A":
            min_dom,

        "Min_RMSD_to_targeted_alternative_A":
            min_alt,

        "Nearest_targeted_family":
            nearest_family,

        "Nearest_targeted_family_RMSD_A":
            nearest_family_rmsd,

        "Centroid_within_6A_of_AP5":
            centroid_distance <= 6.0,

        "At_least_80pct_atoms_within_6A_of_AP5":
            percent_within_6 >= 80.0
    })


# ============================================================
# SAVE CSV
# ============================================================

CSV_FILE = (
    OUT_DIR
    / "AK1_BLIND_seed1001_all20_modes_analysis.csv"
)


headers = list(
    rows[0].keys()
)


with open(
    CSV_FILE,
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
# SAVE EXCEL
# ============================================================

XLSX_FILE = (
    OUT_DIR
    / "AK1_BLIND_seed1001_PILOT_ANALYSIS.xlsx"
)


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


num_fmt = workbook.add_format({
    "border": 1,
    "num_format": "0.000"
})


ws = workbook.add_worksheet(
    "All_20_Modes"
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

        value = row[
            header
        ]


        if isinstance(
            value,
            float
        ):

            ws.write_number(
                row_i,
                col,
                value,
                num_fmt
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
    30
)


workbook.close()


# ============================================================
# PRINT RESULTS
# ============================================================

print()
print("=" * 80)
print("BLIND MODES RANKED BY VINA SCORE")
print("=" * 80)


for row in rows:

    print()

    print(
        f"Mode {row['Mode']:2d}"
        f" | score {row['Vina_score']:.3f}"
        f" | AP5 centroid {row['ISO_AP5_centroid_distance_A']:.3f} A"
        f" | within6 {row['ISO_atoms_within_6A_percent']:.1f}%"
        f" | nearest targeted family {row['Nearest_targeted_family']}"
        f" | RMSD {row['Nearest_targeted_family_RMSD_A']:.3f} A"
    )


# ============================================================
# AP5-LIKE MODES
# ============================================================

ap5_region_modes = [

    row

    for row in rows

    if (
        row[
            "Centroid_within_6A_of_AP5"
        ]
        or
        row[
            "At_least_80pct_atoms_within_6A_of_AP5"
        ]
    )
]


print()
print("=" * 80)
print("AP5-REGION CANDIDATE MODES")
print("=" * 80)


if ap5_region_modes:

    for row in ap5_region_modes:

        print(
            f"Mode {row['Mode']:2d}"
            f" | score {row['Vina_score']:.3f}"
            f" | centroid {row['ISO_AP5_centroid_distance_A']:.3f} A"
            f" | within6 {row['ISO_atoms_within_6A_percent']:.1f}%"
            f" | targeted-family RMSD {row['Nearest_targeted_family_RMSD_A']:.3f} A"
        )

else:

    print(
        "No mode met the broad AP5-region criteria."
    )


print()
print("=" * 80)
print("FILES SAVED")
print("=" * 80)

print(
    CSV_FILE
)

print(
    XLSX_FILE
)

print()
print(
    "Blind pilot analysis complete."
)