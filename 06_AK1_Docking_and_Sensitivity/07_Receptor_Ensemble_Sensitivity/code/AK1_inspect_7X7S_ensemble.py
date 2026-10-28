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

STRUCTURES = (
    PROJECT /
    "structures"
)

INPUT_7X7S = (
    STRUCTURES /
    "AK1_7X7S_NMR.pdb"
)

INPUT_1Z83 = (
    STRUCTURES /
    "AK1_1Z83.pdb"
)

OUT_DIR = (
    PROJECT /
    "results" /
    "AK1_ensemble_sensitivity" /
    "7X7S_inspection"
)

OUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# OUTPUTS
# ============================================================

SUMMARY_CSV = (
    OUT_DIR /
    "AK1_7X7S_model_summary.csv"
)

PAIRWISE_CSV = (
    OUT_DIR /
    "AK1_7X7S_pairwise_CA_RMSD.csv"
)

CANDIDATE_CSV = (
    OUT_DIR /
    "AK1_7X7S_candidate_conformers.csv"
)

XLSX_FILE = (
    OUT_DIR /
    "AK1_7X7S_ENSEMBLE_INSPECTION.xlsx"
)


# ============================================================
# READ MULTI-MODEL PDB
#
# Stores C-alpha coordinates for chain A.
# ============================================================

def read_nmr_models(path):

    if not path.exists():
        raise FileNotFoundError(path)

    models = {}

    current_model = None

    with open(
        path,
        "r",
        encoding="utf-8",
        errors="ignore"
    ) as handle:

        for line in handle:

            if line.startswith("MODEL"):

                fields = line.split()

                if len(fields) >= 2:

                    current_model = int(
                        fields[1]
                    )

                    models[
                        current_model
                    ] = {}

                continue


            if line.startswith("ENDMDL"):

                current_model = None
                continue


            if current_model is None:
                continue


            if not line.startswith("ATOM"):
                continue


            chain = (
                line[21]
                .strip()
            )

            if chain != "A":
                continue


            atom_name = (
                line[12:16]
                .strip()
            )

            if atom_name != "CA":
                continue


            altloc = (
                line[16]
                .strip()
            )

            if altloc not in (
                "",
                "A"
            ):
                continue


            resname = (
                line[17:20]
                .strip()
            )

            resseq = int(
                line[22:26]
            )

            icode = (
                line[26]
                .strip()
            )


            x = float(
                line[30:38]
            )

            y = float(
                line[38:46]
            )

            z = float(
                line[46:54]
            )


            residue_key = (
                resseq,
                icode
            )


            models[
                current_model
            ][
                residue_key
            ] = {

                "resname":
                    resname,

                "coord":
                    np.array(
                        [
                            x,
                            y,
                            z
                        ],
                        dtype=float
                    )
            }


    if not models:

        raise RuntimeError(
            "No MODEL records were read from 7X7S."
        )


    return models


# ============================================================
# READ SINGLE PDB CHAIN A C-ALPHA
# ============================================================

def read_single_structure_ca(
    path
):

    if not path.exists():
        raise FileNotFoundError(path)

    residues = {}

    with open(
        path,
        "r",
        encoding="utf-8",
        errors="ignore"
    ) as handle:

        for line in handle:

            if not line.startswith("ATOM"):
                continue


            chain = (
                line[21]
                .strip()
            )

            if chain != "A":
                continue


            atom_name = (
                line[12:16]
                .strip()
            )

            if atom_name != "CA":
                continue


            altloc = (
                line[16]
                .strip()
            )

            if altloc not in (
                "",
                "A"
            ):
                continue


            resname = (
                line[17:20]
                .strip()
            )

            resseq = int(
                line[22:26]
            )

            icode = (
                line[26]
                .strip()
            )


            x = float(
                line[30:38]
            )

            y = float(
                line[38:46]
            )

            z = float(
                line[46:54]
            )


            key = (
                resseq,
                icode
            )


            # Only first occurrence for chain A
            if key not in residues:

                residues[
                    key
                ] = {

                    "resname":
                        resname,

                    "coord":
                        np.array(
                            [
                                x,
                                y,
                                z
                            ],
                            dtype=float
                        )
                }


    if not residues:

        raise RuntimeError(
            f"No chain-A C-alpha atoms found in {path}"
        )


    return residues


# ============================================================
# RESIDUE HELPER
# ============================================================

def find_residue_by_number(
    model,
    number
):

    matches = [

        value

        for (
            resseq,
            icode
        ), value in model.items()

        if resseq == number
    ]


    if len(matches) != 1:

        raise RuntimeError(
            f"Expected one residue numbered {number}, "
            f"found {len(matches)}."
        )


    return matches[0]


# ============================================================
# DISTANCE
# ============================================================

def ca_distance(
    model,
    residue_a,
    residue_b
):

    a = find_residue_by_number(
        model,
        residue_a
    )["coord"]

    b = find_residue_by_number(
        model,
        residue_b
    )["coord"]


    return float(
        np.linalg.norm(
            a - b
        )
    )


# ============================================================
# KABSCH RMSD
#
# Pairwise C-alpha RMSD after optimal rigid-body
# superposition of the residues common to all 20 models.
# ============================================================

def kabsch_rmsd(
    coords_a,
    coords_b
):

    if coords_a.shape != coords_b.shape:

        raise ValueError(
            "Coordinate arrays have different shapes."
        )


    a = (
        coords_a -
        np.mean(
            coords_a,
            axis=0
        )
    )

    b = (
        coords_b -
        np.mean(
            coords_b,
            axis=0
        )
    )


    covariance = (
        a.T @ b
    )


    u, s, vt = np.linalg.svd(
        covariance
    )


    rotation = (
        u @ vt
    )


    # Avoid improper reflection
    if np.linalg.det(
        rotation
    ) < 0:

        vt[-1, :] *= -1

        rotation = (
            u @ vt
        )


    aligned_a = (
        a @ rotation
    )


    difference = (
        aligned_a -
        b
    )


    return math.sqrt(

        np.mean(

            np.sum(
                difference ** 2,
                axis=1
            )
        )
    )


# ============================================================
# READ STRUCTURES
# ============================================================

models = read_nmr_models(
    INPUT_7X7S
)

closed_1z83 = read_single_structure_ca(
    INPUT_1Z83
)


model_numbers = sorted(
    models.keys()
)


print()
print("=" * 80)
print("7X7S ENSEMBLE INSPECTION")
print("=" * 80)

print()
print(
    "Number of 7X7S models:",
    len(model_numbers)
)


# ============================================================
# RESIDUE COVERAGE
# ============================================================

residue_sets = {

    model_number:
        set(
            models[
                model_number
            ].keys()
        )

    for model_number
    in model_numbers
}


all_residues = set.union(
    *residue_sets.values()
)

common_residues = set.intersection(
    *residue_sets.values()
)


def residue_sort_key(key):

    return (
        key[0],
        key[1]
    )


all_residues_sorted = sorted(
    all_residues,
    key=residue_sort_key
)

common_residues_sorted = sorted(
    common_residues,
    key=residue_sort_key
)


print(
    "Union of modeled C-alpha residues:",
    len(all_residues)
)

print(
    "C-alpha residues common to all 20 models:",
    len(common_residues)
)


# ============================================================
# VERIFY KEY OPENING RESIDUES
#
# These are the same three C-alpha sites used in the
# 7X7S publication to illustrate NMP/LID opening:
#
# V47 - R132
# V47 - T145
# ============================================================

for required_number in (
    47,
    132,
    145
):

    for model_number in model_numbers:

        find_residue_by_number(
            models[
                model_number
            ],
            required_number
        )


# ============================================================
# PRECOMPUTE COMMON CA ARRAYS
# ============================================================

coords_by_model = {}


for model_number in model_numbers:

    coords = []

    for residue_key in (
        common_residues_sorted
    ):

        coords.append(

            models[
                model_number
            ][
                residue_key
            ][
                "coord"
            ]
        )


    coords_by_model[
        model_number
    ] = np.array(
        coords,
        dtype=float
    )


# ============================================================
# PAIRWISE RMSD MATRIX
# ============================================================

n_models = len(
    model_numbers
)


pairwise = np.zeros(
    (
        n_models,
        n_models
    ),
    dtype=float
)


for i, model_i in enumerate(
    model_numbers
):

    for j, model_j in enumerate(
        model_numbers
    ):

        if j < i:

            pairwise[
                i,
                j
            ] = pairwise[
                j,
                i
            ]

            continue


        value = kabsch_rmsd(

            coords_by_model[
                model_i
            ],

            coords_by_model[
                model_j
            ]
        )


        pairwise[
            i,
            j
        ] = value

        pairwise[
            j,
            i
        ] = value


# ============================================================
# MODEL SUMMARY
# ============================================================

summary_rows = []


for i, model_number in enumerate(
    model_numbers
):

    model = models[
        model_number
    ]


    residue_numbers = sorted(
        [
            key[0]
            for key
            in model.keys()
        ]
    )


    missing = sorted(
        all_residues -
        residue_sets[
            model_number
        ],
        key=residue_sort_key
    )


    missing_string = ", ".join(

        [
            (
                f"{number}{icode}"
                if icode
                else str(number)
            )

            for number, icode
            in missing
        ]
    )


    d_47_132 = ca_distance(
        model,
        47,
        132
    )

    d_47_145 = ca_distance(
        model,
        47,
        145
    )


    opening_mean = (
        d_47_132 +
        d_47_145
    ) / 2.0


    other_distances = [

        pairwise[
            i,
            j
        ]

        for j in range(
            n_models
        )

        if j != i
    ]


    summary_rows.append({

        "Model":
            model_number,

        "CA_count":
            len(
                model
            ),

        "First_residue":
            min(
                residue_numbers
            ),

        "Last_residue":
            max(
                residue_numbers
            ),

        "Missing_vs_union":
            missing_string,

        "V47_R132_CA_distance_A":
            d_47_132,

        "V47_T145_CA_distance_A":
            d_47_145,

        "Mean_opening_distance_A":
            opening_mean,

        "Mean_pairwise_CA_RMSD_A":
            float(
                np.mean(
                    other_distances
                )
            ),

        "Max_pairwise_CA_RMSD_A":
            float(
                np.max(
                    other_distances
                )
            )
    })


# ============================================================
# 1Z83 CLOSED-STATE DISTANCES
# ============================================================

closed_d_47_132 = ca_distance(
    closed_1z83,
    47,
    132
)

closed_d_47_145 = ca_distance(
    closed_1z83,
    47,
    145
)

closed_opening_mean = (
    closed_d_47_132 +
    closed_d_47_145
) / 2.0


# ============================================================
# IDENTIFY CANDIDATE REPRESENTATIVE MODELS
# ============================================================

medoid_row = min(
    summary_rows,
    key=lambda row:
        row[
            "Mean_pairwise_CA_RMSD_A"
        ]
)


least_open_row = min(
    summary_rows,
    key=lambda row:
        row[
            "Mean_opening_distance_A"
        ]
)


most_open_row = max(
    summary_rows,
    key=lambda row:
        row[
            "Mean_opening_distance_A"
        ]
)


medoid_index = model_numbers.index(
    medoid_row[
        "Model"
    ]
)


farthest_index = int(
    np.argmax(
        pairwise[
            medoid_index,
            :
        ]
    )
)


farthest_model = (
    model_numbers[
        farthest_index
    ]
)


farthest_row = next(

    row

    for row in summary_rows

    if row[
        "Model"
    ] == farthest_model
)


candidate_rows = [

    {
        "Candidate_type":
            "Ensemble medoid",

        "Model":
            medoid_row[
                "Model"
            ],

        "Mean_opening_distance_A":
            medoid_row[
                "Mean_opening_distance_A"
            ],

        "Mean_pairwise_CA_RMSD_A":
            medoid_row[
                "Mean_pairwise_CA_RMSD_A"
            ],

        "Reason":
            (
                "Lowest mean C-alpha RMSD to the other "
                "7X7S conformers."
            )
    },

    {
        "Candidate_type":
            "Least-open ensemble member",

        "Model":
            least_open_row[
                "Model"
            ],

        "Mean_opening_distance_A":
            least_open_row[
                "Mean_opening_distance_A"
            ],

        "Mean_pairwise_CA_RMSD_A":
            least_open_row[
                "Mean_pairwise_CA_RMSD_A"
            ],

        "Reason":
            (
                "Smallest mean of V47-R132 and "
                "V47-T145 C-alpha distances."
            )
    },

    {
        "Candidate_type":
            "Most-open ensemble member",

        "Model":
            most_open_row[
                "Model"
            ],

        "Mean_opening_distance_A":
            most_open_row[
                "Mean_opening_distance_A"
            ],

        "Mean_pairwise_CA_RMSD_A":
            most_open_row[
                "Mean_pairwise_CA_RMSD_A"
            ],

        "Reason":
            (
                "Largest mean of V47-R132 and "
                "V47-T145 C-alpha distances."
            )
    },

    {
        "Candidate_type":
            "Farthest from ensemble medoid",

        "Model":
            farthest_row[
                "Model"
            ],

        "Mean_opening_distance_A":
            farthest_row[
                "Mean_opening_distance_A"
            ],

        "Mean_pairwise_CA_RMSD_A":
            farthest_row[
                "Mean_pairwise_CA_RMSD_A"
            ],

        "Reason":
            (
                "Largest whole-chain C-alpha RMSD "
                "from the ensemble medoid."
            )
    }
]


# ============================================================
# CONSOLE SUMMARY
# ============================================================

print()
print("=" * 80)
print("OPENING DISTANCES ACROSS 7X7S")
print("=" * 80)

for row in summary_rows:

    print(

        f"Model {row['Model']:>2} | "
        f"V47-R132 = "
        f"{row['V47_R132_CA_distance_A']:.2f} A | "
        f"V47-T145 = "
        f"{row['V47_T145_CA_distance_A']:.2f} A | "
        f"mean opening = "
        f"{row['Mean_opening_distance_A']:.2f} A | "
        f"mean ensemble RMSD = "
        f"{row['Mean_pairwise_CA_RMSD_A']:.3f} A"
    )


print()
print("=" * 80)
print("1Z83 CLOSED REFERENCE")
print("=" * 80)

print(
    f"V47-R132 = "
    f"{closed_d_47_132:.2f} A"
)

print(
    f"V47-T145 = "
    f"{closed_d_47_145:.2f} A"
)

print(
    f"Mean opening distance = "
    f"{closed_opening_mean:.2f} A"
)


print()
print("=" * 80)
print("CANDIDATE 7X7S CONFORMERS")
print("=" * 80)

for row in candidate_rows:

    print()

    print(
        row[
            "Candidate_type"
        ]
    )

    print(
        f"  Model: "
        f"{row['Model']}"
    )

    print(
        f"  Mean opening distance: "
        f"{row['Mean_opening_distance_A']:.2f} A"
    )

    print(
        f"  Mean ensemble RMSD: "
        f"{row['Mean_pairwise_CA_RMSD_A']:.3f} A"
    )

    print(
        f"  Reason: "
        f"{row['Reason']}"
    )


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


write_dict_csv(
    SUMMARY_CSV,
    summary_rows
)

write_dict_csv(
    CANDIDATE_CSV,
    candidate_rows
)


# ============================================================
# PAIRWISE RMSD CSV
# ============================================================

with open(
    PAIRWISE_CSV,
    "w",
    newline="",
    encoding="utf-8"
) as handle:

    writer = csv.writer(
        handle
    )

    writer.writerow(
        ["Model"] +
        model_numbers
    )

    for i, model_number in enumerate(
        model_numbers
    ):

        writer.writerow(

            [
                model_number
            ]
            +
            [
                float(
                    value
                )
                for value in
                pairwise[
                    i,
                    :
                ]
            ]
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


# ============================================================
# MODEL SUMMARY SHEET
# ============================================================

ws = workbook.add_worksheet(
    "Model_Summary"
)

summary_headers = list(
    summary_rows[0].keys()
)

for col, header in enumerate(
    summary_headers
):

    ws.write(
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
    len(summary_headers) - 1,
    25
)


# ============================================================
# PAIRWISE RMSD SHEET
# ============================================================

ws2 = workbook.add_worksheet(
    "Pairwise_CA_RMSD"
)


ws2.write(
    0,
    0,
    "Model",
    header_fmt
)


for col, model_number in enumerate(
    model_numbers,
    start=1
):

    ws2.write(
        0,
        col,
        model_number,
        header_fmt
    )


for row, model_number in enumerate(
    model_numbers,
    start=1
):

    ws2.write(
        row,
        0,
        model_number,
        header_fmt
    )

    for col in range(
        n_models
    ):

        ws2.write_number(
            row,
            col + 1,
            float(
                pairwise[
                    row - 1,
                    col
                ]
            ),
            number_fmt
        )


ws2.set_column(
    0,
    n_models,
    12
)


# ============================================================
# CANDIDATE MODELS
# ============================================================

ws3 = workbook.add_worksheet(
    "Candidate_Models"
)


candidate_headers = list(
    candidate_rows[0].keys()
)


for col, header in enumerate(
    candidate_headers
):

    ws3.write(
        0,
        col,
        header,
        header_fmt
    )


for row_i, row in enumerate(
    candidate_rows,
    start=1
):

    for col, header in enumerate(
        candidate_headers
    ):

        value = row[
            header
        ]

        if isinstance(
            value,
            float
        ):

            ws3.write_number(
                row_i,
                col,
                value,
                number_fmt
            )

        else:

            ws3.write(
                row_i,
                col,
                value,
                cell_fmt
            )


ws3.set_column(
    0,
    len(candidate_headers) - 1,
    35
)


# ============================================================
# RESIDUE COVERAGE
# ============================================================

ws4 = workbook.add_worksheet(
    "Residue_Coverage"
)


coverage_headers = [
    "Model",
    "CA_count",
    "Missing_vs_union"
]


for col, header in enumerate(
    coverage_headers
):

    ws4.write(
        0,
        col,
        header,
        header_fmt
    )


for row_i, row in enumerate(
    summary_rows,
    start=1
):

    ws4.write(
        row_i,
        0,
        row[
            "Model"
        ],
        cell_fmt
    )

    ws4.write(
        row_i,
        1,
        row[
            "CA_count"
        ],
        cell_fmt
    )

    ws4.write(
        row_i,
        2,
        row[
            "Missing_vs_union"
        ],
        cell_fmt
    )


ws4.set_column(
    "A:B",
    16
)

ws4.set_column(
    "C:C",
    90
)


# ============================================================
# CLOSED REFERENCE
# ============================================================

ws5 = workbook.add_worksheet(
    "1Z83_Reference"
)


closed_rows = [

    [
        "V47-R132 CA distance (A)",
        closed_d_47_132
    ],

    [
        "V47-T145 CA distance (A)",
        closed_d_47_145
    ],

    [
        "Mean opening distance (A)",
        closed_opening_mean
    ]
]


for row_i, (
    label,
    value
) in enumerate(
    closed_rows
):

    ws5.write(
        row_i,
        0,
        label,
        header_fmt
    )

    ws5.write_number(
        row_i,
        1,
        float(value),
        number_fmt
    )


ws5.set_column(
    "A:A",
    38
)

ws5.set_column(
    "B:B",
    18
)


# ============================================================
# METHOD NOTES
# ============================================================

ws6 = workbook.add_worksheet(
    "Method_Notes"
)


notes = [

    [
        "7X7S",
        (
            "Wild-type human AK1 apo solution-NMR ensemble. "
            "The PDB contains 20 submitted conformers."
        )
    ],

    [
        "Opening metrics",
        (
            "V47-R132 and V47-T145 C-alpha distances were "
            "calculated because these sites were used in the "
            "7X7S publication to illustrate relative NMP/LID "
            "opening."
        )
    ],

    [
        "Pairwise RMSD",
        (
            "C-alpha RMSD after optimal Kabsch rigid-body "
            "superposition using residues present in every "
            "7X7S model."
        )
    ],

    [
        "Ensemble medoid",
        (
            "The model with the smallest mean pairwise C-alpha "
            "RMSD to all other conformers."
        )
    ],

    [
        "Least/most open",
        (
            "Defined by the minimum/maximum mean of the "
            "V47-R132 and V47-T145 C-alpha distances."
        )
    ],

    [
        "1Z83",
        (
            "Existing AP5-bound closed human AK1 structure used "
            "as the current docking receptor and structural "
            "reference."
        )
    ],

    [
        "Next step",
        (
            "Representative 7X7S conformers will be selected only "
            "after inspecting ensemble coverage, opening metrics "
            "and pairwise RMSD."
        )
    ]
]


for row_i, (
    label,
    description
) in enumerate(
    notes
):

    ws6.write(
        row_i,
        0,
        label,
        header_fmt
    )

    ws6.write(
        row_i,
        1,
        description,
        cell_fmt
    )


ws6.set_column(
    "A:A",
    28
)

ws6.set_column(
    "B:B",
    110
)


workbook.close()


# ============================================================
# FINAL
# ============================================================

print()
print("=" * 80)
print("FILES SAVED")
print("=" * 80)

print()
print(
    SUMMARY_CSV
)

print(
    PAIRWISE_CSV
)

print(
    CANDIDATE_CSV
)

print(
    XLSX_FILE
)

print()
print(
    "7X7S ensemble inspection complete."
)