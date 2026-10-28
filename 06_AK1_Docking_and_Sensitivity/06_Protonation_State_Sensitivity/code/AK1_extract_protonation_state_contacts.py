from pathlib import Path
from collections import defaultdict
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

RECEPTOR_FILE = (
    PROJECT /
    "receptors" /
    "AK1_chainA_FIXED.pdb"
)

DOCKING_DIR = (
    PROJECT /
    "docking"
)

OUT_DIR = (
    PROJECT /
    "results" /
    "AK1_ligand_state_sensitivity" /
    "contact_fingerprint_v1"
)

OUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# CONTACT DEFINITIONS
# ============================================================

ANY_CONTACT_CUTOFF = 4.0

POLAR_CONTACT_CUTOFF = 3.6


# ============================================================
# GROUP DEFINITIONS
# ============================================================

GROUPS = {

    "Neutral dominant": {

        "seeds": [
            1001,
            2002,
            3003,
            4004,
            5005,
            7007,
            9009,
            10101
        ],

        "template":
            "AK1_ISO_seed{seed}.pdbqt"
    },

    "Neutral alternative": {

        "seeds": [
            6006,
            8008
        ],

        "template":
            "AK1_ISO_seed{seed}.pdbqt"
    },

    "7-O- monoanion": {

        "seeds": [
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
        ],

        "template":
            "AK1_ISO_7Ominus_seed{seed}.pdbqt"
    },

    "7-O-/4'-O- dianion": {

        "seeds": [
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
        ],

        "template":
            "AK1_ISO_dianion_seed{seed}.pdbqt"
    }
}


# ============================================================
# BASIC ELEMENT HELPERS
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

    if not letters:
        return ""

    # Standard protein atom naming:
    # CA = carbon alpha, not calcium.
    return letters[0]


def pdbqt_element(
    atom_type,
    atom_name
):

    atom_type = (
        atom_type
        .strip()
        .upper()
    )

    atom_name = (
        atom_name
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

    if atom_type in (
        "C",
        "A"
    ):
        return "C"

    if atom_type.startswith("P"):
        return "P"

    letters = "".join(
        c
        for c in atom_name
        if c.isalpha()
    )

    if letters:
        return letters[0]

    return ""


def is_polar_element(
    element
):

    return element.upper() in {
        "N",
        "O",
        "S"
    }


# ============================================================
# READ AK1 PROTEIN
# ============================================================

def read_protein(
    path
):

    if not path.exists():

        raise FileNotFoundError(
            path
        )

    residues = defaultdict(
        list
    )

    with open(
        path,
        "r",
        encoding="utf-8",
        errors="ignore"
    ) as handle:

        for line in handle:

            if not line.startswith(
                "ATOM"
            ):
                continue

            chain = (
                line[21]
                .strip()
            )

            if chain not in (
                "",
                "A"
            ):
                continue

            atom_name = (
                line[12:16]
                .strip()
            )

            residue_name = (
                line[17:20]
                .strip()
            )

            residue_number = (
                line[22:26]
                .strip()
            )

            insertion_code = (
                line[26]
                .strip()
            )

            element = (
                infer_pdb_element(
                    line
                )
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

            residue_id = (
                f"{residue_name}"
                f"{residue_number}"
            )

            if insertion_code:

                residue_id += (
                    insertion_code
                )

            residues[
                residue_id
            ].append({

                "atom_name":
                    atom_name,

                "element":
                    element,

                "coord":
                    np.array(
                        [
                            x,
                            y,
                            z
                        ],
                        dtype=float
                    )
            })

    if not residues:

        raise RuntimeError(
            "No protein residues were read."
        )

    return residues


# ============================================================
# READ BEST DOCKED POSE FROM PDBQT
# ============================================================

def read_best_ligand_pose(
    path
):

    if not path.exists():

        raise FileNotFoundError(
            path
        )

    atoms = []

    current_model = None
    has_models = False

    lines = path.read_text(
        errors="ignore"
    ).splitlines()

    has_models = any(
        line.startswith(
            "MODEL"
        )
        for line in lines
    )

    for line in lines:

        if line.startswith(
            "MODEL"
        ):

            fields = (
                line.split()
            )

            if len(fields) >= 2:

                try:

                    current_model = int(
                        fields[1]
                    )

                except Exception:

                    current_model = None

            continue

        if line.startswith(
            "ENDMDL"
        ):

            current_model = None
            continue

        if not (
            line.startswith("ATOM")
            or
            line.startswith("HETATM")
        ):
            continue

        if has_models:

            if current_model != 1:
                continue

        atom_name = (
            line[12:16]
            .strip()
        )

        fields = (
            line.split()
        )

        atom_type = (
            fields[-1]
        )

        element = pdbqt_element(
            atom_type,
            atom_name
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

        except Exception as exc:

            raise RuntimeError(
                f"Could not parse coordinates from:\n"
                f"{path}\n{line}"
            ) from exc

        atoms.append({

            "atom_name":
                atom_name,

            "atom_type":
                atom_type,

            "element":
                element,

            "coord":
                np.array(
                    [
                        x,
                        y,
                        z
                    ],
                    dtype=float
                )
        })

    if not atoms:

        raise RuntimeError(
            f"No heavy ligand atoms read from {path}"
        )

    return atoms


# ============================================================
# DISTANCE
# ============================================================

def distance(
    a,
    b
):

    return float(
        np.linalg.norm(
            a - b
        )
    )


# ============================================================
# EXTRACT CONTACTS FOR ONE POSE
# ============================================================

def calculate_pose_contacts(
    protein_residues,
    ligand_atoms
):

    results = []

    ligand_polar = [
        atom
        for atom in ligand_atoms
        if is_polar_element(
            atom["element"]
        )
    ]

    for residue_id, protein_atoms in (
        protein_residues.items()
    ):

        min_any = math.inf
        min_any_pair = None

        min_polar = math.inf
        min_polar_pair = None


        # ----------------------------------------------------
        # ANY HEAVY-ATOM CONTACT
        # ----------------------------------------------------

        for protein_atom in protein_atoms:

            for ligand_atom in ligand_atoms:

                d = distance(
                    protein_atom["coord"],
                    ligand_atom["coord"]
                )

                if d < min_any:

                    min_any = d

                    min_any_pair = (
                        protein_atom,
                        ligand_atom
                    )


        # ----------------------------------------------------
        # POLAR HEAVY-ATOM PROXIMITY
        # ----------------------------------------------------

        protein_polar = [
            atom
            for atom in protein_atoms
            if is_polar_element(
                atom["element"]
            )
        ]

        for protein_atom in protein_polar:

            for ligand_atom in ligand_polar:

                d = distance(
                    protein_atom["coord"],
                    ligand_atom["coord"]
                )

                if d < min_polar:

                    min_polar = d

                    min_polar_pair = (
                        protein_atom,
                        ligand_atom
                    )


        any_contact = (
            min_any <=
            ANY_CONTACT_CUTOFF
        )

        polar_contact = (
            min_polar <=
            POLAR_CONTACT_CUTOFF
        )


        if (
            any_contact
            or
            polar_contact
        ):

            row = {

                "Residue":
                    residue_id,

                "Any_contact":
                    any_contact,

                "Any_min_distance_A":
                    (
                        min_any
                        if math.isfinite(
                            min_any
                        )
                        else None
                    ),

                "Polar_contact":
                    polar_contact,

                "Polar_min_distance_A":
                    (
                        min_polar
                        if math.isfinite(
                            min_polar
                        )
                        else None
                    ),

                "Any_protein_atom":
                    (
                        min_any_pair[0][
                            "atom_name"
                        ]
                        if min_any_pair
                        else ""
                    ),

                "Any_ligand_atom":
                    (
                        min_any_pair[1][
                            "atom_name"
                        ]
                        if min_any_pair
                        else ""
                    ),

                "Polar_protein_atom":
                    (
                        min_polar_pair[0][
                            "atom_name"
                        ]
                        if min_polar_pair
                        else ""
                    ),

                "Polar_ligand_atom":
                    (
                        min_polar_pair[1][
                            "atom_name"
                        ]
                        if min_polar_pair
                        else ""
                    )
            }

            results.append(
                row
            )

    return results


# ============================================================
# LOAD PROTEIN
# ============================================================

protein = read_protein(
    RECEPTOR_FILE
)


print()
print("=" * 80)
print("AK1 PROTONATION-STATE CONTACT ANALYSIS")
print("=" * 80)

print()
print(
    "Protein residues loaded:",
    len(protein)
)

print(
    f"Any-contact cutoff: "
    f"{ANY_CONTACT_CUTOFF:.1f} A"
)

print(
    f"Polar-contact cutoff: "
    f"{POLAR_CONTACT_CUTOFF:.1f} A"
)


# ============================================================
# RUN ALL POSES
# ============================================================

pose_contact_rows = []


for group_name, info in (
    GROUPS.items()
):

    print()
    print(
        "-" * 80
    )

    print(
        group_name
    )

    print(
        "-" * 80
    )

    for seed in info[
        "seeds"
    ]:

        docking_file = (
            DOCKING_DIR /
            info[
                "template"
            ].format(
                seed=seed
            )
        )

        ligand_atoms = (
            read_best_ligand_pose(
                docking_file
            )
        )

        contacts = (
            calculate_pose_contacts(
                protein,
                ligand_atoms
            )
        )

        print(
            f"Seed {seed:<5} : "
            f"{sum(r['Any_contact'] for r in contacts)} "
            f"any contacts, "
            f"{sum(r['Polar_contact'] for r in contacts)} "
            f"polar residues"
        )

        for row in contacts:

            row_out = {

                "Group":
                    group_name,

                "Seed":
                    seed,

                **row
            }

            pose_contact_rows.append(
                row_out
            )


# ============================================================
# GET ALL RESIDUES SEEN IN ANY GROUP
# ============================================================

all_contact_residues = sorted(
    set(
        row["Residue"]
        for row
        in pose_contact_rows
    )
)


# ============================================================
# BUILD RECURRENCE SUMMARY
# ============================================================

recurrence_rows = []


for group_name, info in (
    GROUPS.items()
):

    group_seeds = (
        info["seeds"]
    )

    n_total = len(
        group_seeds
    )

    for residue in (
        all_contact_residues
    ):

        relevant = [
            row
            for row in pose_contact_rows
            if (
                row["Group"] ==
                group_name
                and
                row["Residue"] ==
                residue
            )
        ]


        # ----------------------------------------------------
        # Contact counts
        # ----------------------------------------------------

        any_rows = [
            row
            for row in relevant
            if row[
                "Any_contact"
            ]
        ]

        polar_rows = [
            row
            for row in relevant
            if row[
                "Polar_contact"
            ]
        ]


        any_n = len(
            any_rows
        )

        polar_n = len(
            polar_rows
        )


        # ----------------------------------------------------
        # Distances
        # ----------------------------------------------------

        any_distances = [
            row[
                "Any_min_distance_A"
            ]
            for row in any_rows
            if row[
                "Any_min_distance_A"
            ] is not None
        ]

        polar_distances = [
            row[
                "Polar_min_distance_A"
            ]
            for row in polar_rows
            if row[
                "Polar_min_distance_A"
            ] is not None
        ]


        recurrence_rows.append({

            "Group":
                group_name,

            "Residue":
                residue,

            "n_poses":
                n_total,

            "Any_contact_n":
                any_n,

            "Any_contact_percent":
                (
                    100.0 *
                    any_n /
                    n_total
                ),

            "Any_mean_min_distance_A":
                (
                    float(
                        np.mean(
                            any_distances
                        )
                    )
                    if any_distances
                    else None
                ),

            "Any_shortest_distance_A":
                (
                    float(
                        np.min(
                            any_distances
                        )
                    )
                    if any_distances
                    else None
                ),

            "Polar_contact_n":
                polar_n,

            "Polar_contact_percent":
                (
                    100.0 *
                    polar_n /
                    n_total
                ),

            "Polar_mean_min_distance_A":
                (
                    float(
                        np.mean(
                            polar_distances
                        )
                    )
                    if polar_distances
                    else None
                ),

            "Polar_shortest_distance_A":
                (
                    float(
                        np.min(
                            polar_distances
                        )
                    )
                    if polar_distances
                    else None
                )
        })


# ============================================================
# FILTER RESIDUES WITH AT LEAST ONE POLAR CONTACT
# ============================================================

polar_recurrence_rows = [
    row
    for row in recurrence_rows
    if row[
        "Polar_contact_n"
    ] > 0
]


# ============================================================
# BUILD COMPARISON MATRIX
#
# One row per residue.
# Columns = recurrence counts and percentages per state.
# ============================================================

comparison_rows = []


for residue in (
    all_contact_residues
):

    row_out = {
        "Residue":
            residue
    }

    for group_name, info in (
        GROUPS.items()
    ):

        matching = [
            row
            for row in recurrence_rows
            if (
                row["Group"] ==
                group_name
                and
                row["Residue"] ==
                residue
            )
        ]

        if not matching:
            continue

        row = matching[0]

        safe_group = (
            group_name
            .replace(" ", "_")
            .replace("/", "_")
            .replace("'", "")
            .replace("-", "")
        )

        row_out[
            f"{safe_group}_n"
        ] = row[
            "n_poses"
        ]

        row_out[
            f"{safe_group}_Any_n"
        ] = row[
            "Any_contact_n"
        ]

        row_out[
            f"{safe_group}_Any_percent"
        ] = row[
            "Any_contact_percent"
        ]

        row_out[
            f"{safe_group}_Polar_n"
        ] = row[
            "Polar_contact_n"
        ]

        row_out[
            f"{safe_group}_Polar_percent"
        ] = row[
            "Polar_contact_percent"
        ]

    comparison_rows.append(
        row_out
    )


# ============================================================
# RECURRENT CONTACT SUBSETS
#
# Define "recurrent" as >=80% of poses.
# For neutral dominant: >=7/8
# Neutral alternative: >=2/2
# Charged states: >=8/10
# ============================================================

RECURRENCE_THRESHOLD = 80.0


recurrent_any_rows = [
    row
    for row in recurrence_rows
    if row[
        "Any_contact_percent"
    ] >= RECURRENCE_THRESHOLD
]


recurrent_polar_rows = [
    row
    for row in recurrence_rows
    if row[
        "Polar_contact_percent"
    ] >= RECURRENCE_THRESHOLD
]


# ============================================================
# PRINT GROUP SUMMARIES
# ============================================================

print()
print("=" * 80)
print("RECURRENT CONTACT SUMMARY (>=80%)")
print("=" * 80)


for group_name in GROUPS:

    print()
    print(
        group_name
    )

    print(
        "-" * len(
            group_name
        )
    )

    any_residues = [
        row["Residue"]
        for row in recurrent_any_rows
        if row["Group"] ==
        group_name
    ]

    polar_residues = [
        row["Residue"]
        for row in recurrent_polar_rows
        if row["Group"] ==
        group_name
    ]

    print(
        "Recurrent ANY:"
    )

    print(
        ", ".join(
            any_residues
        )
        if any_residues
        else "None"
    )

    print()

    print(
        "Recurrent POLAR:"
    )

    print(
        ", ".join(
            polar_residues
        )
        if polar_residues
        else "None"
    )


# ============================================================
# VALIDATE NEUTRAL DOMINANT POLAR CONTACTS
# ============================================================

expected_neutral_polar = {
    "LYS21",
    "MET61",
    "ARG97",
    "ARG132",
    "ARG138",
    "ARG149"
}


observed_neutral_polar_100 = {
    row["Residue"]
    for row in recurrence_rows
    if (
        row["Group"] ==
        "Neutral dominant"
        and
        row["Polar_contact_n"] ==
        8
    )
}


print()
print("=" * 80)
print("NEUTRAL-DOMINANT VALIDATION")
print("=" * 80)

print(
    "Previously established 8/8 polar residues:"
)

print(
    ", ".join(
        sorted(
            expected_neutral_polar
        )
    )
)

print()

print(
    "Recovered by this script:"
)

print(
    ", ".join(
        sorted(
            observed_neutral_polar_100
        )
    )
    if observed_neutral_polar_100
    else "None"
)


missing_expected = (
    expected_neutral_polar
    -
    observed_neutral_polar_100
)

unexpected = (
    observed_neutral_polar_100
    -
    expected_neutral_polar
)


if not missing_expected and not unexpected:

    print()
    print(
        "VALIDATION PASSED: exact 8/8 polar set reproduced."
    )

else:

    print()
    print(
        "WARNING: neutral-dominant polar fingerprint "
        "does not exactly reproduce the previous analysis."
    )

    if missing_expected:

        print(
            "Missing:",
            ", ".join(
                sorted(
                    missing_expected
                )
            )
        )

    if unexpected:

        print(
            "Additional:",
            ", ".join(
                sorted(
                    unexpected
                )
            )
        )


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
# SAVE CSV FILES
# ============================================================

POSE_CONTACTS_CSV = (
    OUT_DIR /
    "AK1_protonation_state_all_pose_contacts.csv"
)

RECURRENCE_CSV = (
    OUT_DIR /
    "AK1_protonation_state_contact_recurrence.csv"
)

POLAR_RECURRENCE_CSV = (
    OUT_DIR /
    "AK1_protonation_state_polar_contact_recurrence.csv"
)

COMPARISON_CSV = (
    OUT_DIR /
    "AK1_protonation_state_contact_comparison_matrix.csv"
)

RECURRENT_ANY_CSV = (
    OUT_DIR /
    "AK1_protonation_state_recurrent_any_contacts_80pct.csv"
)

RECURRENT_POLAR_CSV = (
    OUT_DIR /
    "AK1_protonation_state_recurrent_polar_contacts_80pct.csv"
)


write_csv(
    POSE_CONTACTS_CSV,
    pose_contact_rows
)

write_csv(
    RECURRENCE_CSV,
    recurrence_rows
)

write_csv(
    POLAR_RECURRENCE_CSV,
    polar_recurrence_rows
)

write_csv(
    COMPARISON_CSV,
    comparison_rows
)

write_csv(
    RECURRENT_ANY_CSV,
    recurrent_any_rows
)

write_csv(
    RECURRENT_POLAR_CSV,
    recurrent_polar_rows
)


# ============================================================
# EXCEL WORKBOOK
# ============================================================

XLSX_FILE = (
    OUT_DIR /
    "AK1_PROTONATION_STATE_CONTACT_FINGERPRINTS.xlsx"
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


number_fmt = workbook.add_format({
    "border": 1,
    "num_format": "0.000"
})


percent_fmt = workbook.add_format({
    "border": 1,
    "num_format": "0.0"
})


# ============================================================
# GENERIC EXCEL WRITER
# ============================================================

def write_sheet(
    workbook,
    sheet_name,
    rows,
    width=22
):

    worksheet = workbook.add_worksheet(
        sheet_name
    )

    if not rows:

        worksheet.write(
            0,
            0,
            "No rows"
        )

        return

    headers = list(
        rows[0].keys()
    )

    for col, header in enumerate(
        headers
    ):

        worksheet.write(
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

                worksheet.write_blank(
                    row_i,
                    col,
                    None,
                    cell_fmt
                )

            elif isinstance(
                value,
                (
                    bool,
                    np.bool_
                )
            ):

                worksheet.write(
                    row_i,
                    col,
                    str(
                        bool(
                            value
                        )
                    ),
                    cell_fmt
                )

            elif isinstance(
                value,
                (
                    float,
                    np.floating
                )
            ):

                if (
                    "percent"
                    in header.lower()
                ):

                    worksheet.write_number(
                        row_i,
                        col,
                        float(value),
                        percent_fmt
                    )

                else:

                    worksheet.write_number(
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

                worksheet.write_number(
                    row_i,
                    col,
                    int(value),
                    cell_fmt
                )

            else:

                worksheet.write(
                    row_i,
                    col,
                    value,
                    cell_fmt
                )

    worksheet.freeze_panes(
        1,
        0
    )

    worksheet.autofilter(
        0,
        0,
        len(rows),
        len(headers) - 1
    )

    worksheet.set_column(
        0,
        len(headers) - 1,
        width
    )


# ============================================================
# EXCEL SHEETS
# ============================================================

write_sheet(
    workbook,
    "All_Pose_Contacts",
    pose_contact_rows,
    width=22
)

write_sheet(
    workbook,
    "Recurrence_All",
    recurrence_rows,
    width=24
)

write_sheet(
    workbook,
    "Polar_Recurrence",
    polar_recurrence_rows,
    width=24
)

write_sheet(
    workbook,
    "Comparison_Matrix",
    comparison_rows,
    width=22
)

write_sheet(
    workbook,
    "Recurrent_Any_80pct",
    recurrent_any_rows,
    width=24
)

write_sheet(
    workbook,
    "Recurrent_Polar_80pct",
    recurrent_polar_rows,
    width=24
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
            "Compare AK1 residue-contact fingerprints across "
            "isorhamnetin protonation states and pose families."
        )
    ],

    [
        "Neutral dominant",
        (
            "Eight previously established dominant neutral poses: "
            "seeds 1001, 2002, 3003, 4004, 5005, 7007, 9009, 10101."
        )
    ],

    [
        "Neutral alternative",
        (
            "Two previously established alternative neutral poses: "
            "seeds 6006 and 8008."
        )
    ],

    [
        "7-O- monoanion",
        (
            "Ten best poses from the 7-O-deprotonated ligand. "
            "All previously assigned to the neutral alternative "
            "pose family."
        )
    ],

    [
        "Dianion",
        (
            "Ten best poses from the 7-O-/4'-O-deprotonated ligand. "
            "All previously assigned to the neutral dominant "
            "pose family."
        )
    ],

    [
        "Any contact definition",
        (
            f"Minimum heavy-atom distance <= "
            f"{ANY_CONTACT_CUTOFF:.1f} A between ligand and "
            "any heavy atom of an AK1 residue."
        )
    ],

    [
        "Polar proximity definition",
        (
            f"Minimum O/N/S heavy-atom distance <= "
            f"{POLAR_CONTACT_CUTOFF:.1f} A. This is reported as "
            "polar heavy-atom proximity and should not automatically "
            "be interpreted as a confirmed hydrogen bond."
        )
    ],

    [
        "Recurrent contact",
        (
            f"A contact present in >= {RECURRENCE_THRESHOLD:.0f}% "
            "of the poses in that group."
        )
    ],

    [
        "Receptor",
        str(
            RECEPTOR_FILE
        )
    ],

    [
        "Interpretation",
        (
            "Changes in residue-contact recurrence indicate how "
            "ligand protonation and pose-family selection alter "
            "the local interaction environment. Docking contacts "
            "do not establish experimental binding or inhibition."
        )
    ]
]


for row_i, (
    label,
    description
) in enumerate(
    notes
):

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
    31
)

notes_sheet.set_column(
    "B:B",
    115
)


workbook.close()


# ============================================================
# FINAL OUTPUT
# ============================================================

print()
print("=" * 80)
print("FILES SAVED")
print("=" * 80)

print()
print(
    XLSX_FILE
)

print()
print(
    POSE_CONTACTS_CSV
)

print(
    RECURRENCE_CSV
)

print(
    POLAR_RECURRENCE_CSV
)

print(
    COMPARISON_CSV
)

print(
    RECURRENT_ANY_CSV
)

print(
    RECURRENT_POLAR_CSV
)

print()
print(
    "Contact-fingerprint extraction complete."
)