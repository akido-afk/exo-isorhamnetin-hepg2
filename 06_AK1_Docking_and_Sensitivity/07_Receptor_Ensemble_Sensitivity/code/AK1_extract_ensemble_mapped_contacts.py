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

RECEPTORS = PROJECT / "receptors"
DOCKING = PROJECT / "docking"

MAPPING_FILE = (
    PROJECT
    / "results"
    / "AK1_ensemble_sensitivity"
    / "7X7S_selected_models"
    / "AK1_7X7S_to_1Z83_sequence_mapping.csv"
)

OUT_DIR = (
    PROJECT
    / "results"
    / "AK1_ensemble_sensitivity"
    / "mapped_contact_fingerprint_v1"
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

RECURRENCE_THRESHOLD = 80.0


# ============================================================
# GROUP DEFINITIONS
#
# Mapping model:
#   None = already uses 1Z83 numbering
#   20/8/1 = map 7X7S numbering back to 1Z83
# ============================================================

GROUPS = {

    "1Z83 closed dominant": {

        "receptor":
            RECEPTORS /
            "AK1_chainA_FIXED.pdb",

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

        "dock_template":
            "AK1_ISO_seed{seed}.pdbqt",

        "mapping_model":
            None
    },

    "7X7S least-open": {

        "receptor":
            RECEPTORS /
            "AK1_7X7S_model20_least_open_FIXED.pdb",

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

        "dock_template":
            "AK1_7X7S_model20_least_open_ISO_seed{seed}.pdbqt",

        "mapping_model":
            20
    },

    "7X7S medoid": {

        "receptor":
            RECEPTORS /
            "AK1_7X7S_model8_ensemble_medoid_FIXED.pdb",

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

        "dock_template":
            "AK1_7X7S_model8_ensemble_medoid_ISO_seed{seed}.pdbqt",

        "mapping_model":
            8
    },

    "7X7S most-open": {

        "receptor":
            RECEPTORS /
            "AK1_7X7S_model1_most_open_FIXED.pdb",

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

        "dock_template":
            "AK1_7X7S_model1_most_open_ISO_seed{seed}.pdbqt",

        "mapping_model":
            1
    }
}


# ============================================================
# ELEMENT HELPERS
# ============================================================

def infer_pdb_element(line):

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

    # Protein CA = alpha carbon, not calcium.
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

    if atom_type.startswith("P"):
        return "P"

    if atom_type in (
        "C",
        "A"
    ):
        return "C"

    letters = "".join(
        c
        for c in atom_name
        if c.isalpha()
    )

    if letters:
        return letters[0]

    return ""


def is_polar_element(element):

    return element.upper() in {
        "N",
        "O",
        "S"
    }


# ============================================================
# READ 7X7S -> 1Z83 SEQUENCE MAPPING
# ============================================================

if not MAPPING_FILE.exists():

    raise FileNotFoundError(
        MAPPING_FILE
    )


mapping_by_model = {
    20: {},
    8: {},
    1: {}
}


with open(
    MAPPING_FILE,
    "r",
    encoding="utf-8"
) as handle:

    reader = csv.DictReader(
        handle
    )

    for row in reader:

        model = int(
            row["Model"]
        )

        if model not in mapping_by_model:
            continue

        mobile_resnum = int(
            row["Mobile_resnum"]
        )

        mapping_by_model[
            model
        ][
            mobile_resnum
        ] = {

            "Reference_resnum":
                int(
                    row[
                        "Reference_resnum"
                    ]
                ),

            "Reference_resname":
                row[
                    "Reference_resname"
                ],

            "Mobile_resname":
                row[
                    "Mobile_resname"
                ],

            "Identical":
                (
                    row[
                        "Identical"
                    ]
                    .strip()
                    .lower()
                    ==
                    "true"
                )
        }


print()
print("=" * 80)
print("RESIDUE MAPPING")
print("=" * 80)

for model in (
    20,
    8,
    1
):

    n_map = len(
        mapping_by_model[
            model
        ]
    )

    n_identical = sum(
        1
        for item in
        mapping_by_model[
            model
        ].values()
        if item[
            "Identical"
        ]
    )

    print(
        f"7X7S model {model}: "
        f"{n_map} mapped residues; "
        f"{n_identical} sequence-identical"
    )


# ============================================================
# READ PROTEIN RECEPTOR
# ============================================================

def read_protein(path):

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

            atom_name = (
                line[12:16]
                .strip()
            )

            resname = (
                line[17:20]
                .strip()
            )

            resnum = int(
                line[22:26]
            )

            icode = (
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


            residue_key = (
                chain,
                resname,
                resnum,
                icode
            )


            residues[
                residue_key
            ].append({

                "atom_name":
                    atom_name,

                "element":
                    element,

                "coord":
                    np.array(
                        [x, y, z],
                        dtype=float
                    )
            })


    if not residues:

        raise RuntimeError(
            f"No protein residues read from {path}"
        )


    return residues


# ============================================================
# READ BEST DOCKED ISO POSE
# ============================================================

def read_best_ligand_pose(path):

    if not path.exists():

        raise FileNotFoundError(
            path
        )

    atoms = []

    lines = path.read_text(
        errors="ignore"
    ).splitlines()

    has_models = any(
        line.startswith("MODEL")
        for line in lines
    )

    current_model = None


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

            if current_model != 1:
                continue


        atom_name = (
            line[12:16]
            .strip()
        )

        fields = line.split()

        atom_type = (
            fields[-1]
        )

        element = pdbqt_element(
            atom_type,
            atom_name
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
            f"No heavy ligand atoms read from {path}"
        )


    return atoms


# ============================================================
# MAP PROTEIN RESIDUE TO 1Z83 NUMBERING
# ============================================================

def map_residue_to_reference(
    group_name,
    mapping_model,
    residue_key
):

    (
        chain,
        mobile_resname,
        mobile_resnum,
        icode
    ) = residue_key


    # --------------------------------------------------------
    # 1Z83 already uses reference numbering
    # --------------------------------------------------------

    if mapping_model is None:

        return {

            "Mapped":
                True,

            "Reference_resnum":
                mobile_resnum,

            "Reference_resname":
                mobile_resname,

            "Reference_residue":
                f"{mobile_resname}{mobile_resnum}",

            "Mobile_resnum":
                mobile_resnum,

            "Mobile_resname":
                mobile_resname,

            "Sequence_identical":
                True
        }


    # --------------------------------------------------------
    # 7X7S -> 1Z83 mapping
    # --------------------------------------------------------

    mapping = (
        mapping_by_model[
            mapping_model
        ]
    )


    if mobile_resnum not in mapping:

        return {

            "Mapped":
                False,

            "Reference_resnum":
                None,

            "Reference_resname":
                "",

            "Reference_residue":
                "",

            "Mobile_resnum":
                mobile_resnum,

            "Mobile_resname":
                mobile_resname,

            "Sequence_identical":
                False
        }


    item = mapping[
        mobile_resnum
    ]


    return {

        "Mapped":
            True,

        "Reference_resnum":
            item[
                "Reference_resnum"
            ],

        "Reference_resname":
            item[
                "Reference_resname"
            ],

        "Reference_residue":
            (
                f"{item['Reference_resname']}"
                f"{item['Reference_resnum']}"
            ),

        "Mobile_resnum":
            mobile_resnum,

        "Mobile_resname":
            mobile_resname,

        "Sequence_identical":
            item[
                "Identical"
            ]
    }


# ============================================================
# DISTANCE
# ============================================================

def distance(
    coord_a,
    coord_b
):

    return float(
        np.linalg.norm(
            coord_a -
            coord_b
        )
    )


# ============================================================
# CALCULATE CONTACTS FOR ONE POSE
# ============================================================

def calculate_pose_contacts(
    protein_residues,
    ligand_atoms,
    group_name,
    mapping_model
):

    rows = []


    ligand_polar_atoms = [

        atom

        for atom in ligand_atoms

        if is_polar_element(
            atom[
                "element"
            ]
        )
    ]


    for residue_key, protein_atoms in (
        protein_residues.items()
    ):

        mapping = (
            map_residue_to_reference(
                group_name,
                mapping_model,
                residue_key
            )
        )


        min_any = math.inf
        min_any_pair = None

        min_polar = math.inf
        min_polar_pair = None


        # ----------------------------------------------------
        # ANY CONTACT
        # ----------------------------------------------------

        for protein_atom in protein_atoms:

            for ligand_atom in ligand_atoms:

                d = distance(
                    protein_atom[
                        "coord"
                    ],
                    ligand_atom[
                        "coord"
                    ]
                )

                if d < min_any:

                    min_any = d

                    min_any_pair = (
                        protein_atom,
                        ligand_atom
                    )


        # ----------------------------------------------------
        # POLAR PROXIMITY
        # ----------------------------------------------------

        protein_polar_atoms = [

            atom

            for atom in protein_atoms

            if is_polar_element(
                atom[
                    "element"
                ]
            )
        ]


        for protein_atom in protein_polar_atoms:

            for ligand_atom in ligand_polar_atoms:

                d = distance(
                    protein_atom[
                        "coord"
                    ],
                    ligand_atom[
                        "coord"
                    ]
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


        if not (
            any_contact
            or
            polar_contact
        ):
            continue


        rows.append({

            "Reference_residue":
                mapping[
                    "Reference_residue"
                ],

            "Reference_resnum":
                mapping[
                    "Reference_resnum"
                ],

            "Reference_resname":
                mapping[
                    "Reference_resname"
                ],

            "Mobile_resnum":
                mapping[
                    "Mobile_resnum"
                ],

            "Mobile_resname":
                mapping[
                    "Mobile_resname"
                ],

            "Mapped":
                mapping[
                    "Mapped"
                ],

            "Sequence_identical":
                mapping[
                    "Sequence_identical"
                ],

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

            "Any_protein_atom":
                (
                    min_any_pair[
                        0
                    ][
                        "atom_name"
                    ]
                    if min_any_pair
                    else ""
                ),

            "Any_ligand_atom":
                (
                    min_any_pair[
                        1
                    ][
                        "atom_name"
                    ]
                    if min_any_pair
                    else ""
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

            "Polar_protein_atom":
                (
                    min_polar_pair[
                        0
                    ][
                        "atom_name"
                    ]
                    if min_polar_pair
                    else ""
                ),

            "Polar_ligand_atom":
                (
                    min_polar_pair[
                        1
                    ][
                        "atom_name"
                    ]
                    if min_polar_pair
                    else ""
                )
        })


    return rows


# ============================================================
# RUN ALL FOUR RECEPTOR GROUPS
# ============================================================

all_pose_contact_rows = []

unmapped_contact_rows = []


print()
print("=" * 80)
print("AK1 ENSEMBLE MAPPED CONTACT ANALYSIS")
print("=" * 80)

print(
    f"Any-contact cutoff: "
    f"{ANY_CONTACT_CUTOFF:.1f} A"
)

print(
    f"Polar-proximity cutoff: "
    f"{POLAR_CONTACT_CUTOFF:.1f} A"
)


for group_name, info in (
    GROUPS.items()
):

    print()
    print("-" * 80)
    print(group_name)
    print("-" * 80)


    protein = read_protein(
        info[
            "receptor"
        ]
    )


    print(
        "Protein residues:",
        len(
            protein
        )
    )


    for seed in info[
        "seeds"
    ]:

        docking_file = (
            DOCKING
            /
            info[
                "dock_template"
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
                ligand_atoms,
                group_name,
                info[
                    "mapping_model"
                ]
            )
        )


        mapped_contacts = [
            row
            for row in contacts
            if row[
                "Mapped"
            ]
        ]


        unmapped_contacts = [
            row
            for row in contacts
            if not row[
                "Mapped"
            ]
        ]


        any_count = sum(
            row[
                "Any_contact"
            ]
            for row in mapped_contacts
        )


        polar_count = sum(
            row[
                "Polar_contact"
            ]
            for row in mapped_contacts
        )


        print(
            f"Seed {seed:<5} : "
            f"{any_count} mapped any contacts, "
            f"{polar_count} mapped polar residues"
        )


        if unmapped_contacts:

            print(
                f"             WARNING: "
                f"{len(unmapped_contacts)} unmapped contacting residues"
            )


        for row in contacts:

            output_row = {

                "Group":
                    group_name,

                "Seed":
                    seed,

                **row
            }


            all_pose_contact_rows.append(
                output_row
            )


            if not row[
                "Mapped"
            ]:

                unmapped_contact_rows.append(
                    output_row
                )


# ============================================================
# RETAIN ONLY MAPPED RESIDUES FOR CROSS-RECEPTOR COMPARISON
# ============================================================

mapped_pose_rows = [

    row

    for row in all_pose_contact_rows

    if row[
        "Mapped"
    ]
]


all_reference_residues = sorted(

    set(
        row[
            "Reference_residue"
        ]

        for row in mapped_pose_rows

        if row[
            "Reference_residue"
        ]
    ),

    key=lambda x: (
        int(
            "".join(
                c
                for c in x
                if c.isdigit()
            )
        ),
        x
    )
)


# ============================================================
# RECURRENCE SUMMARY
# ============================================================

recurrence_rows = []


for group_name, info in (
    GROUPS.items()
):

    group_seeds = (
        info[
            "seeds"
        ]
    )

    n_total = len(
        group_seeds
    )


    for reference_residue in (
        all_reference_residues
    ):

        relevant = [

            row

            for row in mapped_pose_rows

            if (
                row[
                    "Group"
                ] == group_name
                and
                row[
                    "Reference_residue"
                ] == reference_residue
            )
        ]


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


        # Find reference identity from any matching row
        reference_resnum = None
        reference_resname = ""

        if relevant:

            reference_resnum = relevant[
                0
            ][
                "Reference_resnum"
            ]

            reference_resname = relevant[
                0
            ][
                "Reference_resname"
            ]

        else:

            digits = "".join(
                c
                for c in reference_residue
                if c.isdigit()
            )

            letters = "".join(
                c
                for c in reference_residue
                if c.isalpha()
            )

            if digits:

                reference_resnum = int(
                    digits
                )

            reference_resname = letters


        recurrence_rows.append({

            "Group":
                group_name,

            "Reference_residue":
                reference_residue,

            "Reference_resnum":
                reference_resnum,

            "Reference_resname":
                reference_resname,

            "n_poses":
                n_total,

            "Any_contact_n":
                len(
                    any_rows
                ),

            "Any_contact_percent":
                (
                    100.0 *
                    len(
                        any_rows
                    ) /
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
                len(
                    polar_rows
                ),

            "Polar_contact_percent":
                (
                    100.0 *
                    len(
                        polar_rows
                    ) /
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
# RECURRENT SETS >=80%
# ============================================================

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
# PRINT RECURRENT SETS
# ============================================================

print()
print("=" * 80)
print("RECURRENT MAPPED CONTACTS (>=80%)")
print("=" * 80)


for group_name in GROUPS:

    any_set = [

        row[
            "Reference_residue"
        ]

        for row in recurrent_any_rows

        if row[
            "Group"
        ] == group_name
    ]


    polar_set = [

        row[
            "Reference_residue"
        ]

        for row in recurrent_polar_rows

        if row[
            "Group"
        ] == group_name
    ]


    print()
    print(group_name)
    print("-" * len(group_name))

    print(
        "Recurrent ANY:"
    )

    print(
        ", ".join(
            any_set
        )
        if any_set
        else "None"
    )

    print()

    print(
        "Recurrent POLAR:"
    )

    print(
        ", ".join(
            polar_set
        )
        if polar_set
        else "None"
    )


# ============================================================
# VALIDATE ORIGINAL CLOSED DOMINANT POLAR FINGERPRINT
# ============================================================

expected_closed_polar = {
    "LYS21",
    "MET61",
    "ARG97",
    "ARG132",
    "ARG138",
    "ARG149"
}


observed_closed_polar_100 = {

    row[
        "Reference_residue"
    ]

    for row in recurrence_rows

    if (
        row[
            "Group"
        ] ==
        "1Z83 closed dominant"
        and
        row[
            "Polar_contact_n"
        ] == 8
    )
}


print()
print("=" * 80)
print("1Z83 CLOSED-DOMINANT VALIDATION")
print("=" * 80)

print(
    "Expected exact 8/8 polar set:"
)

print(
    ", ".join(
        sorted(
            expected_closed_polar
        )
    )
)

print()

print(
    "Recovered exact 8/8 polar set:"
)

print(
    ", ".join(
        sorted(
            observed_closed_polar_100
        )
    )
)


if (
    observed_closed_polar_100
    ==
    expected_closed_polar
):

    print()
    print(
        "VALIDATION PASSED."
    )

else:

    print()
    print(
        "WARNING: fingerprint does not exactly match "
        "the previous closed-state analysis."
    )


# ============================================================
# OVERLAP AGAINST CLOSED DOMINANT REFERENCE
# ============================================================

def recurrent_set(
    group_name,
    contact_type
):

    if contact_type == "Any":

        rows = recurrent_any_rows

    elif contact_type == "Polar":

        rows = recurrent_polar_rows

    else:

        raise ValueError(
            contact_type
        )


    return {

        row[
            "Reference_residue"
        ]

        for row in rows

        if row[
            "Group"
        ] == group_name
    }


closed_any = recurrent_set(
    "1Z83 closed dominant",
    "Any"
)

closed_polar = recurrent_set(
    "1Z83 closed dominant",
    "Polar"
)


overlap_rows = []


for group_name in (
    "7X7S least-open",
    "7X7S medoid",
    "7X7S most-open"
):

    for contact_type in (
        "Any",
        "Polar"
    ):

        reference_set = (
            closed_any
            if contact_type == "Any"
            else closed_polar
        )

        test_set = recurrent_set(
            group_name,
            contact_type
        )


        shared = (
            reference_set
            &
            test_set
        )

        lost = (
            reference_set
            -
            test_set
        )

        gained = (
            test_set
            -
            reference_set
        )

        union_set = (
            reference_set
            |
            test_set
        )


        overlap_rows.append({

            "Comparison":
                (
                    f"{group_name} vs "
                    f"1Z83 closed dominant"
                ),

            "Contact_type":
                contact_type,

            "Closed_reference_n":
                len(
                    reference_set
                ),

            "Open_state_n":
                len(
                    test_set
                ),

            "Shared_n":
                len(
                    shared
                ),

            "Closed_reference_retained_percent":
                (
                    100.0 *
                    len(
                        shared
                    ) /
                    len(
                        reference_set
                    )
                    if reference_set
                    else 0.0
                ),

            "Open_state_overlap_percent":
                (
                    100.0 *
                    len(
                        shared
                    ) /
                    len(
                        test_set
                    )
                    if test_set
                    else 0.0
                ),

            "Jaccard":
                (
                    len(
                        shared
                    ) /
                    len(
                        union_set
                    )
                    if union_set
                    else 0.0
                ),

            "Shared_residues":
                ", ".join(
                    sorted(
                        shared
                    )
                ),

            "Lost_from_closed":
                ", ".join(
                    sorted(
                        lost
                    )
                ),

            "Gained_in_open":
                ", ".join(
                    sorted(
                        gained
                    )
                )
        })


# ============================================================
# WIDE COMPARISON MATRIX
# ============================================================

comparison_rows = []


for reference_residue in (
    all_reference_residues
):

    row_out = {
        "Reference_residue":
            reference_residue
    }


    for group_name in GROUPS:

        match = next(

            (
                row

                for row in recurrence_rows

                if (
                    row[
                        "Group"
                    ] == group_name
                    and
                    row[
                        "Reference_residue"
                    ] == reference_residue
                )
            ),

            None
        )


        if match is None:
            continue


        safe_group = (
            group_name
            .replace(" ", "_")
            .replace("-", "_")
        )


        row_out[
            f"{safe_group}_n_poses"
        ] = match[
            "n_poses"
        ]

        row_out[
            f"{safe_group}_Any_n"
        ] = match[
            "Any_contact_n"
        ]

        row_out[
            f"{safe_group}_Any_percent"
        ] = match[
            "Any_contact_percent"
        ]

        row_out[
            f"{safe_group}_Polar_n"
        ] = match[
            "Polar_contact_n"
        ]

        row_out[
            f"{safe_group}_Polar_percent"
        ] = match[
            "Polar_contact_percent"
        ]


    comparison_rows.append(
        row_out
    )


# ============================================================
# CSV WRITER
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
# OUTPUT FILES
# ============================================================

ALL_CONTACTS_CSV = (
    OUT_DIR /
    "AK1_ensemble_all_pose_contacts_mapped_to_1Z83.csv"
)

RECURRENCE_CSV = (
    OUT_DIR /
    "AK1_ensemble_contact_recurrence_mapped_to_1Z83.csv"
)

RECURRENT_ANY_CSV = (
    OUT_DIR /
    "AK1_ensemble_recurrent_ANY_contacts_80pct.csv"
)

RECURRENT_POLAR_CSV = (
    OUT_DIR /
    "AK1_ensemble_recurrent_POLAR_contacts_80pct.csv"
)

OVERLAP_CSV = (
    OUT_DIR /
    "AK1_ensemble_contact_overlap_vs_closed.csv"
)

COMPARISON_CSV = (
    OUT_DIR /
    "AK1_ensemble_contact_comparison_matrix.csv"
)

UNMAPPED_CSV = (
    OUT_DIR /
    "AK1_ensemble_unmapped_contacting_residues.csv"
)


write_csv(
    ALL_CONTACTS_CSV,
    all_pose_contact_rows
)

write_csv(
    RECURRENCE_CSV,
    recurrence_rows
)

write_csv(
    RECURRENT_ANY_CSV,
    recurrent_any_rows
)

write_csv(
    RECURRENT_POLAR_CSV,
    recurrent_polar_rows
)

write_csv(
    OVERLAP_CSV,
    overlap_rows
)

write_csv(
    COMPARISON_CSV,
    comparison_rows
)

if unmapped_contact_rows:

    write_csv(
        UNMAPPED_CSV,
        unmapped_contact_rows
    )


# ============================================================
# EXCEL WORKBOOK
# ============================================================

XLSX_FILE = (
    OUT_DIR /
    "AK1_ENSEMBLE_MAPPED_CONTACT_FINGERPRINTS.xlsx"
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


# ============================================================
# GENERIC EXCEL WRITER
# ============================================================

def write_sheet(
    name,
    rows,
    width=24
):

    ws = workbook.add_worksheet(
        name
    )

    if not rows:

        ws.write(
            0,
            0,
            "No rows"
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


# ============================================================
# WRITE EXCEL DATA
# ============================================================

write_sheet(
    "All_Pose_Contacts",
    all_pose_contact_rows,
    width=23
)

write_sheet(
    "Recurrence",
    recurrence_rows,
    width=25
)

write_sheet(
    "Recurrent_ANY_80pct",
    recurrent_any_rows,
    width=25
)

write_sheet(
    "Recurrent_POLAR_80pct",
    recurrent_polar_rows,
    width=25
)

write_sheet(
    "Overlap_vs_Closed",
    overlap_rows,
    width=32
)

write_sheet(
    "Comparison_Matrix",
    comparison_rows,
    width=23
)

write_sheet(
    "Unmapped_Contacts",
    unmapped_contact_rows,
    width=25
)


# ============================================================
# METHOD NOTES
# ============================================================

ws_notes = workbook.add_worksheet(
    "Method_Notes"
)


notes = [

    [
        "Purpose",
        (
            "Compare residue-contact fingerprints of neutral ISO "
            "across closed 1Z83 and selected open 7X7S receptor "
            "conformations."
        )
    ],

    [
        "Reference numbering",
        (
            "Contacts detected in 7X7S receptors were mapped back "
            "onto the corresponding 1Z83 residue numbers using the "
            "previously generated sequence alignment."
        )
    ],

    [
        "Closed reference",
        (
            "Eight poses belonging to the established dominant "
            "neutral ISO family in 1Z83."
        )
    ],

    [
        "Open-state groups",
        (
            "Ten best poses each for 7X7S least-open model 20, "
            "ensemble-medoid model 8 and most-open model 1."
        )
    ],

    [
        "Any-contact definition",
        (
            f"Minimum ligand-residue heavy-atom distance <= "
            f"{ANY_CONTACT_CUTOFF:.1f} A."
        )
    ],

    [
        "Polar-proximity definition",
        (
            f"Minimum O/N/S heavy-atom distance <= "
            f"{POLAR_CONTACT_CUTOFF:.1f} A. "
            "These are proximity contacts and are not automatically "
            "classified as hydrogen bonds."
        )
    ],

    [
        "Recurrent contact",
        (
            f"Present in >= {RECURRENCE_THRESHOLD:.0f}% "
            "of poses within the receptor-state group."
        )
    ],

    [
        "Sequence mismatches",
        (
            "Residues were mapped by aligned sequence position even "
            "when the residue identity differed between structures. "
            "Sequence identity status is preserved in the per-pose "
            "contact table."
        )
    ],

    [
        "Interpretation",
        (
            "Contact conservation across receptor states supports "
            "persistence of the same functional residue neighborhood "
            "even when the exact ligand orientation changes. "
            "Docking contacts do not establish experimental binding "
            "or inhibition."
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
# PRINT OVERLAP SUMMARY
# ============================================================

print()
print("=" * 80)
print("OPEN-STATE CONTACT OVERLAP VS CLOSED DOMINANT")
print("=" * 80)


for row in overlap_rows:

    print()
    print(
        row[
            "Comparison"
        ],
        "-",
        row[
            "Contact_type"
        ]
    )

    print(
        f"  closed reference residues : "
        f"{row['Closed_reference_n']}"
    )

    print(
        f"  open-state residues       : "
        f"{row['Open_state_n']}"
    )

    print(
        f"  shared                    : "
        f"{row['Shared_n']}"
    )

    print(
        f"  closed fingerprint retained : "
        f"{row['Closed_reference_retained_percent']:.1f}%"
    )

    print(
        f"  Jaccard                   : "
        f"{row['Jaccard']:.3f}"
    )

    print(
        f"  shared residues           : "
        f"{row['Shared_residues'] or 'None'}"
    )

    print(
        f"  lost from closed          : "
        f"{row['Lost_from_closed'] or 'None'}"
    )

    print(
        f"  gained in open            : "
        f"{row['Gained_in_open'] or 'None'}"
    )


# ============================================================
# FINAL
# ============================================================

print()
print("=" * 80)
print("FILES SAVED")
print("=" * 80)

print(
    XLSX_FILE
)

print(
    ALL_CONTACTS_CSV
)

print(
    RECURRENCE_CSV
)

print(
    RECURRENT_ANY_CSV
)

print(
    RECURRENT_POLAR_CSV
)

print(
    OVERLAP_CSV
)

print(
    COMPARISON_CSV
)

if unmapped_contact_rows:

    print(
        UNMAPPED_CSV
    )

else:

    print(
        "No unmapped contacting residues were detected."
    )


print()
print(
    "Mapped ensemble contact analysis complete."
)