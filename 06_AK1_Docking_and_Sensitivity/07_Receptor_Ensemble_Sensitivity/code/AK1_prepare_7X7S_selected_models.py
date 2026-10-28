from pathlib import Path
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

STRUCTURES = (
    PROJECT /
    "structures"
)

NMR_FILE = (
    STRUCTURES /
    "AK1_7X7S_NMR.pdb"
)

REFERENCE_FILE = (
    STRUCTURES /
    "AK1_1Z83.pdb"
)

OUT_DIR = (
    PROJECT /
    "results" /
    "AK1_ensemble_sensitivity" /
    "7X7S_selected_models"
)

OUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# SELECTED 7X7S MODELS
# ============================================================

SELECTED_MODELS = {
    20: "least_open",
    8: "ensemble_medoid",
    1: "most_open"
}


# ============================================================
# OPERATIONAL CORE DEFINITION
#
# Defined in the 1Z83 residue-numbering frame.
#
# Mobile NMP:
#   approximately 30-59
#
# Mobile LID:
#   approximately 122-159
#
# Alignment CORE:
#   1-29
#   60-121
#   160-194
# ============================================================

def is_core_reference_residue(
    resnum
):

    return (
        1 <= resnum <= 29
        or
        60 <= resnum <= 121
        or
        160 <= resnum <= 194
    )


# ============================================================
# IMPORTANT 1Z83 POCKET RESIDUES
# ============================================================

KEY_REFERENCE_RESIDUES = [
    15,
    16,
    17,
    18,
    21,
    22,
    38,
    61,
    64,
    93,
    94,
    97,
    132,
    138,
    140,
    149,
    153
]


# ============================================================
# AMINO-ACID CODE
# ============================================================

AA3_TO_1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "MSE": "M"
}


# ============================================================
# READ SINGLE-STRUCTURE CHAIN A
# ============================================================

def read_reference_structure(
    path
):

    atoms = []
    ca_ordered = []
    seen_ca = set()

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

            if chain != "A":
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

            x = float(
                line[30:38]
            )

            y = float(
                line[38:46]
            )

            z = float(
                line[46:54]
            )

            atom = {
                "line":
                    line.rstrip("\n"),

                "atom_name":
                    atom_name,

                "resname":
                    resname,

                "resnum":
                    resnum,

                "icode":
                    icode,

                "coord":
                    np.array(
                        [x, y, z],
                        dtype=float
                    )
            }

            atoms.append(
                atom
            )

            key = (
                resnum,
                icode
            )

            if (
                atom_name == "CA"
                and
                key not in seen_ca
            ):

                seen_ca.add(
                    key
                )

                ca_ordered.append({

                    "resnum":
                        resnum,

                    "icode":
                        icode,

                    "resname":
                        resname,

                    "coord":
                        atom[
                            "coord"
                        ]
                })

    if not atoms:

        raise RuntimeError(
            "No chain-A atoms read from 1Z83."
        )

    return (
        atoms,
        ca_ordered
    )


# ============================================================
# READ ONE 7X7S MODEL
# ============================================================

def read_nmr_model(
    path,
    wanted_model
):

    atoms = []
    ca_ordered = []
    seen_ca = set()

    current_model = None

    with open(
        path,
        "r",
        encoding="utf-8",
        errors="ignore"
    ) as handle:

        for line in handle:

            if line.startswith(
                "MODEL"
            ):

                fields = (
                    line.split()
                )

                current_model = int(
                    fields[1]
                )

                continue

            if line.startswith(
                "ENDMDL"
            ):

                if current_model == wanted_model:
                    break

                current_model = None
                continue

            if current_model != wanted_model:
                continue

            if not line.startswith(
                "ATOM"
            ):
                continue

            chain = (
                line[21]
                .strip()
            )

            if chain != "A":
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

            x = float(
                line[30:38]
            )

            y = float(
                line[38:46]
            )

            z = float(
                line[46:54]
            )

            atom = {
                "line":
                    line.rstrip("\n"),

                "atom_name":
                    atom_name,

                "resname":
                    resname,

                "resnum":
                    resnum,

                "icode":
                    icode,

                "coord":
                    np.array(
                        [x, y, z],
                        dtype=float
                    )
            }

            atoms.append(
                atom
            )

            key = (
                resnum,
                icode
            )

            if (
                atom_name == "CA"
                and
                key not in seen_ca
            ):

                seen_ca.add(
                    key
                )

                ca_ordered.append({

                    "resnum":
                        resnum,

                    "icode":
                        icode,

                    "resname":
                        resname,

                    "coord":
                        atom[
                            "coord"
                        ]
                })

    if not atoms:

        raise RuntimeError(
            f"No atoms found for 7X7S MODEL {wanted_model}."
        )

    return (
        atoms,
        ca_ordered
    )


# ============================================================
# CONVERT ORDERED RESIDUES TO SEQUENCE
# ============================================================

def sequence_from_residues(
    residues
):

    letters = []

    for residue in residues:

        resname = (
            residue[
                "resname"
            ]
        )

        letters.append(
            AA3_TO_1.get(
                resname,
                "X"
            )
        )

    return "".join(
        letters
    )


# ============================================================
# GLOBAL NEEDLEMAN-WUNSCH ALIGNMENT
#
# Simple sequence alignment used only to establish
# structural residue correspondence.
# ============================================================

def global_align(
    seq_a,
    seq_b,
    match_score=2,
    mismatch_score=-1,
    gap_score=-2
):

    n = len(
        seq_a
    )

    m = len(
        seq_b
    )

    score = np.zeros(
        (
            n + 1,
            m + 1
        ),
        dtype=int
    )

    traceback = np.zeros(
        (
            n + 1,
            m + 1
        ),
        dtype=np.int8
    )


    # Directions:
    # 1 = diagonal
    # 2 = up
    # 3 = left

    for i in range(
        1,
        n + 1
    ):

        score[
            i,
            0
        ] = (
            i *
            gap_score
        )

        traceback[
            i,
            0
        ] = 2


    for j in range(
        1,
        m + 1
    ):

        score[
            0,
            j
        ] = (
            j *
            gap_score
        )

        traceback[
            0,
            j
        ] = 3


    for i in range(
        1,
        n + 1
    ):

        for j in range(
            1,
            m + 1
        ):

            diagonal = (
                score[
                    i - 1,
                    j - 1
                ]
                +
                (
                    match_score
                    if seq_a[
                        i - 1
                    ] == seq_b[
                        j - 1
                    ]
                    else mismatch_score
                )
            )

            up = (
                score[
                    i - 1,
                    j
                ]
                +
                gap_score
            )

            left = (
                score[
                    i,
                    j - 1
                ]
                +
                gap_score
            )

            best = max(
                diagonal,
                up,
                left
            )

            score[
                i,
                j
            ] = best

            if best == diagonal:

                traceback[
                    i,
                    j
                ] = 1

            elif best == up:

                traceback[
                    i,
                    j
                ] = 2

            else:

                traceback[
                    i,
                    j
                ] = 3


    aligned_a = []
    aligned_b = []

    index_pairs = []

    i = n
    j = m


    while (
        i > 0
        or
        j > 0
    ):

        direction = (
            traceback[
                i,
                j
            ]
        )

        if (
            i > 0
            and
            j > 0
            and
            direction == 1
        ):

            aligned_a.append(
                seq_a[
                    i - 1
                ]
            )

            aligned_b.append(
                seq_b[
                    j - 1
                ]
            )

            index_pairs.append(
                (
                    i - 1,
                    j - 1
                )
            )

            i -= 1
            j -= 1

        elif (
            i > 0
            and
            (
                j == 0
                or
                direction == 2
            )
        ):

            aligned_a.append(
                seq_a[
                    i - 1
                ]
            )

            aligned_b.append(
                "-"
            )

            index_pairs.append(
                (
                    i - 1,
                    None
                )
            )

            i -= 1

        else:

            aligned_a.append(
                "-"
            )

            aligned_b.append(
                seq_b[
                    j - 1
                ]
            )

            index_pairs.append(
                (
                    None,
                    j - 1
                )
            )

            j -= 1


    aligned_a.reverse()
    aligned_b.reverse()
    index_pairs.reverse()


    return (
        "".join(
            aligned_a
        ),
        "".join(
            aligned_b
        ),
        index_pairs,
        int(
            score[
                n,
                m
            ]
        )
    )


# ============================================================
# BUILD RESIDUE CORRESPONDENCE
# ============================================================

def build_residue_mapping(
    reference_residues,
    mobile_residues
):

    ref_seq = (
        sequence_from_residues(
            reference_residues
        )
    )

    mob_seq = (
        sequence_from_residues(
            mobile_residues
        )
    )

    (
        aligned_ref,
        aligned_mob,
        index_pairs,
        alignment_score
    ) = global_align(
        ref_seq,
        mob_seq
    )


    mapping_rows = []

    aligned_pairs = []

    identical_count = 0
    paired_count = 0


    for alignment_position, (
        ref_index,
        mob_index
    ) in enumerate(
        index_pairs,
        start=1
    ):

        if (
            ref_index is None
            or
            mob_index is None
        ):
            continue

        ref_residue = (
            reference_residues[
                ref_index
            ]
        )

        mob_residue = (
            mobile_residues[
                mob_index
            ]
        )

        ref_letter = (
            ref_seq[
                ref_index
            ]
        )

        mob_letter = (
            mob_seq[
                mob_index
            ]
        )

        identical = (
            ref_letter ==
            mob_letter
        )

        paired_count += 1

        if identical:
            identical_count += 1


        row = {

            "Alignment_position":
                alignment_position,

            "Reference_index":
                ref_index + 1,

            "Reference_resnum":
                ref_residue[
                    "resnum"
                ],

            "Reference_resname":
                ref_residue[
                    "resname"
                ],

            "Reference_AA":
                ref_letter,

            "Mobile_index":
                mob_index + 1,

            "Mobile_resnum":
                mob_residue[
                    "resnum"
                ],

            "Mobile_resname":
                mob_residue[
                    "resname"
                ],

            "Mobile_AA":
                mob_letter,

            "Identical":
                identical
        }

        mapping_rows.append(
            row
        )

        aligned_pairs.append(
            (
                ref_residue,
                mob_residue,
                identical
            )
        )


    sequence_identity = (
        100.0 *
        identical_count /
        paired_count
        if paired_count > 0
        else 0.0
    )


    return {
        "reference_sequence":
            ref_seq,

        "mobile_sequence":
            mob_seq,

        "aligned_reference":
            aligned_ref,

        "aligned_mobile":
            aligned_mob,

        "alignment_score":
            alignment_score,

        "paired_count":
            paired_count,

        "identical_count":
            identical_count,

        "identity_percent":
            sequence_identity,

        "mapping_rows":
            mapping_rows,

        "aligned_pairs":
            aligned_pairs
    }


# ============================================================
# KABSCH TRANSFORMATION
# ============================================================

def calculate_transform(
    mobile_coords,
    reference_coords
):

    mobile_center = np.mean(
        mobile_coords,
        axis=0
    )

    reference_center = np.mean(
        reference_coords,
        axis=0
    )

    mobile_zero = (
        mobile_coords -
        mobile_center
    )

    reference_zero = (
        reference_coords -
        reference_center
    )

    covariance = (
        mobile_zero.T
        @
        reference_zero
    )

    u, s, vt = np.linalg.svd(
        covariance
    )

    rotation = (
        u @ vt
    )

    if np.linalg.det(
        rotation
    ) < 0:

        vt[
            -1,
            :
        ] *= -1

        rotation = (
            u @ vt
        )


    return (
        rotation,
        mobile_center,
        reference_center
    )


# ============================================================
# APPLY TRANSFORMATION
# ============================================================

def apply_transform(
    coord,
    rotation,
    mobile_center,
    reference_center
):

    return (
        (
            coord -
            mobile_center
        )
        @
        rotation
        +
        reference_center
    )


# ============================================================
# RMSD
# ============================================================

def rmsd(
    coordinates_a,
    coordinates_b
):

    difference = (
        coordinates_a -
        coordinates_b
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
# PDB OUTPUT
# ============================================================

def replace_coordinates(
    line,
    coord
):

    x, y, z = coord

    return (
        line[:30]
        +
        f"{x:8.3f}"
        +
        f"{y:8.3f}"
        +
        f"{z:8.3f}"
        +
        line[54:]
    )


def write_raw_pdb(
    path,
    atoms,
    title
):

    with open(
        path,
        "w",
        encoding="utf-8"
    ) as handle:

        handle.write(
            f"REMARK {title}\n"
        )

        for atom in atoms:

            handle.write(
                atom[
                    "line"
                ]
                +
                "\n"
            )

        handle.write(
            "TER\nEND\n"
        )


def write_aligned_pdb(
    path,
    atoms,
    rotation,
    mobile_center,
    reference_center,
    title
):

    with open(
        path,
        "w",
        encoding="utf-8"
    ) as handle:

        handle.write(
            f"REMARK {title}\n"
        )

        handle.write(
            "REMARK Sequence-mapped CORE alignment to "
            "1Z83 chain A\n"
        )

        for atom in atoms:

            transformed = (
                apply_transform(
                    atom[
                        "coord"
                    ],
                    rotation,
                    mobile_center,
                    reference_center
                )
            )

            new_line = (
                replace_coordinates(
                    atom[
                        "line"
                    ],
                    transformed
                )
            )

            handle.write(
                new_line +
                "\n"
            )

        handle.write(
            "TER\nEND\n"
        )


# ============================================================
# READ REFERENCE 1Z83
# ============================================================

(
    reference_atoms,
    reference_ca
) = read_reference_structure(
    REFERENCE_FILE
)


reference_sequence = (
    sequence_from_residues(
        reference_ca
    )
)


print()
print("=" * 80)
print("1Z83 REFERENCE")
print("=" * 80)

print(
    "Chain A C-alpha residues:",
    len(
        reference_ca
    )
)

print(
    "First residue:",
    reference_ca[
        0
    ][
        "resnum"
    ],
    reference_ca[
        0
    ][
        "resname"
    ]
)

print(
    "Last residue:",
    reference_ca[
        -1
    ][
        "resnum"
    ],
    reference_ca[
        -1
    ][
        "resname"
    ]
)


# ============================================================
# RESULT CONTAINERS
# ============================================================

summary_rows = []
all_mapping_rows = []
pocket_mapping_rows = []


# ============================================================
# PROCESS SELECTED 7X7S MODELS
# ============================================================

for model_number, model_label in (
    SELECTED_MODELS.items()
):

    print()
    print("=" * 80)
    print(
        f"PROCESSING 7X7S MODEL "
        f"{model_number} ({model_label})"
    )
    print("=" * 80)


    (
        model_atoms,
        model_ca
    ) = read_nmr_model(
        NMR_FILE,
        model_number
    )


    # ========================================================
    # SEQUENCE MAPPING
    # ========================================================

    mapping_result = (
        build_residue_mapping(
            reference_ca,
            model_ca
        )
    )


    print(
        "1Z83 sequence length:",
        len(
            mapping_result[
                "reference_sequence"
            ]
        )
    )

    print(
        "7X7S sequence length:",
        len(
            mapping_result[
                "mobile_sequence"
            ]
        )
    )

    print(
        "Aligned residue pairs:",
        mapping_result[
            "paired_count"
        ]
    )

    print(
        "Identical aligned residues:",
        mapping_result[
            "identical_count"
        ]
    )

    print(
        f"Sequence identity: "
        f"{mapping_result['identity_percent']:.2f}%"
    )


    # ========================================================
    # BUILD CORE STRUCTURAL PAIRS
    #
    # Use sequence-corresponding residues.
    # Require same amino-acid identity.
    # Domain selection uses 1Z83 numbering.
    # ========================================================

    core_pairs = []

    all_identical_pairs = []


    for (
        ref_residue,
        mob_residue,
        identical
    ) in mapping_result[
        "aligned_pairs"
    ]:

        if not identical:
            continue

        all_identical_pairs.append(
            (
                ref_residue,
                mob_residue
            )
        )

        if is_core_reference_residue(
            ref_residue[
                "resnum"
            ]
        ):

            core_pairs.append(
                (
                    ref_residue,
                    mob_residue
                )
            )


    if len(
        core_pairs
    ) < 80:

        raise RuntimeError(
            f"Only {len(core_pairs)} sequence-mapped "
            f"CORE residues for model {model_number}; "
            f"alignment is still suspicious."
        )


    mobile_core = np.array(
        [
            mob[
                "coord"
            ]
            for ref, mob
            in core_pairs
        ],
        dtype=float
    )


    reference_core = np.array(
        [
            ref[
                "coord"
            ]
            for ref, mob
            in core_pairs
        ],
        dtype=float
    )


    # ========================================================
    # ALIGN CORE
    # ========================================================

    (
        rotation,
        mobile_center,
        reference_center
    ) = calculate_transform(
        mobile_core,
        reference_core
    )


    aligned_core = np.array(
        [
            apply_transform(
                coord,
                rotation,
                mobile_center,
                reference_center
            )
            for coord
            in mobile_core
        ],
        dtype=float
    )


    core_rmsd = rmsd(
        aligned_core,
        reference_core
    )


    # ========================================================
    # WHOLE CHAIN RMSD AFTER CORE ALIGNMENT
    # ========================================================

    mobile_all = np.array(
        [
            mob[
                "coord"
            ]
            for ref, mob
            in all_identical_pairs
        ],
        dtype=float
    )


    reference_all = np.array(
        [
            ref[
                "coord"
            ]
            for ref, mob
            in all_identical_pairs
        ],
        dtype=float
    )


    aligned_all = np.array(
        [
            apply_transform(
                coord,
                rotation,
                mobile_center,
                reference_center
            )
            for coord
            in mobile_all
        ],
        dtype=float
    )


    whole_rmsd = rmsd(
        aligned_all,
        reference_all
    )


    # ========================================================
    # POCKET RESIDUE MAPPING
    # ========================================================

    reference_to_mobile = {}


    for (
        ref_residue,
        mob_residue,
        identical
    ) in mapping_result[
        "aligned_pairs"
    ]:

        reference_to_mobile[
            ref_residue[
                "resnum"
            ]
        ] = (
            mob_residue,
            identical
        )


    missing_key = []


    for ref_resnum in (
        KEY_REFERENCE_RESIDUES
    ):

        if ref_resnum not in (
            reference_to_mobile
        ):

            missing_key.append(
                ref_resnum
            )

            pocket_mapping_rows.append({

                "Model":
                    model_number,

                "Model_label":
                    model_label,

                "Reference_1Z83_resnum":
                    ref_resnum,

                "Reference_resname":
                    "",

                "7X7S_resnum":
                    "",

                "7X7S_resname":
                    "",

                "Sequence_identical":
                    False,

                "Mapped":
                    False
            })

            continue


        (
            mob_residue,
            identical
        ) = reference_to_mobile[
            ref_resnum
        ]


        ref_match = next(

            residue

            for residue in reference_ca

            if residue[
                "resnum"
            ] == ref_resnum
        )


        pocket_mapping_rows.append({

            "Model":
                model_number,

            "Model_label":
                model_label,

            "Reference_1Z83_resnum":
                ref_resnum,

            "Reference_resname":
                ref_match[
                    "resname"
                ],

            "7X7S_resnum":
                mob_residue[
                    "resnum"
                ],

            "7X7S_resname":
                mob_residue[
                    "resname"
                ],

            "Sequence_identical":
                identical,

            "Mapped":
                True
        })


    # ========================================================
    # SAVE FULL MAPPING
    # ========================================================

    for row in mapping_result[
        "mapping_rows"
    ]:

        row_copy = dict(
            row
        )

        row_copy[
            "Model"
        ] = model_number

        row_copy[
            "Model_label"
        ] = model_label

        all_mapping_rows.append(
            row_copy
        )


    # ========================================================
    # OUTPUT FILES
    # ========================================================

    raw_file = (
        STRUCTURES /
        f"AK1_7X7S_model{model_number}_"
        f"{model_label}_RAW.pdb"
    )


    aligned_file = (
        STRUCTURES /
        f"AK1_7X7S_model{model_number}_"
        f"{model_label}_ALIGNED_1Z83.pdb"
    )


    write_raw_pdb(
        raw_file,
        model_atoms,
        (
            f"7X7S model {model_number} "
            f"({model_label}) raw coordinates"
        )
    )


    write_aligned_pdb(
        aligned_file,
        model_atoms,
        rotation,
        mobile_center,
        reference_center,
        (
            f"7X7S model {model_number} "
            f"({model_label}) aligned to 1Z83"
        )
    )


    # ========================================================
    # SUMMARY
    # ========================================================

    summary_rows.append({

        "Model":
            model_number,

        "Label":
            model_label,

        "1Z83_CA_count":
            len(
                reference_ca
            ),

        "7X7S_CA_count":
            len(
                model_ca
            ),

        "Sequence_aligned_pairs":
            mapping_result[
                "paired_count"
            ],

        "Sequence_identical_pairs":
            mapping_result[
                "identical_count"
            ],

        "Sequence_identity_percent":
            mapping_result[
                "identity_percent"
            ],

        "CORE_CA_pairs":
            len(
                core_pairs
            ),

        "CORE_CA_RMSD_after_alignment_A":
            core_rmsd,

        "Whole_chain_identical_CA_pairs":
            len(
                all_identical_pairs
            ),

        "Whole_chain_CA_RMSD_after_CORE_alignment_A":
            whole_rmsd,

        "Missing_key_reference_residues":
            (
                ", ".join(
                    str(x)
                    for x
                    in missing_key
                )
                if missing_key
                else ""
            ),

        "Raw_PDB":
            str(
                raw_file
            ),

        "Aligned_PDB":
            str(
                aligned_file
            )
    })


    # ========================================================
    # CONSOLE
    # ========================================================

    print()
    print(
        "Sequence-mapped CORE CA pairs:",
        len(
            core_pairs
        )
    )

    print(
        f"CORE RMSD after alignment: "
        f"{core_rmsd:.3f} A"
    )

    print(
        f"Whole-chain RMSD after CORE alignment: "
        f"{whole_rmsd:.3f} A"
    )

    print(
        "Missing key 1Z83 pocket residues:",
        (
            ", ".join(
                str(x)
                for x in missing_key
            )
            if missing_key
            else "None"
        )
    )

    print()
    print(
        "Aligned PDB:"
    )

    print(
        aligned_file
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

SUMMARY_CSV = (
    OUT_DIR /
    "AK1_7X7S_selected_model_alignment_summary.csv"
)

MAPPING_CSV = (
    OUT_DIR /
    "AK1_7X7S_to_1Z83_sequence_mapping.csv"
)

POCKET_CSV = (
    OUT_DIR /
    "AK1_7X7S_pocket_residue_mapping.csv"
)


write_csv(
    SUMMARY_CSV,
    summary_rows
)

write_csv(
    MAPPING_CSV,
    all_mapping_rows
)

write_csv(
    POCKET_CSV,
    pocket_mapping_rows
)


# ============================================================
# EXCEL
# ============================================================

XLSX_FILE = (
    OUT_DIR /
    "AK1_7X7S_SELECTED_MODEL_PREPARATION_V2.xlsx"
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
# GENERIC SHEET WRITER
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
# WRITE DATA SHEETS
# ============================================================

write_sheet(
    "Alignment_Summary",
    summary_rows,
    width=30
)

write_sheet(
    "Sequence_Mapping",
    all_mapping_rows,
    width=22
)

write_sheet(
    "Pocket_Residue_Mapping",
    pocket_mapping_rows,
    width=24
)


# ============================================================
# METHOD NOTES
# ============================================================

ws_notes = workbook.add_worksheet(
    "Method_Notes"
)


notes = [

    [
        "Problem in V1",
        (
            "The first script matched residues using raw PDB "
            "residue numbers plus residue identity. 1Z83 and "
            "7X7S use different residue-numbering schemes, "
            "so only six residues happened to match directly."
        )
    ],

    [
        "Correction",
        (
            "V2 establishes 1Z83-to-7X7S residue correspondence "
            "using global amino-acid sequence alignment before "
            "performing structural alignment."
        )
    ],

    [
        "Selected conformers",
        (
            "7X7S model 8 = ensemble medoid; model 20 = "
            "least-open ensemble member; model 1 = most-open "
            "ensemble member."
        )
    ],

    [
        "CORE definition",
        (
            "CORE selection is defined in the 1Z83 residue-number "
            "frame as residues 1-29, 60-121 and 160-194. "
            "Corresponding 7X7S residues are obtained from the "
            "sequence alignment."
        )
    ],

    [
        "Structural alignment",
        (
            "Optimal Kabsch rigid-body alignment of sequence-"
            "corresponding identical CORE C-alpha atoms."
        )
    ],

    [
        "Why rigid-body alignment",
        (
            "The transformation places each open 7X7S conformer "
            "into the 1Z83 Cartesian frame without changing its "
            "internal open-state geometry."
        )
    ],

    [
        "Docking consequence",
        (
            "Once aligned, the same physical 1Z83 nucleotide-region "
            "search box can be meaningfully applied to the open "
            "7X7S receptors."
        )
    ],

    [
        "Next step",
        (
            "Inspect aligned receptors and residue mapping, then "
            "prepare each model with the same protonation/receptor "
            "preparation pipeline used for the original AK1 docking."
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
# FINAL CONSOLE SUMMARY
# ============================================================

print()
print("=" * 80)
print("7X7S SELECTED-MODEL SEQUENCE-MAPPED ALIGNMENT COMPLETE")
print("=" * 80)


for row in summary_rows:

    print()

    print(
        f"Model {row['Model']} "
        f"({row['Label']})"
    )

    print(
        f"  Sequence identity : "
        f"{row['Sequence_identity_percent']:.2f}%"
    )

    print(
        f"  CORE CA pairs     : "
        f"{row['CORE_CA_pairs']}"
    )

    print(
        f"  CORE RMSD         : "
        f"{row['CORE_CA_RMSD_after_alignment_A']:.3f} A"
    )

    print(
        f"  Whole-chain RMSD  : "
        f"{row['Whole_chain_CA_RMSD_after_CORE_alignment_A']:.3f} A"
    )

    print(
        f"  Missing key pocket residues: "
        f"{row['Missing_key_reference_residues'] or 'None'}"
    )

    print(
        f"  Aligned PDB: "
        f"{row['Aligned_PDB']}"
    )


print()
print("Saved:")
print(SUMMARY_CSV)
print(MAPPING_CSV)
print(POCKET_CSV)
print(XLSX_FILE)

print()
print(
    "Preparation V2 complete."
)