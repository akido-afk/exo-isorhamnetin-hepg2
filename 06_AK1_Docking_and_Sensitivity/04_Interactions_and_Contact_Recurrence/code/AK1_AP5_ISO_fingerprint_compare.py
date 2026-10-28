from pathlib import Path
from collections import defaultdict
import math
import csv
import xlsxwriter


# ============================================================
# PATHS
# ============================================================

PROJECT = Path(r"C:\Users\akidi\Desktop\ISO_AK_Docking")

PROTEIN_FILE = PROJECT / "receptors" / "AK1_chainA_FIXED.pdb"

AP5_FILE = PROJECT / "ligand" / "AK1_native_AP5_chainA.pdb"

ISO_FREQ_FILE = (
    PROJECT
    / "results"
    / "AK1_dominant8_contact_frequency"
    / "AK1_dominant8_contact_frequency.csv"
)

OUTPUT_DIR = (
    PROJECT
    / "results"
    / "AK1_AP5_ISO_fingerprint"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

CSV_OUTPUT = (
    OUTPUT_DIR
    / "AK1_AP5_vs_ISO_residue_fingerprint.csv"
)

XLSX_OUTPUT = (
    OUTPUT_DIR
    / "AK1_AP5_vs_ISO_residue_fingerprint.xlsx"
)


# ============================================================
# DEFINITIONS
# ============================================================

ANY_CUTOFF = 4.0
POLAR_CUTOFF = 3.6


def distance(a, b):

    return math.sqrt(
        (a[0] - b[0]) ** 2
        + (a[1] - b[1]) ** 2
        + (a[2] - b[2]) ** 2
    )


# ============================================================
# READ AK1 PROTEIN ATOMS
# ============================================================

protein_atoms = []

with open(
    PROTEIN_FILE,
    "r"
) as f:

    for line in f:

        if not line.startswith("ATOM"):
            continue

        atom_name = line[12:16].strip()
        resname = line[17:20].strip()
        resid = line[22:26].strip()
        chain = line[21].strip()

        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])

        element = line[76:78].strip()

        if not element:
            element = atom_name[0]

        element = element.upper()

        if element == "H":
            continue

        residue = f"{resname}{resid}"

        protein_atoms.append({
            "residue": residue,
            "atom": atom_name,
            "element": element,
            "chain": chain,
            "coord": (x, y, z)
        })


# ============================================================
# READ NATIVE AP5 ATOMS
# ============================================================

ap5_atoms = []

with open(
    AP5_FILE,
    "r"
) as f:

    for line in f:

        if not line.startswith("HETATM"):
            continue

        atom_name = line[12:16].strip()

        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])

        element = line[76:78].strip()

        if not element:
            element = atom_name[0]

        element = element.upper()

        if element == "H":
            continue

        ap5_atoms.append({
            "atom": atom_name,
            "element": element,
            "coord": (x, y, z)
        })


# ============================================================
# AP5 CONTACT FINGERPRINT
# ============================================================

ap5_any = {}
ap5_polar = {}

for p in protein_atoms:

    residue = p["residue"]

    for ligand_atom in ap5_atoms:

        d = distance(
            p["coord"],
            ligand_atom["coord"]
        )

        # ----------------------------
        # Any heavy-atom contact
        # ----------------------------

        if d <= ANY_CUTOFF:

            if (
                residue not in ap5_any
                or d < ap5_any[residue]
            ):
                ap5_any[residue] = d

        # ----------------------------
        # Polar heavy-atom proximity
        # ----------------------------

        if (
            p["element"] in {"O", "N", "S"}
            and
            ligand_atom["element"] in {"O", "N", "S"}
            and
            d <= POLAR_CUTOFF
        ):

            if (
                residue not in ap5_polar
                or d < ap5_polar[residue]
            ):
                ap5_polar[residue] = d


# ============================================================
# READ ISO DOMINANT-CLUSTER FREQUENCIES
# ============================================================

iso_rows = {}

with open(
    ISO_FREQ_FILE,
    newline="",
    encoding="utf-8"
) as f:

    reader = csv.DictReader(f)

    for row in reader:

        residue = row["Residue"]

        iso_rows[residue] = {
            "any_n": int(row["Any contact n"]),
            "any_pct": float(row["Any contact %"]),
            "polar_n": int(row["Polar contact n"]),
            "polar_pct": float(row["Polar contact %"])
        }


# ============================================================
# UNION OF ALL RESIDUES
# ============================================================

all_residues = sorted(
    set(ap5_any.keys())
    | set(ap5_polar.keys())
    | set(iso_rows.keys()),
    key=lambda r: int(
        "".join(
            c for c in r
            if c.isdigit()
        )
    )
)


# ============================================================
# BUILD COMPARISON TABLE
# ============================================================

rows = []

for residue in all_residues:

    ap5_any_yes = residue in ap5_any
    ap5_polar_yes = residue in ap5_polar

    iso_any_n = (
        iso_rows.get(
            residue,
            {}
        ).get(
            "any_n",
            0
        )
    )

    iso_any_pct = (
        iso_rows.get(
            residue,
            {}
        ).get(
            "any_pct",
            0.0
        )
    )

    iso_polar_n = (
        iso_rows.get(
            residue,
            {}
        ).get(
            "polar_n",
            0
        )
    )

    iso_polar_pct = (
        iso_rows.get(
            residue,
            {}
        ).get(
            "polar_pct",
            0.0
        )
    )

    iso_consensus_any = (
        iso_any_n >= 6
    )

    iso_consensus_polar = (
        iso_polar_n >= 6
    )

    if ap5_any_yes and iso_consensus_any:

        any_class = "Shared"

    elif ap5_any_yes:

        any_class = "AP5 only"

    elif iso_consensus_any:

        any_class = "ISO only"

    else:

        any_class = "Neither consensus"


    if ap5_polar_yes and iso_consensus_polar:

        polar_class = "Shared"

    elif ap5_polar_yes:

        polar_class = "AP5 only"

    elif iso_consensus_polar:

        polar_class = "ISO only"

    else:

        polar_class = "Neither consensus"


    rows.append({
        "Residue": residue,

        "AP5 any contact": int(
            ap5_any_yes
        ),

        "AP5 any min distance (A)": (
            ap5_any.get(residue)
        ),

        "ISO any contact n/8": (
            iso_any_n
        ),

        "ISO any contact %": (
            iso_any_pct
        ),

        "Any-contact class": (
            any_class
        ),

        "AP5 polar contact": int(
            ap5_polar_yes
        ),

        "AP5 polar min distance (A)": (
            ap5_polar.get(residue)
        ),

        "ISO polar contact n/8": (
            iso_polar_n
        ),

        "ISO polar contact %": (
            iso_polar_pct
        ),

        "Polar-contact class": (
            polar_class
        )
    })


# ============================================================
# SUMMARY SETS
# ============================================================

ap5_any_set = {
    r for r in all_residues
    if r in ap5_any
}

iso_consensus_any_set = {
    r for r in all_residues
    if iso_rows.get(
        r,
        {}
    ).get(
        "any_n",
        0
    ) >= 6
}

shared_any = (
    ap5_any_set
    & iso_consensus_any_set
)

ap5_only_any = (
    ap5_any_set
    - iso_consensus_any_set
)

iso_only_any = (
    iso_consensus_any_set
    - ap5_any_set
)


ap5_polar_set = {
    r for r in all_residues
    if r in ap5_polar
}

iso_consensus_polar_set = {
    r for r in all_residues
    if iso_rows.get(
        r,
        {}
    ).get(
        "polar_n",
        0
    ) >= 6
}

shared_polar = (
    ap5_polar_set
    & iso_consensus_polar_set
)

ap5_only_polar = (
    ap5_polar_set
    - iso_consensus_polar_set
)

iso_only_polar = (
    iso_consensus_polar_set
    - ap5_polar_set
)


# ============================================================
# JACCARD OVERLAP
# ============================================================

def jaccard(a, b):

    union = a | b

    if not union:
        return 0.0

    return (
        len(a & b)
        / len(union)
    )


any_jaccard = jaccard(
    ap5_any_set,
    iso_consensus_any_set
)

polar_jaccard = jaccard(
    ap5_polar_set,
    iso_consensus_polar_set
)


# ============================================================
# PRINT RESULTS
# ============================================================

print()
print("=" * 90)
print(
    "AK1 NATIVE AP5 vs ISO DOMINANT-CLUSTER INTERACTION FINGERPRINT"
)
print("=" * 90)

print()
print("ANY HEAVY-ATOM CONTACTS")
print("-----------------------")

print(
    "AP5 residues:",
    sorted(ap5_any_set)
)

print(
    "ISO consensus residues:",
    sorted(iso_consensus_any_set)
)

print(
    "Shared:",
    sorted(shared_any)
)

print(
    "AP5 only:",
    sorted(ap5_only_any)
)

print(
    "ISO only:",
    sorted(iso_only_any)
)

print(
    "Jaccard overlap:",
    f"{any_jaccard:.3f}"
)


print()
print("POLAR CONTACTS")
print("--------------")

print(
    "AP5 residues:",
    sorted(ap5_polar_set)
)

print(
    "ISO consensus residues:",
    sorted(iso_consensus_polar_set)
)

print(
    "Shared:",
    sorted(shared_polar)
)

print(
    "AP5 only:",
    sorted(ap5_only_polar)
)

print(
    "ISO only:",
    sorted(iso_only_polar)
)

print(
    "Jaccard overlap:",
    f"{polar_jaccard:.3f}"
)


# ============================================================
# SAVE CSV
# ============================================================

headers = list(
    rows[0].keys()
)

with open(
    CSV_OUTPUT,
    "w",
    newline="",
    encoding="utf-8"
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=headers
    )

    writer.writeheader()
    writer.writerows(rows)


# ============================================================
# SAVE EXCEL
# ============================================================

workbook = xlsxwriter.Workbook(
    XLSX_OUTPUT
)

title_fmt = workbook.add_format({
    "bold": True,
    "font_size": 15,
    "font_color": "white",
    "bg_color": "#1F4E78"
})

header_fmt = workbook.add_format({
    "bold": True,
    "font_color": "white",
    "bg_color": "#4472C4",
    "border": 1,
    "align": "center",
    "text_wrap": True
})

normal_fmt = workbook.add_format({
    "border": 1
})

distance_fmt = workbook.add_format({
    "border": 1,
    "num_format": "0.000"
})

shared_fmt = workbook.add_format({
    "border": 1,
    "bg_color": "#E2F0D9"
})

ap5_fmt = workbook.add_format({
    "border": 1,
    "bg_color": "#DDEBF7"
})

iso_fmt = workbook.add_format({
    "border": 1,
    "bg_color": "#FFF2CC"
})


# ============================================================
# SHEET 1 — RESIDUE FINGERPRINT
# ============================================================

ws = workbook.add_worksheet(
    "Residue_Fingerprint"
)

ws.merge_range(
    0,
    0,
    0,
    len(headers) - 1,
    "AK1 Native AP5 vs Isorhamnetin Interaction Fingerprint",
    title_fmt
)

for col, header in enumerate(headers):

    ws.write(
        2,
        col,
        header,
        header_fmt
    )


for row_index, row in enumerate(
    rows,
    start=3
):

    class_value = row[
        "Any-contact class"
    ]

    if class_value == "Shared":

        row_fmt = shared_fmt

    elif class_value == "AP5 only":

        row_fmt = ap5_fmt

    elif class_value == "ISO only":

        row_fmt = iso_fmt

    else:

        row_fmt = normal_fmt


    for col, header in enumerate(
        headers
    ):

        value = row[header]

        if value is None:

            ws.write(
                row_index,
                col,
                "",
                row_fmt
            )

        elif isinstance(
            value,
            float
        ) and "distance" in header.lower():

            ws.write_number(
                row_index,
                col,
                value,
                distance_fmt
            )

        else:

            ws.write(
                row_index,
                col,
                value,
                row_fmt
            )


ws.set_column(
    0,
    0,
    15
)

ws.set_column(
    1,
    len(headers) - 1,
    22
)


# ============================================================
# SHEET 2 — SUMMARY
# ============================================================

ws2 = workbook.add_worksheet(
    "Overlap_Summary"
)

summary_rows = [
    [
        "Any-contact AP5 residues",
        ", ".join(
            sorted(ap5_any_set)
        )
    ],
    [
        "Any-contact ISO consensus residues",
        ", ".join(
            sorted(
                iso_consensus_any_set
            )
        )
    ],
    [
        "Any-contact shared residues",
        ", ".join(
            sorted(shared_any)
        )
    ],
    [
        "Any-contact AP5-only residues",
        ", ".join(
            sorted(ap5_only_any)
        )
    ],
    [
        "Any-contact ISO-only residues",
        ", ".join(
            sorted(iso_only_any)
        )
    ],
    [
        "Any-contact Jaccard overlap",
        any_jaccard
    ],

    [
        "Polar AP5 residues",
        ", ".join(
            sorted(ap5_polar_set)
        )
    ],
    [
        "Polar ISO consensus residues",
        ", ".join(
            sorted(
                iso_consensus_polar_set
            )
        )
    ],
    [
        "Polar shared residues",
        ", ".join(
            sorted(shared_polar)
        )
    ],
    [
        "Polar AP5-only residues",
        ", ".join(
            sorted(ap5_only_polar)
        )
    ],
    [
        "Polar ISO-only residues",
        ", ".join(
            sorted(iso_only_polar)
        )
    ],
    [
        "Polar Jaccard overlap",
        polar_jaccard
    ]
]


ws2.write(
    0,
    0,
    "Metric",
    header_fmt
)

ws2.write(
    0,
    1,
    "Value",
    header_fmt
)


for row_index, values in enumerate(
    summary_rows,
    start=1
):

    ws2.write(
        row_index,
        0,
        values[0],
        normal_fmt
    )

    ws2.write(
        row_index,
        1,
        values[1],
        normal_fmt
    )


ws2.set_column(
    0,
    0,
    35
)

ws2.set_column(
    1,
    1,
    100
)


workbook.close()


print()
print("=" * 90)
print("FILES SAVED")
print("=" * 90)

print(CSV_OUTPUT)
print(XLSX_OUTPUT)