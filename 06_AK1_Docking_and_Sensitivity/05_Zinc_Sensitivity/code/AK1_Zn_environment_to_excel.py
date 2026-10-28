from pathlib import Path
import math
import csv
import xlsxwriter


# ============================================================
# PATHS
# ============================================================

PROJECT = Path(r"C:\Users\akidi\Desktop\ISO_AK_Docking")

PDB = (
    PROJECT
    / "structures"
    / "AK1_1Z83.pdb"
)

OUTPUT_DIR = (
    PROJECT
    / "results"
    / "AK1_Zn_sensitivity"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

XLSX_OUTPUT = (
    OUTPUT_DIR
    / "AK1_1Z83_Zn_environment_analysis.xlsx"
)

CSV_OUTPUT = (
    OUTPUT_DIR
    / "AK1_1Z83_Zn_neighbors_within5A.csv"
)


# ============================================================
# DISTANCE
# ============================================================

def dist(a, b):

    return math.sqrt(
        (a[0] - b[0]) ** 2
        + (a[1] - b[1]) ** 2
        + (a[2] - b[2]) ** 2
    )


# ============================================================
# READ PDB
# ============================================================

atoms = []

with open(PDB, "r") as f:

    for line in f:

        if not (
            line.startswith("ATOM")
            or line.startswith("HETATM")
        ):
            continue

        record = line[0:6].strip()
        atom = line[12:16].strip()
        resn = line[17:20].strip()
        chain = line[21].strip()
        resi = line[22:26].strip()

        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])

        element = line[76:78].strip()

        if not element:
            element = atom[0]

        atoms.append({
            "record": record,
            "atom": atom,
            "resn": resn,
            "chain": chain,
            "resi": resi,
            "element": element.upper(),
            "coord": (x, y, z)
        })


# ============================================================
# FIND CHAIN-A Zn
# ============================================================

zn_atoms = [
    a for a in atoms
    if a["resn"] == "ZN"
    and a["chain"] == "A"
]

if not zn_atoms:
    raise RuntimeError(
        "No chain-A Zn atoms found."
    )


# ============================================================
# BUILD 5-A NEIGHBOR TABLE
# ============================================================

neighbors = []

for zn in zn_atoms:

    zn_id = f"ZN{zn['resi']}"

    for a in atoms:

        if a is zn:
            continue

        d = dist(
            zn["coord"],
            a["coord"]
        )

        if d > 5.0:
            continue

        if a["record"] == "ATOM":

            category = "PROTEIN"

        elif a["resn"] == "AP5":

            category = "AP5"

        elif a["resn"] in {
            "HOH",
            "WAT"
        }:

            category = "WATER"

        else:

            category = "OTHER"


        neighbors.append({
            "Zn": zn_id,
            "Zn residue": zn["resi"],
            "Zn X": zn["coord"][0],
            "Zn Y": zn["coord"][1],
            "Zn Z": zn["coord"][2],
            "Distance (A)": d,
            "Category": category,
            "Residue": f"{a['resn']}{a['resi']}",
            "Residue name": a["resn"],
            "Residue number": a["resi"],
            "Chain": a["chain"],
            "Atom": a["atom"],
            "Element": a["element"]
        })


neighbors.sort(
    key=lambda x: (
        int(x["Zn residue"]),
        x["Distance (A)"]
    )
)


# ============================================================
# SUMMARIES FOR EACH Zn
# ============================================================

summary = []

for zn in zn_atoms:

    zn_id = f"ZN{zn['resi']}"

    rows = [
        r for r in neighbors
        if r["Zn"] == zn_id
    ]

    ap5_rows = [
        r for r in rows
        if r["Category"] == "AP5"
    ]

    protein_rows = [
        r for r in rows
        if r["Category"] == "PROTEIN"
    ]

    water_rows = [
        r for r in rows
        if r["Category"] == "WATER"
    ]


    nearest_ap5 = (
        min(
            ap5_rows,
            key=lambda r:
                r["Distance (A)"]
        )
        if ap5_rows
        else None
    )

    nearest_protein = (
        min(
            protein_rows,
            key=lambda r:
                r["Distance (A)"]
        )
        if protein_rows
        else None
    )

    nearest_water = (
        min(
            water_rows,
            key=lambda r:
                r["Distance (A)"]
        )
        if water_rows
        else None
    )


    if zn["resi"] == "601":

        interpretation = (
            "Directly associated with the AP5 phosphate region: "
            "Zn is within ~2.05–2.07 A of AP5 phosphate oxygens "
            "and is additionally coordinated by crystallographic waters."
        )

        relevance = (
            "High relevance for AP5-site docking sensitivity."
        )

    elif zn["resi"] == "641":

        interpretation = (
            "Primarily coordinated by AK1 GLU26/HIS36 and waters; "
            "not directly associated with AP5 within the immediate "
            "coordination shell."
        )

        relevance = (
            "Likely lower direct relevance to AP5-site ligand docking."
        )

    else:

        interpretation = ""
        relevance = ""


    summary.append({
        "Zn": zn_id,
        "X": zn["coord"][0],
        "Y": zn["coord"][1],
        "Z": zn["coord"][2],

        "Nearest AP5 atom": (
            nearest_ap5["Atom"]
            if nearest_ap5 else ""
        ),

        "Nearest AP5 distance (A)": (
            nearest_ap5["Distance (A)"]
            if nearest_ap5 else None
        ),

        "Nearest protein residue": (
            nearest_protein["Residue"]
            if nearest_protein else ""
        ),

        "Nearest protein atom": (
            nearest_protein["Atom"]
            if nearest_protein else ""
        ),

        "Nearest protein distance (A)": (
            nearest_protein["Distance (A)"]
            if nearest_protein else None
        ),

        "Nearest water": (
            nearest_water["Residue"]
            if nearest_water else ""
        ),

        "Nearest water distance (A)": (
            nearest_water["Distance (A)"]
            if nearest_water else None
        ),

        "AP5 atoms within 5 A": len(
            ap5_rows
        ),

        "Protein atoms within 5 A": len(
            protein_rows
        ),

        "Waters within 5 A": len(
            water_rows
        ),

        "Interpretation": interpretation,

        "Docking relevance": relevance
    })


# ============================================================
# SAVE CSV — ALL NEIGHBORS
# ============================================================

neighbor_headers = list(
    neighbors[0].keys()
)

with open(
    CSV_OUTPUT,
    "w",
    newline="",
    encoding="utf-8"
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=neighbor_headers
    )

    writer.writeheader()
    writer.writerows(neighbors)


# ============================================================
# EXCEL FORMATTING
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

dist_fmt = workbook.add_format({
    "border": 1,
    "num_format": "0.000"
})

important_fmt = workbook.add_format({
    "border": 1,
    "bg_color": "#FFF2CC"
})

ap5_fmt = workbook.add_format({
    "border": 1,
    "bg_color": "#DDEBF7"
})

water_fmt = workbook.add_format({
    "border": 1,
    "bg_color": "#E2F0D9"
})


# ============================================================
# SHEET 1 — Zn SUMMARY
# ============================================================

ws = workbook.add_worksheet(
    "Zn_Summary"
)

summary_headers = list(
    summary[0].keys()
)

ws.merge_range(
    0,
    0,
    0,
    len(summary_headers) - 1,
    "AK1 PDB 1Z83 Chain-A Zinc Environment Analysis",
    title_fmt
)

for col, header in enumerate(
    summary_headers
):

    ws.write(
        2,
        col,
        header,
        header_fmt
    )


for row_i, row in enumerate(
    summary,
    start=3
):

    for col, header in enumerate(
        summary_headers
    ):

        value = row[header]

        fmt = (
            important_fmt
            if row["Zn"] == "ZN601"
            else normal_fmt
        )

        if value is None:

            ws.write(
                row_i,
                col,
                "",
                fmt
            )

        elif isinstance(
            value,
            float
        ):

            ws.write_number(
                row_i,
                col,
                value,
                dist_fmt
            )

        else:

            ws.write(
                row_i,
                col,
                value,
                fmt
            )


ws.set_column(
    0,
    3,
    14
)

ws.set_column(
    4,
    13,
    24
)

ws.set_column(
    14,
    15,
    65
)


# ============================================================
# SHEET 2 — ALL NEIGHBORS WITHIN 5 A
# ============================================================

ws2 = workbook.add_worksheet(
    "Neighbors_within_5A"
)

ws2.merge_range(
    0,
    0,
    0,
    len(neighbor_headers) - 1,
    "Atoms Within 5.0 A of Chain-A Zn Ions",
    title_fmt
)


for col, header in enumerate(
    neighbor_headers
):

    ws2.write(
        2,
        col,
        header,
        header_fmt
    )


for row_i, row in enumerate(
    neighbors,
    start=3
):

    category = row["Category"]

    if category == "AP5":

        row_fmt = ap5_fmt

    elif category == "WATER":

        row_fmt = water_fmt

    elif row["Zn"] == "ZN601":

        row_fmt = important_fmt

    else:

        row_fmt = normal_fmt


    for col, header in enumerate(
        neighbor_headers
    ):

        value = row[header]

        if isinstance(
            value,
            float
        ):

            ws2.write_number(
                row_i,
                col,
                value,
                dist_fmt
            )

        else:

            ws2.write(
                row_i,
                col,
                value,
                row_fmt
            )


ws2.set_column(
    0,
    len(neighbor_headers) - 1,
    16
)

ws2.freeze_panes(
    3,
    0
)


# ============================================================
# SHEET 3 — ZN601 AP5-SITE DETAILS
# ============================================================

ws3 = workbook.add_worksheet(
    "ZN601_AP5_Site"
)

zn601_rows = [
    r for r in neighbors
    if r["Zn"] == "ZN601"
]

ws3.merge_range(
    0,
    0,
    0,
    len(neighbor_headers) - 1,
    "ZN601: AP5-Site Zinc Environment",
    title_fmt
)

for col, header in enumerate(
    neighbor_headers
):

    ws3.write(
        2,
        col,
        header,
        header_fmt
    )


for row_i, row in enumerate(
    zn601_rows,
    start=3
):

    if row["Category"] == "AP5":

        row_fmt = ap5_fmt

    elif row["Category"] == "WATER":

        row_fmt = water_fmt

    else:

        row_fmt = normal_fmt


    for col, header in enumerate(
        neighbor_headers
    ):

        value = row[header]

        if isinstance(
            value,
            float
        ):

            ws3.write_number(
                row_i,
                col,
                value,
                dist_fmt
            )

        else:

            ws3.write(
                row_i,
                col,
                value,
                row_fmt
            )


ws3.set_column(
    0,
    len(neighbor_headers) - 1,
    16
)


# ============================================================
# SHEET 4 — METHOD / INTERPRETATION
# ============================================================

ws4 = workbook.add_worksheet(
    "Method_Notes"
)

notes = [
    [
        "Structure",
        "Human AK1 PDB 1Z83, chain A."
    ],

    [
        "Neighbor cutoff",
        "All protein, ligand and water atoms within 5.0 A of chain-A Zn ions were recorded."
    ],

    [
        "ZN601",
        "Located directly adjacent to AP5 phosphate oxygens O2D and O1G at approximately 2.05–2.07 A, with several coordinating crystallographic waters."
    ],

    [
        "ZN641",
        "Nearest direct coordination involves GLU26 OE2 (~1.97 A), HIS36 ND1 (~2.32 A), and waters; AP5 is not in its immediate coordination shell."
    ],

    [
        "Sensitivity interpretation",
        "ZN601 is the priority zinc ion for testing whether inclusion of the crystallographic metal changes the predicted AK1-isorhamnetin docking pose or score."
    ],

    [
        "Important caution",
        "Retention of a crystallographic metal in docking requires careful receptor preparation and appropriate treatment of metal charge/coordination. A metal-containing Vina run should therefore be interpreted as a sensitivity analysis rather than an experimentally measured metal-dependent affinity."
    ]
]


ws4.write(
    0,
    0,
    "Item",
    header_fmt
)

ws4.write(
    0,
    1,
    "Description",
    header_fmt
)


for row_i, values in enumerate(
    notes,
    start=1
):

    ws4.write(
        row_i,
        0,
        values[0],
        normal_fmt
    )

    ws4.write(
        row_i,
        1,
        values[1],
        normal_fmt
    )


ws4.set_column(
    0,
    0,
    28
)

ws4.set_column(
    1,
    1,
    110
)


workbook.close()


# ============================================================
# DONE
# ============================================================

print()
print("=" * 90)
print("Zn ANALYSIS FILES SAVED")
print("=" * 90)

print()
print(CSV_OUTPUT)
print(XLSX_OUTPUT)