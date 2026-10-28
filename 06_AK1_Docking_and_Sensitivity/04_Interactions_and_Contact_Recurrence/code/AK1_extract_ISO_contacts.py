from pathlib import Path
from math import sqrt
from collections import defaultdict

from rdkit import Chem
from rdkit.Chem import ChemicalFeatures
from rdkit import RDConfig

import xlsxwriter


# ============================================================
# PATHS
# ============================================================

PROJECT = Path(r"C:\Users\akidi\Desktop\ISO_AK_Docking")

PROTEIN_FILE = PROJECT / "receptors" / "AK1_chainA_FIXED.pdb"

ISO_FILE = (
    PROJECT
    / "results"
    / "AK1_repr9009"
    / "AK1_ISO_seed9009_best.sdf"
)

RESULTS = PROJECT / "results" / "AK1_panelB"

CSV_OUTPUT = RESULTS / "AK1_ISO_seed9009_polar_contacts.csv"
XLSX_OUTPUT = RESULTS / "AK1_ISO_seed9009_contact_analysis.xlsx"

RESULTS.mkdir(exist_ok=True)


# ============================================================
# SETTINGS
# ============================================================

# Heavy-atom cutoff for screening candidate polar contacts
CUTOFF = 3.6

# More conservative distance usually associated with
# conventional hydrogen-bond-compatible heavy atom contacts
HBOND_CUTOFF = 3.5


# ============================================================
# DISTANCE FUNCTION
# ============================================================

def distance(a, b):
    return sqrt(
        (a[0] - b[0]) ** 2
        + (a[1] - b[1]) ** 2
        + (a[2] - b[2]) ** 2
    )


# ============================================================
# READ PROTEIN ATOMS FROM PDB
# ============================================================

protein_atoms = []

with open(PROTEIN_FILE, "r") as f:

    for line in f:

        if not line.startswith("ATOM"):
            continue

        atom_name = line[12:16].strip()
        resname = line[17:20].strip()
        chain = line[21].strip()
        resid = line[22:26].strip()

        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])

        element = line[76:78].strip()

        if not element:
            element = atom_name[0]

        # Only polar heavy atoms
        if element.upper() not in {"O", "N", "S"}:
            continue

        protein_atoms.append({
            "atom": atom_name,
            "resname": resname,
            "resid": resid,
            "chain": chain,
            "element": element.upper(),
            "coord": (x, y, z)
        })


# ============================================================
# PROTEIN DONOR / ACCEPTOR RULES
# ============================================================

def protein_roles(resname, atom_name):

    roles = set()

    # -------------------------
    # Backbone
    # -------------------------

    if atom_name == "N" and resname != "PRO":
        roles.add("donor")

    if atom_name in {"O", "OXT"}:
        roles.add("acceptor")


    # -------------------------
    # Side-chain donors
    # -------------------------

    donor_atoms = {
        "ARG": {"NE", "NH1", "NH2"},
        "ASN": {"ND2"},
        "GLN": {"NE2"},
        "HIS": {"ND1", "NE2"},
        "LYS": {"NZ"},
        "SER": {"OG"},
        "THR": {"OG1"},
        "TRP": {"NE1"},
        "TYR": {"OH"},
        "CYS": {"SG"},
    }

    if atom_name in donor_atoms.get(resname, set()):
        roles.add("donor")


    # -------------------------
    # Side-chain acceptors
    # -------------------------

    acceptor_atoms = {
        "ASP": {"OD1", "OD2"},
        "GLU": {"OE1", "OE2"},
        "ASN": {"OD1"},
        "GLN": {"OE1"},
        "HIS": {"ND1", "NE2"},
        "SER": {"OG"},
        "THR": {"OG1"},
        "TYR": {"OH"},
        "CYS": {"SG"},
        "MET": {"SD"},
    }

    if atom_name in acceptor_atoms.get(resname, set()):
        roles.add("acceptor")

    return roles


# ============================================================
# READ ISO USING RDKIT
# ============================================================

supplier = Chem.SDMolSupplier(
    str(ISO_FILE),
    removeHs=False
)

iso = None

for mol in supplier:
    if mol is not None:
        iso = mol
        break

if iso is None:
    raise RuntimeError(
        "Could not read the seed 9009 ISO SDF."
    )


conf = iso.GetConformer()


# ============================================================
# IDENTIFY ISO DONORS / ACCEPTORS
# ============================================================

feature_file = (
    Path(RDConfig.RDDataDir)
    / "BaseFeatures.fdef"
)

factory = ChemicalFeatures.BuildFeatureFactory(
    str(feature_file)
)

features = factory.GetFeaturesForMol(iso)


iso_roles = defaultdict(set)

for feature in features:

    family = feature.GetFamily()

    if family == "Donor":

        for atom_id in feature.GetAtomIds():
            iso_roles[atom_id].add("donor")

    elif family == "Acceptor":

        for atom_id in feature.GetAtomIds():
            iso_roles[atom_id].add("acceptor")


# ============================================================
# GIVE ISO ATOMS READABLE NAMES
# O1, O2, N1, etc.
# ============================================================

element_counts = defaultdict(int)
iso_labels = {}

for atom in iso.GetAtoms():

    symbol = atom.GetSymbol()

    element_counts[symbol] += 1

    iso_labels[atom.GetIdx()] = (
        f"{symbol}{element_counts[symbol]}"
    )


# ============================================================
# FIND POLAR CONTACTS
# ============================================================

contacts = []

for iso_atom in iso.GetAtoms():

    iso_idx = iso_atom.GetIdx()
    iso_element = iso_atom.GetSymbol().upper()

    if iso_element not in {"O", "N", "S"}:
        continue

    pos = conf.GetAtomPosition(iso_idx)

    iso_coord = (
        float(pos.x),
        float(pos.y),
        float(pos.z)
    )

    ligand_roles = iso_roles.get(
        iso_idx,
        set()
    )

    for p in protein_atoms:

        d = distance(
            iso_coord,
            p["coord"]
        )

        if d > CUTOFF:
            continue

        p_roles = protein_roles(
            p["resname"],
            p["atom"]
        )

        # ====================================================
        # CONTACT CLASSIFICATION
        # ====================================================

        donor_acceptor_match = (
            (
                "donor" in ligand_roles
                and "acceptor" in p_roles
            )
            or
            (
                "acceptor" in ligand_roles
                and "donor" in p_roles
            )
        )

        if donor_acceptor_match and d <= HBOND_CUTOFF:

            contact_type = "Candidate H-bond"

        else:

            contact_type = "Polar proximity"


        contacts.append({
            "Residue": f"{p['resname']}{p['resid']}",
            "Chain": p["chain"],
            "Protein atom": p["atom"],
            "Protein role": "/".join(sorted(p_roles)) if p_roles else "polar atom",
            "ISO atom": iso_labels[iso_idx],
            "ISO atom index": iso_idx + 1,
            "ISO role": "/".join(sorted(ligand_roles)) if ligand_roles else "polar atom",
            "Distance (A)": d,
            "Contact type": contact_type
        })


# Sort shortest first
contacts.sort(
    key=lambda x: x["Distance (A)"]
)


# ============================================================
# PRINT TO CMD
# ============================================================

print()
print("=" * 100)
print("AK1 - ISORHAMNETIN SEED 9009 POLAR CONTACT ANALYSIS")
print("=" * 100)

print()
print(
    f"{'Residue':<10}"
    f"{'Prot atom':<12}"
    f"{'ISO atom':<10}"
    f"{'Distance':<12}"
    f"{'Contact type'}"
)

print("-" * 100)

for c in contacts:

    print(
        f"{c['Residue']:<10}"
        f"{c['Protein atom']:<12}"
        f"{c['ISO atom']:<10}"
        f"{c['Distance (A)']:<12.3f}"
        f"{c['Contact type']}"
    )


print()
print("Total polar contacts:", len(contacts))

candidate_hbonds = [
    c for c in contacts
    if c["Contact type"] == "Candidate H-bond"
]

print(
    "Candidate H-bonds:",
    len(candidate_hbonds)
)


# ============================================================
# RESIDUE-LEVEL SUMMARY
# ============================================================

residue_summary = {}

for c in contacts:

    residue = c["Residue"]

    if residue not in residue_summary:

        residue_summary[residue] = {
            "Residue": residue,
            "Minimum distance": c["Distance (A)"],
            "Number of polar contacts": 0,
            "Candidate H-bonds": 0
        }

    residue_summary[residue][
        "Number of polar contacts"
    ] += 1

    residue_summary[residue][
        "Minimum distance"
    ] = min(
        residue_summary[residue][
            "Minimum distance"
        ],
        c["Distance (A)"]
    )

    if c["Contact type"] == "Candidate H-bond":

        residue_summary[residue][
            "Candidate H-bonds"
        ] += 1


residue_summary = sorted(
    residue_summary.values(),
    key=lambda x: x["Minimum distance"]
)


# ============================================================
# SAVE CSV
# ============================================================

import csv

columns = [
    "Residue",
    "Chain",
    "Protein atom",
    "Protein role",
    "ISO atom",
    "ISO atom index",
    "ISO role",
    "Distance (A)",
    "Contact type"
]

with open(
    CSV_OUTPUT,
    "w",
    newline="",
    encoding="utf-8"
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=columns
    )

    writer.writeheader()

    for c in contacts:

        row = c.copy()
        row["Distance (A)"] = round(
            row["Distance (A)"],
            3
        )

        writer.writerow(row)


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
    "align": "center"
})

normal_fmt = workbook.add_format({
    "border": 1
})

distance_fmt = workbook.add_format({
    "border": 1,
    "num_format": "0.000",
    "align": "center"
})

hbond_fmt = workbook.add_format({
    "border": 1,
    "bg_color": "#E2F0D9"
})

polar_fmt = workbook.add_format({
    "border": 1,
    "bg_color": "#FFF2CC"
})


# ============================================================
# SHEET 1 — ALL CONTACTS
# ============================================================

ws = workbook.add_worksheet(
    "Atom_Contacts"
)

ws.merge_range(
    "A1:I1",
    "AK1-Isorhamnetin Seed 9009 Polar Contact Analysis",
    title_fmt
)

for col, header in enumerate(columns):

    ws.write(
        2,
        col,
        header,
        header_fmt
    )


for row, c in enumerate(
    contacts,
    start=3
):

    fmt = (
        hbond_fmt
        if c["Contact type"]
        == "Candidate H-bond"
        else polar_fmt
    )

    values = [
        c["Residue"],
        c["Chain"],
        c["Protein atom"],
        c["Protein role"],
        c["ISO atom"],
        c["ISO atom index"],
        c["ISO role"],
        c["Distance (A)"],
        c["Contact type"]
    ]

    for col, value in enumerate(values):

        if col == 7:

            ws.write_number(
                row,
                col,
                value,
                distance_fmt
            )

        else:

            ws.write(
                row,
                col,
                value,
                fmt
            )


ws.set_column("A:A", 13)
ws.set_column("B:B", 9)
ws.set_column("C:C", 15)
ws.set_column("D:D", 20)
ws.set_column("E:E", 12)
ws.set_column("F:F", 15)
ws.set_column("G:G", 20)
ws.set_column("H:H", 16)
ws.set_column("I:I", 24)

ws.freeze_panes(3, 0)


# ============================================================
# SHEET 2 — RESIDUE SUMMARY
# ============================================================

ws2 = workbook.add_worksheet(
    "Residue_Summary"
)

ws2.merge_range(
    "A1:D1",
    "AK1 Residues Contacting Isorhamnetin",
    title_fmt
)

summary_headers = [
    "Residue",
    "Minimum distance (A)",
    "Number of polar contacts",
    "Candidate H-bonds"
]

for col, header in enumerate(
    summary_headers
):

    ws2.write(
        2,
        col,
        header,
        header_fmt
    )


for row, r in enumerate(
    residue_summary,
    start=3
):

    ws2.write(
        row,
        0,
        r["Residue"],
        normal_fmt
    )

    ws2.write_number(
        row,
        1,
        r["Minimum distance"],
        distance_fmt
    )

    ws2.write_number(
        row,
        2,
        r["Number of polar contacts"],
        normal_fmt
    )

    ws2.write_number(
        row,
        3,
        r["Candidate H-bonds"],
        normal_fmt
    )


ws2.set_column("A:A", 15)
ws2.set_column("B:D", 24)


# ============================================================
# SHEET 3 — METHOD NOTES
# ============================================================

ws3 = workbook.add_worksheet(
    "Interpretation"
)

ws3.write(
    "A1",
    "Important interpretation",
    title_fmt
)

notes = [
    [
        "Representative pose",
        "AK1-ISO seed 9009, selected from the dominant 8/10 pose cluster."
    ],
    [
        "Screening cutoff",
        "Protein and ligand polar heavy atoms within 3.6 Å."
    ],
    [
        "Candidate H-bond definition",
        "Donor/acceptor-compatible heavy-atom pair at <=3.5 Å."
    ],
    [
        "Caution",
        "Candidate H-bond classification is based on atom type and heavy-atom distance. "
        "Full hydrogen-bond assignment should additionally consider hydrogen geometry and donor-H-acceptor angle."
    ],
    [
        "Other interactions",
        "This table does not yet quantify hydrophobic, aromatic/pi, cation-pi or other non-polar contacts."
    ]
]

for row, values in enumerate(
    notes,
    start=2
):

    ws3.write(
        row,
        0,
        values[0],
        header_fmt
    )

    ws3.write(
        row,
        1,
        values[1],
        normal_fmt
    )


ws3.set_column("A:A", 28)
ws3.set_column("B:B", 100)


workbook.close()


print()
print("=" * 100)
print("FILES SAVED")
print("=" * 100)

print()
print(CSV_OUTPUT)
print(XLSX_OUTPUT)

print()