from pathlib import Path
from copy import deepcopy
import csv
import xlsxwriter

from rdkit import Chem


# ============================================================
# PATHS
# ============================================================

PROJECT = Path(
    r"C:\Users\akidi\Desktop\ISO_AK_Docking"
)

INPUT_SDF = (
    PROJECT
    / "ligand"
    / "isorhamnetin_3D.sdf"
)

OUTPUT_DIR = (
    PROJECT
    / "ligand"
    / "ISO_protonation_states"
)

RESULT_DIR = (
    PROJECT
    / "results"
    / "AK1_ligand_state_sensitivity"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

RESULT_DIR.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# OUTPUT FILES
# ============================================================

NEUTRAL_SDF = (
    OUTPUT_DIR
    / "isorhamnetin_state0_neutral.sdf"
)

MONO_SDF = (
    OUTPUT_DIR
    / "isorhamnetin_state1_7O_minus.sdf"
)

DI_SDF = (
    OUTPUT_DIR
    / "isorhamnetin_state2_7O_4primeO_dianion.sdf"
)

CSV_OUTPUT = (
    RESULT_DIR
    / "isorhamnetin_protonation_state_mapping.csv"
)

XLSX_OUTPUT = (
    RESULT_DIR
    / "isorhamnetin_protonation_state_mapping.xlsx"
)


# ============================================================
# READ ORIGINAL MOLECULE
# ============================================================

supplier = Chem.SDMolSupplier(
    str(INPUT_SDF),
    removeHs=False
)

mol = None

for candidate in supplier:

    if candidate is not None:
        mol = candidate
        break

if mol is None:

    raise RuntimeError(
        "Could not read the original isorhamnetin SDF."
    )


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def atom_label(atom):

    return (
        f"{atom.GetSymbol()}"
        f"{atom.GetIdx() + 1}"
    )


def total_formal_charge(molecule):

    return sum(
        atom.GetFormalCharge()
        for atom in molecule.GetAtoms()
    )


def find_explicit_h_neighbor(
    molecule,
    oxygen_idx
):

    oxygen = molecule.GetAtomWithIdx(
        oxygen_idx
    )

    h_neighbors = [
        n.GetIdx()
        for n in oxygen.GetNeighbors()
        if n.GetAtomicNum() == 1
    ]

    if len(h_neighbors) != 1:

        raise RuntimeError(
            f"Expected exactly one explicit H attached "
            f"to oxygen {oxygen_idx + 1}; "
            f"found {len(h_neighbors)}."
        )

    return h_neighbors[0]


# ============================================================
# IDENTIFY ALL EXPLICIT PHENOLIC OH OXYGENS
# ============================================================

oh_oxygen_indices = []

for atom in mol.GetAtoms():

    if atom.GetAtomicNum() != 8:
        continue

    h_neighbors = [
        n
        for n in atom.GetNeighbors()
        if n.GetAtomicNum() == 1
    ]

    carbon_neighbors = [
        n
        for n in atom.GetNeighbors()
        if n.GetAtomicNum() == 6
    ]

    if (
        len(h_neighbors) == 1
        and
        len(carbon_neighbors) == 1
    ):

        oh_oxygen_indices.append(
            atom.GetIdx()
        )


if len(oh_oxygen_indices) != 4:

    raise RuntimeError(
        f"Expected 4 OH oxygens, "
        f"found {len(oh_oxygen_indices)}."
    )


# ============================================================
# IDENTIFY METHOXY GROUP
#
# Methoxy oxygen:
# oxygen bound to two carbons,
# one of which is a methyl carbon.
# ============================================================

methoxy_ring_carbon_idx = None

for atom in mol.GetAtoms():

    if atom.GetAtomicNum() != 8:
        continue

    carbon_neighbors = [
        n
        for n in atom.GetNeighbors()
        if n.GetAtomicNum() == 6
    ]

    if len(carbon_neighbors) != 2:
        continue

    methyl_candidates = []

    aromatic_candidates = []

    for carbon in carbon_neighbors:

        if carbon.GetIsAromatic():
            aromatic_candidates.append(
                carbon
            )

        heavy_neighbors = [
            n
            for n in carbon.GetNeighbors()
            if n.GetAtomicNum() != 1
        ]

        if (
            not carbon.GetIsAromatic()
            and
            len(heavy_neighbors) == 1
        ):
            methyl_candidates.append(
                carbon
            )

    if (
        len(methyl_candidates) == 1
        and
        len(aromatic_candidates) == 1
    ):

        methoxy_ring_carbon_idx = (
            aromatic_candidates[0]
            .GetIdx()
        )

        break


if methoxy_ring_carbon_idx is None:

    raise RuntimeError(
        "Could not identify the methoxy-substituted "
        "aromatic carbon."
    )


# ============================================================
# IDENTIFY AROMATIC SIX-MEMBERED CARBON RINGS
# ============================================================

ring_info = mol.GetRingInfo()

six_carbon_aromatic_rings = []

for ring in ring_info.AtomRings():

    if len(ring) != 6:
        continue

    ring_atoms = [
        mol.GetAtomWithIdx(i)
        for i in ring
    ]

    if all(
        a.GetAtomicNum() == 6
        and a.GetIsAromatic()
        for a in ring_atoms
    ):

        six_carbon_aromatic_rings.append(
            set(ring)
        )


if len(six_carbon_aromatic_rings) < 2:

    raise RuntimeError(
        "Could not identify both aromatic carbon rings."
    )


# ============================================================
# IDENTIFY B RING
#
# B ring contains methoxy-substituted carbon.
# ============================================================

b_ring = None

for ring in six_carbon_aromatic_rings:

    if methoxy_ring_carbon_idx in ring:

        b_ring = ring
        break


if b_ring is None:

    raise RuntimeError(
        "Could not identify the B ring."
    )


# ============================================================
# IDENTIFY 4'-OH
#
# The B ring contains:
# - the methoxy substituent
# - one phenolic OH
#
# That OH is 4'-OH.
# ============================================================

four_prime_oxygen_idx = None
four_prime_carbon_idx = None

for oxygen_idx in oh_oxygen_indices:

    oxygen = mol.GetAtomWithIdx(
        oxygen_idx
    )

    carbon = next(
        n
        for n in oxygen.GetNeighbors()
        if n.GetAtomicNum() == 6
    )

    if carbon.GetIdx() in b_ring:

        four_prime_oxygen_idx = (
            oxygen_idx
        )

        four_prime_carbon_idx = (
            carbon.GetIdx()
        )

        break


if four_prime_oxygen_idx is None:

    raise RuntimeError(
        "Could not identify 4'-OH."
    )


# ============================================================
# IDENTIFY A RING
#
# The other all-carbon aromatic ring
# containing two OH substituents.
# ============================================================

a_ring_candidates = []

for ring in six_carbon_aromatic_rings:

    if ring == b_ring:
        continue

    attached_oh = []

    for oxygen_idx in oh_oxygen_indices:

        oxygen = mol.GetAtomWithIdx(
            oxygen_idx
        )

        carbon = next(
            n
            for n in oxygen.GetNeighbors()
            if n.GetAtomicNum() == 6
        )

        if carbon.GetIdx() in ring:

            attached_oh.append(
                (
                    oxygen_idx,
                    carbon.GetIdx()
                )
            )

    if len(attached_oh) == 2:

        a_ring_candidates.append(
            (
                ring,
                attached_oh
            )
        )


if len(a_ring_candidates) != 1:

    raise RuntimeError(
        f"Expected one A-ring candidate, "
        f"found {len(a_ring_candidates)}."
    )


a_ring, a_ring_oh = (
    a_ring_candidates[0]
)


# ============================================================
# IDENTIFY CARBONYL CARBON
# ============================================================

carbonyl_carbon_idx = None
carbonyl_oxygen_idx = None

for bond in mol.GetBonds():

    if bond.GetBondType() != Chem.BondType.DOUBLE:
        continue

    a1 = bond.GetBeginAtom()
    a2 = bond.GetEndAtom()

    if (
        a1.GetAtomicNum() == 6
        and
        a2.GetAtomicNum() == 8
    ):

        carbonyl_carbon_idx = (
            a1.GetIdx()
        )

        carbonyl_oxygen_idx = (
            a2.GetIdx()
        )

        break

    if (
        a2.GetAtomicNum() == 6
        and
        a1.GetAtomicNum() == 8
    ):

        carbonyl_carbon_idx = (
            a2.GetIdx()
        )

        carbonyl_oxygen_idx = (
            a1.GetIdx()
        )

        break


if carbonyl_carbon_idx is None:

    raise RuntimeError(
        "Could not identify the flavone carbonyl."
    )


# ============================================================
# DISTINGUISH 5-OH FROM 7-OH
#
# 5-OH lies closer in the molecular graph
# to the carbonyl carbon.
#
# 7-OH is the other A-ring hydroxyl.
# ============================================================

distance_matrix = Chem.GetDistanceMatrix(
    mol
)

a_ring_oh_distances = []

for oxygen_idx, carbon_idx in a_ring_oh:

    distance_to_carbonyl = int(
        distance_matrix[
            carbon_idx,
            carbonyl_carbon_idx
        ]
    )

    a_ring_oh_distances.append(
        (
            distance_to_carbonyl,
            oxygen_idx,
            carbon_idx
        )
    )


a_ring_oh_distances.sort(
    key=lambda x: x[0]
)


# closer OH = 5-OH
five_distance, five_oxygen_idx, five_carbon_idx = (
    a_ring_oh_distances[0]
)

# farther OH = 7-OH
seven_distance, seven_oxygen_idx, seven_carbon_idx = (
    a_ring_oh_distances[1]
)


if seven_distance <= five_distance:

    raise RuntimeError(
        "Unexpected A-ring OH topology."
    )


# ============================================================
# IDENTIFY REMAINING OH = 3-OH
# ============================================================

assigned_oh = {
    four_prime_oxygen_idx,
    five_oxygen_idx,
    seven_oxygen_idx
}

remaining_oh = [
    idx
    for idx in oh_oxygen_indices
    if idx not in assigned_oh
]


if len(remaining_oh) != 1:

    raise RuntimeError(
        "Could not uniquely identify 3-OH."
    )


three_oxygen_idx = (
    remaining_oh[0]
)

three_oxygen = mol.GetAtomWithIdx(
    three_oxygen_idx
)

three_carbon_idx = next(
    n.GetIdx()
    for n in three_oxygen.GetNeighbors()
    if n.GetAtomicNum() == 6
)


# ============================================================
# REPORT IDENTIFIED OH POSITIONS
# ============================================================

assignments = [
    (
        "3-OH",
        three_oxygen_idx,
        three_carbon_idx
    ),

    (
        "5-OH",
        five_oxygen_idx,
        five_carbon_idx
    ),

    (
        "7-OH",
        seven_oxygen_idx,
        seven_carbon_idx
    ),

    (
        "4'-OH",
        four_prime_oxygen_idx,
        four_prime_carbon_idx
    )
]


print()
print("=" * 80)
print("ISORHAMNETIN OH POSITION ASSIGNMENT")
print("=" * 80)


for name, oxygen_idx, carbon_idx in assignments:

    h_idx = find_explicit_h_neighbor(
        mol,
        oxygen_idx
    )

    print(
        f"{name:<7}  "
        f"O atom = {oxygen_idx + 1:<3}  "
        f"attached C = {carbon_idx + 1:<3}  "
        f"H atom = {h_idx + 1:<3}"
    )


print()
print(
    "Original formal charge:",
    total_formal_charge(mol)
)


# ============================================================
# FUNCTION TO DEPROTONATE SPECIFIC OXYGENS
# ============================================================

def make_deprotonated_state(
    original_mol,
    oxygen_indices
):

    rw = Chem.RWMol(
        original_mol
    )

    # Find explicit H indices BEFORE deleting atoms
    removals = []

    for oxygen_idx in oxygen_indices:

        h_idx = find_explicit_h_neighbor(
            original_mol,
            oxygen_idx
        )

        removals.append(
            (
                oxygen_idx,
                h_idx
            )
        )

    # Set oxygen formal charges first
    for oxygen_idx, h_idx in removals:

        oxygen = rw.GetAtomWithIdx(
            oxygen_idx
        )

        oxygen.SetFormalCharge(
            -1
        )

        oxygen.SetNoImplicit(
            True
        )

    # Remove hydrogens in descending atom index order
    # so earlier indices do not shift unexpectedly.
    h_indices = sorted(
        [
            h_idx
            for oxygen_idx, h_idx
            in removals
        ],
        reverse=True
    )

    for h_idx in h_indices:

        rw.RemoveAtom(
            h_idx
        )

    result = rw.GetMol()

    Chem.SanitizeMol(
        result
    )

    return result


# ============================================================
# CREATE STATES
# ============================================================

neutral = Chem.Mol(
    mol
)

mono_7O = make_deprotonated_state(
    mol,
    [
        seven_oxygen_idx
    ]
)

dianion_7O_4primeO = make_deprotonated_state(
    mol,
    [
        seven_oxygen_idx,
        four_prime_oxygen_idx
    ]
)


# ============================================================
# SET NAMES / METADATA
# ============================================================

neutral.SetProp(
    "_Name",
    "Isorhamnetin neutral"
)

neutral.SetProp(
    "State",
    "Neutral"
)

neutral.SetProp(
    "FormalCharge",
    str(
        total_formal_charge(
            neutral
        )
    )
)


mono_7O.SetProp(
    "_Name",
    "Isorhamnetin 7-O monoanion"
)

mono_7O.SetProp(
    "State",
    "7-O- monoanion"
)

mono_7O.SetProp(
    "FormalCharge",
    str(
        total_formal_charge(
            mono_7O
        )
    )
)


dianion_7O_4primeO.SetProp(
    "_Name",
    "Isorhamnetin 7-O / 4'-O dianion"
)

dianion_7O_4primeO.SetProp(
    "State",
    "7-O- / 4'-O- dianion"
)

dianion_7O_4primeO.SetProp(
    "FormalCharge",
    str(
        total_formal_charge(
            dianion_7O_4primeO
        )
    )
)


# ============================================================
# WRITE SDF FILES
# ============================================================

def write_sdf(
    molecule,
    path
):

    writer = Chem.SDWriter(
        str(path)
    )

    writer.write(
        molecule
    )

    writer.close()


write_sdf(
    neutral,
    NEUTRAL_SDF
)

write_sdf(
    mono_7O,
    MONO_SDF
)

write_sdf(
    dianion_7O_4primeO,
    DI_SDF
)


# ============================================================
# BUILD STATE SUMMARY
# ============================================================

state_rows = [
    {
        "State":
            "Neutral",

        "Deprotonated_sites":
            "None",

        "Formal_charge":
            total_formal_charge(
                neutral
            ),

        "Output_SDF":
            str(
                NEUTRAL_SDF
            )
    },

    {
        "State":
            "7-O- monoanion",

        "Deprotonated_sites":
            "7-OH",

        "Formal_charge":
            total_formal_charge(
                mono_7O
            ),

        "Output_SDF":
            str(
                MONO_SDF
            )
    },

    {
        "State":
            "7-O-/4'-O- dianion",

        "Deprotonated_sites":
            "7-OH; 4'-OH",

        "Formal_charge":
            total_formal_charge(
                dianion_7O_4primeO
            ),

        "Output_SDF":
            str(
                DI_SDF
            )
    }
]


mapping_rows = []

for name, oxygen_idx, carbon_idx in assignments:

    h_idx = find_explicit_h_neighbor(
        mol,
        oxygen_idx
    )

    mapping_rows.append(
        {
            "Functional_group":
                name,

            "O_atom_index_1based":
                oxygen_idx + 1,

            "Attached_C_atom_index_1based":
                carbon_idx + 1,

            "H_atom_index_1based":
                h_idx + 1
        }
    )


# ============================================================
# SAVE CSV
# ============================================================

with open(
    CSV_OUTPUT,
    "w",
    newline="",
    encoding="utf-8"
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=[
            "Functional_group",
            "O_atom_index_1based",
            "Attached_C_atom_index_1based",
            "H_atom_index_1based"
        ]
    )

    writer.writeheader()

    writer.writerows(
        mapping_rows
    )


# ============================================================
# SAVE EXCEL
# ============================================================

workbook = xlsxwriter.Workbook(
    XLSX_OUTPUT
)

header_fmt = workbook.add_format({
    "bold": True,
    "font_color": "white",
    "bg_color": "#4472C4",
    "border": 1
})

normal_fmt = workbook.add_format({
    "border": 1
})


# ------------------------------------------------------------
# OH Mapping
# ------------------------------------------------------------

ws = workbook.add_worksheet(
    "OH_Position_Mapping"
)

mapping_headers = list(
    mapping_rows[0].keys()
)

for col, header in enumerate(
    mapping_headers
):

    ws.write(
        0,
        col,
        header,
        header_fmt
    )


for row_i, row in enumerate(
    mapping_rows,
    start=1
):

    for col, header in enumerate(
        mapping_headers
    ):

        ws.write(
            row_i,
            col,
            row[header],
            normal_fmt
        )


ws.set_column(
    0,
    len(mapping_headers) - 1,
    28
)


# ------------------------------------------------------------
# State summary
# ------------------------------------------------------------

ws2 = workbook.add_worksheet(
    "Protonation_States"
)

state_headers = list(
    state_rows[0].keys()
)

for col, header in enumerate(
    state_headers
):

    ws2.write(
        0,
        col,
        header,
        header_fmt
    )


for row_i, row in enumerate(
    state_rows,
    start=1
):

    for col, header in enumerate(
        state_headers
    ):

        ws2.write(
            row_i,
            col,
            row[header],
            normal_fmt
        )


ws2.set_column(
    "A:C",
    25
)

ws2.set_column(
    "D:D",
    95
)


# ------------------------------------------------------------
# Method notes
# ------------------------------------------------------------

ws3 = workbook.add_worksheet(
    "Method_Notes"
)

notes = [
    [
        "Baseline ligand",
        "Original PubChem-derived 3D SDF used in prior AK1 docking."
    ],
    [
        "Baseline formal charge",
        "0 (neutral isorhamnetin)."
    ],
    [
        "Monoanion",
        "7-OH proton removed; oxygen formal charge set to -1."
    ],
    [
        "Dianion",
        "7-OH and 4'-OH protons removed; both oxygens assigned formal charge -1."
    ],
    [
        "Heavy-atom coordinates",
        "Inherited from the original ligand so protonation-state sensitivity is not confounded by an independently generated starting conformer."
    ],
    [
        "Purpose",
        "Assess whether the AK1 docking pose and recurrent interaction network are robust to plausible isorhamnetin protonation states."
    ]
]


for row_i, values in enumerate(
    notes
):

    ws3.write(
        row_i,
        0,
        values[0],
        header_fmt
    )

    ws3.write(
        row_i,
        1,
        values[1],
        normal_fmt
    )


ws3.set_column(
    "A:A",
    30
)

ws3.set_column(
    "B:B",
    110
)


workbook.close()


# ============================================================
# FINAL VERIFICATION
# ============================================================

print()
print("=" * 80)
print("GENERATED PROTONATION STATES")
print("=" * 80)

print(
    f"Neutral charge:  "
    f"{total_formal_charge(neutral):+d}"
)

print(
    f"7-O- charge:     "
    f"{total_formal_charge(mono_7O):+d}"
)

print(
    f"Dianion charge:  "
    f"{total_formal_charge(dianion_7O_4primeO):+d}"
)

print()

print("Files:")
print(NEUTRAL_SDF)
print(MONO_SDF)
print(DI_SDF)

print()
print("Documentation:")
print(CSV_OUTPUT)
print(XLSX_OUTPUT)

print()
print("=" * 80)