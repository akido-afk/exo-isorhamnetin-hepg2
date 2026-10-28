from pathlib import Path
import csv

import xlsxwriter
from pdbfixer import PDBFixer

try:
    from openmm.app import PDBFile
except ImportError:
    from simtk.openmm.app import PDBFile


# ============================================================
# PROJECT PATHS
# ============================================================

PROJECT = Path(
    r"C:\Users\akidi\Desktop\ISO_AK_Docking"
)

STRUCTURES = PROJECT / "structures"
RECEPTORS = PROJECT / "receptors"

OUT_DIR = (
    PROJECT
    / "results"
    / "AK1_ensemble_sensitivity"
    / "receptor_preparation"
)

RECEPTORS.mkdir(
    parents=True,
    exist_ok=True
)

OUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# INPUT MODELS
# ============================================================

MODELS = {

    "7X7S_model20_least_open": (
        STRUCTURES
        / "AK1_7X7S_model20_least_open_ALIGNED_1Z83.pdb"
    ),

    "7X7S_model8_ensemble_medoid": (
        STRUCTURES
        / "AK1_7X7S_model8_ensemble_medoid_ALIGNED_1Z83.pdb"
    ),

    "7X7S_model1_most_open": (
        STRUCTURES
        / "AK1_7X7S_model1_most_open_ALIGNED_1Z83.pdb"
    )
}


# ============================================================
# HELPERS
# ============================================================

def count_topology(topology):

    chains = list(
        topology.chains()
    )

    residues = list(
        topology.residues()
    )

    atoms = list(
        topology.atoms()
    )

    return (
        len(chains),
        len(residues),
        len(atoms)
    )


# ============================================================
# PREPARE RECEPTORS
# ============================================================

summary_rows = []


for label, input_file in MODELS.items():

    if not input_file.exists():

        raise FileNotFoundError(
            input_file
        )


    print()
    print("=" * 80)
    print(
        f"PROCESSING {label}"
    )
    print("=" * 80)


    fixer = PDBFixer(
        filename=str(
            input_file
        )
    )


    # --------------------------------------------------------
    # Initial topology
    # --------------------------------------------------------

    (
        chains_before,
        residues_before,
        atoms_before
    ) = count_topology(
        fixer.topology
    )


    # --------------------------------------------------------
    # Missing residues
    #
    # We inspect but deliberately do not reconstruct
    # whole missing residues/terminal extensions.
    # The 7X7S models already contain the complete
    # experimentally deposited 194-residue chain.
    # --------------------------------------------------------

    fixer.findMissingResidues()

    missing_residue_positions = (
        len(
            fixer.missingResidues
        )
    )


    # Do NOT build missing whole residues.
    fixer.missingResidues = {}


    # --------------------------------------------------------
    # Nonstandard residues
    # --------------------------------------------------------

    fixer.findNonstandardResidues()

    nonstandard_count = len(
        fixer.nonstandardResidues
    )


    if nonstandard_count > 0:

        fixer.replaceNonstandardResidues()


    # --------------------------------------------------------
    # Remove heterogens / waters
    #
    # keepWater=False means remove heterogens including water.
    # This matches our protein-only receptor sensitivity setup.
    # --------------------------------------------------------

    fixer.removeHeterogens(
        False
    )


    # --------------------------------------------------------
    # Missing heavy atoms
    # --------------------------------------------------------

    fixer.findMissingAtoms()


    missing_atom_count = sum(
        len(atom_list)
        for atom_list
        in fixer.missingAtoms.values()
    )


    missing_terminal_atom_count = sum(
        len(atom_list)
        for atom_list
        in fixer.missingTerminals.values()
    )


    # --------------------------------------------------------
    # Repair present residues
    # --------------------------------------------------------

    fixer.addMissingAtoms()


    # --------------------------------------------------------
    # Add H at pH 7.4
    # --------------------------------------------------------

    fixer.addMissingHydrogens(
        7.4
    )


    # --------------------------------------------------------
    # Final topology
    # --------------------------------------------------------

    (
        chains_after,
        residues_after,
        atoms_after
    ) = count_topology(
        fixer.topology
    )


    output_file = (
        RECEPTORS
        / f"AK1_{label}_FIXED.pdb"
    )


    with open(
        output_file,
        "w"
    ) as handle:

        PDBFile.writeFile(
            fixer.topology,
            fixer.positions,
            handle,
            keepIds=True
        )


    summary_rows.append({

        "Receptor":
            label,

        "Input_PDB":
            str(
                input_file
            ),

        "Output_PDB":
            str(
                output_file
            ),

        "Chains_before":
            chains_before,

        "Residues_before":
            residues_before,

        "Atoms_before":
            atoms_before,

        "Missing_residue_positions_detected":
            missing_residue_positions,

        "Nonstandard_residues_detected":
            nonstandard_count,

        "Missing_heavy_atoms_detected":
            missing_atom_count,

        "Missing_terminal_atoms_detected":
            missing_terminal_atom_count,

        "Chains_after":
            chains_after,

        "Residues_after":
            residues_after,

        "Atoms_after_H_addition":
            atoms_after
    })


    print(
        f"Residues before : {residues_before}"
    )

    print(
        f"Atoms before    : {atoms_before}"
    )

    print(
        f"Missing residue positions detected : "
        f"{missing_residue_positions}"
    )

    print(
        f"Nonstandard residues detected      : "
        f"{nonstandard_count}"
    )

    print(
        f"Missing heavy atoms detected        : "
        f"{missing_atom_count}"
    )

    print(
        f"Missing terminal atoms detected     : "
        f"{missing_terminal_atom_count}"
    )

    print(
        f"Residues after  : {residues_after}"
    )

    print(
        f"Atoms after H   : {atoms_after}"
    )

    print()
    print(
        "Saved:"
    )

    print(
        output_file
    )


# ============================================================
# SAVE CSV
# ============================================================

CSV_FILE = (
    OUT_DIR
    / "AK1_7X7S_receptor_preparation_summary.csv"
)


headers = list(
    summary_rows[0].keys()
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
        summary_rows
    )


# ============================================================
# SAVE EXCEL
# ============================================================

XLSX_FILE = (
    OUT_DIR
    / "AK1_7X7S_RECEPTOR_PREPARATION.xlsx"
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


ws = workbook.add_worksheet(
    "Preparation_Summary"
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
    summary_rows,
    start=1
):

    for col, header in enumerate(
        headers
    ):

        ws.write(
            row_i,
            col,
            row[
                header
            ],
            cell_fmt
        )


ws.freeze_panes(
    1,
    0
)

ws.set_column(
    0,
    len(headers) - 1,
    30
)


# ------------------------------------------------------------
# Method notes
# ------------------------------------------------------------

ws2 = workbook.add_worksheet(
    "Method_Notes"
)


notes = [

    [
        "Purpose",
        (
            "Standardize three aligned 7X7S apo conformers "
            "before receptor-state sensitivity docking."
        )
    ],

    [
        "Whole missing residues",
        (
            "Detected missing-residue positions were inspected "
            "but not reconstructed. The deposited 7X7S models "
            "contain the experimentally represented 194-residue "
            "human AK1 chain."
        )
    ],

    [
        "Missing atoms",
        (
            "Missing atoms within residues already present in "
            "the structure were repaired with PDBFixer."
        )
    ],

    [
        "Hydrogens",
        (
            "Hydrogens added by PDBFixer at pH 7.4."
        )
    ],

    [
        "Heterogens",
        (
            "Heterogens and water molecules removed for the "
            "protein-only ensemble sensitivity experiment."
        )
    ],

    [
        "Alignment",
        (
            "The input conformers had already been rigid-body "
            "aligned into the 1Z83 coordinate frame using "
            "sequence-mapped CORE C-alpha atoms."
        )
    ],

    [
        "Ligand",
        (
            "The subsequent ensemble docking will use the same "
            "neutral isorhamnetin PDBQT used in the original "
            "1Z83 experiment, isolating receptor conformation "
            "as the changing variable."
        )
    ]
]


for row_i, (
    label,
    description
) in enumerate(
    notes
):

    ws2.write(
        row_i,
        0,
        label,
        header_fmt
    )

    ws2.write(
        row_i,
        1,
        description,
        cell_fmt
    )


ws2.set_column(
    "A:A",
    28
)

ws2.set_column(
    "B:B",
    110
)


workbook.close()


# ============================================================
# FINAL
# ============================================================

print()
print("=" * 80)
print("7X7S RECEPTOR PREPARATION COMPLETE")
print("=" * 80)


for row in summary_rows:

    print()

    print(
        row[
            "Receptor"
        ]
    )

    print(
        f"  residues: "
        f"{row['Residues_before']} -> "
        f"{row['Residues_after']}"
    )

    print(
        f"  atoms: "
        f"{row['Atoms_before']} -> "
        f"{row['Atoms_after_H_addition']}"
    )

    print(
        f"  missing heavy atoms repaired: "
        f"{row['Missing_heavy_atoms_detected']}"
    )

    print(
        f"  output: "
        f"{row['Output_PDB']}"
    )


print()
print("Saved:")
print(CSV_FILE)
print(XLSX_FILE)

print()
print(
    "PDBFixer stage complete."
)