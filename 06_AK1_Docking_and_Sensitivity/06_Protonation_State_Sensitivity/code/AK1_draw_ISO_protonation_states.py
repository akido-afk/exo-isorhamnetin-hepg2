from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image


# ============================================================
# PATHS
# ============================================================

PROJECT = Path(
    r"C:\Users\akidi\Desktop\ISO_AK_Docking"
)

LIGAND_DIR = (
    PROJECT /
    "ligand" /
    "ISO_protonation_states"
)

OUT_DIR = (
    PROJECT /
    "results" /
    "AK1_ligand_state_sensitivity" /
    "figures" /
    "chemical_states"
)

OUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# INPUT FILES
# ============================================================

FILES = {

    "Neutral ISO":
        LIGAND_DIR /
        "isorhamnetin_state0_neutral.sdf",

    "7-O\u207B monoanion":
        LIGAND_DIR /
        "isorhamnetin_state1_7O_minus.sdf",

    "7-O\u207B/4\u2032-O\u207B dianion":
        LIGAND_DIR /
        "isorhamnetin_state2_7O_4primeO_dianion.sdf"
}


# ============================================================
# OUTPUT FILES
# ============================================================

PNG_NEUTRAL = (
    OUT_DIR /
    "AK1_ISO_state0_neutral_structure.png"
)

PNG_MONO = (
    OUT_DIR /
    "AK1_ISO_state1_7O_monoanion_structure.png"
)

PNG_DI = (
    OUT_DIR /
    "AK1_ISO_state2_7O_4primeO_dianion_structure.png"
)

PNG_COMBINED = (
    OUT_DIR /
    "AK1_ISO_protonation_states_2D.png"
)


# ============================================================
# LOAD MOLECULE
# ============================================================

def load_molecule(path):

    supplier = Chem.SDMolSupplier(
        str(path),
        removeHs=False
    )

    mol = None

    for candidate in supplier:

        if candidate is not None:
            mol = candidate
            break

    if mol is None:

        raise RuntimeError(
            f"Could not read molecule:\n{path}"
        )

    # Remove explicit H for cleaner publication structure
    mol = Chem.RemoveHs(
        mol
    )

    # Generate clean 2D coordinates
    AllChem.Compute2DCoords(
        mol
    )

    return mol


# ============================================================
# FIND NEGATIVELY CHARGED OXYGENS
# ============================================================

def find_negative_oxygens(mol):

    return [
        atom.GetIdx()
        for atom in mol.GetAtoms()
        if (
            atom.GetAtomicNum() == 8
            and
            atom.GetFormalCharge() == -1
        )
    ]


# ============================================================
# DRAW ONE MOLECULE
# ============================================================

def draw_molecule(
    mol,
    title,
    output_png,
    width=1200,
    height=850
):

    negative_oxygens = (
        find_negative_oxygens(
            mol
        )
    )


    drawer = rdMolDraw2D.MolDraw2DCairo(
        width,
        height
    )


    options = (
        drawer.drawOptions()
    )

    options.addAtomIndices = False

    options.bondLineWidth = 3.0

    options.baseFontSize = 0.60

    options.padding = 0.10

    options.legendFontSize = 36


    # Highlight only formally negative oxygens
    highlight_atoms = (
        negative_oxygens
    )


    # Let RDKit draw the formal charge directly
    drawer.DrawMolecule(
        mol,
        legend=title,
        highlightAtoms=highlight_atoms
    )


    drawer.FinishDrawing()


    output_png.write_bytes(
        drawer.GetDrawingText()
    )


# ============================================================
# LOAD STATES
# ============================================================

neutral = load_molecule(
    FILES["Neutral ISO"]
)

mono = load_molecule(
    FILES["7-O\u207B monoanion"]
)

dianion = load_molecule(
    FILES[
        "7-O\u207B/4\u2032-O\u207B dianion"
    ]
)


# ============================================================
# DRAW INDIVIDUAL STRUCTURES
# ============================================================

draw_molecule(
    neutral,
    "Neutral isorhamnetin",
    PNG_NEUTRAL
)

draw_molecule(
    mono,
    "7-O\u207B isorhamnetin",
    PNG_MONO
)

draw_molecule(
    dianion,
    "7-O\u207B / 4\u2032-O\u207B isorhamnetin",
    PNG_DI
)


# ============================================================
# BUILD COMBINED IMAGE
# ============================================================

img_neutral = Image.open(
    PNG_NEUTRAL
).convert(
    "RGB"
)

img_mono = Image.open(
    PNG_MONO
).convert(
    "RGB"
)

img_di = Image.open(
    PNG_DI
).convert(
    "RGB"
)


width = max(
    img_neutral.width,
    img_mono.width,
    img_di.width
)

height = max(
    img_neutral.height,
    img_mono.height,
    img_di.height
)


canvas = Image.new(
    "RGB",
    (
        width * 3,
        height
    ),
    "white"
)


canvas.paste(
    img_neutral,
    (
        0,
        0
    )
)

canvas.paste(
    img_mono,
    (
        width,
        0
    )
)

canvas.paste(
    img_di,
    (
        width * 2,
        0
    )
)


canvas.save(
    PNG_COMBINED,
    dpi=(
        600,
        600
    )
)


# ============================================================
# VERIFY FORMAL CHARGES
# ============================================================

def total_charge(mol):

    return sum(
        atom.GetFormalCharge()
        for atom in mol.GetAtoms()
    )


print()
print("=" * 80)
print("ISORHAMNETIN PROTONATION-STATE STRUCTURES")
print("=" * 80)

print()
print(
    "Neutral formal charge:",
    total_charge(neutral)
)

print(
    "Monoanion formal charge:",
    total_charge(mono)
)

print(
    "Dianion formal charge:",
    total_charge(dianion)
)

print()

print(
    "Neutral negative oxygens:",
    find_negative_oxygens(
        neutral
    )
)

print(
    "Monoanion negative oxygens:",
    find_negative_oxygens(
        mono
    )
)

print(
    "Dianion negative oxygens:",
    find_negative_oxygens(
        dianion
    )
)

print()
print("Saved:")
print(PNG_NEUTRAL)
print(PNG_MONO)
print(PNG_DI)
print(PNG_COMBINED)

print()
print("=" * 80)