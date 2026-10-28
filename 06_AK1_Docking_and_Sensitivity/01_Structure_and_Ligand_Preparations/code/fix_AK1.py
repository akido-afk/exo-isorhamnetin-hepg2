from pdbfixer import PDBFixer
from openmm.app import PDBFile

fixer = PDBFixer(filename=r"..\receptors\AK1_chainA.pdb")

fixer.findMissingResidues()
fixer.findMissingAtoms()

# Do not create entire missing loops
fixer.missingResidues = {}

# Repair missing atoms in residues that are already present
fixer.addMissingAtoms()

# Add hydrogens at approximately physiological pH
fixer.addMissingHydrogens(7.4)

with open(r"..\receptors\AK1_chainA_FIXED.pdb", "w") as f:
    PDBFile.writeFile(
        fixer.topology,
        fixer.positions,
        f,
        keepIds=True
    )

print("AK1 repair finished")