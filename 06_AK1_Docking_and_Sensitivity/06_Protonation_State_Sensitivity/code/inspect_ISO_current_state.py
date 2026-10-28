from pathlib import Path
from rdkit import Chem

PROJECT = Path(r"C:\Users\akidi\Desktop\ISO_AK_Docking")

SDF = PROJECT / "ligand" / "isorhamnetin_3D.sdf"

if not SDF.exists():
    raise FileNotFoundError(SDF)

supplier = Chem.SDMolSupplier(
    str(SDF),
    removeHs=False
)

mol = None

for m in supplier:
    if m is not None:
        mol = m
        break

if mol is None:
    raise RuntimeError("Could not read isorhamnetin SDF.")

print("=" * 75)
print("CURRENT ISORHAMNETIN LIGAND STATE")
print("=" * 75)

print()
print("Input:")
print(SDF)

print()
print("Molecular formula-like atom counts:")

counts = {}

for atom in mol.GetAtoms():
    symbol = atom.GetSymbol()
    counts[symbol] = counts.get(symbol, 0) + 1

for symbol in sorted(counts):
    print(f"  {symbol}: {counts[symbol]}")

formal_charge = sum(
    atom.GetFormalCharge()
    for atom in mol.GetAtoms()
)

print()
print("Total formal charge:", formal_charge)

print()
print("OXYGEN ATOMS")
print("-" * 75)

for atom in mol.GetAtoms():

    if atom.GetSymbol() != "O":
        continue

    idx = atom.GetIdx()

    neighbors = []

    for nbr in atom.GetNeighbors():

        neighbors.append(
            f"{nbr.GetSymbol()}{nbr.GetIdx()+1}"
        )

    print(
        f"O atom index {idx+1:>2} | "
        f"formal charge {atom.GetFormalCharge():>2} | "
        f"explicit H {atom.GetNumExplicitHs()} | "
        f"implicit H {atom.GetNumImplicitHs()} | "
        f"neighbors: {', '.join(neighbors)}"
    )

print()
print("Canonical SMILES:")
print(
    Chem.MolToSmiles(
        mol,
        canonical=True
    )
)

print()
print("=" * 75)