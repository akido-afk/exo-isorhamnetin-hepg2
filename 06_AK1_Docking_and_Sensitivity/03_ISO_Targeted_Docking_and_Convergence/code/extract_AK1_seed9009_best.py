from rdkit import Chem
from pathlib import Path

inp = Path(r"..\results\AK1_repr9009\AK1_ISO_seed9009_allposes.sdf")
out = Path(r"..\results\AK1_repr9009\AK1_ISO_seed9009_best.sdf")

supplier = Chem.SDMolSupplier(str(inp), removeHs=False)

best = None
for mol in supplier:
    if mol is not None:
        best = mol
        break

if best is None:
    raise RuntimeError("No valid molecule found in SDF.")

writer = Chem.SDWriter(str(out))
writer.write(best)
writer.close()

print("Best pose saved to:")
print(out)