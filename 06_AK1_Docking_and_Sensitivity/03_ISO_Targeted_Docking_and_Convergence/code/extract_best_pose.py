from rdkit import Chem

jobs = [
    (r"..\docking\AK1_ISO_seed1001.sdf", r"..\docking\AK1_ISO_bestpose.sdf"),
    (r"..\docking\AK1_AP5_redock_seed1001.sdf", r"..\docking\AK1_AP5_redock_bestpose.sdf"),
]

for infile, outfile in jobs:
    suppl = Chem.SDMolSupplier(infile, removeHs=False)
    mol = None
    for m in suppl:
        if m is not None:
            mol = m
            break

    if mol is None:
        raise ValueError(f"No valid molecule found in {infile}")

    writer = Chem.SDWriter(outfile)
    writer.write(mol)
    writer.close()

    print("Saved:", outfile)