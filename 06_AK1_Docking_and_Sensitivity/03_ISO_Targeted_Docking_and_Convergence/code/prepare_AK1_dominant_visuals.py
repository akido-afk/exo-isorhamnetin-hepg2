from pathlib import Path
import subprocess
import sys

from rdkit import Chem


PROJECT = Path(r"C:\Users\akidi\Desktop\ISO_AK_Docking")

DOCKING = PROJECT / "docking"
RESULTS = PROJECT / "results"

VISUALS = RESULTS / "AK1_final_visuals"
POSES = VISUALS / "poses"
FRAMES = VISUALS / "rotation_frames"

VISUALS.mkdir(exist_ok=True)
POSES.mkdir(exist_ok=True)
FRAMES.mkdir(exist_ok=True)


# Dominant 8-pose cluster
SEEDS = [
    1001,
    2002,
    3003,
    4004,
    5005,
    7007,
    9009,
    10101
]


meeko_export = (
    Path(sys.prefix)
    / "Lib"
    / "site-packages"
    / "meeko"
    / "cli"
    / "mk_export.py"
)


for seed in SEEDS:

    pdbqt = DOCKING / f"AK1_ISO_seed{seed}.pdbqt"

    all_sdf = POSES / f"AK1_ISO_seed{seed}_ALL.sdf"

    best_sdf = POSES / f"AK1_ISO_seed{seed}_best.sdf"

    print()
    print("=" * 60)
    print(f"Preparing seed {seed}")
    print("=" * 60)

    subprocess.run(
        [
            sys.executable,
            str(meeko_export),
            str(pdbqt),
            "-s",
            str(all_sdf)
        ],
        check=True
    )

    supplier = Chem.SDMolSupplier(
        str(all_sdf),
        removeHs=False
    )

    first_molecule = None

    for mol in supplier:

        if mol is not None:
            first_molecule = mol
            break

    if first_molecule is None:
        raise RuntimeError(
            f"Could not extract best pose from seed {seed}"
        )

    writer = Chem.SDWriter(
        str(best_sdf)
    )

    writer.write(first_molecule)
    writer.close()

    print("Best pose saved:")
    print(best_sdf)


print()
print("=" * 60)
print("ALL 8 DOMINANT POSES PREPARED")
print("=" * 60)
print()
print("Representative pose = seed 9009")