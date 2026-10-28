import subprocess
from pathlib import Path

PROJECT = Path(r"C:\Users\akidi\Desktop\ISO_AK_Docking")

VINA = Path(r"C:\Users\akidi\Desktop\Docking tools\vina_win.exe")
RECEPTOR = PROJECT / "receptors" / "AK1.pdbqt"
LIGAND = PROJECT / "ligand" / "isorhamnetin.pdbqt"
OUTPUT_DIR = PROJECT / "docking"

# Same experimentally defined AP5 pocket used for seed 1001
CENTER_X = 25.044
CENTER_Y = 44.571
CENTER_Z = -3.708

SIZE_X = 23.687
SIZE_Y = 23.364
SIZE_Z = 30.062

# Ten independent seeds
SEEDS = [
    1001,
    2002,
    3003,
    4004,
    5005,
    6006,
    7007,
    8008,
    9009,
    10101
]

for seed in SEEDS:

    output_file = OUTPUT_DIR / f"AK1_ISO_seed{seed}.pdbqt"

    # Protect completed runs from accidental overwrite
    if output_file.exists():
        print()
        print("=" * 60)
        print(f"Seed {seed} already exists - SKIPPING")
        print("=" * 60)
        continue

    print()
    print("=" * 60)
    print(f"STARTING AK1 + ISO : SEED {seed}")
    print("=" * 60)

    command = [
        str(VINA),

        "--receptor", str(RECEPTOR),
        "--ligand", str(LIGAND),

        "--center_x", str(CENTER_X),
        "--center_y", str(CENTER_Y),
        "--center_z", str(CENTER_Z),

        "--size_x", str(SIZE_X),
        "--size_y", str(SIZE_Y),
        "--size_z", str(SIZE_Z),

        "--exhaustiveness", "16",
        "--num_modes", "10",

        "--seed", str(seed),

        "--out", str(output_file)
    ]

    subprocess.run(command, check=True)

    print()
    print(f"FINISHED SEED {seed}")
    print(f"Saved: {output_file}")


print()
print("=" * 60)
print("ALL AVAILABLE AK1 MULTI-SEED RUNS FINISHED")
print("=" * 60)