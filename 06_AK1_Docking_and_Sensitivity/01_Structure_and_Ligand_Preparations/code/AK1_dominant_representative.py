from pathlib import Path
import pandas as pd

PROJECT = Path(r"C:\Users\akidi\Desktop\ISO_AK_Docking")

matrix_file = (
    PROJECT
    / "results"
    / "AK1_ISO_pose_RMSD_matrix.csv"
)

df = pd.read_csv(matrix_file, index_col=0)

df.index = df.index.astype(int)
df.columns = df.columns.astype(int)

dominant = [
    1001,
    2002,
    3003,
    4004,
    5005,
    7007,
    9009,
    10101
]

sub = df.loc[dominant, dominant]

means = {}

for seed in dominant:

    values = [
        sub.loc[seed, other]
        for other in dominant
        if other != seed
    ]

    means[seed] = sum(values) / len(values)


print()
print("=" * 60)
print("AK1 DOMINANT-CLUSTER REPRESENTATIVE POSE")
print("=" * 60)

for seed in sorted(means):
    print(
        f"Seed {seed}: "
        f"mean RMSD to dominant cluster = "
        f"{means[seed]:.3f} A"
    )

representative = min(
    means,
    key=means.get
)

print()
print("-------------------------------------------")
print(
    f"Representative dominant pose: seed {representative}"
)
print(
    f"Mean RMSD within dominant cluster: "
    f"{means[representative]:.3f} A"
)
print("-------------------------------------------")