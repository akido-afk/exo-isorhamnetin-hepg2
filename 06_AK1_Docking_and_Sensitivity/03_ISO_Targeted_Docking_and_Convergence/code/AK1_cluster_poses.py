from pathlib import Path
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

PROJECT = Path(r"C:\Users\akidi\Desktop\ISO_AK_Docking")

matrix_file = (
    PROJECT
    / "results"
    / "AK1_ISO_pose_RMSD_matrix.csv"
)

df = pd.read_csv(matrix_file, index_col=0)

# Convert row/column labels to seed numbers
df.index = df.index.astype(int)
df.columns = df.columns.astype(int)

distance_matrix = df.values

# Convert full RMSD matrix into condensed format
condensed = squareform(
    distance_matrix,
    checks=False
)

# Average-linkage clustering
Z = linkage(
    condensed,
    method="average"
)

print()
print("=" * 60)
print("AK1 ISO POSE CLUSTERING")
print("=" * 60)

# Examine both 2 A and 3 A thresholds
for threshold in [2.0, 3.0]:

    labels = fcluster(
        Z,
        t=threshold,
        criterion="distance"
    )

    clusters = {}

    for seed, cluster in zip(df.index, labels):
        clusters.setdefault(cluster, []).append(seed)

    print()
    print(f"RMSD threshold = {threshold:.1f} A")
    print("-" * 40)

    ordered = sorted(
        clusters.items(),
        key=lambda x: len(x[1]),
        reverse=True
    )

    for number, (cluster, seeds) in enumerate(
        ordered,
        start=1
    ):

        print(
            f"Cluster {number}: "
            f"{len(seeds)} poses -> "
            f"seeds {seeds}"
        )

print()