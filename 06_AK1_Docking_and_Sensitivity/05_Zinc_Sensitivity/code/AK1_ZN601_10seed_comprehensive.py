from pathlib import Path
import subprocess
import sys
import math
import csv
import statistics
from collections import Counter

from rdkit import Chem
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
import numpy as np
import xlsxwriter


# ============================================================
# PROJECT PATHS
# ============================================================

PROJECT = Path(r"C:\Users\akidi\Desktop\ISO_AK_Docking")

DOCKING = PROJECT / "docking"
RECEPTORS = PROJECT / "receptors"

OUTPUT_DIR = (
    PROJECT
    / "results"
    / "AK1_ZN_sensitivity"
    / "ZN601_10seed"
)

POSE_DIR = OUTPUT_DIR / "poses"

ORIGINAL_POSE_DIR = (
    OUTPUT_DIR
    / "original_Znfree_reference_poses"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

POSE_DIR.mkdir(
    parents=True,
    exist_ok=True
)

ORIGINAL_POSE_DIR.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# SEEDS
# ============================================================

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


# ============================================================
# ORIGINAL Zn-FREE CLUSTERS
# ============================================================

ORIGINAL_DOMINANT = {
    1001,
    2002,
    3003,
    4004,
    5005,
    7007,
    9009,
    10101
}

ORIGINAL_ALTERNATIVE = {
    6006,
    8008
}


# ============================================================
# MEEKO EXPORT
# ============================================================

MEEKO_EXPORT = (
    Path(sys.prefix)
    / "Lib"
    / "site-packages"
    / "meeko"
    / "cli"
    / "mk_export.py"
)


# ============================================================
# ZN601 COORDINATE
# ============================================================

ZN_RECEPTOR = (
    RECEPTORS
    / "AK1_chainA_FIXED_ZN601.pdb"
)

zn_coord = None

with open(ZN_RECEPTOR, "r") as f:

    for line in f:

        if (
            line.startswith("HETATM")
            and line[17:20].strip() == "ZN"
            and line[22:26].strip() == "601"
        ):

            zn_coord = (
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54])
            )

            break


if zn_coord is None:

    raise RuntimeError(
        "Could not locate ZN601."
    )


# ============================================================
# READ BEST VINA SCORE
# ============================================================

def read_best_score(pdbqt):

    with open(
        pdbqt,
        "r",
        encoding="utf-8"
    ) as f:

        for line in f:

            if line.startswith(
                "REMARK VINA RESULT:"
            ):

                return float(
                    line.split()[3]
                )

    raise RuntimeError(
        f"No score found in {pdbqt}"
    )


# ============================================================
# EXPORT BEST POSE
# ============================================================

def export_best_pose(
    pdbqt,
    output_dir,
    prefix
):

    all_sdf = (
        output_dir
        / f"{prefix}_allposes.sdf"
    )

    best_sdf = (
        output_dir
        / f"{prefix}_best.sdf"
    )


    if not all_sdf.exists():

        subprocess.run(
            [
                sys.executable,
                str(MEEKO_EXPORT),
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


    mol = None

    for candidate in supplier:

        if candidate is not None:

            mol = candidate
            break


    if mol is None:

        raise RuntimeError(
            f"Could not extract pose from {pdbqt}"
        )


    if not best_sdf.exists():

        writer = Chem.SDWriter(
            str(best_sdf)
        )

        writer.write(mol)
        writer.close()


    return mol, best_sdf


# ============================================================
# POSITIONAL RMSD
#
# No superposition:
# all poses are already in the common AK1 coordinate frame.
# ============================================================

def positional_rmsd(
    mol1,
    mol2
):

    if (
        mol1.GetNumAtoms()
        != mol2.GetNumAtoms()
    ):

        raise RuntimeError(
            "Different ligand atom counts."
        )


    c1 = mol1.GetConformer()
    c2 = mol2.GetConformer()

    values = []


    for i in range(
        mol1.GetNumAtoms()
    ):

        a1 = mol1.GetAtomWithIdx(i)
        a2 = mol2.GetAtomWithIdx(i)

        if (
            a1.GetAtomicNum()
            != a2.GetAtomicNum()
        ):

            raise RuntimeError(
                f"Atom-order mismatch at index {i}"
            )


        if a1.GetAtomicNum() == 1:
            continue


        p1 = c1.GetAtomPosition(i)
        p2 = c2.GetAtomPosition(i)


        values.append(
            (p1.x - p2.x) ** 2
            +
            (p1.y - p2.y) ** 2
            +
            (p1.z - p2.z) ** 2
        )


    return math.sqrt(
        sum(values)
        /
        len(values)
    )


# ============================================================
# ZN601 -> ISO DISTANCES
# ============================================================

def zn_distances(mol):

    conf = mol.GetConformer()

    minimum_heavy = None
    heavy_atom = None

    minimum_polar = None
    polar_atom = None


    for atom in mol.GetAtoms():

        if atom.GetAtomicNum() == 1:
            continue


        i = atom.GetIdx()

        p = conf.GetAtomPosition(i)


        d = math.sqrt(
            (p.x - zn_coord[0]) ** 2
            +
            (p.y - zn_coord[1]) ** 2
            +
            (p.z - zn_coord[2]) ** 2
        )


        atom_name = (
            f"{atom.GetSymbol()}{i + 1}"
        )


        if (
            minimum_heavy is None
            or d < minimum_heavy
        ):

            minimum_heavy = d
            heavy_atom = atom_name


        if (
            atom.GetSymbol().upper()
            in {"O", "N", "S"}
        ):

            if (
                minimum_polar is None
                or d < minimum_polar
            ):

                minimum_polar = d
                polar_atom = atom_name


    return (
        minimum_heavy,
        heavy_atom,
        minimum_polar,
        polar_atom
    )


# ============================================================
# LOAD ALL ZN601 POSES + SCORES
# ============================================================

zn_poses = {}
zn_scores = {}

zn_rows = []


for seed in SEEDS:

    pdbqt = (
        DOCKING
        / f"AK1_ZN601_ISO_seed{seed}.pdbqt"
    )


    if not pdbqt.exists():

        raise FileNotFoundError(
            pdbqt
        )


    score = read_best_score(
        pdbqt
    )


    mol, best_sdf = export_best_pose(
        pdbqt,
        POSE_DIR,
        f"AK1_ZN601_ISO_seed{seed}"
    )


    zn_poses[seed] = mol
    zn_scores[seed] = score


    (
        min_heavy,
        heavy_atom,
        min_polar,
        polar_atom
    ) = zn_distances(mol)


    zn_rows.append({
        "Seed": seed,
        "Best Vina score": score,
        "Closest ISO heavy atom": heavy_atom,
        "Closest heavy-atom distance to ZN601 (A)": min_heavy,
        "Closest ISO polar atom": polar_atom,
        "Closest polar distance to ZN601 (A)": min_polar
    })


# ============================================================
# ORIGINAL Zn-FREE POSES + SCORES
# ============================================================

original_poses = {}
original_scores = {}


for seed in SEEDS:

    pdbqt = (
        DOCKING
        / f"AK1_ISO_seed{seed}.pdbqt"
    )


    if not pdbqt.exists():

        raise FileNotFoundError(
            pdbqt
        )


    original_scores[seed] = (
        read_best_score(pdbqt)
    )


    mol, _ = export_best_pose(
        pdbqt,
        ORIGINAL_POSE_DIR,
        f"AK1_ISO_seed{seed}"
    )


    original_poses[seed] = mol


# ============================================================
# ZN601 PAIRWISE RMSD MATRIX
# ============================================================

n = len(SEEDS)

rmsd_matrix = np.zeros(
    (n, n),
    dtype=float
)


for i, seed_i in enumerate(SEEDS):

    for j, seed_j in enumerate(SEEDS):

        if j < i:

            rmsd_matrix[
                i,
                j
            ] = rmsd_matrix[
                j,
                i
            ]

            continue


        if i == j:

            value = 0.0

        else:

            value = positional_rmsd(
                zn_poses[seed_i],
                zn_poses[seed_j]
            )


        rmsd_matrix[
            i,
            j
        ] = value

        rmsd_matrix[
            j,
            i
        ] = value


# ============================================================
# CLUSTERING FUNCTION
# Complete-link hierarchical clustering
# ============================================================

def cluster_at_threshold(
    matrix,
    threshold
):

    condensed = squareform(
        matrix,
        checks=False
    )


    Z = linkage(
        condensed,
        method="complete"
    )


    raw_labels = fcluster(
        Z,
        t=threshold,
        criterion="distance"
    )


    # Remap labels by cluster size
    counts = Counter(
        raw_labels
    )


    ordered_old_labels = sorted(
        counts.keys(),
        key=lambda lab:
            (
                -counts[lab],
                lab
            )
    )


    mapping = {
        old:
        new_index + 1

        for new_index, old
        in enumerate(
            ordered_old_labels
        )
    }


    labels = [
        mapping[x]
        for x in raw_labels
    ]


    return labels


clusters_2A = cluster_at_threshold(
    rmsd_matrix,
    2.0
)

clusters_3A = cluster_at_threshold(
    rmsd_matrix,
    3.0
)


# ============================================================
# DOMINANT ZN601 CLUSTER + REPRESENTATIVE POSE
# ============================================================

counts_2A = Counter(
    clusters_2A
)

dominant_cluster_id = max(
    counts_2A,
    key=counts_2A.get
)


dominant_indices = [
    i for i, c
    in enumerate(clusters_2A)
    if c == dominant_cluster_id
]


dominant_seeds = [
    SEEDS[i]
    for i in dominant_indices
]


representative_candidates = []


for i in dominant_indices:

    other_indices = [
        j for j
        in dominant_indices
        if j != i
    ]


    if not other_indices:

        mean_rmsd = 0.0

    else:

        mean_rmsd = statistics.mean(
            rmsd_matrix[
                i,
                j
            ]
            for j in other_indices
        )


    representative_candidates.append(
        (
            mean_rmsd,
            SEEDS[i]
        )
    )


representative_candidates.sort()

representative_mean_rmsd = (
    representative_candidates[0][0]
)

representative_seed = (
    representative_candidates[0][1]
)


# ============================================================
# COMPARE ZN601 POSES WITH ORIGINAL Zn-FREE FAMILIES
# ============================================================

comparison_rows = []


for zn_seed in SEEDS:

    rmsd_to_original = {}


    for original_seed in SEEDS:

        rmsd_to_original[
            original_seed
        ] = positional_rmsd(
            zn_poses[zn_seed],
            original_poses[
                original_seed
            ]
        )


    dominant_values = [
        rmsd_to_original[s]
        for s in ORIGINAL_DOMINANT
    ]


    alternative_values = [
        rmsd_to_original[s]
        for s in ORIGINAL_ALTERNATIVE
    ]


    mean_dom = statistics.mean(
        dominant_values
    )

    mean_alt = statistics.mean(
        alternative_values
    )


    nearest_seed = min(
        rmsd_to_original,
        key=rmsd_to_original.get
    )


    nearest_rmsd = (
        rmsd_to_original[
            nearest_seed
        ]
    )


    nearest_family = (
        "Dominant"
        if nearest_seed
        in ORIGINAL_DOMINANT
        else "Alternative"
    )


    closer_family = (
        "Dominant"
        if mean_dom < mean_alt
        else "Alternative"
    )


    # Conservative match classification
    if nearest_rmsd <= 2.0:

        family_match = (
            f"Matches original {nearest_family} family"
        )

    else:

        family_match = (
            "Distinct from original Zn-free families"
        )


    comparison_rows.append({
        "ZN601 seed": zn_seed,
        "Mean RMSD to original dominant (A)": mean_dom,
        "Mean RMSD to original alternative (A)": mean_alt,
        "Closer original family": closer_family,
        "Nearest original seed": nearest_seed,
        "Nearest original family": nearest_family,
        "Nearest original RMSD (A)": nearest_rmsd,
        "Family match interpretation": family_match
    })


# ============================================================
# SCORE COMPARISON
# ============================================================

score_rows = []


for seed in SEEDS:

    diff = (
        zn_scores[seed]
        -
        original_scores[seed]
    )


    score_rows.append({
        "Seed": seed,
        "Zn-free score": original_scores[seed],
        "ZN601-retained score": zn_scores[seed],
        "Score difference ZN601-minus-Znfree": diff
    })


zn_score_values = [
    zn_scores[s]
    for s in SEEDS
]

original_score_values = [
    original_scores[s]
    for s in SEEDS
]

score_differences = [
    r[
        "Score difference ZN601-minus-Znfree"
    ]
    for r in score_rows
]


# ============================================================
# PRINT SUMMARY
# ============================================================

print()
print("=" * 90)
print("AK1 ZN601 — 10-SEED COMPREHENSIVE ANALYSIS")
print("=" * 90)

print()

print("ZN601-retained scores")
print("---------------------")

for row in zn_rows:

    print(
        f"Seed {row['Seed']}: "
        f"{row['Best Vina score']:.3f} kcal/mol   "
        f"Zn-polar distance = "
        f"{row['Closest polar distance to ZN601 (A)']:.3f} A"
    )


print()
print("Score statistics")
print("----------------")

print(
    "ZN601 mean:",
    f"{statistics.mean(zn_score_values):.3f}"
)

print(
    "ZN601 median:",
    f"{statistics.median(zn_score_values):.3f}"
)

print(
    "ZN601 SD:",
    f"{statistics.stdev(zn_score_values):.3f}"
)

print(
    "ZN601 range:",
    f"{max(zn_score_values)-min(zn_score_values):.3f}"
)

print()

print(
    "Zn-free mean:",
    f"{statistics.mean(original_score_values):.3f}"
)

print(
    "Mean paired score shift:",
    f"{statistics.mean(score_differences):+.3f} kcal/mol"
)


print()
print("2.0 A clustering")
print("----------------")

for cluster_id in sorted(
    set(clusters_2A)
):

    seeds = [
        SEEDS[i]
        for i, c
        in enumerate(clusters_2A)
        if c == cluster_id
    ]

    print(
        f"Cluster {cluster_id}: "
        f"{len(seeds)} poses -> {seeds}"
    )


print()
print("3.0 A clustering")
print("----------------")

for cluster_id in sorted(
    set(clusters_3A)
):

    seeds = [
        SEEDS[i]
        for i, c
        in enumerate(clusters_3A)
        if c == cluster_id
    ]

    print(
        f"Cluster {cluster_id}: "
        f"{len(seeds)} poses -> {seeds}"
    )


print()
print(
    "Dominant ZN601 cluster:",
    dominant_seeds
)

print(
    "Representative ZN601 pose:",
    representative_seed
)

print(
    "Representative mean RMSD:",
    f"{representative_mean_rmsd:.3f} A"
)


print()
print(
    "Comparison with original Zn-free families"
)

print(
    "-----------------------------------------"
)


for row in comparison_rows:

    print(
        f"Seed {row['ZN601 seed']}: "
        f"mean-to-dominant "
        f"{row['Mean RMSD to original dominant (A)']:.3f} A, "
        f"mean-to-alternative "
        f"{row['Mean RMSD to original alternative (A)']:.3f} A, "
        f"nearest original "
        f"{row['Nearest original seed']} "
        f"({row['Nearest original family']}) "
        f"{row['Nearest original RMSD (A)']:.3f} A"
    )


# ============================================================
# SAVE CSV FILES
# ============================================================

def save_csv(
    path,
    rows
):

    with open(
        path,
        "w",
        newline="",
        encoding="utf-8"
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=list(
                rows[0].keys()
            )
        )

        writer.writeheader()
        writer.writerows(rows)


scores_csv = (
    OUTPUT_DIR
    / "AK1_ZN601_scores.csv"
)

zn_csv = (
    OUTPUT_DIR
    / "AK1_ZN601_Zn_distances.csv"
)

comparison_csv = (
    OUTPUT_DIR
    / "AK1_ZN601_vs_Znfree_families.csv"
)

rmsd_csv = (
    OUTPUT_DIR
    / "AK1_ZN601_pose_RMSD_matrix.csv"
)


save_csv(
    scores_csv,
    score_rows
)

save_csv(
    zn_csv,
    zn_rows
)

save_csv(
    comparison_csv,
    comparison_rows
)


with open(
    rmsd_csv,
    "w",
    newline="",
    encoding="utf-8"
) as f:

    writer = csv.writer(f)

    writer.writerow(
        ["Seed"] + SEEDS
    )

    for i, seed in enumerate(
        SEEDS
    ):

        writer.writerow(
            [seed]
            +
            [
                f"{x:.6f}"
                for x in rmsd_matrix[i]
            ]
        )


# ============================================================
# EXCEL WORKBOOK
# ============================================================

XLSX = (
    OUTPUT_DIR
    / "AK1_ZN601_10seed_comprehensive.xlsx"
)

workbook = xlsxwriter.Workbook(
    XLSX
)


title_fmt = workbook.add_format({
    "bold": True,
    "font_size": 15,
    "font_color": "white",
    "bg_color": "#1F4E78"
})

header_fmt = workbook.add_format({
    "bold": True,
    "font_color": "white",
    "bg_color": "#4472C4",
    "border": 1,
    "align": "center",
    "text_wrap": True
})

normal_fmt = workbook.add_format({
    "border": 1
})

number_fmt = workbook.add_format({
    "border": 1,
    "num_format": "0.000"
})

dominant_fmt = workbook.add_format({
    "border": 1,
    "bg_color": "#FFF2CC"
})


# ============================================================
# SHEET 1 — SUMMARY
# ============================================================

ws = workbook.add_worksheet(
    "Summary"
)


summary_rows = [

    [
        "ZN601 mean Vina score",
        statistics.mean(
            zn_score_values
        )
    ],

    [
        "ZN601 median Vina score",
        statistics.median(
            zn_score_values
        )
    ],

    [
        "ZN601 score SD",
        statistics.stdev(
            zn_score_values
        )
    ],

    [
        "ZN601 score range",
        max(zn_score_values)
        -
        min(zn_score_values)
    ],

    [
        "Zn-free mean Vina score",
        statistics.mean(
            original_score_values
        )
    ],

    [
        "Mean paired score shift",
        statistics.mean(
            score_differences
        )
    ],

    [
        "Dominant ZN601 cluster size",
        len(
            dominant_seeds
        )
    ],

    [
        "Dominant ZN601 cluster seeds",
        ", ".join(
            str(x)
            for x in dominant_seeds
        )
    ],

    [
        "Representative ZN601 seed",
        representative_seed
    ],

    [
        "Representative mean RMSD",
        representative_mean_rmsd
    ],

    [
        "Mean closest ZN601-polar distance",
        statistics.mean(
            row[
                "Closest polar distance to ZN601 (A)"
            ]
            for row in zn_rows
        )
    ]
]


ws.write(
    0,
    0,
    "Metric",
    header_fmt
)

ws.write(
    0,
    1,
    "Value",
    header_fmt
)


for row_i, values in enumerate(
    summary_rows,
    start=1
):

    ws.write(
        row_i,
        0,
        values[0],
        normal_fmt
    )

    if isinstance(
        values[1],
        float
    ):

        ws.write_number(
            row_i,
            1,
            values[1],
            number_fmt
        )

    else:

        ws.write(
            row_i,
            1,
            values[1],
            normal_fmt
        )


ws.set_column(
    "A:A",
    42
)

ws.set_column(
    "B:B",
    55
)


# ============================================================
# SHEET 2 — SCORES
# ============================================================

ws2 = workbook.add_worksheet(
    "Score_Comparison"
)

headers = list(
    score_rows[0].keys()
)


for col, h in enumerate(
    headers
):

    ws2.write(
        0,
        col,
        h,
        header_fmt
    )


for row_i, row in enumerate(
    score_rows,
    start=1
):

    for col, h in enumerate(
        headers
    ):

        value = row[h]

        if isinstance(
            value,
            float
        ):

            ws2.write_number(
                row_i,
                col,
                value,
                number_fmt
            )

        else:

            ws2.write(
                row_i,
                col,
                value,
                normal_fmt
            )


ws2.set_column(
    0,
    len(headers)-1,
    24
)


# ============================================================
# SHEET 3 — RMSD MATRIX
# ============================================================

ws3 = workbook.add_worksheet(
    "ZN601_RMSD_Matrix"
)


ws3.write(
    0,
    0,
    "Seed",
    header_fmt
)


for j, seed in enumerate(
    SEEDS,
    start=1
):

    ws3.write_number(
        0,
        j,
        seed,
        header_fmt
    )


for i, seed in enumerate(
    SEEDS,
    start=1
):

    ws3.write_number(
        i,
        0,
        seed,
        header_fmt
    )


    for j, value in enumerate(
        rmsd_matrix[
            i-1
        ],
        start=1
    ):

        ws3.write_number(
            i,
            j,
            float(value),
            number_fmt
        )


ws3.set_column(
    0,
    len(SEEDS),
    13
)


# ============================================================
# SHEET 4 — CLUSTERS
# ============================================================

ws4 = workbook.add_worksheet(
    "Clusters"
)

cluster_headers = [
    "Seed",
    "Cluster at 2 A",
    "Cluster at 3 A",
    "Representative dominant pose"
]


for col, h in enumerate(
    cluster_headers
):

    ws4.write(
        0,
        col,
        h,
        header_fmt
    )


for i, seed in enumerate(
    SEEDS,
    start=1
):

    is_rep = (
        "YES"
        if seed == representative_seed
        else ""
    )


    fmt = (
        dominant_fmt
        if clusters_2A[
            i-1
        ] == dominant_cluster_id
        else normal_fmt
    )


    ws4.write_number(
        i,
        0,
        seed,
        fmt
    )

    ws4.write_number(
        i,
        1,
        clusters_2A[
            i-1
        ],
        fmt
    )

    ws4.write_number(
        i,
        2,
        clusters_3A[
            i-1
        ],
        fmt
    )

    ws4.write(
        i,
        3,
        is_rep,
        fmt
    )


ws4.set_column(
    "A:D",
    25
)


# ============================================================
# SHEET 5 — Zn DISTANCES
# ============================================================

ws5 = workbook.add_worksheet(
    "ZN601_Distances"
)

headers = list(
    zn_rows[0].keys()
)


for col, h in enumerate(
    headers
):

    ws5.write(
        0,
        col,
        h,
        header_fmt
    )


for row_i, row in enumerate(
    zn_rows,
    start=1
):

    for col, h in enumerate(
        headers
    ):

        value = row[h]

        if isinstance(
            value,
            float
        ):

            ws5.write_number(
                row_i,
                col,
                value,
                number_fmt
            )

        else:

            ws5.write(
                row_i,
                col,
                value,
                normal_fmt
            )


ws5.set_column(
    0,
    len(headers)-1,
    30
)


# ============================================================
# SHEET 6 — FAMILY COMPARISON
# ============================================================

ws6 = workbook.add_worksheet(
    "Vs_ZnFree_Families"
)

headers = list(
    comparison_rows[0].keys()
)


for col, h in enumerate(
    headers
):

    ws6.write(
        0,
        col,
        h,
        header_fmt
    )


for row_i, row in enumerate(
    comparison_rows,
    start=1
):

    for col, h in enumerate(
        headers
    ):

        value = row[h]

        if isinstance(
            value,
            float
        ):

            ws6.write_number(
                row_i,
                col,
                value,
                number_fmt
            )

        else:

            ws6.write(
                row_i,
                col,
                value,
                normal_fmt
            )


ws6.set_column(
    0,
    len(headers)-1,
    31
)


# ============================================================
# SHEET 7 — METHOD
# ============================================================

ws7 = workbook.add_worksheet(
    "Method"
)

method_rows = [

    [
        "Receptor",
        "AK1 fixed chain A with crystallographic ZN601 retained."
    ],

    [
        "Scoring",
        "AutoDock Vina 1.2.7 standard Vina scoring."
    ],

    [
        "Seeds",
        ", ".join(
            str(s)
            for s in SEEDS
        )
    ],

    [
        "Exhaustiveness",
        "16"
    ],

    [
        "Modes",
        "10 per seed"
    ],

    [
        "RMSD",
        "Direct heavy-atom positional RMSD in the common AK1 coordinate frame; no ligand superposition."
    ],

    [
        "Clustering",
        "Complete-link hierarchical clustering of pairwise positional RMSD; evaluated at 2.0 A and 3.0 A thresholds."
    ],

    [
        "Representative",
        "Medoid-like pose: member of the largest 2.0 A cluster with the lowest mean RMSD to other members."
    ],

    [
        "Metal interpretation",
        "Standard-Vina Zn-retained analysis is a receptor-state sensitivity test and not a dedicated metallodocking/free-energy calculation."
    ]
]


for row_i, values in enumerate(
    method_rows
):

    ws7.write(
        row_i,
        0,
        values[0],
        header_fmt
    )

    ws7.write(
        row_i,
        1,
        values[1],
        normal_fmt
    )


ws7.set_column(
    "A:A",
    30
)

ws7.set_column(
    "B:B",
    110
)


workbook.close()


# ============================================================
# FINISHED
# ============================================================

print()
print("=" * 90)
print("FILES SAVED")
print("=" * 90)

print()
print(XLSX)

print()
print(scores_csv)
print(rmsd_csv)
print(comparison_csv)
print(zn_csv)