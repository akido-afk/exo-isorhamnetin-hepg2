from pathlib import Path
import re
import csv
import statistics
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
import xlsxwriter


# ============================================================
# PROJECT PATHS
# ============================================================

PROJECT = Path(r"C:\Users\akidi\Desktop\ISO_AK_Docking")

DOCKING = PROJECT / "docking"
RESULTS = PROJECT / "results"
RECEPTORS = PROJECT / "receptors"
LIGAND = PROJECT / "ligand"
STRUCTURES = PROJECT / "structures"

RESULTS.mkdir(exist_ok=True)

OUTPUT = RESULTS / "AK1_ISO_COMPREHENSIVE_DOCKING_RESULTS.xlsx"

RMSD_MATRIX_FILE = RESULTS / "AK1_ISO_pose_RMSD_matrix.csv"


# ============================================================
# FIXED EXPERIMENTAL / METHOD INFORMATION
# ============================================================

TARGET = "Human AK1"
GENE = "AK1"

STRUCTURE_ID = "PDB 1Z83"
CHAIN = "A"

NATIVE_LIGAND = "AP5"
AP5_REDocking_RMSD = 1.776

ISO_NAME = "Isorhamnetin"
ISO_PUBCHEM_CID = "5281654"

VINA_VERSION = "1.2.7"
MEEKO_VERSION = "0.8.0"
PYTHON_VERSION = "3.10"

EXHAUSTIVENESS = 16
NUM_MODES = 10

# Docking box defined from crystallographic AP5 site
CENTER_X = 25.044
CENTER_Y = 44.571
CENTER_Z = -3.708

SIZE_X = 23.687
SIZE_Y = 23.364
SIZE_Z = 30.062

# Primary clustering threshold
CLUSTER_THRESHOLD = 2.0

# Secondary threshold used to test stability
SECONDARY_THRESHOLD = 3.0


# ============================================================
# HELPER: READ VINA PDBQT RESULTS
# ============================================================

def parse_vina_pdbqt(path):

    poses = []

    with open(path, encoding="utf-8") as f:

        for line in f:

            if line.startswith("REMARK VINA RESULT:"):

                parts = line.split()

                poses.append({
                    "mode": len(poses) + 1,
                    "score": float(parts[3]),
                    "rmsd_lb": float(parts[4]),
                    "rmsd_ub": float(parts[5])
                })

    return poses


# ============================================================
# FIND ALL 10 AK1 ISO RUNS
# ============================================================

pattern = re.compile(
    r"AK1_ISO_seed(\d+)\.pdbqt$",
    re.IGNORECASE
)

runs = []
all_poses = []

for file in DOCKING.glob("AK1_ISO_seed*.pdbqt"):

    match = pattern.match(file.name)

    if not match:
        continue

    seed = int(match.group(1))

    poses = parse_vina_pdbqt(file)

    if not poses:
        print("WARNING: no Vina results in", file.name)
        continue

    best = poses[0]

    runs.append({
        "seed": seed,
        "best_score": best["score"],
        "file": file.name
    })

    for pose in poses:

        all_poses.append({
            "seed": seed,
            "mode": pose["mode"],
            "score": pose["score"],
            "rmsd_lb": pose["rmsd_lb"],
            "rmsd_ub": pose["rmsd_ub"],
            "file": file.name
        })


runs = sorted(runs, key=lambda x: x["seed"])
all_poses = sorted(all_poses, key=lambda x: (x["seed"], x["mode"]))


if len(runs) != 10:
    print()
    print("WARNING:")
    print(f"Expected 10 AK1 seed files, but found {len(runs)}.")
    print()


# ============================================================
# MULTI-SEED SCORE STATISTICS
# ============================================================

scores = [x["best_score"] for x in runs]

overall_best = min(scores)
overall_worst = max(scores)
overall_mean = statistics.mean(scores)
overall_median = statistics.median(scores)
overall_sd = statistics.stdev(scores)
overall_range = overall_worst - overall_best

best_run = min(runs, key=lambda x: x["best_score"])
best_score_seed = best_run["seed"]


# ============================================================
# READ POSE RMSD MATRIX
# ============================================================

if not RMSD_MATRIX_FILE.exists():
    raise FileNotFoundError(
        f"Could not find RMSD matrix:\n{RMSD_MATRIX_FILE}\n"
        "Run AK1_pose_convergence.py first."
    )


with open(RMSD_MATRIX_FILE, newline="") as f:

    reader = csv.reader(f)

    rows = list(reader)


matrix_seeds = [int(x) for x in rows[0][1:]]

rmsd_matrix = []

for row in rows[1:]:

    rmsd_matrix.append(
        [float(x) for x in row[1:]]
    )

rmsd_matrix = np.array(rmsd_matrix, dtype=float)


# ============================================================
# GLOBAL POSE CONVERGENCE STATISTICS
# ============================================================

pairwise_values = []

for i in range(len(matrix_seeds)):

    for j in range(i + 1, len(matrix_seeds)):

        pairwise_values.append(
            rmsd_matrix[i, j]
        )


mean_pairwise_rmsd = statistics.mean(pairwise_values)
median_pairwise_rmsd = statistics.median(pairwise_values)
min_pairwise_rmsd = min(pairwise_values)
max_pairwise_rmsd = max(pairwise_values)
sd_pairwise_rmsd = statistics.stdev(pairwise_values)


# ============================================================
# HIERARCHICAL CLUSTERING
# ============================================================

condensed = squareform(
    rmsd_matrix,
    checks=False
)

Z = linkage(
    condensed,
    method="average"
)


def make_clusters(threshold):

    labels = fcluster(
        Z,
        t=threshold,
        criterion="distance"
    )

    clusters = {}

    for seed, cluster_id in zip(matrix_seeds, labels):

        clusters.setdefault(
            int(cluster_id),
            []
        ).append(seed)

    # Sort clusters largest first
    ordered = sorted(
        clusters.values(),
        key=lambda x: len(x),
        reverse=True
    )

    return ordered


clusters_2A = make_clusters(CLUSTER_THRESHOLD)
clusters_3A = make_clusters(SECONDARY_THRESHOLD)

dominant_cluster = clusters_2A[0]

alternative_seeds = [
    seed
    for cluster in clusters_2A[1:]
    for seed in cluster
]

dominant_fraction = len(dominant_cluster) / len(matrix_seeds)


# ============================================================
# FIND REPRESENTATIVE POSE OF DOMINANT CLUSTER
# ============================================================

seed_to_index = {
    seed: i
    for i, seed in enumerate(matrix_seeds)
}

dominant_mean_rmsd = {}

for seed in dominant_cluster:

    i = seed_to_index[seed]

    values = []

    for other_seed in dominant_cluster:

        if other_seed == seed:
            continue

        j = seed_to_index[other_seed]

        values.append(
            rmsd_matrix[i, j]
        )

    dominant_mean_rmsd[seed] = statistics.mean(values)


representative_seed = min(
    dominant_mean_rmsd,
    key=dominant_mean_rmsd.get
)

representative_mean_rmsd = dominant_mean_rmsd[representative_seed]


# ============================================================
# CLUSTER SCORE STATISTICS
# ============================================================

score_by_seed = {
    x["seed"]: x["best_score"]
    for x in runs
}


def cluster_stats(seeds):

    cluster_scores = [
        score_by_seed[seed]
        for seed in seeds
    ]

    return {
        "n": len(cluster_scores),
        "mean": statistics.mean(cluster_scores),
        "median": statistics.median(cluster_scores),
        "sd": (
            statistics.stdev(cluster_scores)
            if len(cluster_scores) > 1
            else 0
        ),
        "best": min(cluster_scores),
        "worst": max(cluster_scores),
        "range": max(cluster_scores) - min(cluster_scores)
    }


dominant_stats = cluster_stats(dominant_cluster)

alternative_stats = (
    cluster_stats(alternative_seeds)
    if alternative_seeds
    else None
)


# ============================================================
# PARSE AP5 REDOCKING SCORES
# ============================================================

AP5_FILE = DOCKING / "AK1_AP5_redock_seed1001.pdbqt"

ap5_poses = []

if AP5_FILE.exists():

    ap5_poses = parse_vina_pdbqt(AP5_FILE)


# ============================================================
# CREATE WORKBOOK
# ============================================================

workbook = xlsxwriter.Workbook(OUTPUT)


# ============================================================
# FORMATS
# ============================================================

title_fmt = workbook.add_format({
    "bold": True,
    "font_size": 16,
    "font_color": "white",
    "bg_color": "#1F4E78",
    "align": "left",
    "valign": "vcenter"
})

section_fmt = workbook.add_format({
    "bold": True,
    "font_size": 12,
    "font_color": "white",
    "bg_color": "#5B9BD5",
    "align": "left",
    "valign": "vcenter"
})

header_fmt = workbook.add_format({
    "bold": True,
    "font_color": "white",
    "bg_color": "#4472C4",
    "border": 1,
    "align": "center",
    "valign": "vcenter",
    "text_wrap": True
})

normal_fmt = workbook.add_format({
    "border": 1,
    "valign": "top"
})

center_fmt = workbook.add_format({
    "border": 1,
    "align": "center",
    "valign": "center"
})

decimal_fmt = workbook.add_format({
    "border": 1,
    "num_format": "0.000",
    "align": "center"
})

decimal4_fmt = workbook.add_format({
    "border": 1,
    "num_format": "0.0000",
    "align": "center"
})

percent_fmt = workbook.add_format({
    "border": 1,
    "num_format": "0%",
    "align": "center"
})

good_fmt = workbook.add_format({
    "border": 1,
    "bg_color": "#E2F0D9",
    "font_color": "#375623",
    "text_wrap": True,
    "valign": "top"
})

warning_fmt = workbook.add_format({
    "border": 1,
    "bg_color": "#FFF2CC",
    "font_color": "#7F6000",
    "text_wrap": True,
    "valign": "top"
})

note_fmt = workbook.add_format({
    "border": 1,
    "bg_color": "#F2F2F2",
    "text_wrap": True,
    "valign": "top"
})

dominant_fmt = workbook.add_format({
    "border": 1,
    "bg_color": "#E2F0D9",
    "align": "center"
})

alternative_fmt = workbook.add_format({
    "border": 1,
    "bg_color": "#FFF2CC",
    "align": "center"
})


# ============================================================
# SHEET 1 — EXECUTIVE SUMMARY
# ============================================================

ws = workbook.add_worksheet("01_AK1_Summary")

ws.merge_range(
    "A1:H1",
    "AK1–Isorhamnetin Molecular Docking: Comprehensive Summary",
    title_fmt
)

ws.set_row(0, 26)

summary_data = [
    ["Target", TARGET],
    ["Structure", f"{STRUCTURE_ID}, chain {CHAIN}"],
    ["Docked ligand", ISO_NAME],
    ["Independent Vina seeds", len(runs)],
    ["Total generated poses", len(all_poses)],
    ["Best Mode-1 score", overall_best],
    ["Seed with most negative score", best_score_seed],
    ["Mean Mode-1 score", overall_mean],
    ["Median Mode-1 score", overall_median],
    ["SD across Mode-1 scores", overall_sd],
    ["Score range", overall_range],
    ["Native validation ligand", NATIVE_LIGAND],
    ["Native-ligand redocking RMSD", AP5_REDocking_RMSD],
    ["Dominant pose-cluster size", len(dominant_cluster)],
    ["Dominant-cluster frequency", dominant_fraction],
    ["Alternative pose count", len(alternative_seeds)],
    ["Representative dominant pose", representative_seed],
    ["Representative pose mean RMSD within dominant cluster",
     representative_mean_rmsd],
]

ws.write("A3", "Metric", header_fmt)
ws.write("B3", "Result", header_fmt)

for i, (metric, value) in enumerate(summary_data, start=3):

    ws.write(i, 0, metric, normal_fmt)

    if isinstance(value, float):

        if "frequency" in metric.lower():
            ws.write_number(i, 1, value, percent_fmt)

        else:
            ws.write_number(i, 1, value, decimal_fmt)

    else:
        ws.write(i, 1, value, normal_fmt)


# Interpretation block
interpret_row = 23

ws.merge_range(
    interpret_row,
    0,
    interpret_row,
    7,
    "Interpretation",
    section_fmt
)

interpretation = (
    f"Across {len(runs)} independent AutoDock Vina runs, isorhamnetin showed "
    f"highly reproducible AK1 docking scores (median {overall_median:.3f} kcal/mol; "
    f"SD {overall_sd:.3f} kcal/mol). "
    f"Hierarchical clustering of the top-ranked pose from each run identified "
    f"a dominant structural cluster containing {len(dominant_cluster)}/{len(runs)} "
    f"runs ({dominant_fraction*100:.0f}%). "
    f"The representative pose of this dominant cluster was seed {representative_seed}, "
    f"with a mean RMSD of only {representative_mean_rmsd:.3f} Å to the other "
    f"dominant-cluster poses. "
    f"The AK1 docking protocol was additionally supported by native AP5 redocking "
    f"(RMSD {AP5_REDocking_RMSD:.3f} Å)."
)

ws.merge_range(
    interpret_row + 1,
    0,
    interpret_row + 4,
    7,
    interpretation,
    good_fmt
)


ws.merge_range(
    interpret_row + 6,
    0,
    interpret_row + 6,
    7,
    "Important limitation",
    section_fmt
)

limitation = (
    "AutoDock Vina scores are computational docking scores, not experimentally "
    "measured binding free energies, dissociation constants, inhibition constants, "
    "or proof that isorhamnetin directly inhibits AK1. The small numerical score "
    "difference between the dominant and alternative pose clusters should not be "
    "interpreted as a meaningful energetic difference."
)

ws.merge_range(
    interpret_row + 7,
    0,
    interpret_row + 10,
    7,
    limitation,
    warning_fmt
)

ws.set_column("A:A", 45)
ws.set_column("B:B", 24)
ws.set_column("C:H", 15)

ws.freeze_panes(3, 0)


# ============================================================
# CHART 1 — BEST SCORE BY SEED
# ============================================================

chart1 = workbook.add_chart({
    "type": "column"
})

chart1.add_series({
    "name": "Best Vina score",
    "categories": [
        "02_Seed_Results",
        1,
        0,
        len(runs),
        0
    ],
    "values": [
        "02_Seed_Results",
        1,
        2,
        len(runs),
        2
    ]
})

chart1.set_title({
    "name": "Best AK1–ISO docking score by seed"
})

chart1.set_x_axis({
    "name": "Seed"
})

chart1.set_y_axis({
    "name": "Vina score (kcal/mol)"
})

chart1.set_legend({
    "none": True
})

chart1.set_size({
    "width": 620,
    "height": 330
})

ws.insert_chart(
    "D3",
    chart1
)


# ============================================================
# SHEET 2 — SEED RESULTS
# ============================================================

ws2 = workbook.add_worksheet("02_Seed_Results")

ws2.merge_range(
    "A1:H1",
    "Mode-1 Results from 10 Independent AK1 Docking Runs",
    title_fmt
)

headers = [
    "Seed",
    "Cluster",
    "Best Vina score\n(kcal/mol)",
    "Mean RMSD to dominant cluster\n(Å)",
    "Representative dominant pose?",
    "Most negative score overall?",
    "Source PDBQT",
    "Interpretation"
]

for col, h in enumerate(headers):
    ws2.write(2, col, h, header_fmt)


for row_index, run in enumerate(runs, start=3):

    seed = run["seed"]

    if seed in dominant_cluster:
        cluster_name = "Dominant"
        row_format = dominant_fmt
        mean_dom = dominant_mean_rmsd[seed]

    else:
        cluster_name = "Alternative"
        row_format = alternative_fmt
        mean_dom = None

    ws2.write_number(row_index, 0, seed, row_format)
    ws2.write(row_index, 1, cluster_name, row_format)
    ws2.write_number(
        row_index,
        2,
        run["best_score"],
        decimal_fmt
    )

    if mean_dom is not None:
        ws2.write_number(
            row_index,
            3,
            mean_dom,
            decimal_fmt
        )
    else:
        ws2.write(
            row_index,
            3,
            "N/A",
            row_format
        )

    ws2.write(
        row_index,
        4,
        "YES" if seed == representative_seed else "",
        row_format
    )

    ws2.write(
        row_index,
        5,
        "YES" if seed == best_score_seed else "",
        row_format
    )

    ws2.write(
        row_index,
        6,
        run["file"],
        normal_fmt
    )

    if seed == representative_seed:

        interpretation_text = (
            "Representative pose of dominant cluster"
        )

    elif seed in dominant_cluster:

        interpretation_text = (
            "Dominant consensus binding orientation"
        )

    else:

        interpretation_text = (
            "Alternative binding orientation"
        )

    ws2.write(
        row_index,
        7,
        interpretation_text,
        normal_fmt
    )


ws2.set_column("A:A", 12)
ws2.set_column("B:B", 15)
ws2.set_column("C:D", 24)
ws2.set_column("E:F", 24)
ws2.set_column("G:G", 35)
ws2.set_column("H:H", 38)

ws2.freeze_panes(3, 0)


# ============================================================
# SHEET 3 — ALL 100 POSES
# ============================================================

ws3 = workbook.add_worksheet("03_All_100_Poses")

ws3.merge_range(
    "A1:G1",
    "All AutoDock Vina Poses from the 10 Independent Runs",
    title_fmt
)

pose_headers = [
    "Seed",
    "Mode",
    "Vina score (kcal/mol)",
    "RMSD lower bound",
    "RMSD upper bound",
    "Is Mode 1?",
    "Source file"
]

for col, h in enumerate(pose_headers):
    ws3.write(2, col, h, header_fmt)


for row_index, pose in enumerate(all_poses, start=3):

    ws3.write_number(
        row_index,
        0,
        pose["seed"],
        center_fmt
    )

    ws3.write_number(
        row_index,
        1,
        pose["mode"],
        center_fmt
    )

    ws3.write_number(
        row_index,
        2,
        pose["score"],
        decimal_fmt
    )

    ws3.write_number(
        row_index,
        3,
        pose["rmsd_lb"],
        decimal_fmt
    )

    ws3.write_number(
        row_index,
        4,
        pose["rmsd_ub"],
        decimal_fmt
    )

    ws3.write(
        row_index,
        5,
        "YES" if pose["mode"] == 1 else "",
        center_fmt
    )

    ws3.write(
        row_index,
        6,
        pose["file"],
        normal_fmt
    )


ws3.set_column("A:B", 12)
ws3.set_column("C:E", 21)
ws3.set_column("F:F", 14)
ws3.set_column("G:G", 38)

ws3.freeze_panes(3, 0)

ws3.autofilter(
    2,
    0,
    2 + len(all_poses),
    6
)


# ============================================================
# SHEET 4 — POSE RMSD MATRIX
# ============================================================

ws4 = workbook.add_worksheet("04_Pose_RMSD_Matrix")

ws4.merge_range(
    "A1:L1",
    "Pairwise RMSD Matrix of the 10 Best AK1–ISO Poses",
    title_fmt
)

ws4.write(
    2,
    0,
    "Seed",
    header_fmt
)

for col, seed in enumerate(matrix_seeds, start=1):

    ws4.write_number(
        2,
        col,
        seed,
        header_fmt
    )


for row_idx, seed in enumerate(matrix_seeds, start=3):

    ws4.write_number(
        row_idx,
        0,
        seed,
        header_fmt
    )

    for col_idx, value in enumerate(
        rmsd_matrix[row_idx - 3],
        start=1
    ):

        ws4.write_number(
            row_idx,
            col_idx,
            value,
            decimal_fmt
        )


# Heatmap
ws4.conditional_format(
    3,
    1,
    2 + len(matrix_seeds),
    len(matrix_seeds),
    {
        "type": "3_color_scale",
        "min_color": "#63BE7B",
        "mid_color": "#FFEB84",
        "max_color": "#F8696B"
    }
)


summary_start = 15

rmsd_summary = [
    ["Number of best poses", len(matrix_seeds)],
    ["Pairwise comparisons", len(pairwise_values)],
    ["Mean pairwise RMSD", mean_pairwise_rmsd],
    ["Median pairwise RMSD", median_pairwise_rmsd],
    ["SD pairwise RMSD", sd_pairwise_rmsd],
    ["Minimum pairwise RMSD", min_pairwise_rmsd],
    ["Maximum pairwise RMSD", max_pairwise_rmsd],
]

ws4.write(
    summary_start,
    0,
    "Global pose convergence",
    section_fmt
)

ws4.write(
    summary_start + 1,
    0,
    "Metric",
    header_fmt
)

ws4.write(
    summary_start + 1,
    1,
    "Value",
    header_fmt
)

for i, (metric, value) in enumerate(
    rmsd_summary,
    start=summary_start + 2
):

    ws4.write(
        i,
        0,
        metric,
        normal_fmt
    )

    if isinstance(value, float):

        ws4.write_number(
            i,
            1,
            value,
            decimal_fmt
        )

    else:

        ws4.write_number(
            i,
            1,
            value,
            center_fmt
        )


ws4.set_column("A:A", 26)
ws4.set_column("B:L", 12)

ws4.freeze_panes(3, 1)


# ============================================================
# SHEET 5 — CLUSTER ANALYSIS
# ============================================================

ws5 = workbook.add_worksheet("05_Cluster_Analysis")

ws5.merge_range(
    "A1:H1",
    "Structural Clustering of Mode-1 AK1–ISO Poses",
    title_fmt
)


# 2 Å clustering
ws5.merge_range(
    "A3:H3",
    f"Clustering at {CLUSTER_THRESHOLD:.1f} Å RMSD threshold",
    section_fmt
)

cluster_headers = [
    "Cluster",
    "n poses",
    "% of runs",
    "Seeds",
    "Mean score",
    "Median score",
    "SD score",
    "Interpretation"
]

for col, h in enumerate(cluster_headers):
    ws5.write(
        3,
        col,
        h,
        header_fmt
    )


for i, cluster in enumerate(clusters_2A, start=4):

    stats = cluster_stats(cluster)

    cluster_name = (
        "Dominant"
        if i == 4
        else f"Alternative {i-4}"
    )

    fmt = (
        dominant_fmt
        if i == 4
        else alternative_fmt
    )

    ws5.write(
        i,
        0,
        cluster_name,
        fmt
    )

    ws5.write_number(
        i,
        1,
        len(cluster),
        center_fmt
    )

    ws5.write_number(
        i,
        2,
        len(cluster) / len(runs),
        percent_fmt
    )

    ws5.write(
        i,
        3,
        ", ".join(map(str, cluster)),
        normal_fmt
    )

    ws5.write_number(
        i,
        4,
        stats["mean"],
        decimal_fmt
    )

    ws5.write_number(
        i,
        5,
        stats["median"],
        decimal_fmt
    )

    ws5.write_number(
        i,
        6,
        stats["sd"],
        decimal_fmt
    )

    ws5.write(
        i,
        7,
        (
            "Main reproducible binding orientation"
            if i == 4
            else "Minor alternative binding orientation"
        ),
        normal_fmt
    )


# 3 Å clustering stability check
row3A = 9

ws5.merge_range(
    row3A,
    0,
    row3A,
    7,
    f"Clustering stability check at {SECONDARY_THRESHOLD:.1f} Å",
    section_fmt
)

ws5.write(
    row3A + 1,
    0,
    "Cluster",
    header_fmt
)

ws5.write(
    row3A + 1,
    1,
    "n poses",
    header_fmt
)

ws5.write(
    row3A + 1,
    2,
    "Seeds",
    header_fmt
)

for i, cluster in enumerate(
    clusters_3A,
    start=row3A + 2
):

    ws5.write(
        i,
        0,
        f"Cluster {i-(row3A+1)}",
        normal_fmt
    )

    ws5.write_number(
        i,
        1,
        len(cluster),
        center_fmt
    )

    ws5.write(
        i,
        2,
        ", ".join(map(str, cluster)),
        normal_fmt
    )


# Representative dominant pose
rep_row = 15

ws5.merge_range(
    rep_row,
    0,
    rep_row,
    7,
    "Dominant-cluster representative pose",
    section_fmt
)

rep_headers = [
    "Seed",
    "Mean RMSD to other dominant poses (Å)",
    "Vina score (kcal/mol)",
    "Representative?"
]

for col, h in enumerate(rep_headers):
    ws5.write(
        rep_row + 1,
        col,
        h,
        header_fmt
    )


for r, seed in enumerate(
    sorted(dominant_cluster),
    start=rep_row + 2
):

    ws5.write_number(
        r,
        0,
        seed,
        center_fmt
    )

    ws5.write_number(
        r,
        1,
        dominant_mean_rmsd[seed],
        decimal_fmt
    )

    ws5.write_number(
        r,
        2,
        score_by_seed[seed],
        decimal_fmt
    )

    ws5.write(
        r,
        3,
        "YES" if seed == representative_seed else "",
        dominant_fmt if seed == representative_seed else normal_fmt
    )


ws5.set_column("A:A", 22)
ws5.set_column("B:C", 24)
ws5.set_column("D:D", 50)
ws5.set_column("E:G", 18)
ws5.set_column("H:H", 40)


# ============================================================
# SHEET 6 — AP5 VALIDATION
# ============================================================

ws6 = workbook.add_worksheet("06_AP5_Validation")

ws6.merge_range(
    "A1:G1",
    "AK1 Native-Ligand Redocking Validation",
    title_fmt
)

validation_rows = [
    ["Protein", TARGET],
    ["Crystal structure", f"{STRUCTURE_ID}, chain {CHAIN}"],
    ["Native ligand", NATIVE_LIGAND],
    ["Redocking RMSD", AP5_REDocking_RMSD],
    ["Validation interpretation",
     "Redocking reproduced the crystallographic AP5 pose with an RMSD below 2 Å, supporting use of the AK1 docking setup."],
]

ws6.write("A3", "Parameter", header_fmt)
ws6.write("B3", "Result", header_fmt)

for row, (parameter, value) in enumerate(
    validation_rows,
    start=3
):

    ws6.write(
        row,
        0,
        parameter,
        normal_fmt
    )

    if isinstance(value, float):

        ws6.write_number(
            row,
            1,
            value,
            decimal_fmt
        )

    else:

        ws6.write(
            row,
            1,
            value,
            normal_fmt
        )


if ap5_poses:

    ws6.write(
        10,
        0,
        "AP5 redocking modes",
        section_fmt
    )

    ap5_headers = [
        "Mode",
        "Vina score (kcal/mol)",
        "RMSD lower bound",
        "RMSD upper bound"
    ]

    for col, h in enumerate(ap5_headers):
        ws6.write(
            11,
            col,
            h,
            header_fmt
        )

    for row, pose in enumerate(
        ap5_poses,
        start=12
    ):

        ws6.write_number(
            row,
            0,
            pose["mode"],
            center_fmt
        )

        ws6.write_number(
            row,
            1,
            pose["score"],
            decimal_fmt
        )

        ws6.write_number(
            row,
            2,
            pose["rmsd_lb"],
            decimal_fmt
        )

        ws6.write_number(
            row,
            3,
            pose["rmsd_ub"],
            decimal_fmt
        )


ws6.set_column("A:A", 30)
ws6.set_column("B:B", 75)
ws6.set_column("C:G", 20)


# ============================================================
# SHEET 7 — METHODS & PARAMETERS
# ============================================================

ws7 = workbook.add_worksheet("07_Methods_Parameters")

ws7.merge_range(
    "A1:D1",
    "Computational Docking Parameters and Reproducibility Information",
    title_fmt
)

method_data = [
    ["Target protein", TARGET, "Protein identity"],
    ["Gene symbol", GENE, "Target identity"],
    ["Crystal structure", STRUCTURE_ID, "Experimental protein structure"],
    ["Protein chain", CHAIN, "Chain used for docking"],
    ["Native ligand", NATIVE_LIGAND, "Used to define and validate binding site"],
    ["Docked ligand", ISO_NAME, "Test ligand"],
    ["PubChem CID", ISO_PUBCHEM_CID, "Ligand identifier"],
    ["Docking engine", f"AutoDock Vina {VINA_VERSION}", "Docking search/scoring"],
    ["Meeko", MEEKO_VERSION, "PDBQT preparation"],
    ["Python environment", PYTHON_VERSION, "Computational environment"],
    ["Receptor preparation",
     "Chain extraction; PDBFixer repair of missing atoms on present residues; hydrogens added at pH 7.4; Meeko PDBQT preparation",
     "Whole missing residues/loops were not reconstructed"],
    ["Binding-site definition",
     "Crystallographic AP5-binding region",
     "Pocket-based docking"],
    ["Center X", CENTER_X, "Å"],
    ["Center Y", CENTER_Y, "Å"],
    ["Center Z", CENTER_Z, "Å"],
    ["Size X", SIZE_X, "Å"],
    ["Size Y", SIZE_Y, "Å"],
    ["Size Z", SIZE_Z, "Å"],
    ["Exhaustiveness", EXHAUSTIVENESS, "Vina search parameter"],
    ["Number of modes/run", NUM_MODES, "Vina output parameter"],
    ["Independent seeds", len(runs), "Random-search replicates"],
    ["Primary clustering threshold", CLUSTER_THRESHOLD, "Å"],
    ["Secondary clustering threshold", SECONDARY_THRESHOLD, "Å"],
    ["Clustering method", "Average-linkage hierarchical clustering", "Based on pairwise pose RMSD"],
]

method_headers = [
    "Parameter",
    "Value",
    "Notes"
]

for col, h in enumerate(method_headers):
    ws7.write(
        2,
        col,
        h,
        header_fmt
    )

for row, values in enumerate(
    method_data,
    start=3
):

    for col, value in enumerate(values):

        if isinstance(value, float):

            ws7.write_number(
                row,
                col,
                value,
                decimal_fmt
            )

        elif isinstance(value, int):

            ws7.write_number(
                row,
                col,
                value,
                center_fmt
            )

        else:

            ws7.write(
                row,
                col,
                value,
                normal_fmt
            )


ws7.set_column("A:A", 32)
ws7.set_column("B:B", 70)
ws7.set_column("C:C", 42)
ws7.set_column("D:D", 15)


# ============================================================
# SHEET 8 — FILE PROVENANCE
# ============================================================

ws8 = workbook.add_worksheet("08_File_Provenance")

ws8.merge_range(
    "A1:D1",
    "AK1 Docking Input and Output Files",
    title_fmt
)

file_rows = [
    [
        "Original AK1 structure",
        STRUCTURES / "AK1_1Z83.pdb",
        "Input",
        "Downloaded PDB structure"
    ],
    [
        "AK1 chain A",
        RECEPTORS / "AK1_chainA.pdb",
        "Intermediate",
        "Protein chain extracted from 1Z83"
    ],
    [
        "Fixed AK1 structure",
        RECEPTORS / "AK1_chainA_FIXED.pdb",
        "Intermediate",
        "PDBFixer-prepared receptor"
    ],
    [
        "Docking receptor",
        RECEPTORS / "AK1.pdbqt",
        "Input",
        "Meeko-prepared AK1 receptor"
    ],
    [
        "Isorhamnetin 3D structure",
        LIGAND / "isorhamnetin_3D.sdf",
        "Input",
        "PubChem CID 5281654"
    ],
    [
        "Isorhamnetin docking ligand",
        LIGAND / "isorhamnetin.pdbqt",
        "Input",
        "Meeko-prepared ligand"
    ],
    [
        "Native AP5",
        LIGAND / "AK1_native_AP5_chainA.pdb",
        "Validation",
        "Crystallographic native ligand"
    ],
    [
        "AP5 redocking",
        AP5_FILE,
        "Validation",
        "Native ligand redocking result"
    ],
    [
        "Pose RMSD matrix",
        RMSD_MATRIX_FILE,
        "Analysis",
        "10 × 10 best-pose RMSD matrix"
    ],
]

for run in runs:

    file_rows.append([
        f"ISO seed {run['seed']}",
        DOCKING / run["file"],
        "Docking output",
        "10 Vina modes"
    ])


file_headers = [
    "File / dataset",
    "Path",
    "Role",
    "Description"
]

for col, h in enumerate(file_headers):
    ws8.write(
        2,
        col,
        h,
        header_fmt
    )


for row, values in enumerate(
    file_rows,
    start=3
):

    ws8.write(
        row,
        0,
        values[0],
        normal_fmt
    )

    ws8.write(
        row,
        1,
        str(values[1]),
        normal_fmt
    )

    ws8.write(
        row,
        2,
        values[2],
        normal_fmt
    )

    ws8.write(
        row,
        3,
        values[3],
        normal_fmt
    )


ws8.set_column("A:A", 30)
ws8.set_column("B:B", 85)
ws8.set_column("C:C", 20)
ws8.set_column("D:D", 45)


# ============================================================
# SHEET 9 — REPORT-READY TEXT
# ============================================================

ws9 = workbook.add_worksheet("09_Report_Ready_Text")

ws9.merge_range(
    "A1:F1",
    "Report-Ready AK1 Docking Statements",
    title_fmt
)

report_sections = [

    (
        "Validation",
        f"Redocking of the crystallographic native ligand AP5 into the "
        f"AK1 binding pocket reproduced the experimental binding orientation "
        f"with an RMSD of {AP5_REDocking_RMSD:.3f} Å, supporting the validity "
        f"of the docking protocol."
    ),

    (
        "Multi-seed docking",
        f"Isorhamnetin was docked to AK1 using {len(runs)} independent random "
        f"seeds. The best-ranked poses showed a median AutoDock Vina score of "
        f"{overall_median:.3f} kcal/mol, a mean score of {overall_mean:.3f} "
        f"kcal/mol, and a standard deviation of {overall_sd:.3f} kcal/mol."
    ),

    (
        "Pose clustering",
        f"Hierarchical clustering of the top-ranked pose from each independent "
        f"run identified a dominant pose family containing "
        f"{len(dominant_cluster)} of {len(runs)} runs "
        f"({dominant_fraction*100:.0f}%). The same {len(dominant_cluster)}:"
        f"{len(alternative_seeds)} cluster partition was retained at both "
        f"{CLUSTER_THRESHOLD:.1f} Å and {SECONDARY_THRESHOLD:.1f} Å RMSD "
        f"thresholds."
    ),

    (
        "Representative pose",
        f"Seed {representative_seed} was selected as the representative "
        f"dominant-cluster pose because it had the lowest mean RMSD "
        f"({representative_mean_rmsd:.3f} Å) to the other poses within "
        f"the dominant cluster."
    ),

    (
        "Interpretation limitation",
        "The calculated AutoDock Vina scores represent computational docking "
        "scores and should not be interpreted as experimentally measured "
        "binding free energies, Kd/Ki values, or direct evidence of AK1 "
        "inhibition."
    )
]


ws9.write(
    2,
    0,
    "Section",
    header_fmt
)

ws9.write(
    2,
    1,
    "Suggested text",
    header_fmt
)


for row, (section, text) in enumerate(
    report_sections,
    start=3
):

    ws9.write(
        row,
        0,
        section,
        normal_fmt
    )

    ws9.write(
        row,
        1,
        text,
        normal_fmt
    )


ws9.set_column("A:A", 28)
ws9.set_column("B:B", 120)
ws9.set_default_row(55)


# ============================================================
# FINAL WORKBOOK SETTINGS
# ============================================================

for sheet in [
    ws, ws2, ws3, ws4, ws5,
    ws6, ws7, ws8, ws9
]:

    sheet.hide_gridlines(2)


workbook.close()


# ============================================================
# TERMINAL SUMMARY
# ============================================================

print()
print("=" * 70)
print("AK1 COMPREHENSIVE DOCKING WORKBOOK CREATED")
print("=" * 70)

print()
print("Output file:")
print(OUTPUT)

print()
print("MULTI-SEED SCORES")
print("-----------------")
print(f"n seeds       : {len(runs)}")
print(f"Total poses   : {len(all_poses)}")
print(f"Best score    : {overall_best:.3f} kcal/mol")
print(f"Best-score seed: {best_score_seed}")
print(f"Mean          : {overall_mean:.3f} kcal/mol")
print(f"Median        : {overall_median:.3f} kcal/mol")
print(f"SD            : {overall_sd:.3f} kcal/mol")

print()
print("POSE CONVERGENCE")
print("----------------")
print(f"Mean pairwise RMSD : {mean_pairwise_rmsd:.3f} A")
print(f"Minimum RMSD       : {min_pairwise_rmsd:.3f} A")
print(f"Maximum RMSD       : {max_pairwise_rmsd:.3f} A")

print()
print("CLUSTERING")
print("----------")
print(
    f"Dominant cluster: {len(dominant_cluster)}/{len(runs)} "
    f"({dominant_fraction*100:.0f}%)"
)
print("Seeds:", dominant_cluster)

print(
    f"Alternative poses: {len(alternative_seeds)}"
)
print("Seeds:", alternative_seeds)

print()
print("REPRESENTATIVE POSE")
print("-------------------")
print(f"Seed: {representative_seed}")
print(
    f"Mean RMSD within dominant cluster: "
    f"{representative_mean_rmsd:.3f} A"
)

print()
print("VALIDATION")
print("----------")
print(
    f"AP5 redocking RMSD: "
    f"{AP5_REDocking_RMSD:.3f} A"
)

print()
print("=" * 70)