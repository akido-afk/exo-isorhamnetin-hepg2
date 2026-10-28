from pathlib import Path
import re
import statistics
import xlsxwriter

# ============================================================
# PROJECT PATHS
# ============================================================

PROJECT = Path(r"C:\Users\akidi\Desktop\ISO_AK_Docking")

DOCKING = PROJECT / "docking"
RESULTS = PROJECT / "results"

RESULTS.mkdir(exist_ok=True)

OUTPUT = RESULTS / "ISO_AK_Docking_Master_Results.xlsx"


# ============================================================
# CURRENT AK1 VALIDATION INFORMATION
# ============================================================

validation = {
    "AK1": {
        "structure": "PDB 1Z83 chain A",
        "validation_ligand": "AP5",
        "redocking_rmsd_A": 1.776,
        "box_center_x": 25.044,
        "box_center_y": 44.571,
        "box_center_z": -3.708,
        "box_size_x": 23.687,
        "box_size_y": 23.364,
        "box_size_z": 30.062,
    }
}


# ============================================================
# READ VINA PDBQT FILES
# ============================================================

def parse_vina_pdbqt(path):

    results = []

    with open(path, encoding="utf-8") as f:

        for line in f:

            if line.startswith("REMARK VINA RESULT:"):

                parts = line.split()

                score = float(parts[3])
                rmsd_lb = float(parts[4])
                rmsd_ub = float(parts[5])

                results.append({
                    "mode": len(results) + 1,
                    "score": score,
                    "rmsd_lb": rmsd_lb,
                    "rmsd_ub": rmsd_ub
                })

    return results


# ============================================================
# FIND ALL ISO MULTI-SEED RESULTS AUTOMATICALLY
# Example:
# AK1_ISO_seed1001.pdbqt
# ============================================================

pattern = re.compile(
    r"(.+)_ISO_seed(\d+)\.pdbqt$",
    re.IGNORECASE
)

all_runs = []
all_modes = []

for file in sorted(DOCKING.glob("AK1_ISO_seed*.pdbqt")):

    match = pattern.match(file.name)

    if not match:
        continue

    target = match.group(1)
    seed = int(match.group(2))

    poses = parse_vina_pdbqt(file)

    if not poses:
        print("WARNING: no Vina result found in", file.name)
        continue

    best = poses[0]

    all_runs.append({
        "target": target,
        "seed": seed,
        "best_score": best["score"],
        "file": file.name
    })

    for pose in poses:

        all_modes.append({
            "target": target,
            "seed": seed,
            "mode": pose["mode"],
            "score": pose["score"],
            "rmsd_lb": pose["rmsd_lb"],
            "rmsd_ub": pose["rmsd_ub"],
            "file": file.name
        })


# ============================================================
# SUMMARY BY TARGET
# ============================================================

targets = sorted(set(x["target"] for x in all_runs))

summary = []

for target in targets:

    scores = [
        x["best_score"]
        for x in all_runs
        if x["target"] == target
    ]

    info = validation.get(target, {})

    summary.append({
        "target": target,
        "structure": info.get("structure", ""),
        "n_runs": len(scores),
        "best_score": min(scores),
        "mean_score": statistics.mean(scores),
        "median_score": statistics.median(scores),
        "sd_score": statistics.stdev(scores) if len(scores) > 1 else None,
        "min_score": min(scores),
        "max_score": max(scores),
        "score_range": max(scores) - min(scores),
        "validation_ligand": info.get("validation_ligand", ""),
        "redocking_rmsd": info.get("redocking_rmsd_A", None)
    })


# More negative median score = stronger predicted docking
summary.sort(key=lambda x: x["median_score"])

for rank, row in enumerate(summary, start=1):
    row["rank"] = rank


# ============================================================
# CREATE EXCEL WORKBOOK
# ============================================================

workbook = xlsxwriter.Workbook(OUTPUT)


# ----------------------------
# FORMATS
# ----------------------------

title = workbook.add_format({
    "bold": True,
    "font_size": 15,
    "font_color": "white",
    "bg_color": "#1F4E78",
    "align": "left",
    "valign": "vcenter"
})

header = workbook.add_format({
    "bold": True,
    "font_color": "white",
    "bg_color": "#4472C4",
    "border": 1,
    "align": "center",
    "valign": "vcenter",
    "text_wrap": True
})

normal = workbook.add_format({
    "border": 1
})

decimal = workbook.add_format({
    "border": 1,
    "num_format": "0.000"
})

warning = workbook.add_format({
    "bg_color": "#FFF2CC",
    "border": 1,
    "text_wrap": True,
    "valign": "top"
})


# ============================================================
# SHEET 1 — SUMMARY
# ============================================================

ws = workbook.add_worksheet("Summary")

ws.merge_range(
    "A1:L1",
    "ISO–AK Family Docking: Automatically Parsed Results",
    title
)

headers = [
    "Rank",
    "Target / domain",
    "Structure",
    "n independent seeds",
    "Best score",
    "Mean best score",
    "Median best score",
    "SD",
    "Worst best score",
    "Range",
    "Validation ligand",
    "Redocking RMSD (Å)"
]

for col, h in enumerate(headers):
    ws.write(2, col, h, header)

for row_i, x in enumerate(summary, start=3):

    values = [
        x["rank"],
        x["target"],
        x["structure"],
        x["n_runs"],
        x["best_score"],
        x["mean_score"],
        x["median_score"],
        x["sd_score"],
        x["max_score"],
        x["score_range"],
        x["validation_ligand"],
        x["redocking_rmsd"]
    ]

    for col, value in enumerate(values):

        if col in [4, 5, 6, 7, 8, 9, 11] and value is not None:
            ws.write_number(row_i, col, value, decimal)
        else:
            ws.write(row_i, col, value, normal)


note_row = len(summary) + 5

ws.write(note_row, 0, "Interpretation", header)

ws.merge_range(
    note_row,
    1,
    note_row + 2,
    11,
    "More negative Vina scores indicate more favorable predicted docking "
    "within the same standardized protocol. Final target ranking should use "
    "the distribution of best scores across independent seeds, especially "
    "the median and variability, together with pose convergence and validation. "
    "Docking scores are not experimental binding affinities.",
    warning
)

ws.set_column("A:A", 10)
ws.set_column("B:B", 20)
ws.set_column("C:C", 26)
ws.set_column("D:D", 20)
ws.set_column("E:J", 18)
ws.set_column("K:L", 20)

ws.freeze_panes(3, 0)


# ============================================================
# SHEET 2 — SEED-BY-SEED RESULTS
# ============================================================

ws2 = workbook.add_worksheet("ISO_Runs")

headers2 = [
    "Target",
    "Seed",
    "Best Vina score (kcal/mol)",
    "Source file"
]

for col, h in enumerate(headers2):
    ws2.write(0, col, h, header)

sorted_runs = sorted(
    all_runs,
    key=lambda x: (x["target"], x["seed"])
)

for row_i, x in enumerate(sorted_runs, start=1):

    ws2.write(row_i, 0, x["target"], normal)
    ws2.write_number(row_i, 1, x["seed"], normal)
    ws2.write_number(row_i, 2, x["best_score"], decimal)
    ws2.write(row_i, 3, x["file"], normal)

ws2.set_column("A:A", 20)
ws2.set_column("B:B", 14)
ws2.set_column("C:C", 26)
ws2.set_column("D:D", 42)

ws2.freeze_panes(1, 0)


# ============================================================
# SHEET 3 — ALL 100 POSES
# ============================================================

ws3 = workbook.add_worksheet("All_Poses")

headers3 = [
    "Target",
    "Seed",
    "Mode",
    "Vina score (kcal/mol)",
    "RMSD lower bound",
    "RMSD upper bound",
    "Source file"
]

for col, h in enumerate(headers3):
    ws3.write(0, col, h, header)

sorted_modes = sorted(
    all_modes,
    key=lambda x: (
        x["target"],
        x["seed"],
        x["mode"]
    )
)

for row_i, x in enumerate(sorted_modes, start=1):

    ws3.write(row_i, 0, x["target"], normal)
    ws3.write_number(row_i, 1, x["seed"], normal)
    ws3.write_number(row_i, 2, x["mode"], normal)
    ws3.write_number(row_i, 3, x["score"], decimal)
    ws3.write_number(row_i, 4, x["rmsd_lb"], decimal)
    ws3.write_number(row_i, 5, x["rmsd_ub"], decimal)
    ws3.write(row_i, 6, x["file"], normal)

ws3.set_column("A:A", 18)
ws3.set_column("B:C", 13)
ws3.set_column("D:F", 22)
ws3.set_column("G:G", 42)

ws3.freeze_panes(1, 0)


# ============================================================
# SHEET 4 — AK1 VALIDATION + BOX
# ============================================================

ws4 = workbook.add_worksheet("Validation_and_Box")

headers4 = [
    "Target",
    "Structure",
    "Validation ligand",
    "Redocking RMSD (Å)",
    "Center X",
    "Center Y",
    "Center Z",
    "Size X",
    "Size Y",
    "Size Z"
]

for col, h in enumerate(headers4):
    ws4.write(0, col, h, header)

for row_i, (target, x) in enumerate(validation.items(), start=1):

    vals = [
        target,
        x.get("structure"),
        x.get("validation_ligand"),
        x.get("redocking_rmsd_A"),
        x.get("box_center_x"),
        x.get("box_center_y"),
        x.get("box_center_z"),
        x.get("box_size_x"),
        x.get("box_size_y"),
        x.get("box_size_z")
    ]

    for col, value in enumerate(vals):

        if col >= 3 and value is not None:
            ws4.write_number(row_i, col, value, decimal)
        else:
            ws4.write(row_i, col, value, normal)

ws4.set_column("A:C", 24)
ws4.set_column("D:J", 18)


# ============================================================
# SHEET 5 — METHODS / COMPUTATIONAL SETTINGS
# ============================================================

ws5 = workbook.add_worksheet("Methods")

ws5.merge_range(
    "A1:B1",
    "AK1 Docking Computational Parameters",
    title
)

method_rows = [
    ["Target", "Human adenylate kinase 1 (AK1)"],
    ["Structure", "PDB 1Z83, chain A"],
    ["Docking ligand", "Isorhamnetin"],
    ["Isorhamnetin source", "PubChem CID 5281654"],
    ["Native ligand", "AP5"],
    ["Docking engine", "AutoDock Vina 1.2.7"],
    ["Receptor preparation", "PDBFixer followed by Meeko"],
    ["Meeko version", "0.8.0"],
    ["Docking exhaustiveness", 16],
    ["Number of modes per seed", 10],
    ["Independent seeds", len([x for x in all_runs if x["target"] == "AK1"])],
    ["Binding-site definition", "Box centered on crystallographic AP5 binding region"],
    ["AP5 redocking RMSD", "1.776 Å"],
    ["Interpretation caution",
     "Vina scores are predicted docking scores and should not be interpreted as experimental binding free energies."]
]

for row, item in enumerate(method_rows, start=2):
    ws5.write(row, 0, item[0], header)
    ws5.write(row, 1, item[1], normal)

ws5.set_column("A:A", 28)
ws5.set_column("B:B", 65)


# ============================================================
# SCORE DISTRIBUTION CHART
# ============================================================

if all_runs:

    chart = workbook.add_chart({
        "type": "column"
    })

    nrows = len(sorted_runs)

    chart.add_series({
        "name": "Best Vina score",
        "categories": [
            "ISO_Runs",
            1,
            1,
            nrows,
            1
        ],
        "values": [
            "ISO_Runs",
            1,
            2,
            nrows,
            2
        ]
    })

    chart.set_title({
        "name": "AK1 ISO docking score across independent seeds"
    })

    chart.set_x_axis({
        "name": "Seed"
    })

    chart.set_y_axis({
        "name": "Vina score (kcal/mol)"
    })

    chart.set_legend({
        "none": True
    })

    ws.insert_chart(
        "N3",
        chart,
        {
            "x_scale": 1.4,
            "y_scale": 1.3
        }
    )


# ============================================================
# SAVE
# ============================================================

workbook.close()


print()
print("=" * 60)
print("DOCKING WORKBOOK SUCCESSFULLY CREATED")
print("=" * 60)

print("Output:")
print(OUTPUT)

print()

print("Targets found:", targets)
print("ISO docking runs found:", len(all_runs))
print("Individual poses found:", len(all_modes))

print()

if summary:

    for x in summary:

        print(
            x["target"],
            "| n =", x["n_runs"],
            "| best =", round(x["best_score"], 3),
            "| mean =", round(x["mean_score"], 3),
            "| median =", round(x["median_score"], 3),
            "| SD =",
            round(x["sd_score"], 3)
            if x["sd_score"] is not None
            else "N/A"
        )

print()