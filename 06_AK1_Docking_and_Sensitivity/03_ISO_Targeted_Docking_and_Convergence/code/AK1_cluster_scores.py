from pathlib import Path
import re
import statistics

PROJECT = Path(r"C:\Users\akidi\Desktop\ISO_AK_Docking")
DOCKING = PROJECT / "docking"

dominant = {
    1001, 2002, 3003, 4004,
    5005, 7007, 9009, 10101
}

alternative = {
    6006, 8008
}


def get_best_score(path):

    with open(path, "r") as f:

        for line in f:

            if line.startswith("REMARK VINA RESULT:"):

                return float(
                    line.split()[3]
                )

    raise ValueError(
        f"No Vina result found in {path}"
    )


groups = {
    "Dominant cluster": dominant,
    "Alternative cluster": alternative
}


for name, seeds in groups.items():

    scores = []

    print()
    print("=" * 55)
    print(name)
    print("=" * 55)

    for seed in sorted(seeds):

        file = (
            DOCKING /
            f"AK1_ISO_seed{seed}.pdbqt"
        )

        score = get_best_score(file)

        scores.append(score)

        print(
            f"Seed {seed}: {score:.3f} kcal/mol"
        )

    print()
    print("n =", len(scores))
    print(
        "Mean =",
        f"{statistics.mean(scores):.3f}"
    )
    print(
        "Median =",
        f"{statistics.median(scores):.3f}"
    )

    if len(scores) > 1:
        print(
            "SD =",
            f"{statistics.stdev(scores):.3f}"
        )

    print(
        "Best =",
        f"{min(scores):.3f}"
    )

    print(
        "Worst =",
        f"{max(scores):.3f}"
    )