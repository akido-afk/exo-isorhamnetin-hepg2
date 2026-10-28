from pathlib import Path
import math

PROJECT = Path(r"C:\Users\akidi\Desktop\ISO_AK_Docking")
PDB = PROJECT / "structures" / "AK1_1Z83.pdb"


def dist(a, b):
    return math.sqrt(
        (a[0]-b[0])**2 +
        (a[1]-b[1])**2 +
        (a[2]-b[2])**2
    )


atoms = []

with open(PDB, "r") as f:
    for line in f:

        if not (
            line.startswith("ATOM")
            or line.startswith("HETATM")
        ):
            continue

        record = line[0:6].strip()
        atom = line[12:16].strip()
        resn = line[17:20].strip()
        chain = line[21].strip()
        resi = line[22:26].strip()

        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])

        element = line[76:78].strip()

        if not element:
            element = atom[0]

        atoms.append({
            "record": record,
            "atom": atom,
            "resn": resn,
            "chain": chain,
            "resi": resi,
            "element": element.upper(),
            "coord": (x, y, z)
        })


# Find chain-A zinc
zn_atoms = [
    a for a in atoms
    if a["resn"] == "ZN"
    and a["chain"] == "A"
]

if not zn_atoms:
    raise RuntimeError("No chain-A Zn found.")

print()
print("=" * 80)
print("CHAIN-A Zn ENVIRONMENT IN PDB 1Z83")
print("=" * 80)

for zn in zn_atoms:

    print()
    print(
        f"Zn: chain {zn['chain']} residue {zn['resi']} "
        f"at {zn['coord']}"
    )

    neighbours = []

    for a in atoms:

        if a is zn:
            continue

        d = dist(
            zn["coord"],
            a["coord"]
        )

        if d <= 5.0:
            neighbours.append((d, a))

    neighbours.sort(
        key=lambda x: x[0]
    )

    print()
    print("Atoms within 5.0 A:")
    print("-" * 80)

    for d, a in neighbours:

        category = "OTHER"

        if a["record"] == "ATOM":
            category = "PROTEIN"

        elif a["resn"] == "AP5":
            category = "AP5"

        elif a["resn"] in {"HOH", "WAT"}:
            category = "WATER"

        print(
            f"{d:6.3f} A   "
            f"{category:<8} "
            f"{a['resn']:>4} "
            f"{a['chain']} "
            f"{a['resi']:>4} "
            f"{a['atom']:<5}"
        )


# Specifically summarize nearest AP5 / protein / water distances
def nearest(category_filter):

    hits = []

    for zn in zn_atoms:
        for a in atoms:

            if a is zn:
                continue

            if not category_filter(a):
                continue

            d = dist(
                zn["coord"],
                a["coord"]
            )

            hits.append(
                (d, a)
            )

    hits.sort(
        key=lambda x: x[0]
    )

    return hits[:10]


print()
print("=" * 80)
print("NEAREST AP5 ATOMS")
print("=" * 80)

for d, a in nearest(
    lambda x:
        x["resn"] == "AP5"
        and x["chain"] == "A"
):
    print(
        f"{d:6.3f} A   "
        f"AP5 {a['atom']}"
    )


print()
print("=" * 80)
print("NEAREST AK1 PROTEIN ATOMS")
print("=" * 80)

for d, a in nearest(
    lambda x:
        x["record"] == "ATOM"
        and x["chain"] == "A"
):
    print(
        f"{d:6.3f} A   "
        f"{a['resn']}{a['resi']} "
        f"{a['atom']}"
    )


print()
print("=" * 80)
print("NEAREST WATERS")
print("=" * 80)

for d, a in nearest(
    lambda x:
        x["resn"] in {"HOH", "WAT"}
        and x["chain"] == "A"
):
    print(
        f"{d:6.3f} A   "
        f"{a['resn']}{a['resi']} "
        f"{a['atom']}"
    )