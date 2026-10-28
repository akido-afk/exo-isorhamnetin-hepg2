from pathlib import Path

PROJECT = Path(r"C:\Users\akidi\Desktop\ISO_AK_Docking")

FIXED_PROTEIN = PROJECT / "receptors" / "AK1_chainA_FIXED.pdb"
ORIGINAL_PDB  = PROJECT / "structures" / "AK1_1Z83.pdb"

OUTPUT = PROJECT / "receptors" / "AK1_chainA_FIXED_ZN601.pdb"

if not FIXED_PROTEIN.exists():
    raise FileNotFoundError(FIXED_PROTEIN)

if not ORIGINAL_PDB.exists():
    raise FileNotFoundError(ORIGINAL_PDB)

# ------------------------------------------------------------
# Read repaired AK1 protein
# ------------------------------------------------------------

protein_lines = []

with open(FIXED_PROTEIN, "r") as f:
    for line in f:

        if line.startswith("ATOM"):
            protein_lines.append(line)


# ------------------------------------------------------------
# Extract only chain-A ZN601 from original 1Z83
# ------------------------------------------------------------

zn_lines = []

with open(ORIGINAL_PDB, "r") as f:
    for line in f:

        if not line.startswith("HETATM"):
            continue

        resn = line[17:20].strip()
        chain = line[21].strip()
        resi = line[22:26].strip()

        if (
            resn == "ZN"
            and chain == "A"
            and resi == "601"
        ):
            zn_lines.append(line)


if len(zn_lines) != 1:
    raise RuntimeError(
        f"Expected exactly 1 ZN601 atom, found {len(zn_lines)}"
    )


# ------------------------------------------------------------
# Write combined receptor
# ------------------------------------------------------------

with open(OUTPUT, "w") as f:

    for line in protein_lines:
        f.write(line)

    f.write(zn_lines[0])

    f.write("TER\n")
    f.write("END\n")


print("=" * 70)
print("AK1 FIXED + ZN601 RECEPTOR CREATED")
print("=" * 70)

print("Protein ATOM records:", len(protein_lines))
print("ZN601 records:", len(zn_lines))
print("Output:")
print(OUTPUT)