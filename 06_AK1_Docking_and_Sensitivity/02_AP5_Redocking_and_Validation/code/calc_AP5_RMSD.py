from rdkit import Chem
from rdkit.Geometry import Point3D
import math

native_pdb = r"..\ligand\AK1_native_AP5_chainA.pdb"
ideal_sdf = r"..\ligand\AP5_ideal.sdf"
docked_sdf = r"..\docking\AK1_AP5_redock_seed1001.sdf"

# --------------------------------------------------
# Read native crystallographic AP5 coordinates
# --------------------------------------------------

native_coords = []
native_elements = []

with open(native_pdb) as f:
    for line in f:
        if line.startswith("HETATM"):

            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            element = line[76:78].strip()

            if not element:
                element = ''.join(
                    c for c in line[12:16] if c.isalpha()
                )[0]

            native_coords.append((x, y, z))
            native_elements.append(element.capitalize())

# --------------------------------------------------
# Use the RCSB ideal AP5 molecule as connectivity template
# --------------------------------------------------

ideal = Chem.SDMolSupplier(
    ideal_sdf,
    removeHs=False
)[0]

ideal = Chem.RemoveHs(ideal)

ideal_elements = [
    atom.GetSymbol()
    for atom in ideal.GetAtoms()
]

print("Native heavy atoms:", len(native_coords))
print("Template heavy atoms:", ideal.GetNumAtoms())

if len(native_coords) != ideal.GetNumAtoms():
    raise ValueError(
        "Native and template AP5 atom counts do not match."
    )

if native_elements != ideal_elements:
    print("WARNING: Element ordering differs.")
    print("Native :", native_elements)
    print("Template:", ideal_elements)
    raise ValueError(
        "Stop here because atom order must be checked."
    )

# Replace ideal coordinates with the actual crystal coordinates

reference = Chem.Mol(ideal)

reference.RemoveAllConformers()

conf = Chem.Conformer(reference.GetNumAtoms())

for i, (x, y, z) in enumerate(native_coords):
    conf.SetAtomPosition(i, Point3D(x, y, z))

reference.AddConformer(conf)

# --------------------------------------------------
# Read all docked AP5 poses
# --------------------------------------------------

supplier = Chem.SDMolSupplier(
    docked_sdf,
    removeHs=False
)

poses = [
    Chem.RemoveHs(m)
    for m in supplier
    if m is not None
]

# --------------------------------------------------
# Direct coordinate RMSD with symmetry-aware atom mapping
#
# IMPORTANT:
# We do NOT superimpose the ligand here.
# Both poses already share the same AK1 coordinate system.
# --------------------------------------------------

def direct_rmsd(ref, pose, mapping):

    ref_conf = ref.GetConformer()
    pose_conf = pose.GetConformer()

    squared = []

    for ref_idx, pose_idx in enumerate(mapping):

        a = ref_conf.GetAtomPosition(ref_idx)
        b = pose_conf.GetAtomPosition(pose_idx)

        squared.append(
            (a.x - b.x)**2 +
            (a.y - b.y)**2 +
            (a.z - b.z)**2
        )

    return math.sqrt(sum(squared) / len(squared))


print()
print("AP5 REDOCKING RMSD")
print("------------------")

results = []

for pose_number, pose in enumerate(poses, start=1):

    if pose.GetNumAtoms() != reference.GetNumAtoms():
        print(
            f"Pose {pose_number}: atom-count mismatch"
        )
        continue

    # Find chemically equivalent atom mappings
    mappings = pose.GetSubstructMatches(
        reference,
        uniquify=False,
        maxMatches=10000
    )

    if not mappings:
        raise ValueError(
            f"No atom mapping found for pose {pose_number}"
        )

    rmsds = [
        direct_rmsd(reference, pose, mapping)
        for mapping in mappings
    ]

    best = min(rmsds)

    results.append(best)

    print(
        f"Pose {pose_number}: {best:.3f} A"
    )

print()
print(
    "Best redocking RMSD:",
    f"{min(results):.3f} A"
)