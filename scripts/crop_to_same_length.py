"""
This script crops protein structures given in OpenEye's OEB format to the same common sequence, adds caps at N and C
terminus, removes structures with gaps and saves the processed structures in PDB format.
"""
from pathlib import Path

from openeye import oechem, oespruce
from tqdm import tqdm

from kinoml.modeling.OEModeling import (
    read_molecules,
    assign_caps,
    write_molecules,
    update_residue_identifiers
)


print("Analyzing structure sequences ...")
oeb_files = list(Path(".").glob("*.oeb"))
structure_dict = dict()
for oeb_file in oeb_files:
    structure = read_molecules(oeb_file)[0]
    hier_view = oechem.OEHierView(structure)
    resids = [
        residue.GetResidueNumber()
        for residue in hier_view.GetResidues()
        if residue.GetResidueName().strip() not in ["NME", "ACE"]
        and not residue.GetOEResidue().IsHetAtom()
    ]
    structure_dict[oeb_file] = resids
common_start = max([min(x) for x in structure_dict.values()])
common_end = min([max(x) for x in structure_dict.values()])
common_residue_numbers = list(range(common_start, common_end + 1))
print(f"Common sequence: {common_start}-{common_end}")

Path("prep").mkdir(exist_ok=True)
excluded = []
for oeb_file in tqdm(oeb_files, total=len(oeb_files)):
    resids = structure_dict[oeb_file]
    if all(resid in resids for resid in common_residue_numbers):
        structure = read_molecules(oeb_file)[0]
        hier_view = oechem.OEHierView(structure)
        for residue in hier_view.GetResidues():
            if residue.GetResidueNumber() not in common_residue_numbers:
                for atom in residue.GetAtoms():
                    structure.DeleteAtom(atom)
        structure = assign_caps(structure)
        structure = update_residue_identifiers(structure)
        pdb_file_prep = Path("prep") / (oeb_file.stem + ".pdb")
        write_molecules([structure], pdb_file_prep)
    else:
        excluded.append(str(oeb_file))
print(
    f"Cropped {len(oeb_files) - len(excluded)} "
    "structures to same length."
)
if len(excluded) > 0:
    excluded = "\n" + "\n".join(excluded)
    print(f"Structures with missing residues: {excluded}")
