"""
This script generates a dictionary mapping the canonical UniProt residue IDs of a kinase to the KLIFS pocket residue
IDs and saves it in JSON format.
"""
import json
import tempfile

from openeye import oechem

from kinoml.core.sequences import AminoAcidSequence
from kinoml.features.complexes import OEKLIFSKinaseApoFeaturizer
from kinoml.modeling.OEModeling import (
    read_molecules,
    select_chain,
    remove_non_protein,
    get_structure_sequence_alignment
)
from opencadd.databases.klifs import setup_remote


klifs_kinase_id = 392  # ABL1
# klifs_kinase_id = 406  # EGFR

# get the highest quality structure of the kinase of interest with complete KLIFS pocket
remote = setup_remote()
structures = remote.structures.by_kinase_klifs_id(klifs_kinase_id)
structures["structure.pocket"] = structures["structure.pocket"].str.replace("_", "")
structures = structures[structures["structure.pocket"].str.len() == 85]
structures = structures.sort_values(
    by=["structure.qualityscore", "structure.resolution", "structure.chain", "structure.alternate_model"],
    ascending=[False, True, True, True]
)
klifs_structure = structures.iloc[0]

# read structure and select chain and protein
with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb") as file:
    file.write(remote.coordinates.to_text(klifs_structure["structure.klifs_id"], extension="pdb"))
    structure = read_molecules(file.name)[0]
    structure = select_chain(structure, klifs_structure["structure.chain"])
    structure = remove_non_protein(structure, remove_water=True)

# retrieve sequence from UniProt
uniprot_id = remote.kinases.by_kinase_klifs_id(klifs_kinase_id)["kinase.uniprot"].iloc[0]
sequence = AminoAcidSequence.from_uniprot(uniprot_id)
uniprot_resids = OEKLIFSKinaseApoFeaturizer._get_kinase_residue_numbers(structure, sequence)

# map residue IDs to KLIFS residue IDs
uniprot_resids_to_klifs_resids = dict()
pocket = remote.pockets.by_structure_klifs_id(klifs_structure["structure.klifs_id"])
hier_view = oechem.OEHierView(structure)
structure_resids = [residue.GetResidueNumber() for residue in hier_view.GetResidues()]
for structure_resid, uniprot_resid in zip(structure_resids, uniprot_resids):
    if str(structure_resid) in pocket["residue.id"].unique():
        klifs_residue_id = str(
            pocket[pocket["residue.id"] == str(structure_resid)]["residue.klifs_id"].iloc[0]
        )
        uniprot_resids_to_klifs_resids[str(uniprot_resid)] = klifs_residue_id

# save dictionary
with open("uniprot_resids_to_klifs_resids.json", "w") as file:
    json.dump(uniprot_resids_to_klifs_resids, file, indent=4, sort_keys=True)
