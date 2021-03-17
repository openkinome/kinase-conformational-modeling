"""
This script generates a dictionary containing the KLIFS quality score for each prepared kinase structure and saves it
in JSON format.
"""
import json
from pathlib import Path

from tqdm import tqdm

from opencadd.databases.klifs import setup_remote


pdb_files = list(Path("../osf/kinase_structures/human_abl1").glob("*.pdb"))
# pdb_files = list(Path("../osf/kinase_structures/human_egfr").glob("*.pdb"))

remote = setup_remote()
structure_quality_dictionary = dict()
# get KLIFS quality score for each prepared structure
for pdb_file in tqdm(pdb_files, total=len(pdb_files)):
    # identify structure
    name_content = pdb_file.stem.split("_")
    pdb_id = name_content[3]
    chain_id = name_content[4][5:]
    alternate_location = name_content[5][6:]
    if alternate_location == "None":
        alternate_location = None
    klifs_structure = remote.structures.by_structure_pdb_id(
        pdb_id,
        structure_chain=chain_id,
        structure_alternate_model=alternate_location
    )
    # check if KLIFS query results in exactly one match
    if len(klifs_structure) != 1:
        raise LookupError(
            "Could not assign a quality score. "
            "Querying KLIFS for "
            f"PDB id {pdb_id}, "
            f"chain {chain_id} and "
            f"alternate location {alternate_location} "
            f"resulted in retrieval of {len(klifs_structure)} structures. "
            )
    quality_score = round(float(klifs_structure["structure.qualityscore"].iloc[0]), 1)
    structure_quality_dictionary[pdb_file.stem] = quality_score

# save dictionary
with open("structure_quality_dictionary.json", "w") as file:
    json.dump(structure_quality_dictionary, file, indent=4, sort_keys=True)
