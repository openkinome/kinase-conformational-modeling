"""
This script prepares apo structures of kinases for a given KLIFS kinase ID. NMR structures are ignored.
It uses functionalities from KinoML and opencadd:
https://github.com/openkinome/kinoml/commit/550e902037a62c120f71a83dd04f95e2ad408666
https://github.com/volkamerlab/opencadd/commit/3fe042334c17368deaa17794e5368a632f964e3a

The specified opencadd commit contains a fix for wrongly assigned alternate locations for human ABL1 and EGFR
structures in KLIFS (March 2021). If you run this script for a different kinase, you may experience errors associated
to alternate locations.
"""
import logging

from kinoml.core.proteins import PDBProtein
from kinoml.core.systems import Protein
from kinoml.features.complexes import OEKLIFSKinaseApoFeaturizer
from opencadd.databases.klifs import setup_remote


klifs_kinase_id = 392  # ABL1
# klifs_kinase_id = 406  # EGFR
loop_db = "~/.OpenEye/rcsb_spruce.loop_db"

logging.basicConfig(
   filename=f"apo_kinases_392.log",
   filemode="a",
   level=logging.DEBUG
)
remote = setup_remote()
structures = remote.structures.by_kinase_klifs_id(klifs_kinase_id)
structures = structures[structures['structure.resolution'].notnull()]  # remove NMR structures
docking_featurizer = OEKLIFSKinaseApoFeaturizer(loop_db=loop_db)
print(f"Generating {len(structures)} apo structures ...")
for i, (index, row) in enumerate(structures.iterrows()):
    pdb_id = row["structure.pdb_id"]
    chain_id = row["structure.chain"]
    alternate_location = row["structure.alternate_model"]
    print(
        " ".join([
            f"Preparing structure",
            f"PDB {pdb_id}",
            f"chain {chain_id}",
            f"altloc {alternate_location}",
            f"({i + 1}/{len(structures)})"
        ])
    )
    pdb_protein = PDBProtein(pdb_id=pdb_id)
    pdb_protein.chain_id = chain_id
    pdb_protein.alternate_location = alternate_location
    protein = Protein(components=[pdb_protein])
    system = docking_featurizer.featurize(protein)
print("Finished")
