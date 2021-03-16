"""
This script fetches all files from the OSF project 6b3jr (kinase-conformational-modeling) and stores them
in a directory named osf.
"""
from pathlib import Path

import osfclient


DATA = Path("../osf")

# connect to osf project
osf = osfclient.OSF()
project = osf.project("6b3jr")
storage = project.storage()

# download files
files = list(storage.files)
for i, file in enumerate(files):
    print(f"Downloading file {i + 1}/{len(files)}")
    local_path = DATA / file.path[1:]
    local_path.parent.mkdir(parents=True, exist_ok=True)
    file.write_to(open(local_path, "wb"))
