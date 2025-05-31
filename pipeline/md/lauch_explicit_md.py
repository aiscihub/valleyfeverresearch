# launch_implicit_replicas.py

import os
import glob
import os
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from run_openmm_md import run_multiple_replicas

protein = "CIMG_00533"
pocket_id = "pocket6"
base_dir = f""
#example
base_dir = f"/home/user/valleyfever/dataset/protein_db/md/{protein}/simulation_explicit/{pocket_id}"

# === Ligand Configuration ===
# Set to [] to run all ligands in the directory
selected_ligands = ["Tacrolimus", "Verapamil", "Beauvericin", "Disulfiram", "Milbemycin"]  # Example: run only these ligands
selected_ligands = ["Milbemycin"]
# === Discover all ligands from SDF files
# filessdf_files = glob.glob(os.path.join(base_dir, "*.sdf"))
# ligand_names = [os.path.splitext(os.path.basename(f))[0] for f in sdf_files]

# === Apply filter if configured
# Example call:
for ligand_name in selected_ligands:

    complex_files = glob.glob(os.path.join(base_dir, f"{protein}_prepared_{ligand_name}_{pocket_id}_complex.pdb"))
    n_steps = 1000000
    for complex_path in complex_files:
        print(f"{ligand_name}---{complex_path}---- {pocket_id}")
        prefix = os.path.basename(complex_path).replace("", "")
        print(f"\n Running MD for: {prefix}, n_steps = {n_steps}")
        ligand_name = os.path.basename(complex_path).replace(".pdb", "_aligned_ligand.sdf")

        # run_protein_ligand_md(
        #     base_dir=base_dir,
        #     sdf_filename=ligand_name,
        #     pdb_filename=complex_path,
        #     n_steps=5000000
        # )
        run_multiple_replicas(
            parent_dir=base_dir,
            base_dir=base_dir,
            sdf_filename=ligand_name,
            pdb_filename=complex_path,
            n_steps=n_steps
        )
