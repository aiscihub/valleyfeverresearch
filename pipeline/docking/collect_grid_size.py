import os
from pathlib import Path
import re
import pandas as pd

# Define the root directory containing docking configs
docking_root = Path("/home/user/git/valleyfever/dataset/protein_db/docking")

# Prepare a list to collect config data
config_data = []

# Loop through docking_dynamic folders
for protein_dir in docking_root.glob("*/docking_dynamic/*"):
    protein = protein_dir.parts[-3]
    ligand = protein_dir.parts[-1].replace(f"{protein}_prepared_", "")
    config_files = list(protein_dir.glob("vina_config_pocket*.txt"))

    for config_file in config_files:
        pocket_match = re.search(r"vina_config_pocket(\d+)\.txt", config_file.name)
        pocket = int(pocket_match.group(1)) if pocket_match else None

        # Initialize grid and center values
        size_x = size_y = size_z = center_x = center_y = center_z = None

        with open(config_file) as f:
            for line in f:
                if line.startswith("size_x"):
                    size_x = float(line.split("=")[1].strip())
                elif line.startswith("size_y"):
                    size_y = float(line.split("=")[1].strip())
                elif line.startswith("size_z"):
                    size_z = float(line.split("=")[1].strip())
                elif line.startswith("center_x"):
                    center_x = float(line.split("=")[1].strip())
                elif line.startswith("center_y"):
                    center_y = float(line.split("=")[1].strip())
                elif line.startswith("center_z"):
                    center_z = float(line.split("=")[1].strip())

        config_data.append({
            "protein": protein,
            "ligand": ligand,
            "pocket": pocket,
            "center_x": center_x,
            "center_y": center_y,
            "center_z": center_z,
            "size_x": size_x,
            "size_y": size_y,
            "size_z": size_z,
            "config_path": str(config_file)
        })

# Save as CSV
config_df = pd.DataFrame(config_data)
output_path = ("/home/user/git/valleyfever/dataset/protein_db/results/compositescore/RawDatasetSummary"
               "/vina_config_grid_summary.csv")
config_df.to_csv(output_path, index=False)
output_path
