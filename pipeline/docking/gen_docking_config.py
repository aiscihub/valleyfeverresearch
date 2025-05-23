import os
import csv
import math
from pathlib import Path
from pocket_util import grid_edge_from_sas

def generate_vina_configs(protein, ligand, receptor_file, pocket_csv, output_root, ligand_dir, project_root):
    ligand_file = ligand_dir / f"{ligand}.pdbqt"
    output_folder = output_root / protein / "docking_dynamic" / f"{protein}_prepared_{ligand}"
    output_folder.mkdir(parents=True, exist_ok=True)

    # Define large ligands list
    large_ligands = {"Tacrolimus", "Beauvericin"}

    with open(pocket_csv, newline='') as f:
        reader = csv.DictReader(f)
        reader.fieldnames = [h.strip() for h in reader.fieldnames]
        for i, row in enumerate(reader, 1):
            try:
                x = float(row['center_x'].strip())
                y = float(row['center_y'].strip())
                z = float(row['center_z'].strip())
                sas_points = int(row['sas_points'].strip())
                grid_size = grid_edge_from_sas(sas_points)

                # === Large ligand adjustment ===
                if ligand in large_ligands:
                    grid_size = max(grid_size, 28)  # force at least 28 Å grid if ligand is large
                else:
                    continue
                config_path = output_folder / f"vina_config_pocket{i}.txt"
                with open(config_path, "w") as config:
                    config.write(f"receptor = {receptor_file.relative_to(project_root)}\n")
                    config.write(f"ligand = {ligand_file.relative_to(project_root)}\n\n")
                    config.write(f"center_x = {x:.4f}\ncenter_y = {y:.4f}\ncenter_z = {z:.4f}\n\n")
                    config.write(f"size_x = {grid_size}\nsize_y = {grid_size}\nsize_z = {grid_size}\n\n")
                    config.write("exhaustiveness = 32\nnum_modes = 10\nenergy_range = 4\n")
            except Exception as e:
                print(f"Skipping pocket {i} due to error: {e}")


def main():
    project_root = Path(__file__).resolve().parents[2]
    db_root = project_root / "dataset" / "protein_db"
    docking_dir = db_root / "docking"

    proteins = ["CIMG_01418", "CIMG_00533", "CIMG_00780", "CIMG_06197", "CIMG_09093"]
    #proteins = ["CIMG_00533"]
    for protein in proteins:
        prankweb_dir = docking_dir / f"{protein}" / f"prankweb-{protein}_relaxed"
        protein_dir = docking_dir / protein / "protein"
        receptor_file = protein_dir / f"{protein}_prepared.pdbqt"
        pocket_csv = prankweb_dir / "structure.pdb_predictions.csv"
        ligand_dir = docking_dir / f"{protein}" /"ligands"
        ligands = [l.stem for l in ligand_dir.glob("*.pdbqt")]
        for ligand in ligands:
            generate_vina_configs(protein, ligand, receptor_file, pocket_csv, docking_dir, ligand_dir, project_root)


if __name__ == "__main__":
    main()
