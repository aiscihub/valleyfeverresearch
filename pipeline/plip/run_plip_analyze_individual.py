import os
from plip.structure.preparation import PDBComplex
import pandas as pd
import os
import glob

# Customize this path to the *top-level* protein folder you're working with
protein_folder = "CIMG_00533"
base_path = f"/home/user/valleyfever/dataset/protein_db/md/{protein_folder}/simulation_explicit"

# Output results folder
output_dir = os.path.join(base_path, "plip_results")
os.makedirs(output_dir, exist_ok=True)

summary_data = []
pdb_files = glob.glob(f"{base_path}/**/*explicit_stripped*.pdb", recursive=True)


# Loop through all complex pdbs
for filename in pdb_files:
    if "explicit_stripped" in filename and filename.endswith(".pdb"):
        complex_path = os.path.join(base_path, filename)
        complex_name = filename.replace(".pdb", "")

        print(f"Processing: {filename}")
        complex = PDBComplex()
        complex.load_pdb(complex_path)
        complex.analyze()

        for lig in complex.ligands:
            key = f"{lig.hetid}:{lig.chain}:{lig.position}"
            print(f"{key}")
            if key not in complex.interaction_sets:
                print(f"No interaction set found for: {key} in {filename}")
                continue

            interactions = complex.interaction_sets[key]
            rows = []

            for hb in interactions.hbonds_ldon + interactions.hbonds_pdon:
                rows.append(["HydrogenBond", f"{hb.resnr}{hb.restype}", round(hb.distance_ad, 2), f"{hb.atype}->{hb.dtype}"])

            for hyd in interactions.hydrophobic_contacts:
                rows.append(["Hydrophobic", f"{hyd.resnr}{hyd.restype}", round(hyd.distance, 2), ""])

            for sb in interactions.saltbridge_lneg + interactions.saltbridge_pneg:
                rows.append(["SaltBridge", f"{sb.resnr}{sb.restype}", round(sb.distance, 2), f"{sb.reschain}↔{sb.reschain_l}"])

            for pi in interactions.pistacking:
                rows.append(["PiStacking", f"{pi.resnr}{pi.restype}", round(pi.distance, 2), f"type={pi.type}"])

            for wb in interactions.water_bridges:
                rows.append(["WaterBridge", f"{wb.resnr}{wb.restype}", round(wb.dist_h2o, 2), f"{wb.donor}->{wb.acceptor}"])

            for metal in interactions.metal_complexes:
                rows.append(["MetalComplex", f"{metal.resnr}{metal.restype}", round(metal.distance, 2), f"{metal.metaltype}"])

                if rows:
                    df = pd.DataFrame(rows, columns=["Type", "Residue", "Distance", "Details"])
                    out_path = os.path.join(output_dir, f"{complex_name}_plip.csv")
                    df.to_csv(out_path, index=False)
                    print(f"Saved: {out_path}")

                    for row in rows:
                        summary_data.append([complex_name] + row)
                else:
                    print(f"No interactions found for {key} in {filename}")

            # Add to global summary
            for row in rows:
                summary_data.append([complex_name] + row)

# Write summary file
summary_df = pd.DataFrame(summary_data, columns=["Complex", "Type", "Residue", "Distance", "Details"])
summary_df.to_csv(os.path.join(output_dir, "all_plip_interactions_summary.csv"), index=False)
print("PLIP profiling completed for:", protein_folder)
