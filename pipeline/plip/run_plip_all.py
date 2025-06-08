import os
from plip.structure.preparation import PDBComplex
import pandas as pd

def run_plip_on_protein(protein_id, base_root):
    print(f"Starting PLIP analysis for: {protein_id}")
    base_path = os.path.join(base_root, protein_id, "complex")
    output_dir = os.path.join(base_path, "plip_results")
    os.makedirs(output_dir, exist_ok=True)

    summary_data = []

    for filename in os.listdir(base_path):
        if filename.endswith("_complex.pdb"):
            complex_path = os.path.join(base_path, filename)
            complex_name = filename.replace("_complex.pdb", "")

            print(f"🔬 Processing: {complex_path}/{filename}")
            complex = PDBComplex()
            complex.load_pdb(complex_path)
            complex.analyze()

            for lig in complex.ligands:
                key = f"{lig.hetid}:{lig.chain}:{lig.position}"
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
                   # df.to_csv(out_path, index=False)
                   # print(f"Saved: {out_path}")

                    for row in rows:
                        summary_data.append([complex_name] + row)
                else:
                    print(f"No interactions found for {key} in {filename}")

    # Save summary for this protein
    summary_df = pd.DataFrame(summary_data, columns=["Complex", "Type", "Residue", "Distance", "Details"])
    summary_csv = os.path.join(output_dir, "all_plip_interactions_summary.csv")
    summary_df.to_csv(summary_csv, index=False)
    print(f"Summary saved to: {summary_csv}\n")


# Run on all proteins
protein_list = ["CIMG_00533", "CIMG_00780", "CIMG_01418", "CIMG_06197", "CIMG_09093"]

project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
base_path = os.path.join(project_root, "dataset/protein_db/docking")

for protein_id in protein_list:
    run_plip_on_protein(protein_id, base_path)

print("🎉 All PLIP analyses completed.")
