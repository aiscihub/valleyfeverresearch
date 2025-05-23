import os
import re
import numpy as np
import mdtraj as md
import pandas as pd

def measure_drift_from_pdbs(pdb_initial, pdb_final, ligand_resnames=["UNK", "MOL", "LIG", "UNL", "BEA"]):
    traj_initial = md.load(pdb_initial)
    traj_final = md.load(pdb_final)

    if traj_initial.n_atoms != traj_final.n_atoms:
        raise ValueError("Atom count mismatch between initial and final frame.")

    ligand_atoms = None
    for resname in ligand_resnames:
        atoms = traj_initial.topology.select(f"resname {resname}")
        if len(atoms) > 0:
            ligand_atoms = atoms
            break

    if ligand_atoms is None:
        raise ValueError("No ligand atoms found.")

    com_initial = traj_initial.xyz[0][ligand_atoms].mean(axis=0)
    com_final = traj_final.xyz[0][ligand_atoms].mean(axis=0)
    drift_vector = com_final - com_initial
    drift_distance = np.linalg.norm(drift_vector)

    return drift_distance, drift_vector, len(ligand_atoms)

# === Batch Runner ===
def run_batch_drift_analysis_from_pdbs(base_root_dir, protein, pocket_id):
    pocket_dir = os.path.join(base_root_dir, f"pocket{pocket_id}")
    all_results = []

    for replica in sorted(os.listdir(pocket_dir)):
        if not replica.startswith("replica_"):
            continue

        replica_path = os.path.join(pocket_dir, replica)

        for fname in os.listdir(replica_path):
            if fname.endswith("_explicit_stripped_initial_frame.pdb"):
                prefix = fname.replace("_explicit_stripped_initial_frame.pdb", "")
                initial_pdb = os.path.join(replica_path, fname)
                final_pdb = os.path.join(replica_path, f"{prefix}_explicit_stripped_final_frame.pdb")

                if not os.path.exists(final_pdb):
                    print(f"⚠️ Missing final frame for {prefix}")
                    continue

                try:
                    drift_distance, drift_vector, ligand_atoms = measure_drift_from_pdbs(initial_pdb, final_pdb)
                    status = "Stayed" if drift_distance <= 0.5 else "Drifted"

                    ligand_match = re.search(rf"{protein}_prepared_(.+?)_pocket\d+_complex", prefix)
                    ligand = ligand_match.group(1) if ligand_match else "Unknown"

                    all_results.append({
                        'protein': protein,
                        'pocket': f"pocket{pocket_id}",
                        'replica': replica,
                        'ligand': ligand,
                        'prefix': prefix,
                        'drift_distance': drift_distance,
                        'drift_x': drift_vector[0],
                        'drift_y': drift_vector[1],
                        'drift_z': drift_vector[2],
                        'status': status,
                        'ligand_atoms': ligand_atoms
                    })
                except Exception as e:
                    print(f"❌ Failed {prefix}: {str(e)}")
                    continue

    df = pd.DataFrame(all_results)
    out_path = os.path.join(pocket_dir, f"{protein}_drift_summary_from_pdb.csv")
    df.to_csv(out_path, index=False)
    print(f"✅ Drift summary written to {out_path}")

# Example usage:
run_batch_drift_analysis_from_pdbs(
    base_root_dir="/home/user/valleyfever/dataset/protein_db/md/CIMG_00533/simulation_explicit",
    protein="CIMG_00533",
    pocket_id=3
)
