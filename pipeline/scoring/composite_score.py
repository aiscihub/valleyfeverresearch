#!/usr/bin/env python3
"""
Two-layer composite scoring across all ABC transporter proteins
-----------------------------------------------------------------
- Layer 1: Pocket composite (mean docking, probability, plip count, ABC adjustment)
- Layer 2: Ligand-pocket composite (energy gate, plip per ligand, rescue pockets)
"""

from pathlib import Path
from scipy.spatial.distance import cdist
import os

from pocket_util import grid_edge_from_sas, load_pdb_coords, \
    filter_tmd_domains, filter_nbd_domains, compute_domain_centers, min_dist

# Root project folder is two levels up from this script
PROJECT_ROOT = Path(__file__).resolve().parents[2]

# Input Paths
DOCKING_CSV = PROJECT_ROOT / "dataset" / "protein_db" / "results" / "docking_summary.csv"
INTERPRO_CSV = PROJECT_ROOT / "dataset" / "protein_db" / "results" / "parsed_interproscan_domains.csv"
BASE_DOCK_PATH = PROJECT_ROOT / "dataset" / "protein_db" / "docking"

# Output Path
OUT_BASE = PROJECT_ROOT / "dataset" / "protein_db" / "results" / "compositescore"
TOP_K_POCKETS = 5
LP_THRES = 0.7
TOP_N_LIGANDS = 3
ENERGY_GATE = -8.0
# ───────────────────────────────────────────────────

# Global result holders
all_pocket_composite = []
all_ligand_pocket_composite = []
all_lp_for_md = []
all_pockets_for_md = []

import numpy as np
import pandas as pd


def is_pocket_near_tmd(pocket_center, tmd_coords, threshold=10.0) -> bool:
    """
    Check if pocket center is within threshold distance to any TMD residue atom.
    """
    distances = cdist([pocket_center], tmd_coords)
    min_distance = distances.min()
    return min_distance <= threshold


def normalize(series: pd.Series) -> pd.Series:
    """Min-max normalize a pandas Series safely."""
    if series.max() == series.min():
        return pd.Series(0, index=series.index)  # avoid division by zero
    return (series - series.min()) / (series.max() - series.min())


def calculate_composite_score_updated(protein_id: str, master_df: pd.DataFrame) -> pd.DataFrame:
    """
    Updated: Calculate a probability-weighted average Z-score per pocket.

    Args:
        protein_id (str): Protein (gene) ID.
        master_df (pd.DataFrame): Master docking dataset (already loaded, corrected).

    Returns:
        pd.DataFrame: DataFrame with Pocket, AvgZScore, probability, and WeightedScore.
    """
    # Filter for the protein
    df_protein = master_df[master_df["ProteinID"] == protein_id].copy()
    if df_protein.empty:
        raise ValueError(f"No docking data found for protein {protein_id}")

    pocket_prob_df = df_protein[["PocketID", "PocketProbability"]].drop_duplicates()
    #pocket_prob_df = pocket_prob_df.rename(columns={"PocketID": "Pocket", "PocketProbability": "probability"})
    # Parse ligand name if needed
    df_protein["Ligand"] = df_protein["Ligand"].apply(lambda x: x.split("_prepared_")[1])

    # Remove invalid docking scores (e.g., 0.0 or > 0)
    df_protein = df_protein[df_protein["DockingScore"] < 0.0]

    # Only compute Z-score if enough data points remain
    if df_protein.empty:
        raise ValueError(f"All docking scores for {protein_id} are invalid!")

    df_protein["ZScore"] = df_protein.groupby("Ligand")["DockingScore"].transform(
        lambda x: (x - x.mean()) / x.std(ddof=0) if len(x) > 1 and x.std(ddof=0) > 0 else 0
    )
    # # Z-score normalize per ligand
    # df_protein["ZScore"] = df_protein.groupby("Ligand")["DockingScore"].transform(
    #     lambda x: (x - x.mean()) / x.std()
    # )

    excluded = master_df[
        (master_df["ProteinID"] == protein_id) & (master_df["DockingScore"] >= 0.0)
        ]
    print(f"Excluded {len(excluded)} bad docking results for {protein_id}")
    # Compute average Z-score per pocket
    avg_z = df_protein.groupby("PocketID")["ZScore"].mean().reset_index(name="AvgZScore")
    #avg_z["Pocket"] = avg_z["PocketID"].astype(str)

    # Load probability info
    #pocket_prob_df["Pocket"] = pocket_prob_df["PocketID"].astype(str)
    merged = avg_z.merge(pocket_prob_df[["PocketID", "PocketProbability"]], on="PocketID", how="left")

    # Compute weighted score: sum(Z * P) / sum(P)
    merged["WeightedZ"] = merged["AvgZScore"] * merged["PocketProbability"]
    total_weighted = merged["WeightedZ"].sum()
    total_prob = merged["PocketProbability"].sum()

    return merged.sort_values("WeightedZ", ascending=False)


def calculate_composite_score(protein_id: str, master_df: pd.DataFrame) -> tuple:
    """
    Calculate CompositeScore for a specific protein using Corrected_PLIP_contacts.

    Args:
        protein_id (str): Protein (gene) ID.
        master_df (pd.DataFrame): Master docking dataset (already loaded, corrected).

    Returns:
        tuple: (PocketComposite_df, PocketsForMD_df, LigandPocketComposite_df, LPForMD_df)
    """

    # === 1. Filter MasterDockingDataset for this protein ===
    df_protein = master_df[master_df["ProteinID"] == protein_id]

    if df_protein.empty:
        print(f"No docking data for protein {protein_id}")
        return None, None, None, None
    df_protein["Ligand"] = df_protein["Ligand"].apply(lambda x: x.split("_prepared_")[1])
    # === 2. Prepare tables ===
    df_poc = df_protein[["PocketID", "PocketProbability", "CenterX", "CenterY", "CenterZ", "ABCadj"]].drop_duplicates()
    df_poc = df_poc.rename(columns={"PocketID": "Pocket", "PocketProbability": "probability"})

    # Step 1: Z-score normalize docking scores for each ligand
    df_dock = df_protein[["Ligand", "PocketID", "DockingScore"]].rename(
        columns={"Ligand": "LigandID", "PocketID": "Pocket", "DockingScore": "BestScore"}
    )

    # Use Corrected_PLIP_contacts here instead of old PLIP_contacts
    plip_counts = df_protein[["Ligand", "PocketID", "Corrected_PLIP_contacts"]].rename(
        columns={"Ligand": "LigandID", "PocketID": "Pocket", "Corrected_PLIP_contacts": "PLIPcnt"}
    )

    # === 3. Pocket-level composite (layer 1) ===
    dock_mean = df_dock.groupby("Pocket")["ZScore"].mean().reset_index(name="MeanDock")
    plip_pocket = plip_counts.groupby("Pocket")["PLIPcnt"].sum().reset_index()

    pockets = (
        dock_mean
        .merge(df_poc, on="Pocket", how="left")
        .merge(plip_pocket, on="Pocket", how="left")
    )

    pockets["NormDock"] = 1 - normalize(pockets["MeanDock"])
    pockets["NormPLIP"] = normalize(pockets["PLIPcnt"])

    pockets["PocketScore"] = (
            0.5 * pockets["NormDock"] +
            0.3 * pockets["probability"] +
            0.1 * pockets["NormPLIP"] +
            pockets["ABCadj"]
    )

    pockets["Protein"] = protein_id

    # === 4. Ligand-Pocket-level composite (layer 2) ===
    lp = (
        df_dock
        .groupby(["LigandID", "Pocket"])["BestScore"]
        .mean()
        .reset_index(name="DockLP")
        .merge(df_poc[["Pocket", "probability", "ABCadj"]], on="Pocket", how="left")
        .merge(plip_counts, on=["LigandID", "Pocket"], how="left")
    )

    lp["NormDock"] = 1 - normalize(lp["DockLP"])
    lp["NormPLIP"] = normalize(lp["PLIPcnt"])

    lp["LPscore"] = (
            0.5 * lp["NormDock"] +
            0.3 * lp["probability"] +
            0.2 * lp["NormPLIP"] +
            lp["ABCadj"]
    )

    lp["Protein"] = protein_id

    # === 5. Select pockets for MD ===
    topK_pockets = set(pockets.nlargest(TOP_K_POCKETS, "PocketScore")["Pocket"])
    rescued_pockets = set(lp.query("LPscore >= @LP_THRES")["Pocket"])
    keep_pockets = topK_pockets | rescued_pockets

    # === 6. Shortlist ligand-pocket pairs ===
    short = (
        lp.query("Pocket in @keep_pockets")
        .sort_values(["Protein", "Pocket", "LPscore"], ascending=[True, True, False])
        .groupby(["Protein", "Pocket"])
        .head(TOP_N_LIGANDS)
    )

    # === 7. Output ===
    pockets = pockets.sort_values("PocketScore", ascending=False)
    lp = lp.sort_values("LPscore", ascending=False)

    return pockets, pd.DataFrame({"Protein": protein_id, "Pocket": sorted(keep_pockets)}), lp, short


import shutil


def collect_raw_files(gene_ids=None):
    if gene_ids is None:
        gene_ids = []
    raw_out_dir = OUT_BASE / "RawDatasetSummary"
    raw_out_dir.mkdir(exist_ok=True)

    # Collect main input files
    shutil.copy(DOCKING_CSV, raw_out_dir / "docking_summary.csv")
    shutil.copy(INTERPRO_CSV, raw_out_dir / "parsed_interproscan_domains.csv")

    # Collect all individual gene-level files
    for gene_id in gene_ids:
        # Docking pocket prediction
        pocket_csv = BASE_DOCK_PATH / gene_id / f"prankweb-{gene_id}_relaxed" / "structure.pdb_predictions.csv"
        if pocket_csv.exists():
            shutil.copy(pocket_csv, raw_out_dir / f"{gene_id}_PocketPredictions.csv")

        # PLIP interaction results
        plip_csv = BASE_DOCK_PATH / gene_id / "complex" / "plip_results" / "all_plip_interactions_summary.csv"
        if plip_csv.exists():
            shutil.copy(plip_csv, raw_out_dir / f"{gene_id}_PLIPSummary.csv")

        # Protein PDB (optional)
        pdb_file = BASE_DOCK_PATH / gene_id / "protein" / f"{gene_id}_relaxed.pdb"
        if pdb_file.exists():
            shutil.copy(pdb_file, raw_out_dir / f"{gene_id}_relaxed.pdb")

    print(f"\n All raw datasets collected under {raw_out_dir}")





def gen_master_dataset():
    RAW_DIR = OUT_BASE / "RawDatasetSummary"
    config_path = RAW_DIR / "vina_config_grid_summary.csv"  # Adjust path if needed

    # 1. Load docking summary
    docking_df = pd.read_csv(RAW_DIR / "docking_summary.csv")
    docking_df = docking_df.rename(columns={'Pocket': 'PocketID', 'Ligand': 'Ligand'})
    docking_df['ProteinID'] = docking_df['Ligand'].apply(lambda x: x.split('_prepared_')[0])

    # 2. Load InterPro domain data (TMD/NBD)
    interpro_df = pd.read_csv(RAW_DIR / "parsed_interproscan_domains.csv")

    config_df = pd.read_csv(config_path)
    config_df["LigandID"] = config_df["ligand"]
    config_df["ProteinID"] = config_df["protein"]
    config_df["PocketID"] = config_df["pocket"]

    # 3. Load Full PLIP interaction data
    full_plip_path = RAW_DIR / "AllPLIPInteractions.csv"
    full_plip_df = pd.read_csv(full_plip_path)
    full_plip_df["PocketID"] = full_plip_df["PocketID"].astype(int)  # safety

    # 4. Prepare empty list to collect master records
    master_records = []

    # 5. Loop through docking_df
    for idx, row in docking_df.iterrows():
        protein_id = row['ProteinID']
        pocket_id = row['PocketID']
        ligand = row['Ligand']
        docking_score = row['Best Score (kcal/mol)'] if 'Best Score (kcal/mol)' in row else row['DockingScore']

        # Load PocketPrediction for this protein
        pocket_pred_file = RAW_DIR / f"{protein_id}_PocketPredictions.csv"
        if pocket_pred_file.exists():
            pocket_df = pd.read_csv(pocket_pred_file)
            pocket_df.columns = pocket_df.columns.str.strip()
            pocket_df['name'] = pocket_df['name'].str.strip()
            pocket_df['PocketID'] = pocket_df['name'].str.extract(r'pocket(\d+)').astype(int)
        else:
            pocket_df = pd.DataFrame()

        # Merge pocket information
        pocket_info = pocket_df[pocket_df['PocketID'] == pocket_id].iloc[0] if not pocket_df.empty and pocket_id in \
                                                                               pocket_df['PocketID'].values else {}

        # Use FullPlipData for PLIP information
        ligand_clean = ligand.replace(f"{protein_id}_prepared_", "")

        plip_subset = full_plip_df[
            (full_plip_df["ProteinID"] == protein_id) &
            (full_plip_df["LigandID"] == ligand_clean) &
            (full_plip_df["PocketID"] == pocket_id)
            ]

        if not plip_subset.empty:
            plip_contact_count = plip_subset["ResidueNumber"].nunique()
            plip_interaction_types = ";".join(sorted(plip_subset["InteractionType"].dropna().unique()))
        else:
            plip_contact_count = 0
            plip_interaction_types = ""

        # Inject real grid size
        grid_row = config_df[
            (config_df["ProteinID"] == protein_id) &
            (config_df["LigandID"] == ligand_clean) &
            (config_df["PocketID"] == pocket_id)
            ]

        if not grid_row.empty:
            grid_edge = grid_row.iloc[0]["size_x"]  # Assume cube
        else:
            grid_edge = None

        # New: Compute ABCadj based on TMD and NBD residue proximity
        abc_adj = 0.0
        domain = 'NA'
        pdb_file = BASE_DOCK_PATH / protein_id / "protein" / f"{protein_id}_relaxed.pdb"
        if pdb_file.exists():
            try:
                ca_coords = load_pdb_coords(pdb_file)
                center_xyz = np.array([
                    pocket_info.get('center_x', 0),
                    pocket_info.get('center_y', 0),
                    pocket_info.get('center_z', 0)
                ])

                # Filter TMD and NBD domains using new helper functions
                tmd_rows = filter_tmd_domains(interpro_df, protein_id)
                nbd_rows = filter_nbd_domains(interpro_df, protein_id)

                # Compute spatial centers
                tmd_centers = compute_domain_centers(tmd_rows, ca_coords)
                nbd_centers = compute_domain_centers(nbd_rows, ca_coords)

                # Apply cutoff-based adjustments
                tmd_adj = 0.05 if min_dist(center_xyz, tmd_centers) <= 10.0 else 0.0
                nbd_adj = 0.05 if min_dist(center_xyz, nbd_centers) <= 10.0 else 0.0

                # Determine domain label
                if tmd_adj > 0:
                    domain = "TMD"
                elif nbd_adj > 0:
                    domain = "NBD"
                else:
                    domain = "NONE"

                abc_adj = tmd_adj + nbd_adj
                # Final adjustment
            except Exception as e:
                print(f"⚠️ Skipping domain proximity check for {protein_id} pocket {pocket_id}: {e}")
                abc_adj = 0.0
            # Final record
            record = {
                'ProteinID': protein_id,
                'PocketID': pocket_id,
                'Ligand': ligand,
                'DockingScore': docking_score,
                'PocketProbability': pocket_info.get('probability', None),
                'CenterX': pocket_info.get('center_x', None),
                'CenterY': pocket_info.get('center_y', None),
                'CenterZ': pocket_info.get('center_z', None),
                'SASpoints': pocket_info.get('sas_points', None),
                'GridEdgeLength': grid_edge,
                'ABCadj': abc_adj,
                'Corrected_PLIP_contacts': plip_contact_count,
                'InteractionTypes': plip_interaction_types,
                'Domain': domain
            }

            master_records.append(record)

    # 6. Create full master DataFrame
    master_df = pd.DataFrame(master_records)

    # 7. Save
    master_df.to_csv(OUT_BASE / "MasterDockingDataset.csv", index=False)

    print(
        "MasterDockingDataset.csv created successfully with corrected PLIP contacts, GridEdgeLength, TMDdistance, and TMDclass!")


def build_full_plip_records(raw_dir: str, output_path: str):
    """
    Build full raw PLIP interaction table from all *_PLIPSummary.csv files.

    Args:
        raw_dir (str): Directory containing PLIPSummary CSVs
        output_path (str): Path to save combined full interaction table
    """

    all_records = []

    for fname in os.listdir(raw_dir):
        if not fname.endswith("_PLIPSummary.csv"):
            continue

        df_plip = pd.read_csv(os.path.join(raw_dir, fname))
        if df_plip.empty:
            continue

        df_plip["Complex"] = df_plip["Complex"].astype(str).str.strip()

        for idx, row in df_plip.iterrows():
            complex_full = row["Complex"]
            type_interaction = row["Type"]
            residue = str(row["Residue"]) if "Residue" in row else None
            distance = row["Distance"] if "Distance" in row else None
            details = row.get("Details", None)

            # Parse complex into ProteinID, LigandID, PocketID
            parts = complex_full.split("_prepared_")
            if len(parts) != 2:
                continue

            protein = parts[0]
            ligand_and_pocket = parts[1]

            if "_pocket" not in ligand_and_pocket:
                continue

            ligand, pocket_id = ligand_and_pocket.split("_pocket")
            pocket_id = int(pocket_id)

            # Generate unique ID
            interaction_id = f"{protein}_{ligand}_{pocket_id:02d}_{idx:04d}"

            all_records.append({
                "InteractionID": interaction_id,
                "ProteinID": protein,
                "LigandID": ligand,
                "PocketID": pocket_id,
                "InteractionType": type_interaction,
                "ResidueName": ''.join([c for c in residue if not c.isdigit()]) if residue else None,
                "ResidueNumber": ''.join([c for c in residue if c.isdigit()]) if residue else None,
                "Distance": distance,
                "AtomDetail": details
            })

    full_df = pd.DataFrame(all_records)
    full_df.to_csv(output_path, index=False)
    print(f"Full raw PLIP interaction table saved to: {output_path}")

    return full_df

def pareto_front(group):
    dominated = []
    for i, row_i in group.iterrows():
        dom = False
        for j, row_j in group.iterrows():
            if (row_j['DockingScore'] <= row_i['DockingScore']) and (row_j['PocketProbability'] >= row_i['PocketProbability']):
                if (row_j['DockingScore'] < row_i['DockingScore']) or (row_j['PocketProbability'] > row_i['PocketProbability']):
                    dom = True
                    break
        dominated.append(dom)
    group['Pareto'] = ~np.array(dominated)
    return group

def layer1_shortlist(out_dir):
    # Load the master dataset
    df = pd.read_csv(out_dir / "Layer1_filtered_df_20.csv")
    df = df.groupby('ProteinID', group_keys=False).apply(pareto_front)

    # select top 3 per protein based on non-dominated; if fewer than 3 pareto, fill with best docking
    shortlists = []
    for prot, sub in df.groupby('ProteinID'):
        pareto = sub[sub['Pareto']]
        if len(pareto) >= 3:
            shortlist = pareto.sort_values(['DockingScore', 'PocketProbability']).head(3)
        else:
            remaining = sub[~sub['Pareto']].sort_values(['DockingScore', 'PocketProbability'])
            shortlist = pd.concat([pareto, remaining]).head(3)
        shortlists.append(shortlist)

    short_df = pd.concat(shortlists).reset_index(drop=True)
    short_df.to_csv(out_dir/"layer1_shortlist_final.csv", index=False)

def layer1_output(out_dir):
    # Load the master dataset
    df = pd.read_csv(out_dir / "MasterDockingDataset.csv")

    # Apply Layer 1 filter:
    # Keep only rows with strong binding and reliable pockets
    filtered_df_and = df[
        (df["DockingScore"] <= -8.0) &
        (df["PocketProbability"] >= 0.5)
        ].copy()

    filtered_df_20 = df[
        (df["DockingScore"] <= -8.0)
       & (df["PocketProbability"] >= 0.2)
        ].copy()

    # Optional: generate a unique identifier for each pocket–ligand–protein combination
    filtered_df_and["Protein_Pocket_Ligand"] = (
            filtered_df_and["ProteinID"].astype(str) + "_" +
            filtered_df_and["PocketID"].astype(str) + "_" +
            filtered_df_and["Ligand"]
    )

    # Optional: generate a unique identifier for each pocket–ligand–protein combination
    filtered_df_20["Protein_Pocket_Ligand"] = (
            filtered_df_20["ProteinID"].astype(str) + "_" +
            filtered_df_20["PocketID"].astype(str) + "_" +
            filtered_df_20["Ligand"]
    )

    # Drop duplicates if needed
    shortlist_df = filtered_df_and.drop_duplicates(subset=["Protein_Pocket_Ligand"])
    filtered_df_20 = filtered_df_20.drop_duplicates(subset=["Protein_Pocket_Ligand"])
    # Save the shortlist for Layer 2 validation
    shortlist_df.to_csv(out_dir/"Layer1_Shortlist_filtered_df_and.csv", index=False)
    filtered_df_20.to_csv(out_dir/"Layer1_filtered_df_20.csv", index=False)

def main():
    out_dir = OUT_BASE
    out_dir.mkdir(exist_ok=True)
    protein_list = ["CIMG_00533", "CIMG_00780", "CIMG_01418", "CIMG_06197", "CIMG_09093"]
    #collect_raw_files(gene_ids=protein_list)
    #build_full_plip_records(raw_dir = f"{OUT_BASE}/RawDatasetSummary", output_path= f"{OUT_BASE}/RawDatasetSummary/AllPLIPInteractions.csv")
    # gen_master_dataset()


   # layer1_output(out_dir=out_dir)
    layer1_shortlist(out_dir=out_dir)
    # Suppose you have a list of proteins to process

    # Load master docking dataset once
    # master_df = pd.read_csv(out_dir / "MasterDockingDataset.csv")

    # all_pocket_composite = []
    # all_pockets_for_md = []
    # for protein_id in protein_list:
    #     try:
    #         # Merge probability into composite calculation
    #         composite_df = calculate_composite_score_updated(protein_id, master_df)
    #
    #         # Select top pockets by weighted Z
    #         top_pockets = composite_df.nlargest(TOP_K_POCKETS, "WeightedZ")[["PocketID"]].copy()
    #         top_pockets["Protein"] = protein_id
    #
    #         # Record all results
    #         composite_df["Protein"] = protein_id
    #         all_pocket_composite.append(composite_df)
    #         all_pockets_for_md.append(top_pockets)
    #
    #     except Exception as e:
    #         print(f"Skipped {protein_id}: {e}")
    #         continue
    #
    # # Save all combined pocket scores
    # pockets_combined = pd.concat(all_pocket_composite).sort_values(by="WeightedZ", ascending=True)
    # pockets_combined.to_csv(out_dir / "PocketComposite_WeightedZ.csv", index=False)

    # Save selected top pockets per protein
    # pockets_for_md_df = pd.concat(all_pockets_for_md).sort_values(["Protein", "PocketID"])
    # pockets_for_md_df.to_csv(out_dir / "PocketComposite_TopK.csv", index=False)
    #
    #
    # pockets_for_md_combined = pockets_combined.sort_values(["Protein", "PocketScore"], ascending=False)
    # pockets_for_md_combined.to_csv(out_dir / "PocketsForMD.csv", index=False)

    # Generate top 5 pockets per protein
#    top5_df = pockets_combined.sort_values(["Protein", "WeightedZ"], ascending=True).groupby("Protein").head(5)

    # Save to new CSV
    # top5_df.to_csv(out_dir / "Top5PocketsPerProtein.csv", index=False)

    # Save LPForMD (sort by Protein, Pocket, LPscore descending)
    # lp_for_md_combined = pd.concat(all_lp_for_md)
    # lp_for_md_combined = lp_for_md_combined.sort_values(["Protein", "Pocket", "LPscore"], ascending=[True, True, False])
    # lp_for_md_combined.to_csv(out_dir / "LPForMD.csv", index=False)
    #
    # # Merge top5 pockets with LP table to get ligand list
    # merged = pd.merge(top5_df[["Protein", "Pocket"]], lp_for_md_combined, on=["Protein", "Pocket"], how="left")
    #
    # # Sort and group
    # merged_sorted = merged.sort_values(["Protein", "LPscore", "Pocket", ], ascending=[True, True, False])
    #
    # merged_sorted.to_csv(out_dir / "Top5PocketsWithLigands.csv", index=False)

    print(f"\n All composite score results sorted and saved under {out_dir}")


if __name__ == "__main__":
    main()
