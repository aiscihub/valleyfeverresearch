import os
import glob
import pandas as pd
from pathlib import Path

def find_docking_files(base_path):
    return glob.glob(str(base_path / "**" / "docking_pocket*.pdbqt"), recursive=True)

def extract_scores(file_path):
    best_score = None
    all_scores = []
    with open(file_path, "r") as f:
        for line in f:
            if line.startswith("REMARK VINA RESULT:"):
                scores = line.split()
                if len(scores) > 3:
                    current_scores = [float(s) for s in scores[3:]]
                    all_scores.extend(current_scores)
    if all_scores:
        best_score = min(all_scores)
    return best_score, all_scores

def collect_docking_results(docking_files):
    results = []
    for file in docking_files:
        try:
            ligand = Path(file).parent.name
            pocket = file.split("_pocket")[-1].split(".")[0]
            best_score, all_scores = extract_scores(file)
            results.append([ligand, pocket, best_score, all_scores])
        except Exception as e:
            print(f"Error processing {file}: {e}")
    return pd.DataFrame(results, columns=["Ligand", "Pocket", "Best Score (kcal/mol)", "All Scores"])

def filter_top_poses(df, score_threshold=-8.0):
    return df[df["Best Score (kcal/mol)"] <= score_threshold].copy()

def enrich_and_filter_subset(df):
    df["Protein"] = df["Ligand"].str.extract(r'^(CIMG_\d+)', expand=False)
    df["Compound"] = df["Ligand"].str.extract(r'_(Verapamil|Tacrolimus|Beauvericin|Disulfiram|Milbemycin)', expand=False)
    df_sorted = df.sort_values(by=["Protein", "Compound", "Best Score (kcal/mol)"])
    return df_sorted.groupby(["Protein", "Compound"], group_keys=False).apply(
        lambda group: group[group["Best Score (kcal/mol)"] <= group["Best Score (kcal/mol)"].min() + 1.0]
    ).reset_index(drop=True)

def main():
    project_root = Path(__file__).resolve().parents[2]
    docking_dir = project_root / "dataset" / "protein_db" / "docking"
    results_dir = project_root / "dataset" / "protein_db" / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    docking_files = find_docking_files(docking_dir)
    df_all = collect_docking_results(docking_files)

    summary_path = results_dir / "docking_summary.csv"
    filtered_path = results_dir / "docking_summary_filtered_lessthan_8.csv"
    subset_path = results_dir / "plip_candidates_within_1kcal.csv"

    df_all.to_csv(summary_path, index=False)

    df_filtered = filter_top_poses(df_all)
    df_filtered.to_csv(filtered_path, index=False)

    df_subset = enrich_and_filter_subset(df_filtered)
    df_subset.to_csv(subset_path, index=False)

    print("Docking summary written to:")
    print(f"  - All results: {summary_path}")
    print(f"  - Filtered (<= -8.0 kcal/mol): {filtered_path}")
    print(f"  - Top poses (within 1 kcal/mol): {subset_path}")

if __name__ == "__main__":
    main()