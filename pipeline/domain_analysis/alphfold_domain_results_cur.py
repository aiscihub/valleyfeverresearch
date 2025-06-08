import os
import re
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from itertools import groupby

# === CONFIG ===
GENE_TO_PDB = {
    "CIMG_00533": "CIMG_00533_relaxed.pdb",
    "CIMG_00780": "CIMG_00780_relaxed.pdb",
    "CIMG_01418": "CIMG_01418_relaxed.pdb",
    "CIMG_06197": "CIMG_06197_relaxed.pdb",
    "CIMG_09093": "CIMG_09093_relaxed.pdb",
}

DOMAIN_MAP = {
    "IPR003439": "Nucleotide-Binding Domain (ABC_NBD)",
    "IPR017871": "Conserved Motif Site (ABC_CS)",
    "IPR027417": "P-loop ATPase Domain (P_LOOP_NTPASE)",
    "IPR010929": "Drug Resistance Domain (PDR/CDR)",
    "IPR029481": "N-terminal Exporter Domain (ABC_NTERM)",
    "IPR013525": "Transmembrane Domain (ABC_TM)",
    "IPR036640": "Transmembrane Domain (ABC_TM)",
    "IPR011527": "Transmembrane Domain (ABC_TM)",
}

def load_domain_annotations(domain_csv_path):
    """Load InterProScan domain annotations."""
    if not os.path.exists(domain_csv_path):
        raise FileNotFoundError(f"Domain CSV not found: {domain_csv_path}")
    df = pd.read_csv(domain_csv_path)
    df["Start"] = df["Start"].astype(int)
    df["End"] = df["End"].astype(int)
    return df

def extract_plddt(pdb_path):
    """Extract per-residue pLDDT scores from PDB."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("model", pdb_path)
    plddt_list = []
    for model in structure:
        for chain in model:
            for res in chain:
                if "CA" in res:
                    res_id = res.get_id()[1]
                    plddt = res["CA"].get_bfactor()
                    plddt_list.append((res_id, plddt))
    return np.array(plddt_list, dtype=[("res_id", int), ("plddt", float)])

def compute_domain_plddt(df_domains, gene_to_pdb, pdb_folder):
    """Compute mean pLDDT per domain."""
    records = []
    for gene_id, pdb_file in gene_to_pdb.items():
        pdb_path = os.path.join(pdb_folder, pdb_file)
        if not os.path.exists(pdb_path):
            print(f"Skipping missing PDB: {pdb_path}")
            continue

        print(f"Processing {gene_id}")
        plddt_array = extract_plddt(pdb_path)
        gene_domains = df_domains[df_domains["Gene_ID"] == gene_id]

        for _, row in gene_domains.iterrows():
            start, end = row["Start"], row["End"]
            region = plddt_array[(plddt_array["res_id"] >= start) & (plddt_array["res_id"] <= end)]
            mean_plddt = round(np.mean(region["plddt"]), 2) if len(region) else None

            records.append({
                "Gene": gene_id,
                "Protein Accession": row["Protein_Accession"],
                "Domain Name": row["Domain_Name"],
                "InterPro ID": row["InterPro_ID"],
                "Start": start,
                "End": end,
                "Length": end - start + 1,
                "Mean pLDDT": mean_plddt
            })
    return pd.DataFrame(records)

def merge_overlapping_ranges(domain_string):
    """Merge overlapping or adjacent domains within ±5 residues."""
    if pd.isna(domain_string) or domain_string.strip() == "":
        return ""

    matches = re.findall(r"([\d.]+) \((\d+)-(\d+)\)", domain_string)
    if not matches:
        return domain_string

    domains = [(float(score), int(start), int(end)) for score, start, end in matches]
    domains.sort(key=lambda x: x[1])

    merged = []
    for k, group in groupby(domains, key=lambda x: x[0]):
        group = list(group)
        group.sort(key=lambda x: x[1])
        current = group[0]
        for nxt in group[1:]:
            if nxt[1] <= current[2] + 5:
                current = (current[0], current[1], max(current[2], nxt[2]))
            else:
                merged.append(current)
                current = nxt
        merged.append(current)

    return ", ".join([f"{score:.1f} ({start}-{end})" for score, start, end in merged])

def save_csv(df, output_path):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    df.to_csv(output_path, index=False)
    print(f"Saved: {output_path}")

def main(domain_csv_path, pdb_dir, output_dir):
    print("Starting domain-level pLDDT processing...")

    # Step 1: Load domain data
    df_domains = load_domain_annotations(domain_csv_path)

    # Step 2: Compute mean pLDDT per domain
    df_output = compute_domain_plddt(df_domains, GENE_TO_PDB, pdb_dir)
    df_output.sort_values(by=["Gene", "Start"], inplace=True)

    save_csv(df_output, os.path.join(output_dir, "domain_plddt_summary.csv"))

    # Step 3: Focus on core domains
    df_core = df_output[df_output["InterPro ID"].isin(DOMAIN_MAP)].copy()
    df_core["Domain_Label"] = df_core["InterPro ID"].map(DOMAIN_MAP)
    df_core["Formatted"] = df_core.apply(
        lambda row: "–" if pd.isna(row["Mean pLDDT"]) else f"{row['Mean pLDDT']:.1f} ({row['Start']}-{row['End']})",
        axis=1
    )

    df_core["pLDDT_Numeric"] = df_core["Mean pLDDT"].fillna(0)
    df_core.sort_values(by=["Gene", "Domain_Label", "Length", "pLDDT_Numeric"],
                        ascending=[True, True, False, False], inplace=True)

    save_csv(df_core, os.path.join(output_dir, "core_domains_long.csv"))

    # Step 4: Group by domain type
    df_summary = (
        df_core.groupby(["Gene", "Domain_Label"])["Formatted"]
        .apply(lambda x: ", ".join(x))
        .unstack()
        .reset_index()
    )

    # Reorder columns
    columns_order = [
        "Gene",
        "Nucleotide-Binding Domain (ABC_NBD)",
        "Conserved Motif Site (ABC_CS)",
        "P-loop ATPase Domain (P_LOOP_NTPASE)",
        "Transmembrane Domain (ABC_TM)",
        "Drug Resistance Domain (PDR/CDR)",
        "N-terminal Exporter Domain (ABC_NTERM)"
    ]
    df_summary = df_summary.reindex(columns=columns_order)

    # Step 5: Merge overlapping ranges
    for col in df_summary.columns[1:]:
        df_summary[col] = df_summary[col].apply(merge_overlapping_ranges)

    save_csv(df_summary, os.path.join(output_dir, "core_domains_plddt_table.csv"))

    print("Completed domain-level pLDDT processing.")

if __name__ == "__main__":
    # Example usage (paths can be passed differently if turned into CLI later)
    DOMAIN_CSV_PATH = os.path.join(BASE_PATH, "results", "parsed_interproscan_domains.csv")
    PDB_DIR = os.path.join(BASE_PATH, "proteins", "pdb")
    OUTPUT_DIR = os.path.join(BASE_PATH, "results")
    main(DOMAIN_CSV_PATH, PDB_DIR, OUTPUT_DIR)
