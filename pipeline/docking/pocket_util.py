import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from Bio.PDB import PDBParser



def grid_edge_from_sas(sas_points, alpha=3.8, beta=9.0, min_edge=20, max_edge=38, step=2):
    """
    Compute grid edge length from SAS points using empirical formula.
    """
    if sas_points <= 0 or pd.isna(sas_points):
        return min_edge
    cube_root = sas_points ** (1/3)
    raw = alpha * cube_root + beta
    snapped = round(raw / step) * step
    return max(min_edge, min(snapped, max_edge))

def load_pdb_coords(pdb_path):
    parser = PDBParser(QUIET=True)
    coords = {}
    structure = parser.get_structure("protein", pdb_path)
    for atom in structure.get_atoms():
        if atom.get_id() == "CA":
            resnum = atom.get_parent().get_id()[1]
            coords[resnum] = atom.get_coord()
    return coords



def filter_tmd_domains(interpro_df, protein_id):
    """
    Filter TMD-related domain annotations from InterProScan results.

    Logic:
    - TMD domains in ABC transporters are typically labeled as 'transmembrane' or 'membrane'
    - This function filters rows where Domain_Name contains 'membrane' or 'transmembrane' (case-insensitive)
    - Also restricts by Gene_ID (protein_id)

    Args:
        interpro_df (pd.DataFrame): Parsed InterProScan annotation table
        protein_id (str): Target protein ID

    Returns:
        pd.DataFrame: Filtered rows corresponding to TMD domains for the given protein
    """
    return interpro_df.query("Gene_ID == @protein_id").copy().loc[
        lambda df: df['Domain_Name'].str.contains("membrane|transmembrane", case=False, na=False)
    ]

def filter_nbd_domains(interpro_df, protein_id):
    """
    Filter NBD-related domain annotations from InterProScan results.

    Logic:
    - ABC transporters have conserved ATPase domains annotated as 'ABC_tran' or 'ABC_trans_N'
    - This function filters rows where Domain_Name contains 'ABC_tran' (covers both common cases)
    - Also restricts by Gene_ID (protein_id)

    Args:
        interpro_df (pd.DataFrame): Parsed InterProScan annotation table
        protein_id (str): Target protein ID

    Returns:
        pd.DataFrame: Filtered rows corresponding to NBD domains for the given protein
    """
    return interpro_df.query("Gene_ID == @protein_id").copy().loc[
        lambda df: df['Domain_Name'].str.contains("ABC_tran", case=False, na=False)
    ]

def compute_domain_centers(domain_rows, coords):
    centers = []
    for start, end in zip(domain_rows["Start"], domain_rows["End"]):
        segment = [coords[r] for r in range(start, end + 1) if r in coords]
        if segment:
            centers.append(np.mean(segment, axis=0))
    return np.array(centers)

def min_dist(pocket_xyz, domain_centers):
    if domain_centers.size == 0:
        return np.inf
    return cdist([pocket_xyz], domain_centers).min()

def compute_abc_adjustment(protein_id, pocket_info, pdb_file_path, interpro_df):
    try:
        ca_coords = load_pdb_coords(pdb_file_path)
        center_xyz = np.array([
            pocket_info.get('center_x', 0),
            pocket_info.get('center_y', 0),
            pocket_info.get('center_z', 0)
        ])

        # Filter domain rows
        tmd_rows = interpro_df.query("Gene_ID == @protein_id and Domain_ID in ['TMhelix','TRANSMEMBRANE']")
        nbd_rows = interpro_df.query("Gene_ID == @protein_id and Domain_ID in ['ABC_NBD','ABC_transporter']")

        # Compute centers
        tmd_centers = compute_domain_centers(tmd_rows, ca_coords)
        nbd_centers = compute_domain_centers(nbd_rows, ca_coords)

        # Apply cutoffs
        tmd_adj = 0.1 if min_dist(center_xyz, tmd_centers) <= 10.0 else 0.0
        nbd_adj = 0.1 if min_dist(center_xyz, nbd_centers) <= 10.0 else 0.0

        abc_adj = tmd_adj + nbd_adj
        return abc_adj

    except Exception as e:
        print(f"Failed to compute domain adjustment for {protein_id}: {e}")
        return 0.0