"""
Protein-Ligand Preparation Script for Molecular Dynamics (MD)

This script prepares protein-ligand complexes for MD simulations by:
1. Cleaning and repairing the protein using PDBFixer.
2. Extracting the docked ligand pose from the complex PDB.
3. Reconstructing the full ligand structure from the docked fragment using SMILES.
4. Aligning the full ligand to the docked pose using heavy atom substructure matching.
5. Grafting the exact docked coordinates onto the corresponding atoms of the full ligand.
6. Saving the aligned, chemically complete ligand in SDF format.

This approach ensures:
- The ligand is fully parameterizable for MD (includes all atoms and hydrogens).
- The aligned core maintains the correct binding pose as predicted by docking.
- Extra unmatched atoms (not present in the original docking result) are safely retained
  with RDKit-generated 3D coordinates, and do not interfere with MD as long as they are chemically valid.

No original ligand SDF file is required: the ligand is rebuilt and aligned directly from the docked complex.
"""

import os
import re
import glob
from pdbfixer import PDBFixer
from openmm.app import PDBFile

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from openff.toolkit.topology import Molecule
from openff.toolkit.typing.engines.smirnoff import ForceField

def align_full_ligand_to_pose(
        complex_pdb = "complex.pdb",
        aligned_sdf = "Milbemycin_aligned.sdf",
        resname     = "UNL"
):
    # Step 1: Extract docked pose from PDB
    with open(complex_pdb) as fh:
        lig_lines = [
            l for l in fh if l.startswith(("HETATM", "ATOM"))
                             and l[17:20].strip() == resname
        ]
    if not lig_lines:
        raise ValueError("No ligand atoms found in complex")

    pdb_block = "MODEL\n" + "".join(lig_lines) + "ENDMDL\n"
    pose = Chem.MolFromPDBBlock(pdb_block, sanitize=True, removeHs=True)
    print(f"Extracted {pose.GetNumAtoms()} atoms from docked pose")

    # Step 2: Regenerate full ligand from correct SMILES
    smiles = Chem.MolToSmiles(pose, isomericSmiles=True)
    full = Chem.MolFromSmiles(smiles)
    full = Chem.AddHs(full)
    AllChem.EmbedMolecule(full, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(full)
    print("🛠️ Regenerated 3D conformer from pose SMILES")

    # Step 3: Find heavy atom map
    heavy_pose = Chem.RemoveHs(pose)
    heavy_full = Chem.RemoveHs(full)
    match = heavy_full.GetSubstructMatch(heavy_pose, useChirality=True)
    if not match:
        match = heavy_full.GetSubstructMatch(heavy_pose, useChirality=False)
    if not match:
        raise RuntimeError("Docked fragment not found in regenerated full ligand")

    hf_to_full = [i for i, a in enumerate(full.GetAtoms()) if a.GetAtomicNum() > 1]
    match_full = [hf_to_full[i] for i in match]
    print(f"Mapping: {len(match_full)} atoms aligned")

    # Step 4: Align regenerated ligand to docked pose
    rdMolAlign.AlignMol(full, pose, atomMap=list(zip(match_full, range(len(match_full)))))

    # Step 5: Graft exact docked coordinates onto aligned full ligand
    conf_full = full.GetConformer()
    conf_pose = pose.GetConformer()
    for f_idx, p_idx in zip(match_full, range(len(match_full))):
        conf_full.SetAtomPosition(f_idx, conf_pose.GetAtomPosition(p_idx))

    # Step 6: Save final aligned ligand
    Chem.SDWriter(aligned_sdf).write(full)
    print(f"Aligned ligand saved ➜ {aligned_sdf}")


# Customize here
protein = "CIMG_00533"
ligands = ["Verapamil", "Milbemycin", "Disulfiram", "Tacrolimus", "Beauvericin"]
ligands = [ "Verapamil"]
# Regex pattern: Match complex files with a pocket
pattern = re.compile(rf"^{protein}_prepared_[A-Za-z0-9]+_pocket\d+_complex\.pdb$")

# Paths
base_path = f"/home/user/Projects/valleyfever/dataset/protein_db/md/{protein}/simulation_explicit/pocket3"
outbase = f"/home/user/Projects/valleyfever/dataset/protein_db/md/{protein}/simulation_explicit/pocket3"

# Process each ligand
for ligand in ligands:
    docking_files = glob.glob(
        os.path.join(base_path, f"{protein}_prepared_{ligand}_pocket*_complex.pdb"),
        recursive=True
    )

    for file in docking_files:
        filename = os.path.basename(file)
        if pattern.match(filename):
            print(f"Processing: {filename}")

            # Define output name
            base_name = filename.replace(".pdb", "")
            output_file = f"{base_name}_clean.pdb"
            output_path = os.path.join(outbase, output_file)

            if os.path.exists(output_path):
                print(f"Skipped (already exists): {output_file}")
                continue

            # Run PDBFixer
            fixer = PDBFixer(filename=file)
            residues_to_delete = [res for res in fixer.topology.residues() if res.name == "UNL"]
            chains_to_remove = set(res.chain.index for res in residues_to_delete)
            fixer.removeChains(list(chains_to_remove))

            fixer.findMissingResidues()
            fixer.findMissingAtoms()
            fixer.addMissingAtoms()
            fixer.addMissingHydrogens(pH=7.0)

            # Write output
            with open(output_path, "w") as f:
                PDBFile.writeFile(fixer.topology, fixer.positions, f)

            print(f"Saved: {output_file}")
            ligand_pose_file = os.path.join(outbase, f"{base_name}_aligned_ligand.sdf")

            align_full_ligand_to_pose(complex_pdb=file, aligned_sdf = ligand_pose_file)
            print(f"Saved: {ligand_pose_file}")
