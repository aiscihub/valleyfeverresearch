# Protein and Ligand Preparation for Docking

This doc outlines the standard operating procedure for preparing protein and ligand structures for molecular docking using AutoDock Vina. The final outputs are `.pdbqt` files compatible with docking, and docking results are further processed for interaction analysis.

https://pubchem.ncbi.nlm.nih.gov/compound/milbemycin
---

## Protein Preparation

### Tools
- **PyMOL** for structure cleanup
- **AutoDockTools** (MGLTools) for hydrogen addition, charge assignment, and PDBQT conversion

### Steps
1. **Load and clean protein structure in PyMOL**
```pymol
load protein/CIMG_00780_relaxed.pdb
remove solvent
remove resn UNK
remove not (polymer.protein)
save protein/CIMG_00780_cleaned.pdb
```

2. **Convert to .pdbqt using AutoDockTools GUI**
- Open cleaned PDB file
- Add hydrogens: `Edit → Hydrogens → Add → All Hydrogens`
- Compute Gasteiger charges: `Edit → Charges → Compute Gasteiger`
- Assign AutoDock 4 atom types: `Edit → Atom Types → Assign AD4 Types`
- Save as `.pdbqt`: `File → Save → Write PDBQT` → `CIMG_00780_prepared.pdbqt`

---

## Ligand Preparation

### Tools
- **Open Babel** for file conversion and 3D structure generation
- **MGLTools** for `.pdbqt` conversion and rotatable bond detection

### Steps
1. **Download ligand from PubChem**
- Example: [Verapamil (CID: 2520)](https://pubchem.ncbi.nlm.nih.gov/compound/2520)
- Save as `.sdf` (e.g., `verapamil.sdf`)

2. **Convert to 3D PDB format**
```bash
obabel verapamil.sdf -O verapamil_3d.pdb --gen3d
```

3. **Convert to PDBQT using MGLTools**
```bash
~/Downloads/mgltools_1.5.7*/bin/pythonsh \
~/Downloads/mgltools_1.5.7*/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py \
-l verapamil_3d.pdb \
-o verapamil.pdbqt
```

⚠️ If you see errors like:
```
Sorry, there are no Gasteiger parameters available for atom UNL1:O
```
Try converting via `.mol2`:
```bash
obabel verapamil.sdf -O verapamil_3d.mol2 --gen3d
```

---

## Docking and Postprocessing

- Docking is performed using Vina with pocket-specific config files.
- Outputs: `docking_pocket*.pdbqt`

### Summary & Filtering Pipeline
- Parses docking scores from each result
- Filters poses with scores ≤ –8.0 kcal/mol
- Retains poses within 1.0 kcal/mol of the best per protein–compound
- Saves final candidates for PLIP analysis

Script: `protein_db_util/scripts/summarize_docking_scores.py`

---

##  Complex & Interaction Analysis

Use PLIP to analyze protein–ligand complexes post-docking.

Script: `protein_db_util/data_processing.py`

---

For questions or reproducibility issues, please ensure all dependencies (PyMOL, Open Babel, MGLTools) are correctly installed and configured.