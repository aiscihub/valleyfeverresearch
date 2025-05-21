# 🧬 ValleyFever ABC Transporter Inhibitor Discovery

This repository contains a reproducible, layered framework for identifying druggable binding pockets in fungal ABC transporters involved in antifungal resistance. The pipeline integrates AlphaFold2-based protein modeling, geometry-aware docking, and molecular dynamics simulation filtering.

## 🔬 Step 1: Protein Structure Prediction with AlphaFold2

The first step of this pipeline involves generating high-confidence 3D structures for five ABC transporter proteins from *Coccidioides immitis* using AlphaFold2. These structures serve as the foundation for downstream pocket detection, docking, and dynamic validation.

### 📁 FASTA Files Provided

We provide FASTA sequences for five transporter genes in:
`dataset/protein_db/proteins/fasta/`

Each file is named by its NCBI protein accession (e.g., `XP_001247647.1.fasta`).

#### Option A: Run AlphaFold2 Locally (GPU needed)
1. Set up AlphaFold2 locally: https://github.com/deepmind/alphafold
2. Run for each FASTA file:
```bash
bash run_alphafold.sh -d /path/to/databases \
  -o ./dataset/protein_db/proteins \
  -m monomer \
  --fasta_paths=./dataset/protein_db/proteins/fasta/XP_001247647.1.fasta
  [TO VERIFY]
  ```

#### Option B: Use ColabFold (Google Colab)
Clone the Colab notebook: https://github.com/sokrypton/ColabFold, upload a FASTA file, run predictions, and download 
the .pdb. The run will hit issues on free tier accounts, the paper describes a two-stage run to get around. 

## 🔬 Step 2: Pocket Detection and Docking
1. Predict pockets using Prankweb, then generate docking grid sizes from predicted pocket geometry: