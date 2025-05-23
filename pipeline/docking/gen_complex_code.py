import os
import pandas as pd

# Base project root
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
base_path = os.path.join(project_root, "dataset/protein_db")

# Load docking summary
input_csv = os.path.join(base_path, "results", "docking_summary.csv")
df = pd.read_csv(input_csv)
df["Protein"] = df["Ligand"].str.extract(r'^(CIMG_\d+)', expand=False)
df["Compound"] = df["Ligand"].str.extract(r'_(Verapamil|Milbemycin|Disulfiram|Tacrolimus|Beauvericin)', expand=False)

# Set to record all needed complex directories
complex_dirs = set()

def format_complex_script(row):
    protein = row["Protein"]
    compound = row["Compound"]
    pocket = int(row["Pocket"])

    # Paths
    protein_file = f"{base_path}/docking/{protein}/protein/{protein}_cleaned.pdb"
    ligand_file = (f"{base_path}/docking/{protein}/docking_dynamic/"
                   f"{protein}_prepared_{compound}/docking_pocket{pocket}.pdbqt")
    output_dir = f"{base_path}/docking/{protein}/complex"
    output_file = f"{output_dir}/{protein}_prepared_{compound}_pocket{pocket}_complex.pdb"

    # Record needed complex folder
    complex_dirs.add(output_dir)

    return f"""load {protein_file}, protein
load {ligand_file}, ligand
create complex, protein or ligand
save {output_file}, complex
delete all
"""

# Generate PyMOL script content
script_blocks = df.apply(format_complex_script, axis=1).tolist()
pml_script = "\n".join(script_blocks)

# Create all missing complex directories before running docking
for dir_path in complex_dirs:
    os.makedirs(dir_path, exist_ok=True)

# Output path
output_script = os.path.join(".", "generate_all_complexes.pml")
os.makedirs(os.path.dirname(output_script), exist_ok=True)

with open(output_script, "w") as f:
    f.write(pml_script)

print(f"✅ PyMOL script generated at: {output_script}")
