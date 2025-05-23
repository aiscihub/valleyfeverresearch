import logging
import os
import glob
from typing import Tuple

import numpy as np
from sys import stdout
import shutil
from openmm.app import *
from openmm import *
from openmm.unit import *
from openff.toolkit.topology import Molecule
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from openmmforcefields.generators import GAFFTemplateGenerator
import mdtraj as md
from rdkit import Chem

from openmm import CustomExternalForce
from openmm import MonteCarloBarostat


# Wrapper function
def run_multiple_replicas(parent_dir, base_dir, sdf_filename, pdb_filename, n_steps=50000, n_replicas=3):
    for i in range(1, n_replicas+1):
        print(f"🚀 Running Replica {i}...")

        # Create a new subfolder for each replica (optional but nice)
        replica_base_dir = os.path.join(base_dir, f"replica_{i}")
        os.makedirs(replica_base_dir, exist_ok=True)

        shutil.copy(os.path.join(base_dir, sdf_filename), replica_base_dir)
        shutil.copy(os.path.join(base_dir, pdb_filename), replica_base_dir)
        # Call your main MD function with different seeds
        # run_explicit_protein_ligand_md(
        #     parent_dir = base_dir,
        #     base_dir=replica_base_dir,
        #     sdf_filename= sdf_filename,
        #     pdb_filename= pdb_filename,
        #     production_steps=n_steps
        # )
        run_explicit_protein_ligand_md_v2(
            parent_dir=parent_dir,
            base_dir=replica_base_dir,
            sdf_filename=sdf_filename,
            pdb_filename=pdb_filename,
            production_steps=n_steps,
            random_seed=42 *i,
            restraint_k = 2.5
        )

def save_frame_pdb(dcd_path, pdb_path, output_path, keep_ligand=True, frame_index=-1):
    import mdtraj as md

    traj = md.load(dcd_path, top=pdb_path)
    if traj.unitcell_lengths is not None:
        traj.image_molecules(inplace=True)

    protein_atoms = traj.topology.select("protein and backbone")
    if len(protein_atoms) >= 3:
        traj.superpose(traj[0], atom_indices=protein_atoms)

    if keep_ligand:
        selection = traj.topology.select("protein or resname UNK or resname MOL or resname LIG or resname BEA")
    else:
        selection = traj.topology.select("protein")

    stripped = traj.atom_slice(selection)
    stripped[frame_index].save_pdb(output_path)
    print(f"✅ Saved frame {frame_index} to {output_path}")


def run_explicit_protein_ligand_md_v2(
        parent_dir: str,
        base_dir: str,
        sdf_filename: str,
        pdb_filename: str,
        production_steps: int = 50000,
        use_solvated_box: bool = False,
        random_seed: int = 42,
        restraint_k: float = 20.0,
        apply_restrain: bool = False# kcal/mol/Å²
) -> Tuple[str, str]:
    """Run stable explicit solvent MD with NaN protection."""
    print(f"Using random_seed={random_seed}")
    # ===== 1. Initialize Paths =====
    os.makedirs(base_dir, exist_ok=True)
    prefix = os.path.basename(pdb_filename).replace(".pdb", "")

    output_files = {
        'final': os.path.join(base_dir, f"{prefix}_explicit_final_frame.pdb"),
        'initial': os.path.join(base_dir, f"{prefix}_explicit_initial_frame.pdb"),
        'final_stripped': os.path.join(base_dir, f"{prefix}_explicit_stripped_final_frame.pdb"),
        'initial_stripped': os.path.join(base_dir, f"{prefix}_explicit_stripped_initial_frame.pdb"),
        'dcd': os.path.join(base_dir, f"{prefix}_explicit_trajectory.dcd"),
        'solvated_box': os.path.join(parent_dir, f"{prefix}_solvated_box.pdb"),
        'log': os.path.join(base_dir, f"{prefix}_simulation.log"),
        'checkpoint': os.path.join(base_dir, f"{prefix}_checkpoint.chk")
    }

    # ===== 2. Enhanced System Setup =====
    try:
        # Load and validate structures
        orig_complex = PDBFile(os.path.join(base_dir, pdb_filename))
        sdf_path = os.path.join(base_dir, sdf_filename)
        rdmol = Chem.SDMolSupplier(sdf_path, removeHs=False)[0]
        ligand = Molecule.from_rdkit(rdmol, allow_undefined_stereo=True)
        ligand.generate_conformers(n_conformers=1)
        ligand.assign_partial_charges('am1bcc')

        # Build system with proper alignment
        modeller = Modeller(orig_complex.topology, orig_complex.positions)
        modeller.delete([res for res in modeller.topology.residues()
                         if res.name in ["UNL", "MOL", "LIG", "UNK"]])

        # Align ligand to original pocket position
        lig_positions = ligand.conformers[0].to_openmm()
        original_lig_pos = np.mean([orig_complex.positions[i].value_in_unit(nanometers)
                                    for i, a in enumerate(orig_complex.topology.atoms())
                                    if a.residue.name in ["UNL", "MOL", "LIG", "UNK"]],
                                   axis=0) * nanometers
        current_center = np.mean([pos.value_in_unit(nanometers) for pos in lig_positions], axis=0) * nanometers
        shift = original_lig_pos - current_center
        lig_positions = [(pos + shift).value_in_unit(nanometers) * nanometers for pos in lig_positions]

        modeller.add(ligand.to_topology().to_openmm(), lig_positions)

        # Force field with GAFF for ligand
        ff = ForceField('amber14/protein.ff14SB.xml', 'amber14/tip3p.xml')
        gaff = GAFFTemplateGenerator(molecules=[ligand], forcefield='gaff-2.11')
        ff.registerTemplateGenerator(gaff.generator)

        # Add hydrogens and solvate
        modeller.addHydrogens(ff)
        print("Adding solvent...")
        modeller.addSolvent(ff, padding=1.0*nanometer, ionicStrength=0.15*molar)  # Slightly smaller padding

        # Save full solvated topology that exactly matches the DCD
        with open(output_files['solvated_box'], 'w') as f:
            PDBFile.writeFile(modeller.topology, modeller.positions, f)

        # Create system with PME
        system = ff.createSystem(
            modeller.topology,
            nonbondedMethod=PME,
            nonbondedCutoff=1*nanometer,
            constraints=HBonds,
            hydrogenMass=4*amu  # Helps stability
        )
        system.addForce(MonteCarloBarostat(1*atmosphere, 300*kelvin, 25))

        # Add soft-core potentials for stability
        for force in system.getForces():
            if isinstance(force, NonbondedForce):
                force.setUseDispersionCorrection(True)
                force.setUseSwitchingFunction(True)
                force.setSwitchingDistance(0.8*nanometer)

        # Add restraints
        restraint_force = CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2")
        restraint_force.addPerParticleParameter("x0")
        restraint_force.addPerParticleParameter("y0")
        restraint_force.addPerParticleParameter("z0")
        restraint_force.addGlobalParameter("k", restraint_k*kilocalories_per_mole/nanometer**2)

        for atom in modeller.topology.atoms():
            if atom.residue.name in ["UNL", "MOL", "LIG", "UNK"]:
                pos = modeller.positions[atom.index]
                restraint_force.addParticle(atom.index, [pos.x, pos.y, pos.z])
        if apply_restrain:
            print(f"addForce....{restraint_k} = restraint_k")
            system.addForce(restraint_force)
        else:
            print("addForce Skipped.")

    except Exception as e:
        raise RuntimeError(f"System setup failed: {str(e)}")

    # ===== 3. Simulation with Stability Checks =====
    try:
        # More stable integrator settings
        integrator = LangevinMiddleIntegrator(
            300*kelvin,
            1/picosecond,
            0.002*picoseconds
        )
        integrator.setRandomNumberSeed(random_seed)
        integrator.setConstraintTolerance(0.00001)

        platform = Platform.getPlatformByName("CUDA")
        properties = {'CudaPrecision': 'mixed', 'DeterministicForces': 'true'}
        simulation = Simulation(
            modeller.topology,
            system,
            integrator,
            platform,
            properties
        )
        simulation.context.setPositions(modeller.positions)

        # Check for NaN positions
        state = simulation.context.getState(getPositions=True)
        if np.any(np.isnan(state.getPositions(asNumpy=True))):
            raise ValueError("NaN positions detected after setup")

    except Exception as e:
        raise RuntimeError(f"Simulation initialization failed: {str(e)}")

    # ===== 4. Run Simulation with Safety Nets =====
    try:
        # Configure reporters with proper totalSteps
        simulation.reporters = [
            DCDReporter(output_files['dcd'], 10000),  # Save every 10 ps
            CheckpointReporter(output_files['checkpoint'], 10000),
            StateDataReporter(
                output_files['log'],
                1000,
                step=True,
                time=True,
                potentialEnergy=True,
                temperature=True,
                density=True,
                totalSteps=production_steps + 50000  # production + equilibration
            ),
            StateDataReporter(
                stdout,
                5000,
                step=True,
                time=True,
                potentialEnergy=True,
                temperature=True,
                progress=True,
                remainingTime=True,
                totalSteps=production_steps + 50000
            )
        ]

        # Save initial frame
        with open(output_files['initial'], 'w') as f:
            PDBFile.writeFile(modeller.topology, modeller.positions, f)

        # Minimization
        print("Minimizing energy...")
        simulation.minimizeEnergy(maxIterations=1000)

        # Equilibration (50 ps)
        print("Equilibrating...")
        simulation.step(25000)  # 50 ps

        # Production MD
        print(f"Starting production MD ({production_steps} steps ~{production_steps*0.002} ps)")
        simulation.step(production_steps)

        # Save final frame
        state = simulation.context.getState(getPositions=True)
        with open(output_files['final'], 'w') as f:
            PDBFile.writeFile(modeller.topology, state.getPositions(), f)

        # Save stripped final frame
        save_frame_pdb(
            dcd_path=output_files['dcd'],
            pdb_path=output_files['final'],
            output_path=output_files['final_stripped'],
            keep_ligand=True,
            frame_index=-1
        )

        # Save stripped initial frame
        save_frame_pdb(
            dcd_path=output_files['dcd'],
            pdb_path=output_files['initial'],  # still valid topology
            output_path=output_files['initial_stripped'],
            keep_ligand=True,
            frame_index=0
        )


        print(f"Simulation completed successfully. Results in {base_dir}")

    except Exception as e:
        raise RuntimeError(f"Simulation failed at step {simulation.currentStep if 'simulation' in locals() else 'N/A'}: {str(e)}")

    return output_files['final_stripped'], output_files['dcd']
import argparse
def cli_main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--base_dir", required=True)
    parser.add_argument("--sdf_filename", required=True)
    parser.add_argument("--pdb_filename", required=True)
    parser.add_argument("--n_steps", type=int, default=50000)
    args = parser.parse_args()

    run_implicit_protein_ligand_md_v2(
        base_dir=args.base_dir,
        sdf_filename=args.sdf_filename,
        pdb_filename=args.pdb_filename,
        n_steps=args.n_steps
    )

if __name__ == "__main__":
   cli_main()

