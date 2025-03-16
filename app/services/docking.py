# app/services/docking.py
"""
Complete molecular docking implementation with integrated functionality
"""

import logging
import os
import subprocess
import tempfile
import shutil
from pathlib import Path
from typing import Dict, List
from django.conf import settings
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom, rdForceFieldHelpers
from MDAnalysis import Universe
from app.models import PredefinedProteinTarget, ScreeningJob

logger = logging.getLogger(__name__)

def execute_docking(target: PredefinedProteinTarget, compound_data: str, input_type: str) -> dict:
    """
    Run molecular docking using AutoDock Vina.
    
    This function prepares the receptor and ligand files based on the provided
    ScreeningJob and compound_data (which may be SMILES or a PDB string), configures
    docking parameters, executes Vina, and parses the results.
    
    Returns a dict with keys:
      - success: bool
      - message: str
      - best_affinity: float (if docking is successful)
      - output_file: str (path to the docking output file)
      - stdout, stderr: logs from the docking run
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdDistGeom, rdForceFieldHelpers
    # Importing MDAnalysis for docking configuration
    try:
        from MDAnalysis import Universe
    except ImportError:
        logger.error("MDAnalysis is required for docking configuration.")
        return {"success": False, "message": "MDAnalysis package not found."}

    # Create a temporary working directory under MEDIA_ROOT
    work_dir = Path(tempfile.mkdtemp(dir=settings.MEDIA_ROOT))
    try:
        # Validate target structure file existence
        target_pdb_path = work_dir / "pdb_files" / target.pdb_file.path
        if not target_pdb_path.exists():
            raise ValueError(f"Target PDB file not found: {target.name}")

        # --- Prepare Receptor ---
        receptor_pdbqt = work_dir / "receptor.pdbqt"
        try:
            subprocess.run([
                "prepare_receptor4.py",
                "-r", str(target_pdb_path),
                "-o", str(receptor_pdbqt),
                "-A", "checkhydrogens"
            ], check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            logger.error(f"Receptor preparation failed: {e.stdout.decode()}")
            raise ValueError("Receptor preparation error") from e


        # --- Prepare Ligand ---
        ligand_pdbqt = work_dir / "ligand.pdbqt"
        if input_type.upper() == "SMILES":
            try:
                # Generate 3D structure from SMILES using RDKit
                mol = Chem.MolFromSmiles(compound_data)
                if not mol:
                    raise ValueError("Invalid SMILES string provided.")
                mol = Chem.AddHs(mol)
                rdDistGeom.EmbedMolecule(mol)
                rdForceFieldHelpers.MMFFOptimizeMolecule(mol)

                # Write temporary PDB file
                temp_pdb = work_dir / "ligand.pdb"
                Chem.MolToPDBFile(mol, str(temp_pdb))

                # Convert PDB to PDBQT using prepare_ligand4.py
                subprocess.run([
                    "prepare_ligand4.py",
                    "-l", str(temp_pdb),
                    "-o", str(ligand_pdbqt),
                    "-A", "hydrogens"
                ], check=True)
            except Exception as e:
                logger.error(f"Ligand SMILES processing failed: {str(e)}")
                raise ValueError("Ligand preparation error (SMILES)") from e

        elif input_type.upper() == "PDB":
            try:
                subprocess.run([
                    "prepare_ligand4.py",
                    "-l", compound_data,
                    "-o", str(ligand_pdbqt),
                    "-A", "hydrogens"
                ], check=True)
            except subprocess.CalledProcessError as e:
                logger.error(f"Ligand PDB conversion failed: {e.stdout.decode()}")
                raise ValueError("Ligand preparation error (PDB)") from e
        else:
            raise ValueError(f"Unsupported input type: {input_type}")

        # --- Configure Docking ---
        try:
            # Use MDAnalysis to compute the centroid of the protein
            u = Universe(str(target_pdb_path))
            protein = u.select_atoms("protein")
            centroid = protein.centroid()
            config = {
                "center": centroid.tolist(),
                "box_size": settings.VINA_BOX_SIZE,
                "exhaustiveness": settings.VINA_EXHAUSTIVENESS,
                "cpu": settings.VINA_CPU_CORES,
                "num_poses": settings.VINA_NUM_POSES
            }
        except Exception as e:
            logger.error(f"Docking configuration failed: {str(e)}")
            raise ValueError("Docking configuration error") from e

        # --- Execute Docking with AutoDock Vina ---
        output_pdbqt = work_dir / "docking_out.pdbqt"
        try:
            vina_cmd = [
                "vina",
                "--receptor", str(receptor_pdbqt),
                "--ligand", str(ligand_pdbqt),
                "--center_x", str(config["center"][0]),
                "--center_y", str(config["center"][1]),
                "--center_z", str(config["center"][2]),
                "--size_x", str(config["box_size"][0]),
                "--size_y", str(config["box_size"][1]),
                "--size_z", str(config["box_size"][2]),
                "--exhaustiveness", str(config["exhaustiveness"]),
                "--cpu", str(config["cpu"]),
                "--num_modes", str(config["num_poses"]),
                "--out", str(output_pdbqt)
            ]
            result = subprocess.run(vina_cmd, capture_output=True, text=True, check=True)
        except subprocess.CalledProcessError as e:
            logger.error(f"Docking execution failed: {e.stderr}")
            raise ValueError("Docking execution error") from e

        # --- Parse Vina Output ---
        best_affinity = None
        for line in result.stdout.splitlines():
            if line.startswith("REMARK VINA RESULT:"):
                parts = line.split()
                if len(parts) >= 4:
                    current_affinity = float(parts[3])
                    if best_affinity is None or current_affinity < best_affinity:
                        best_affinity = current_affinity

        if best_affinity is None:
            raise ValueError("No valid docking results found in Vina output.")

        # Optionally, copy the output file to a permanent location
        final_output = os.path.join(settings.MEDIA_ROOT, f"docking_result_{os.path.basename(str(target.name))}.pdbqt")
        subprocess.run(f"cp {output_pdbqt} {final_output}", shell=True, check=True)

        return {
            "success": True,
            "message": "Docking completed successfully",
            "best_affinity": best_affinity,
            "output_file": final_output,
            "stdout": result.stdout,
        }
    except Exception as exc:
        return {"success": False, "message": str(exc)}
    finally:
        # Clean up working directory unless in DEBUG mode
        if not settings.DEBUG:
            shutil.rmtree(work_dir, ignore_errors=True)