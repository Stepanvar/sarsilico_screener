# app/services/conversions.py
"""Molecular Conversion Engine - Core Implementation"""

from pathlib import Path
from xml.dom import ValidationErr
from rdkit import Chem
from rdkit.Chem import MolFromSmiles, MolFromSequence
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO
from io import StringIO
import subprocess
import os
import tempfile
import time

from django.conf import settings
from django.core.exceptions import ImproperlyConfigured

from inSilicoScreening import settings

def jsme_to_smiles(jsme_data: str) -> str:
    """Convert JSME data to canonical SMILES"""
    mol = Chem.MolFromSmiles(jsme_data) or Chem.MolFromMolBlock(jsme_data)
    if not mol:
        raise ValidationErr("Invalid JSME input format")
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)

def peptide_to_smiles(sequence: str) -> str:
    """Convert peptide sequence to SMILES using RDKit's amino acid handling"""
    mol = MolFromSequence(sequence)
    if not mol:
        raise ValidationErr("Invalid peptide sequence")
    return Chem.MolToSmiles(mol)


import logging
logger = logging.getLogger(__name__)


def sequence_to_cif(sequence: str) -> dict:
    """
    Run AlphaFold via Docker for structure prediction.
    This function writes the sequence to a temporary FASTA file,
    calls the Dockerized AlphaFold prediction, and attempts to retrieve
    a generated CIF file.
    
    Returns a dict with keys:
      - success: bool
      - message: str
      - output_dir: str (if successful)
      - cif_data: str (content of first CIF file found)
      - stdout, stderr: logs from the prediction run
    """
    # Check required settings
    data_dir = getattr(settings, "ALPHAFOLD_DATA_DIR", None)
    if not data_dir:
        raise ImproperlyConfigured("ALPHAFOLD_DATA_DIR must be set in settings.")

    # Create temporary FASTA file
    try:
        fasta_fd, fasta_path = tempfile.mkstemp(suffix=".fasta")
        with os.fdopen(fasta_fd, "w") as f:
            f.write(">query\n" + sequence + "\n")
    except Exception as e:
        return {"success": False, "message": f"Failed to write FASTA file: {str(e)}"}

    # Create a unique output directory under MEDIA_ROOT
    output_dir = os.path.join(settings.MEDIA_ROOT, f"alphafold_{int(time.time())}")
    os.makedirs(output_dir, exist_ok=True)

    # Construct the Docker command.
    # (Adjust the command as necessary to match your Docker image/entrypoint.)
    cmd = [
        "python3", "docker/run_docker.py",
        f"--fasta_paths={fasta_path}",
        "--max_template_date=2022-01-01",
        f"--data_dir={data_dir}",
        f"--output_dir={output_dir}"
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        # Look for a CIF file in the output directory (recursively)
        cif_files = list(Path(output_dir).rglob("*.cif"))
        if cif_files:
            cif_file_path = str(cif_files[0])
            with open(cif_file_path, "r") as cf:
                cif_data = cf.read()
            return {
                "success": True,
                "message": "AlphaFold prediction completed successfully",
                "output_dir": output_dir,
                "cif_data": cif_data,
                "stdout": result.stdout,
                "stderr": result.stderr,
            }
        else:
            return {
                "success": False,
                "message": "AlphaFold prediction ran but no CIF file was generated",
                "stdout": result.stdout,
                "stderr": result.stderr,
            }
    except subprocess.CalledProcessError as e:
        logger.error(f"AlphaFold prediction failed: {e.stderr}")
        return {
            "success": False,
            "message": f"AlphaFold prediction failed: {str(e)}",
            "stdout": e.stdout,
            "stderr": e.stderr,
        }
    finally:
        # Clean up the temporary FASTA file
        if os.path.exists(fasta_path):
            os.unlink(fasta_path)

def cif_to_pdb(cif_data: str) -> str:
    """Convert mmCIF to PDB using BioPython's structure parser"""
    parser = MMCIFParser()
    with StringIO(cif_data) as cif_file:
        structure = parser.get_structure('cif_struct', cif_file)
    pdb_io = StringIO()
    PDBIO().set_structure(structure)
    PDBIO().save(pdb_io)
    return pdb_io.getvalue()

def is_short_sequence(sequence: str, threshold=30) -> bool:
    """Determine appropriate conversion path"""
    return len(sequence) <= threshold

def validate_input(data: str, input_type: str) -> bool:
    """Basic format validation"""
    try:
        if input_type == 'SMILES':
            return bool(MolFromSmiles(data))
        elif input_type == 'PDB':
            # Simple header check
            return data.startswith(('HEADER', 'ATOM'))
        elif input_type == 'sequence':
            return data.isalpha()
        return False
    except:
        return False
