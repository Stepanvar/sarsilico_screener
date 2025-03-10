import subprocess
import logging
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import List, Dict, Tuple, Optional
import re

# RDKit imports with version check
try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem, DataStructs, rdMolDescriptors
    RDLogger.DisableLog('rdApp.*')  # Disable RDKit warnings
except ImportError:
    raise RuntimeError("RDKit not installed. Install with: conda install -c conda-forge rdkit")

# Configure logging
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

# Constants
VINA_RESULT_PATTERN = re.compile(r"REMARK VINA RESULT\s+([-\d.]+)")

def run_pepsmi(sequence: str) -> str:
    """Convert peptide sequence to SMILES using PepSMI with validation."""
    if not isinstance(sequence, str) or not sequence.isupper():
        raise ValueError("Invalid sequence format. Use uppercase amino acid letters")
    
    try:
        result = subprocess.run(
            ["pepsmi", sequence],
            capture_output=True,
            text=True,
            check=True,
            timeout=30  # Add timeout
        )
        return result.stdout.strip()
    except subprocess.CalledProcessError as e:
        logger.error(f"PepSMI failed: {e.stderr}")
        raise RuntimeError(f"PepSMI conversion error: {e.stderr}") from e
    except FileNotFoundError:
        raise RuntimeError("PepSMI not found. Install from: [repository_url]")

def validate_smiles(smiles: str) -> Tuple[bool, Optional[str]]:
    """Validate SMILES string with detailed error reporting."""
    if not isinstance(smiles, str) or len(smiles) == 0:
        return False, "Empty or invalid input"
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES syntax"
    
    try:
        Chem.SanitizeMol(mol)
        return True, None
    except Exception as e:
        return False, f"Sanitization error: {str(e)}"

def validate_pdb_file(pdb_file: Path) -> Tuple[bool, Optional[str]]:
    """Validate PDB file structure with efficient parsing."""
    required_sections = {'ATOM', 'END'}
    found_sections = set()
    
    try:
        with pdb_file.open('r') as f:
            for line in f:
                record = line[:6].strip()
                if record in required_sections:
                    found_sections.add(record)
                if found_sections == required_sections:
                    break
                    
        if not required_sections.issubset(found_sections):
            missing = required_sections - found_sections
            return False, f"Missing required sections: {', '.join(missing)}"
            
        return True, None
    except Exception as e:
        logger.error(f"PDB validation failed: {str(e)}")
        return False, f"File processing error: {str(e)}"

def convert_pdb_to_pdbqt(input_pdb: Path, output_pdbqt: Path) -> Path:
    """Convert PDB to PDBQT format with error handling."""
    if not input_pdb.exists():
        raise FileNotFoundError(f"Input PDB file not found: {input_pdb}")
        
    try:
        subprocess.run(
            ["prepare_receptor4.py", "-r", str(input_pdb), "-o", str(output_pdbqt)],
            check=True,
            timeout=120,
            capture_output=True
        )
        return output_pdbqt
    except subprocess.TimeoutExpired:
        raise RuntimeError("PDBQT conversion timed out after 120 seconds")
    except subprocess.CalledProcessError as e:
        logger.error(f"Conversion failed: {e.stderr.decode()}")
        raise

def convert_smiles_to_pdbqt(smiles: str, output_pdbqt: Path) -> Path:
    """Convert SMILES to 3D structure and PDBQT format."""
    valid, msg = validate_smiles(smiles)
    if not valid:
        raise ValueError(f"Invalid SMILES: {msg}")
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        
        # Generate conformation with ETKDGv3
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        AllChem.MMFFOptimizeMolecule(mol)  # More accurate than UFF
        
        with NamedTemporaryFile(suffix=".pdb", delete=True) as tmp:
            Chem.MolToPDBFile(mol, tmp.name)
            subprocess.run(
                ["prepare_ligand4.py", "-l", tmp.name, "-o", str(output_pdbqt)],
                check=True,
                timeout=60,
                capture_output=True
            )
        return output_pdbqt
    except Exception as e:
        logger.error(f"SMILES to PDBQT conversion failed: {str(e)}")
        raise

def parse_vina_output(output_file: Path) -> List[float]:
    """Parse AutoDock Vina results with regex pattern matching."""
    scores = []
    try:
        with output_file.open('r') as f:
            for line in f:
                match = VINA_RESULT_PATTERN.search(line)
                if match:
                    scores.append(float(match.group(1)))
        return sorted(scores)
    except Exception as e:
        logger.error(f"Vina output parsing failed: {str(e)}")
        raise

def calculate_tanimoto_similarity(input_smiles: str, known_drugs: List[Dict]) -> List[Dict]:
    """Calculate Tanimoto similarity scores with batch processing."""
    valid, msg = validate_smiles(input_smiles)
    if not valid:
        raise ValueError(f"Invalid input SMILES: {msg}")
    
    input_mol = Chem.MolFromSmiles(input_smiles)
    input_fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(input_mol, 2)
    
    results = []
    for drug in known_drugs:
        try:
            drug_mol = Chem.MolFromSmiles(drug['smiles'])
            if drug_mol is None:
                logger.warning(f"Invalid drug SMILES: {drug['name']}")
                continue
                
            drug_fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(drug_mol, 2)
            similarity = DataStructs.TanimotoSimilarity(input_fp, drug_fp)
            
            results.append({
                'name': drug['name'],
                'score': round(similarity, 3),
                'smiles': drug['smiles'],
                'description': drug.get('description', '')
            })
        except Exception as e:
            logger.error(f"Error processing {drug['name']}: {str(e)}")
    
    return sorted(results, key=lambda x: x['score'], reverse=True)

def run_alphafold_convert(input_sequence: str, output_dir: Path) -> Path:
    """Run AlphaFold structure prediction with resource monitoring."""
    output_dir.mkdir(exist_ok=True, parents=True)
    
    try:
        # Validate input sequence
        if not re.match("^[ACDEFGHIKLMNPQRSTVWY]+$", input_sequence):
            raise ValueError("Invalid amino acid sequence")
            
        with NamedTemporaryFile(dir=output_dir, suffix=".fasta", delete=False) as tmp:
            tmp.write(f">input_sequence\n{input_sequence}".encode())
            fasta_path = Path(tmp.name)
            
        cmd = [
            "singularity", "run", "--nv",
            "--bind", f"{output_dir}:/data",
            "/path/to/alphafold.sif",
            "--fasta_paths=/data/input.fasta",
            "--output_dir=/data/output",
            "--max_template_date=2025-03-10",
            "--model_preset=monomer",
            "--db_preset=reduced_dbs"
        ]
        
        result = subprocess.run(
            cmd,
            check=True,
            timeout=3600,  # 1 hour timeout
            capture_output=True,
            text=True
        )
        
        # Find best ranked PDB file
        pdb_files = list(output_dir.glob("output/ranked_*.pdb"))
        if not pdb_files:
            raise FileNotFoundError("No AlphaFold output files generated")
            
        return max(pdb_files, key=lambda f: f.stat().st_size)
        
    except subprocess.TimeoutExpired:
        logger.error("AlphaFold prediction timed out after 1 hour")
        raise
    except Exception as e:
        logger.error(f"AlphaFold failed: {str(e)}")
        raise
