# app/services/conversions.py
"""
Molecular input conversions implementing insilico.txt section 2 requirements
Handles: SMILES, PDB, peptide, and biomolecule inputs
"""

import subprocess
import re
from typing import Dict
from django.conf import settings
from rdkit import Chem
from rdkit.Chem import AllChem, rdDepictor
from rdkit.Chem import Draw
from chemistry_adapters import AminoAcidAdapter
import logging

logger = logging.getLogger(__name__)
rdDepictor.SetPreferCoordGen(True)  # Better 2D coordinate generation

class ConversionEngine:
    """
    Main conversion router for section 2 inputs
    Implements requirements for:
    - MarvinJS → SMILES (2a)
    - Peptide → SMILES (2b)
    - Biomolecule → PDB via AlphaFold (2c)
    """
    
    def __init__(self):
        self.aa_adapter = AminoAcidAdapter()

    def marvin_to_smiles(self, mrv_data: str) -> Dict:
        """
        Convert MarvinJS sketches to validated SMILES
        Args:
            mrv_data: MRV string from MarvinJS
        Returns: Dict with SMILES and validation data
        """
        try:
            mol = Chem.MolFromMrvBlock(mrv_data)
            if not mol:
                return {'error': 'Invalid MRV structure'}
            
            Chem.SanitizeMol(mol)
            smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
            
            if not ValidationTools.validate_smiles(smiles):
                return {'error': 'Invalid SMILES generated'}
            
            return {
                'smiles': smiles,
                'type': 'small_molecule',
                'structure_2d': self._generate_2d_svg(mol),
                'validation': {'valid': True}
            }
            
        except Exception as e:
            logger.error(f"MarvinJS conversion failed: {str(e)}")
            return {'error': str(e)}

    def peptide_to_smiles(self, sequence: str) -> Dict:
        """
        Convert short peptides (<20aa) to SMILES using PepSMI
        Args:
            sequence: Amino acid sequence
        Returns: Dict with SMILES and validation data
        """
        try:
            clean_seq = re.sub(r'\s+', '', sequence).upper()
            validation = ValidationTools.validate_peptide(clean_seq)
            if not validation['valid']:
                return validation
                
            smiles = self.aa_adapter.convert_amino_acid_sequence_to_smiles(clean_seq)
            
            if not ValidationTools.validate_smiles(smiles):
                return {'error': 'Invalid SMILES from peptide'}
                
            return {
                'smiles': smiles,
                'type': 'peptide',
                'structure_2d': self._generate_2d_svg(Chem.MolFromSmiles(smiles)),
                'validation': validation
            }
            
        except Exception as e:
            logger.error(f"Peptide conversion failed: {str(e)}")
            return {'error': str(e)}

    def run_alphafold(self, sequence: str) -> Dict:
        """
        Generate 3D structure via AlphaFold Docker
        Args:
            sequence: Long biomolecule sequence
        Returns: Dict with CIF data and processing info
        """
        try:
            cmd = [
                'docker', 'run', '--rm',
                '-e', f'AF_KEY={settings.ALPHAFOLD_API_KEY}',
                settings.ALPHAFOLD_DOCKER_IMAGE,
                '--sequence', sequence,
                '--output', '/dev/stdout'
            ]
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=settings.SUBPROCESS_TIMEOUT  
            )
            
            if result.returncode != 0:
                return {'error': f"AlphaFold error: {result.stderr}"}
                
            return {
                'cif_data': result.stdout,
                'job_id': f"AF-{hash(sequence)}",
                'warnings': result.stderr.splitlines()
            }
            
        except subprocess.TimeoutExpired:
            return {'error': 'AlphaFold processing timed out'}
        except Exception as e:
            logger.error(f"AlphaFold failed: {str(e)}")
            return {'error': str(e)}

    def cif_to_pdb(self, cif_data: str) -> Dict:
        """
        Convert CIF to PDB using mmcif2pdb
        Args:
            cif_data: CIF format string
        Returns: Dict with PDB data and validation
        """
        try:
            result = subprocess.run(
                ['mmcif2pdb', '-'],
                input=cif_data,
                capture_output=True,
                text=True,
                check=True,
                timeout=settings.SUBPROCESS_TIMEOUT
            )
            
            validation = ValidationTools.validate_pdb(result.stdout)
            if not validation['valid']:
                return {'error': validation.get('error', 'Invalid PDB')}
            
            return {
                'pdb_data': result.stdout,
                'validation': validation,
                'warnings': result.stderr.splitlines()
            }
            
        except subprocess.CalledProcessError as e:
            return {'error': f"Conversion failed: {e.stderr}"}
        except Exception as e:
            logger.error(f"CIF conversion failed: {str(e)}")
            return {'error': str(e)}

    def _generate_2d_svg(self, mol: Chem.Mol) -> str:
        """Generate 2D structure visualization for UI"""
        Draw.PrepareMolForDrawing(mol)
        return Draw.MolToSVG(mol)

class ValidationTools:
    """
    Input validation methods for section 2 requirements
    Implements checks for:
    - SMILES validity
    - PDB quality
    - Sequence validity
    """
    
    @staticmethod
    def validate_smiles(smiles: str) -> bool:
        """Comprehensive SMILES validation using RDKit"""
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None and mol.GetNumAtoms() > 0

    @staticmethod
    def validate_pdb(pdb_content: str) -> Dict:
        """Validate PDB structure and resolution"""
        validation = {'valid': True, 'warnings': []}
        
        # Resolution check
        for line in pdb_content.split('\n'):
            if line.startswith('REMARK   2 RESOLUTION.'):
                try:
                    res = float(line.split()[-1])
                    if res > settings.MAX_PDB_RESOLUTION:
                        validation.update({
                            'valid': False,
                            'error': f"Resolution {res}Å exceeds limit"
                        })
                except (IndexError, ValueError):
                    validation['warnings'].append("Invalid resolution format")
        
        # Structural validation
        mol = Chem.MolFromPDBBlock(pdb_content)
        if not mol or mol.GetNumAtoms() == 0:
            validation.update({
                'valid': False,
                'error': "Invalid PDB structure"
            })
            
        return validation

    @staticmethod
    def validate_peptide(sequence: str) -> Dict:
        """Validate amino acid sequence composition"""
        if len(sequence) == 0:
            return {'valid': False, 'error': 'Empty sequence'}
            
        if not re.fullmatch(r'^[ACDEFGHIKLMNPQRSTVWY]+$', sequence):
            invalid = set(sequence) - set("ACDEFGHIKLMNPQRSTVWY")
            return {
                'valid': False,
                'error': f"Invalid amino acids: {', '.join(invalid)}"
            }
            
        return {'valid': True, 'sequence': sequence}
