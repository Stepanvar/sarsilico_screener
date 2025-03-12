# app/utils.py
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from MDAnalysis import Universe
import logging
import json
from django.conf import settings

logger = logging.getLogger(__name__)

def generate_fingerprint(smiles: str) -> dict:
    """Generate Morgan fingerprint with validation (Section 3)"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return {'error': 'Invalid SMILES structure'}
            
        fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius=2)
        return {
            'fingerprint': fp.ToBitString(),
            'smiles': smiles,
            'num_bits': fp.GetNumBits()
        }
    except Exception as e:
        logger.error(f"Fingerprint generation failed: {str(e)}")
        return {'error': str(e)}

def validate_pdb_resolution(pdb_content: str) -> dict:
    """PDB validation with resolution check (Section 2)"""
    try:
        # Check resolution from REMARK records
        resolution = None
        for line in pdb_content.split('\n'):
            if line.startswith('REMARK   2 RESOLUTION.'):
                resolution = float(line.split()[-1])
                break
        
        validation_result = {
            'valid': True,
            'resolution': resolution,
            'warnings': []
        }
        
        # Check against maximum allowed resolution
        if resolution and resolution > settings.MAX_PDB_RESOLUTION:
            validation_result.update({
                'valid': False,
                'error': f"Resolution {resolution}Å exceeds limit"
            })
        
        # Basic structure validation
        mol = Chem.MolFromPDBBlock(pdb_content)
        if not mol or mol.GetNumAtoms() == 0:
            validation_result.update({
                'valid': False,
                'error': "Invalid PDB structure"
            })
            
        return validation_result
        
    except Exception as e:
        logger.error(f"PDB validation failed: {str(e)}")
        return {'valid': False, 'error': str(e)}

def prepare_visualization_data(results: dict) -> dict:
    """Generate NGL.js configuration (Section 5)"""
    try:
        # Generate visualization script
        script = f"""
        var stage = new NGL.Stage("viewport");
        Promise.all([
            stage.loadFile("{results['target_pdb']}"),
            stage.loadFile("{results['ligand_pdbqt']}")
        ]).then(function(components) {{
            components[0].addRepresentation("cartoon", {{ 
                color: "residueindex" 
            }});
            components[1].addRepresentation("ball+stick");
            stage.autoView();
        }});
        """
        
        return {
            'ngl_script': script,
            'target_pdb': results['target_pdb'],
            'ligand_structure': results['ligand_pdbqt'],
            'interaction_data': results.get('interaction_map', {})
        }
        
    except Exception as e:
        logger.error(f"Visualization prep failed: {str(e)}")
        return {'error': str(e)}

class ChemistryJSONEncoder(json.JSONEncoder):
    """Custom JSON encoder for chemistry data types"""
    def default(self, obj):
        if isinstance(obj, Chem.Mol):
            return Chem.MolToPDBBlock(obj)
        return super().default(obj)
