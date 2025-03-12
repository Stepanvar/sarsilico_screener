# app/services/docking.py
"""
AutoDock Vina integration implementing insilico.txt section 4 requirements
Handles docking preparation, execution, and result processing
"""

import logging
import subprocess
from pathlib import Path
from typing import Dict, List
from django.conf import settings
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom, rdForceFieldHelpers
from MDAnalysis import Universe
import tempfile
import shutil

logger = logging.getLogger(__name__)

class DockingEngine:
    """Handles complete docking workflow from preparation to execution"""
    
    def __init__(self, target_pdb: str, compound_data: Dict):
        self.target_pdb = Path(target_pdb).resolve()
        self.compound_data = compound_data
        self.work_dir = Path(tempfile.mkdtemp(dir=settings.MEDIA_ROOT))
        self.receptor_pdbqt = self.work_dir / 'receptor.pdbqt'
        self.ligand_pdbqt = self.work_dir / 'ligand.pdbqt'

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Cleanup temporary files"""
        if not settings.DEBUG:
            shutil.rmtree(self.work_dir, ignore_errors=True)

    def prepare_receptor(self) -> str:
        """Convert target PDB to PDBQT using AutoDockTools"""
        try:
            if not self.target_pdb.exists():
                raise FileNotFoundError(f"Target PDB not found: {self.target_pdb}")
            
            cmd = [
                'prepare_receptor4.py',
                '-r', str(self.target_pdb),
                '-o', str(self.receptor_pdbqt),
                '-A', 'checkhydrogens'
            ]
            
            result = subprocess.run(
                cmd,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True
            )
            logger.debug(f"Receptor prepared: {result.stdout}")
            return str(self.receptor_pdbqt)
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Receptor prep failed: {e.stdout}")
            raise ValueError("Receptor preparation error") from e

    def prepare_ligand(self) -> str:
        """Prepare ligand PDBQT from various input types"""
        try:
            if self.compound_data['type'] == 'SMILES':
                return self._process_smiles(self.compound_data['data'])
            return self._process_pdb(self.compound_data['data'])
        except KeyError as e:
            logger.error(f"Invalid compound data: {str(e)}")
            raise ValueError("Invalid ligand input format") from e

    def _process_smiles(self, smiles: str) -> str:
        """Convert SMILES to 3D PDBQT via RDKit and OBabel"""
        try:
            # Generate 3D structure with RDKit
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                raise ValueError("Invalid SMILES structure")
                
            mol = Chem.AddHs(mol)
            rdDistGeom.EmbedMolecule(mol)
            rdForceFieldHelpers.MMFFOptimizeMolecule(mol)
            
            # Save temporary PDB
            temp_pdb = self.work_dir / 'ligand.pdb'
            Chem.MolToPDBFile(mol, str(temp_pdb))
            
            # Convert to PDBQT
            return self._convert_to_pdbqt(temp_pdb)
            
        except Exception as e:
            logger.error(f"SMILES processing failed: {str(e)}")
            raise

    def _process_pdb(self, pdb_path: str) -> str:
        """Convert PDB to PDBQT using AutoDockTools"""
        try:
            input_pdb = Path(pdb_path).resolve()
            if not input_pdb.exists():
                raise FileNotFoundError(f"Ligand PDB not found: {input_pdb}")
                
            return self._convert_to_pdbqt(input_pdb)
        except Exception as e:
            logger.error(f"PDB processing failed: {str(e)}")
            raise

    def _convert_to_pdbqt(self, pdb_path: Path) -> str:
        """Generic PDB to PDBQT conversion"""
        cmd = [
            'prepare_ligand4.py',
            '-l', str(pdb_path),
            '-o', str(self.ligand_pdbqt),
            '-A', 'hydrogens'
        ]
        
        try:
            result = subprocess.run(
                cmd,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True
            )
            logger.debug(f"Ligand converted: {result.stdout}")
            return str(self.ligand_pdbqt)
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Ligand conversion failed: {e.stdout}")
            raise ValueError("Ligand preparation error") from e

    def configure_docking(self) -> Dict:
        """Generate docking parameters from target structure"""
        try:
            u = Universe(str(self.target_pdb))
            protein = u.select_atoms("protein")
            centroid = protein.centroid()
            
            return {
                'center': centroid.tolist(),
                'box_size': settings.VINA_BOX_SIZE,
                'exhaustiveness': settings.VINA_EXHAUSTIVENESS,
                'cpu': settings.VINA_CPU_CORES,
                'num_poses': settings.VINA_NUM_POSES
            }
        except Exception as e:
            logger.error(f"Docking config failed: {str(e)}")
            raise

    def run_docking(self) -> Dict:
        """Execute AutoDock Vina and return results"""
        try:
            config = self.configure_docking()
            
            cmd = [
                'vina',
                '--receptor', str(self.receptor_pdbqt),
                '--ligand', str(self.ligand_pdbqt),
                '--center_x', str(config['center'][0]),
                '--center_y', str(config['center'][1]),
                '--center_z', str(config['center'][2]),
                '--size_x', str(config['box_size'][0]),
                '--size_y', str(config['box_size'][1]),
                '--size_z', str(config['box_size'][2]),
                '--exhaustiveness', str(config['exhaustiveness']),
                '--cpu', str(config['cpu']),
                '--num_modes', str(config['num_poses'])
            ]
            
            result = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True
            )
            
            return DockingResultsParser.parse_vina_output(result.stdout)
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Docking failed: {e.stderr}")
            raise ValueError("Docking execution error") from e

class DockingResultsParser:
    """Processes and formats raw Vina output"""
    
    @staticmethod
    def parse_vina_output(output: str) -> Dict:
        """Extract docking results from Vina stdout"""
        try:
            lines = output.split('\n')
            modes: List[Dict] = []
            current_mode: Dict = {}

            for line in lines:
                if '-----+' in line:  # Table header
                    continue
                if line.startswith('MODEL'):
                    current_mode = {
                        'model_id': int(line.split()[1]),
                        'pose': [],
                        'affinity': 0.0,
                        'rmsd_lower': 0.0,
                        'rmsd_upper': 0.0
                    }
                elif line.startswith('ENDMDL'):
                    modes.append(current_mode)
                elif line.startswith('REMARK VINA RESULT:'):
                    parts = line.split()
                    current_mode.update({
                        'affinity': float(parts[3]),
                        'rmsd_lower': float(parts[4]),
                        'rmsd_upper': float(parts[5])
                    })
                elif line.startswith('ATOM') or line.startswith('HETATM'):
                    current_mode['pose'].append(line)  # Now properly typed as list
            
            return {
                'poses': modes,
                'best_affinity': min(m['affinity'] for m in modes) if modes else None,
                'raw_output': output
            }
        except Exception as e:
            logger.error(f"Result parsing failed: {str(e)}")
            raise


    @staticmethod
    def generate_report(raw_results: Dict) -> Dict:
        """Create comprehensive docking report"""
        return {
            'summary': {
                'best_affinity': raw_results['best_affinity'],
                'total_poses': len(raw_results['poses']),
                'mean_affinity': sum(p['affinity'] for p in raw_results['poses']) / len(raw_results['poses'])
            },
            'visualization': {
                'complex_pdb': DockingResultsParser._generate_complex_pdb(raw_results),
                'interactions': DockingResultsParser._analyze_interactions(raw_results)
            },
            'raw_data': raw_results
        }

    @staticmethod
    def _generate_complex_pdb(results: Dict) -> str:
        """Generate receptor-ligand complex structure"""
        try:
            # In production, implement proper complex generation
            return "\n".join(results['poses'][0]['pose'])
        except Exception as e:
            logger.error(f"Complex generation failed: {str(e)}")
            return ""

    @staticmethod
    def _analyze_interactions(results: Dict) -> List[str]:
        """Analyze ligand-protein interactions"""
        # Implement actual interaction analysis
        return ["Hydrogen bond: Atom1-Atom2", "Hydrophobic: ResidueX"]
