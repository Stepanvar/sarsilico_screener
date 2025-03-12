# app/tasks.py
import uuid
from celery import shared_task
from django.core.files.base import ContentFile
from jsonschema import ValidationError
from regex import T
from .models import ScreeningJob, Compound
from .services.conversions import (
    ConversionEngine,
    ValidationTools,
)
from .services.docking import DockingEngine, DockingResultsParser
from .services.similarity_checker import SimilarityChecker
import logging
import subprocess
import json
from django.db import transaction

logger = logging.getLogger(__name__)

# app/tasks.py
@shared_task(bind=True, autoretry_for=(Exception,), retry_kwargs={'max_retries': 3})
def validate_input_task(self, compound_id: str) -> dict:
    """
    Async input validation (section 2)
    Args:
        compound_id: UUID string of Compound instance
    """
    try:
        compound = Compound.objects.get(id=compound_id)
        converter = ConversionEngine()
        result = {}
        isValid = True
        with transaction.atomic():
            if compound.type == 'small_molecule':
                if compound.input_data.startswith('MRV'):
                    result = converter.marvin_to_smiles(compound.input_data)
                else:
                    isValid = ValidationTools.validate_smiles(compound.input_data)
                    if isValid is False:
                        raise ValidationError("Compound is not smiles")
                    
            elif compound.type == 'peptide':
                result = converter.peptide_to_smiles(compound.input_data)
                
            elif compound.type == 'biomolecule':
                result = handle_alphafold_task(str(compound.id))
                return result  # Async chain

            if result['error'] is None:
                compound.validation_status = False
                compound.validation_errors = result
            else:
                compound.validation_status = True
                if compound.type in ['small_molecule', 'peptide']:
                    compound.converted_smiles = result.get('smiles', '')
                
            compound.save()

            if compound.validation_status:
                job = ScreeningJob.objects.get(compound=compound)
                similarity_check_task.delay(int(job.id))

        return result

    except Exception as e:
        logger.error(f"Validation failed: {str(e)}")
        compound.validation_errors = {'error': str(e)}
        compound.save()
        raise self.retry(exc=e)


@shared_task(bind=True)
def similarity_check_task(self, job_id: int) -> dict:
    """
    Drug repurposing workflow (section 3)
    Args:
        job_id: ScreeningJob ID
    """
    job = ScreeningJob.objects.get(id=job_id)
    try:
        job.status = 'SIM_CHECK'
        job.save()
        
        compound = job.compound
        checker = SimilarityChecker()
        
        # Get SMILES based on compound type
        smiles: str = compound.converted_smiles or ''
        if not smiles:
            raise ValueError("Invalid SMILES for similarity check")
        result = checker.check_similarity(smiles)
        
        if 'error' not in result:
            job.similarity_data = result
            job.status = 'DOCKING' if not result['warning'] else 'SIM_CHECK'
        
        job.save()
        return result
        
    except Exception as e:
        logger.error(f"Similarity check failed: {str(e)}")
        job.status = 'FAILED'
        job.save()
        raise

@shared_task(bind=True)
def run_docking_task(self, job_id: int) -> dict:
    """
    AutoDock execution (section 4)
    Args:
        job_id: ScreeningJob ID
    """
    job = ScreeningJob.objects.get(id=job_id)
    try:
        job.status = 'DOCKING'
        job.save()
        
        # Prepare inputs
        target_path = job.target.pdb_file.path
        compound_data = {
            'type': job.compound.type.lower(),
            'data': job.compound.converted_smiles or job.compound.converted_pdb.path
        }
        
        # Run docking
        with DockingEngine(target_path, compound_data) as dock:
            dock.prepare_receptor()
            dock.prepare_ligand()
            raw_results = dock.run_docking()
            
            # Process results
            parser = DockingResultsParser()
            report = parser.generate_report(raw_results)
            
            # Save results
            job.docking_score = raw_results['best_affinity']
            job.ligand_pdbqt.save(
                f'ligand_{job.id}.pdbqt',
                ContentFile(raw_results['poses'][0]['pose'])
            )
            job.complex_structure.save(
                f'complex_{job.id}.pdb',
                ContentFile(report['visualization']['complex_pdb'])
            )
            job.visualization_config = report['visualization']
            job.status = 'COMPLETED'
            job.save()
            
            return report
        
    except Exception as e:
        logger.error(f"Docking failed: {str(e)}")
        job.status = 'FAILED'
        job.save()
        raise

@shared_task(bind=True)
def handle_alphafold_task(self, sequence: str) -> dict:
    """
    AlphaFold prediction (section 2c)
    Args:
        sequence: Biomolecule sequence
    """
    try:
        converter = ConversionEngine()
        result = converter.run_alphafold(sequence)
        
        if 'error' in result:
            return result
            
        # Convert CIF to PDB
        pdb_result = converter.cif_to_pdb(result['cif_data'])
        return {
            'pdb_data': pdb_result['pdb_data'],
            'validation': pdb_result['validation']
        }
        
    except subprocess.TimeoutExpired:
        return {'error': 'AlphaFold processing timed out'}
    except Exception as e:
        logger.error(f"AlphaFold failed: {str(e)}")
        return {'error': str(e)}
