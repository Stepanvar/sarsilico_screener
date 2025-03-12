# app/views.py
import json
import logging
from django.views import View
from django.shortcuts import render, redirect
from django.http import JsonResponse, HttpResponseBadRequest
from django.urls import reverse
from django.db import transaction
from django.core.exceptions import ValidationError

from .models import Target, Compound, ScreeningJob
from django.contrib.auth.forms import UserCreationForm
from .forms import SMILESForm, PDBUploadForm, SequenceForm
from .services.conversions import (
    ConversionEngine, 
    ValidationTools,
)
from .services.similarity_checker import SimilarityChecker
from .services.docking import DockingEngine
from .tasks import (
    validate_input_task,
    similarity_check_task,
    run_docking_task,
    handle_alphafold_task
)

logger = logging.getLogger(__name__)

class IndexView(View):
    template_name = 'app/index.html'
    
    def get(self, request):
        return render(request, self.template_name)

class TargetSelectionView(View):
    template_name = 'app/target_selection.html'

    def get(self, request):
        targets = Target.objects.all()
        return render(request, self.template_name, {'targets': targets})

    def post(self, request):
        try:
            target = Target.objects.get(pk=request.POST['target_id'])
            request.session['selected_target'] = str(target.id)
            return redirect('small_molecule')
        except (KeyError, Target.DoesNotExist):
            return render(request, self.template_name, {'error': 'Invalid target selection'})

class SmallMoleculeInputView(View):
    template_name = 'app/input/small_molecule.html'
    form_class = SMILESForm

    def get(self, request):
        if 'selected_target' not in request.session:
            return redirect('index')
        return render(request, self.template_name, {'form': self.form_class()})

    def post(self, request):
        form = self.form_class(request.POST)
        if not form.is_valid():
            return render(request, self.template_name, {'form': form})

        try:
            with transaction.atomic():
                # Handle conversions
                converter = ConversionEngine()
                if form.cleaned_data.get('marvinjs'):
                    result = converter.marvin_to_smiles(form.cleaned_data['marvinjs'])
                    if 'error' in result:
                        raise ValidationError(result['error'])
                    smiles = result['smiles']
                else:
                    smiles = form.cleaned_data['smiles']

                # Create compound and job
                compound = Compound.objects.create(
                    type='small_molecule',
                    data=smiles,
                    validation_result={'status': 'pending'}
                )
                job = ScreeningJob.objects.create(
                    target_id=request.session['selected_target'],
                    compound=compound,
                    status='validation_pending'
                )
                compound_data = {
                    'type': compound.type,
                    'data': compound.input_data,
                    'compound_id': compound.id
                }
                # Start validation task
                validate_input_task.delay(compound_data)

            return redirect('validation_feedback', job_id=job.id)

        except Exception as e:
            logger.error(f"Small molecule input error: {str(e)}")
            return render(request, self.template_name, {'form': form, 'error': str(e)})

class PeptideInputView(View):
    template_name = 'app/input/peptide.html'
    form_class = SequenceForm

    def get(self, request):
        if 'selected_target' not in request.session:
            return redirect('index')
        return render(request, self.template_name, {'form': self.form_class()})

    def post(self, request):
        form = self.form_class(request.POST)
        if not form.is_valid():
            return render(request, self.template_name, {'form': form})

        try:
            converter = ConversionEngine()
            result = converter.peptide_to_smiles(form.cleaned_data['sequence'])
            
            if 'error' in result:
                raise ValidationError(result['error'])

            with transaction.atomic():
                compound = Compound.objects.create(
                    type='peptide',
                    data=form.cleaned_data['sequence'],
                    converted_data=result['smiles'],
                    validation_result={'status': 'pending'}
                )
                job = ScreeningJob.objects.create(
                    target_id=request.session['selected_target'],
                    compound=compound,
                    status='validation_pending'
                )
                validate_input_task.delay(str(compound.id))

            return redirect('validation_feedback', job_id=job.id)

        except Exception as e:
            logger.error(f"Peptide input error: {str(e)}")
            return render(request, self.template_name, {'form': form, 'error': str(e)})

class BiomoleculeInputView(View):
    template_name = 'app/input/biomolecule.html'
    form_class = SequenceForm

    def get(self, request):
        if 'selected_target' not in request.session:
            return redirect('index')
        return render(request, self.template_name, {'form': self.form_class()})

    def post(self, request):
        form = self.form_class(request.POST)
        if not form.is_valid():
            return render(request, self.template_name, {'form': form})

        try:
            converter = ConversionEngine()
            af_result = converter.run_alphafold(form.cleaned_data['sequence'])
            
            if 'error' in af_result:
                raise ValidationError(af_result['error'])

            with transaction.atomic():
                compound = Compound.objects.create(
                    type='biomolecule',
                    data=form.cleaned_data['sequence'],
                    converted_data=af_result['cif_data'],
                    validation_result={'status': 'pending'}
                )
                job = ScreeningJob.objects.create(
                    target_id=request.session['selected_target'],
                    compound=compound,
                    status='alphafold_processing'
                )
                handle_alphafold_task.delay(str(compound.id))

            return redirect('docking_progress', job_id=job.id)

        except Exception as e:
            logger.error(f"Biomolecule input error: {str(e)}")
            return render(request, self.template_name, {'form': form, 'error': str(e)})

class ValidationFeedbackView(View):
    template_name = 'app/validation_feedback.html'

    def get(self, request, job_id):
        try:
            job = ScreeningJob.objects.get(id=job_id)
            return render(request, self.template_name, {
                'job': job,
                'compound': job.compound,
                'target': job.target
            })
        except ScreeningJob.DoesNotExist:
            return render(request, self.template_name, {'error': 'Invalid job ID'})

class SimilarityCheckView(View):
    template_name = 'app/similarity_check.html'

    def get(self, request, job_id):
        try:
            job = ScreeningJob.objects.get(id=job_id)
            if job.status != 'similarity_check':
                return redirect('validation_feedback', job_id=job_id)
            
            return render(request, self.template_name, {
                'similarity_data': job.similarity_data,
                'compound': job.compound
            })
        except ScreeningJob.DoesNotExist:
            return render(request, self.template_name, {'error': 'Invalid job ID'})

class DockingProgressView(View):
    template_name = 'app/docking_progress.html'

    def get(self, request, job_id):
        try:
            job = ScreeningJob.objects.get(id=job_id)
            return render(request, self.template_name, {
                'job': job,
                'target': job.target
            })
        except ScreeningJob.DoesNotExist:
            return render(request, self.template_name, {'error': 'Invalid job ID'})

class ResultsVisualizationView(View):
    template_name = 'app/results/visualization.html'

    def get(self, request, job_id):
        try:
            job = ScreeningJob.objects.get(id=job_id)
            return render(request, self.template_name, {
                'job': job,
                'visualization_config': self._prepare_ngl_config(job)
            })
        except ScreeningJob.DoesNotExist:
            return render(request, self.template_name, {'error': 'Invalid job ID'})

    def _prepare_ngl_config(self, job):
        return {
            'target_pdb': job.target.pdb_file.url,
            'ligand_pdb': job.compound.converted_pdb.url,
            'interactions': job.docking_results.get('interactions', [])
        }

class ResultsSummaryView(View):
    template_name = 'app/results/summary.html'

    def get(self, request, job_id):
        try:
            job = ScreeningJob.objects.get(id=job_id)
            return render(request, self.template_name, {
                'job': job,
                'summary_data': self._prepare_summary(job)
            })
        except ScreeningJob.DoesNotExist:
            return render(request, self.template_name, {'error': 'Invalid job ID'})

    def _prepare_summary(self, job):
        return {
            'affinity_scores': job.docking_results.get('scores', []),
            'similarity_warnings': job.similarity_data,
            'validation_report': job.compound.validation_result
        }

# API Views
class ValidatePDBAPI(View):
    def post(self, request):
        form = PDBUploadForm(request.POST, request.FILES)
        if not form.is_valid():
            return HttpResponseBadRequest(json.dumps(form.errors))
        
        try:
            validator = ValidationTools()
            result = validator.validate_pdb(request.FILES['pdb_file'].read().decode())
            return JsonResponse(result)
        except Exception as e:
            logger.error(f"PDB validation error: {str(e)}")
            return HttpResponseBadRequest(json.dumps({'error': str(e)}))

class ConvertSequenceAPI(View):
    def post(self, request):
        data = json.loads(request.body)
        sequence = data.get('sequence', '')
        seq_type = data.get('type', 'peptide')
        
        try:
            converter = ConversionEngine()
            if seq_type == 'peptide':
                result = converter.peptide_to_smiles(sequence)
            else:
                result = converter.run_alphafold(sequence)
            
            if 'error' in result:
                return JsonResponse({'error': result['error']}, status=400)
                
            return JsonResponse(result)
        except Exception as e:
            logger.error(f"Sequence conversion error: {str(e)}")
            return JsonResponse({'error': str(e)}, status=500)

def register(request):
    if request.method == 'POST':
        form = UserCreationForm(request.POST)
        if form.is_valid():
            form.save()
        else:
            form = UserCreationForm()
    return render(request, 'app/index.html')


# views.py
from django.shortcuts import render
from django.http import HttpResponse
from .forms import ScreeningForm

# Import view model functions (conversion, validation, docking, etc.)
from .view_models import (
    marvin_to_smiles, peptide_to_smiles, sequence_to_cif, cif_to_pdb,
    validate_input, check_similarity, lookup_medicine_name, execute_docking,
    # Plus any helper functions like is_short_sequence, receive_input, etc.
)

# Import models (PredefinedProteinTarget, Compound, ScreeningJob)
from .models import PredefinedProteinTarget, Compound, ScreeningJob

# Main screening view
def screening_view(request):
    if request.method == "POST":
        form = ScreeningForm(request.POST)
        if form.is_valid():
            # 1. Input Reception
            target_choice = form.cleaned_data['target']
            raw_compound_input = form.cleaned_data['compound_input']
            input_type = form.cleaned_data['input_type']

            # For demo purposes, look up target object by choice id
            target = PredefinedProteinTarget(name=target_choice, pdb_structure="pdb_data_placeholder")

            # 2. Pre-Processing & Conversion
            if input_type == 'marvin_drawn':
                compound_data = marvin_to_smiles(raw_compound_input)
            elif input_type == 'SMILES':
                compound_data = raw_compound_input
            elif input_type == 'PDB':
                compound_data = raw_compound_input
            elif input_type == 'sequence':
                # Assume short sequence threshold of 30 characters.
                if len(raw_compound_input) <= 30:
                    compound_data = peptide_to_smiles(raw_compound_input)
                else:
                    cif_result = sequence_to_cif(raw_compound_input)
                    if 'error' in cif_result:
                        return HttpResponse("Error in CIF conversion: " + cif_result['error'])
                    compound_data = cif_to_pdb(cif_result['cif_data'])
            else:
                return HttpResponse("Unsupported input type.")

            # Create Compound instance
            compound = Compound(input_type=input_type, data=compound_data)

            # 3. Validation
            if not validate_input(compound.input_data, compound.type):
                return HttpResponse("Invalid input provided.")

            # 4. Drug Repurposing / Similarity Check
            similarity_results, repurposing_info = check_similarity(compound.input_data, compound.type)
            if similarity_results.exceeds_threshold():
                medicine_name = lookup_medicine_name(similarity_results)
                repurposing_info['medicine_name'] = medicine_name
                # Optionally, add a message to alert the user here

            # 5. Docking Execution
            docking_score = execute_docking(target, compound.input_data, compound.type)

            # 6. Output Generation (using a helper function to build output data)
            output = {
                'target': target.name,
                'compound': compound.input_data,
                'docking_score': docking_score,
                'repurposing_info': repurposing_info,
            }

            # 7. Save Screening Job
            screening_job = ScreeningJob(target, compound, docking_score, repurposing_info)
            screening_job.save()

            return render(request, "index.html", {'output': output, 'form': form})
    else:
        form = ScreeningForm()

    return render(request, "app/index.html", {'form': form})
