# views.py (improved)
from tkinter import NO
from django.shortcuts import render, redirect, get_object_or_404
from django.contrib import messages
from django.db import transaction
from sqlalchemy import null
from .models import PredefinedProteinTarget, Compound, ScreeningJob
from .services.conversions import (
    jsme_to_smiles, peptide_to_smiles,
    sequence_to_cif, cif_to_pdb, is_short_sequence, validate_input
)
from .services.docking import execute_docking
from .services.similarity_check import check_similarity, lookup_medicine_name
from django.contrib.auth.views import LoginView, LogoutView
from django.contrib.auth.forms import UserCreationForm, AuthenticationForm
from django.urls import reverse_lazy
from django.views.generic import CreateView

class CustomLoginView(LoginView):
    template_name = 'app/login.html'
    form_class = AuthenticationForm

class CustomRegisterView(CreateView):
    form_class = UserCreationForm
    success_url = reverse_lazy('login')
    template_name = 'app/register.html'

class CustomLogoutView(LogoutView):
    template_name = 'app/logout.html'

def index(request):
    """Main view with form handling and pipeline orchestration"""
    if request.method == 'POST':
        try:
            with transaction.atomic():  # Ensure database integrity
                # Data extraction and validation
                target_id = request.POST.get('target')
                input_type = request.POST.get('input_type')
                raw_compound = request.POST.get('compound_input', '').strip()

                if not all([target_id, input_type, raw_compound]):
                    messages.error(request, "Missing required fields")
                    return redirect('index')

                # Retrieve and validate target
                target = get_object_or_404(PredefinedProteinTarget, pk=target_id)

                # Compound processing pipeline
                compound_data = process_compound_input(input_type, raw_compound)
                validate_input(compound_data, input_type)

                # Create compound record
                compound = Compound.objects.create(
                    input_type=input_type,
                    data=compound_data
                )
                check_similarity_flag = request.POST.get('check_similarity') == 'true'
                similarity_results = {}
                if check_similarity_flag:
                    if input_type == "sequence":
                        if is_short_sequence(raw_compound) is False:
                            similarity_results = check_similarity(raw_compound, False)
                        else:
                            similarity_results = check_similarity(compound_data, True)
                    else:
                        similarity_results = check_similarity(compound_data, True)

                # Create screening job record
                screening_job = ScreeningJob.objects.create(
                    target=target,
                    compound=compound,
                    docking_info={"info": ""},
                    repurposing_info={"info": ""}
                )

                # Correct condition: check if medicine_name is a non-empty string.
                medicine_name = similarity_results.get('medicine_name', '')
                if medicine_name and medicine_name.strip():
                    messages.warning(
                        request,
                        'Similar compound found in drug database!: ' + medicine_name.strip()
                    )
                    messages.warning(
                        request,
                        "Suggesting to find drug name on this website: https://covirus.cc/drugs/"
                    )
                    screening_job.repurposing_info = similarity_results
                    screening_job.save()
                    return redirect('results', job_id=screening_job.id)

                # Docking execution stage
                docking_info = execute_docking(target, compound_data, input_type)
                screening_job.docking_info = docking_info
                screening_job.save()
                return redirect('results', job_id=screening_job.id)

        except Exception as e:
            transaction.rollback()
            messages.error(request, f"Processing error: {str(e)}")
            return redirect('index')

    # GET request - show empty form
    return render(request, 'app/index.html', {
        'targets': PredefinedProteinTarget.objects.all()
    })


def process_compound_input(input_type: str, raw_input: str):
    """Centralized compound input processing"""
    if input_type == 'jsme_drawn':
        return jsme_to_smiles(raw_input)
    elif input_type == 'SMILES':
        return raw_input
    elif input_type == 'PDB':
        return raw_input
    elif input_type == 'sequence':
        if is_short_sequence(raw_input):
            return peptide_to_smiles(raw_input)
        cif_result = sequence_to_cif(raw_input)
        if 'error' in cif_result:
            raise ValueError(f"CIF conversion failed: {cif_result['error']}")
        return cif_to_pdb(cif_result['cif_data'])
    raise ValueError(f"Unsupported input type: {input_type}")

def results_view(request, job_id):
    """Results display view with proper context handling"""
    try:
        job = ScreeningJob.objects.select_related('target', 'compound').get(id=job_id)
        return render(request, 'app/results.html', {
            'output': {
                'target': job.target.name,
                'docking_info': job.docking_info,  # Formatted score
                'repurposing_info': job.repurposing_info or {},
                'compound_data': job.compound.data
            }
        })
    except ScreeningJob.DoesNotExist:
        messages.error(request, "Requested results not found")
        return redirect('index')
