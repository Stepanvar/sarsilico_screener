from django.shortcuts import render, redirect, get_object_or_404
from django.contrib.auth.decorators import login_required
from django.urls import reverse
from .models import Job, Target, KnownDrug
from .forms import *
from .tasks import perform_docking_task, perform_similarity_analysis_task
from django.shortcuts import render, redirect
from django.contrib import messages
from django.contrib.auth import login
from django.contrib.auth.forms import UserCreationForm
from django.contrib.auth.decorators import login_required
from .models import Job
from .utils import run_pepsmi, run_alphafold_convert

def home(request):
    return render(request, 'app/home.html')

@login_required
def smiles_input(request):
    if request.method == 'POST':
        form = SMILESInputForm(request.POST, request.FILES)
        if form.is_valid():
            return redirect(reverse('target_selection'))
    else:
        form = SMILESInputForm()
    return render(request, 'app/smiles_input.html', {'form': form})

@login_required
def pdb_upload(request):
    if request.method == 'POST':
        form = PDBUploadForm(request.POST, request.FILES)
        if form.is_valid():
            return redirect(reverse('target_selection'))
    else:
        form = PDBUploadForm()
    return render(request, 'app/pdb_upload.html', {'form': form})

@login_required
def target_selection(request):
    if request.method == 'POST':
        form = DockingForm(request.POST)
        if form.is_valid():
            # Process target selection
            return redirect('docking') 
    else:
        form = DockingForm()
    
    targets = Target.objects.all()
    return render(request, 'app/target_selection.html', {
        'form': form,
        'targets': targets
    })


@login_required
def similarity_analysis(request):
    if request.method == 'POST':
        form = SimilarityForm(request.POST)
        if form.is_valid():
            job = Job.objects.create(user=request.user)
            perform_similarity_analysis_task.delay(
                form.cleaned_data['smiles'],
                form.cleaned_data['threshold'],
                str(job.job_id),
                request.user.email
            )
            return redirect('similarity_results', job_id=job.job_id)
    else:
        form = SimilarityForm()
    return render(request, 'app/similarity_analysis.html', {'form': form})

@login_required
def docking(request):
    if request.method == 'POST':
        form = DockingForm(request.POST)
        if form.is_valid():
            job = Job.objects.create(user=request.user)
            perform_docking_task.delay(
                form.cleaned_data['target'].id,
                str(job.job_id),
                request.user.email
            )
            return redirect('docking_results', job_id=job.job_id)
    else:
        form = DockingForm()
    return render(request, 'app/docking.html', {'form': form})

@login_required
def docking_results(request, job_id):
    job = get_object_or_404(Job, job_id=job_id, user=request.user)
    scores = []
    if job.status == 'COMPLETED' and "Affinity scores:" in str(job.results):
        result_str = str(job.results).split("Affinity scores:")[-1].strip().strip("[]")
        scores = [s.strip("' ") for s in result_str.split(",")] if result_str else []
    return render(request, 'app/docking_results.html', {'job': job, 'scores': scores})

@login_required
def similarity_results(request, job_id):
    """
    Display results of similarity analysis for a specific job
    """
    job = Job.objects.get(job_id=job_id, user=request.user)
    
    context = {
        'job': job,
        'results': eval(job.results) if job.results else None
    }
    return render(request, 'app/similarity_results.html', context)

@login_required
def sequence_conversion(request):
    """
    Handle sequence conversion from amino acid sequences to structures
    """
    converted_structure = None
    error = None
    
    if request.method == 'POST':
        sequence = request.POST.get('sequence', '')
        sequence_file = request.FILES.get('sequence_file')
        
        try:
            if sequence_file:
                # Handle file upload
                sequence = sequence_file.read().decode('utf-8')
                
            # Perform conversion using utility functions
            if len(sequence) <= 50:  # Short peptide
                converted_structure = run_pepsmi(sequence)
            else:  # Long sequence
                converted_structure = run_alphafold_convert(sequence)
                
            messages.success(request, 'Sequence converted successfully!')
            
        except Exception as e:
            error = str(e)
            messages.error(request, f'Conversion failed: {error}')

    return render(request, 'app/sequence_conversion.html', {
        'converted_structure': converted_structure,
        'error': error
    })

def register(request):
    """
    Handle user registration with built-in Django form
    """
    if request.method == 'POST':
        form = UserCreationForm(request.POST)
        if form.is_valid():
            user = form.save()
            login(request, user)
            messages.success(request, 'Registration successful!')
            return redirect('home')
        else:
            messages.error(request, 'Registration failed. Please correct the errors below.')
    else:
        form = UserCreationForm()

    return render(request, 'app/register.html', {'form': form})