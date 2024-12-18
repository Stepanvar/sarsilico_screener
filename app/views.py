from django.shortcuts import render, redirect
from django.contrib.auth.decorators import login_required
from django.http import HttpResponse, JsonResponse
import subprocess
import uuid
from .forms import SMILESInputForm, PDBUploadForm, SequenceForm, TargetSelectionForm, SimilarityForm
from .models import Job, Target, KnownDrug
from .tasks import process_sequence_task, perform_docking_task, similarity_analysis_task
from .utils import validate_pdb_file, validate_smiles
from django.contrib.auth import authenticate, login, logout

def home_view(request):
    return render(request, 'app/home.html')

def about_view(request):
    return render(request, 'app/about.html')

def contact_view(request):
    return render(request, 'app/contact.html')

def login_view(request):
    if request.method == 'POST':
        username = request.POST.get('username')
        password = request.POST.get('password')
        user = authenticate(request, username=username, password=password)
        if user:
            login(request, user)
            return redirect('home')
        else:
            return render(request, 'app/login.html', {'error': 'Invalid credentials'})
    return render(request, 'app/login.html')

def logout_view(request):
    logout(request)
    return redirect('home')

def convert_to_smiles(request):
    if request.method == "POST":
        sequence = request.POST.get("sequence", "").strip()
        if not sequence:
            return JsonResponse({"error": "No sequence provided."}, status=400)

        try:
            result = subprocess.run(
                ["pepsmi", sequence], capture_output=True, text=True, check=True
            )
            smiles = result.stdout.strip()
            return JsonResponse({"smiles": smiles})
        except subprocess.CalledProcessError as e:
            return JsonResponse({"error": f"Conversion failed: {e.stderr}"}, status=500)
    return JsonResponse({"error": "Invalid request method."}, status=405)

@login_required
def smiles_input_view(request):
    if request.method == 'POST':
        form = SMILESInputForm(request.POST, request.FILES)
        if form.is_valid():
            smiles = form.cleaned_data.get('smiles')
            smiles_file = form.cleaned_data.get('smiles_file')
            if smiles_file:
                content = smiles_file.read().decode('utf-8').strip()
                smiles = content
            if smiles and not validate_smiles(smiles):
                form.add_error('smiles', 'Invalid SMILES')
            else:
                # Proceed with similarity check or docking
                return redirect('similarity_analysis')
    else:
        form = SMILESInputForm()
    return render(request, 'app/smiles_input.html', {'form': form})

@login_required
def pdb_upload_view(request):
    if request.method == 'POST':
        form = PDBUploadForm(request.POST, request.FILES)
        if form.is_valid():
            pdb_file = form.cleaned_data.get('pdb_file')
            valid, error = validate_pdb_file(pdb_file)
            if not valid:
                form.add_error('pdb_file', error)
            else:
                # Save and redirect to target selection or docking
                return redirect('target_selection')
    else:
        form = PDBUploadForm()
    return render(request, 'app/pdb_upload.html', {'form': form})

@login_required
def sequence_conversion_view(request):
    if request.method == 'POST':
        form = SequenceForm(request.POST)
        if form.is_valid():
            sequence = form.cleaned_data['sequence']
            job_id = str(uuid.uuid4())
            job = Job.objects.create(user=request.user, job_id=job_id, status='PENDING')
            process_sequence_task.delay(sequence, job_id)
            return redirect('docking_results', job_id=job_id)
    else:
        form = SequenceForm()
    return render(request, 'app/sequence_input.html', {'form': form})

@login_required
def target_selection_view(request):
    targets = Target.objects.all()
    if request.method == 'POST':
        form = TargetSelectionForm(request.POST, targets=targets)
        if form.is_valid():
            target_id = form.cleaned_data['target_id']
            job_id = str(uuid.uuid4())
            job = Job.objects.create(user=request.user, job_id=job_id, status='PENDING')
            # Example ligand path and config path
            ligand_path = "/path/to/ligand.pdbqt"
            config_path = "/path/to/config.txt"
            perform_docking_task.delay(target_id, ligand_path, config_path, job_id, email=request.user.email)
            return redirect('docking_results', job_id=job_id)
    else:
        form = TargetSelectionForm(targets=targets)
    return render(request, 'app/target_selection.html', {'form': form, 'targets': targets})

@login_required
def docking_results_view(request, job_id):
    try:
        job = Job.objects.get(job_id=job_id, user=request.user)
    except Job.DoesNotExist:
        return HttpResponse("Job not found.", status=404)

    context = {'job': job, 'scores': []}
    if job.status == 'COMPLETED' and "Affinity scores" in job.results:
        # parse scores from job.results
        # This is a simplification; in a real scenario, store results more structurally
        scores_str = job.results.split("Affinity scores: ")[-1]
        # Expected scores_str like ['-7.5', '-6.8']
        scores_str = scores_str.strip().strip("[]").replace("'", "")
        context['scores'] = scores_str.split(", ") if scores_str else []

    return render(request, 'app/docking_results.html', context)

@login_required
def similarity_analysis_view(request):
    if request.method == 'POST':
        form = SimilarityForm(request.POST)
        if form.is_valid():
            smiles = form.cleaned_data['smiles']
            threshold = form.cleaned_data['threshold']
            job_id = str(uuid.uuid4())
            job = Job.objects.create(user=request.user, job_id=job_id, status='PENDING')
            similarity_analysis_task.delay(smiles, threshold, job_id, email=request.user.email)
            return redirect('docking_results', job_id=job_id)
    else:
        form = SimilarityForm()
    return render(request, 'app/similarity_analysis.html', {'form': form})
