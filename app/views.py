from django.shortcuts import render, redirect, get_object_or_404
from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.core.files.uploadedfile import SimpleUploadedFile
from .forms import SMILESInputForm, PDBUploadForm, SimilarityForm, DockingForm, SequenceForm, TargetSelectionForm, UserRegistrationForm, TargetSelectionForm
from .models import Job, Target, KnownDrug
from .tasks import perform_docking_task, perform_similarity_analysis_task
from .utils import validate_smiles, validate_pdb_file, run_pepsmi, run_alphafold_convert
from django.contrib.auth import authenticate, login, logout
from django.contrib.auth.decorators import login_required
import uuid

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

def smiles_input_view(request):
    if request.method == 'POST':
        form = SMILESInputForm(request.POST, request.FILES)
        if form.is_valid():
            smiles = form.cleaned_data.get('smiles')
            smiles_file = form.cleaned_data.get('smiles_file')

            if smiles:
                is_valid, error = validate_smiles(smiles)
                if not is_valid:
                    form.add_error('smiles', str(error))
                    return render(request, 'app/smiles_input.html', {'form': form})

            if smiles_file:
                file_content = smiles_file.read().decode('utf-8').strip()
                is_valid, error = validate_smiles(file_content)
                if not is_valid:
                    form.add_error('smiles_file', str(error))
                    return render(request, 'app/smiles_input.html', {'form': form})

            messages.success(request, 'SMILES input successfully validated and processed.')
            return redirect('docking_results', job_id='some-job-id')  # Adjust job_id as needed
    else:
        form = SMILESInputForm()
    return render(request, 'app/smiles_input.html', {'form': form})

def pdb_upload_view(request):
    if request.method == 'POST':
        form = PDBUploadForm(request.POST, request.FILES)
        if form.is_valid():
            pdb_file = form.cleaned_data.get('pdb_file')
            is_valid, error = validate_pdb_file(pdb_file)
            if not is_valid:
                form.add_error('pdb_file', str(error))
                return render(request, 'app/pdb_upload.html', {'form': form})

            messages.success(request, 'PDB file successfully validated and processed.')
            return redirect('docking_results', job_id='some-job-id')  # Adjust job_id as needed
    else:
        form = PDBUploadForm()
    return render(request, 'app/pdb_upload.html', {'form': form})

@login_required
def docking_view(request):
    if request.method == 'POST':
        form = DockingForm(request.POST)
        if form.is_valid():
            target = form.cleaned_data.get('target')
            smiles = form.cleaned_data.get('smiles')
            job = Job.objects.create(user=request.user, status='PENDING')
            perform_docking_task.delay(target.id, smiles, job.job_id, request.user.email)
            messages.success(request, 'Docking simulation initiated.')
            return redirect('docking_results', job_id=job.job_id)
    else:
        form = DockingForm()
    return render(request, 'app/docking.html', {'form': form})

@login_required
def docking_results_view(request, job_id):
    job = get_object_or_404(Job, job_id=job_id, user=request.user)
    # Parse affinity scores if completed
    scores = []
    if job.status == 'COMPLETED' and "Affinity scores:" in str(job.results):
        result_str = str(job.results).split("Affinity scores:")[-1].strip().strip("[]")
        scores = [s.strip("' ") for s in result_str.split(",")] if result_str else []
    context = {'job': job, 'scores': scores}
    return render(request, 'app/docking_results.html', context)

@login_required
def similarity_analysis_view(request):
    if request.method == 'POST':
        form = SimilarityForm(request.POST)
        if form.is_valid():
            smiles = form.cleaned_data['smiles']
            threshold = form.cleaned_data['threshold']
            job = Job.objects.create(user=request.user, status='PENDING')
            perform_similarity_analysis_task.delay(smiles, threshold, job.job_id, request.user.email)
            messages.success(request, 'Similarity analysis initiated.')
            return redirect('similarity_results', job_id=job.job_id)
    else:
        form = SimilarityForm()
    return render(request, 'app/similarity_analysis.html', {'form': form})

@login_required
def similarity_results_view(request, job_id):
    job = get_object_or_404(Job, job_id=job_id, user=request.user)
    results = None
    if job.status == 'COMPLETED' and "Drugs above threshold:" in str(job.results):
        # Extract the JSON-like structure from job.results
        # Expected format: "Drugs above threshold: [{'name': '...', 'smiles': '...', 'score': ..., 'description': '...'}, ...]"
        # You can safely use eval or json.loads after some manipulation if carefully sanitized
        data_str = str(job.results).split("Drugs above threshold:")[-1].strip()
        # Since we trust internal data, we can eval it. For production, parse more safely.
        drugs = eval(data_str)
        results = {'drugs': drugs}
    context = {'job': job, 'results': results}
    return render(request, 'app/similarity_results.html', context)

@login_required
def sequence_conversion_view(request):
    if request.method == 'POST':
        form = SequenceForm(request.POST, request.FILES)
        if form.is_valid():
            sequence = form.cleaned_data.get('sequence')
            sequence_file = form.cleaned_data.get('sequence_file')
            if sequence:
                # Short sequence: run PepSMI
                try:
                    smiles = run_pepsmi(sequence)
                    messages.success(request, f'Sequence converted to SMILES: {smiles}')
                    return redirect('docking_results', job_id='some-job-id')
                except ValueError as e:
                    form.add_error(None, str(e))
            elif sequence_file:
                # Longer sequence: run AlphaFold convert
                input_cif = handle_uploaded_file(sequence_file) # Implement file saving as needed
                output_pdb = f"/path/to/output/{sequence_file.name}.pdb"
                try:
                    pdb = run_alphafold_convert(input_cif, output_pdb)
                    messages.success(request, f'Sequence converted to PDB: {pdb}')
                    return redirect('docking_results', job_id='some-job-id')
                except ValueError as e:
                    form.add_error(None, str(e))
    else:
        form = SequenceForm()
    return render(request, 'app/sequence_conversion.html', {'form': form})

def handle_uploaded_file(f):
    file_path = f"/path/to/uploads/{f.name}"
    with open(file_path, 'wb+') as destination:
        for chunk in f.chunks():
            destination.write(chunk)
    return file_path

def register_view(request):
    if request.method == 'POST':
        form = UserRegistrationForm(request.POST)
        if form.is_valid():
            form.save()
            username = form.cleaned_data.get('username')
            messages.success(request, f'Account created for {username}! You can now log in.')
            return redirect('login')
        else:
            messages.error(request, 'Please correct the error below.')
    else:
        form = UserRegistrationForm()
    return render(request, 'app/register.html', {'form': form})

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
