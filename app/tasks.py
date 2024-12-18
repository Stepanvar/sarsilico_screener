from celery import shared_task
from django.core.mail import send_mail
from .models import Job, KnownDrug, Target
from .utils import (
    run_pepsmi, run_alphafold_convert, convert_to_pdbqt, run_vina, parse_vina_output, calculate_tanimoto_similarity
)
import uuid
import os

@shared_task
def process_sequence_task(sequence, job_id):
    job = Job.objects.get(job_id=job_id)
    try:
        job.status = 'RUNNING'
        job.save()
        smiles = run_pepsmi(sequence)
        # For longer sequences, run AlphaFold and convert CIF to PDB as needed
        # Store SMILES or structural data in job.results for reference
        job.results = f"SMILES: {smiles}"
        job.status = 'COMPLETED'
        job.save()
        send_mail("Sequence Conversion Completed",
                  f"Your sequence conversion (Job ID: {job_id}) has completed successfully.\nSMILES: {smiles}",
                  "noreply@example.com",
                  [job.user.email],
                  fail_silently=True)
    except Exception as e:
        job.status = 'FAILED'
        job.results = str(e)
        job.save()
        send_mail("Sequence Conversion Failed",
                  f"Your sequence conversion (Job ID: {job_id}) failed.\nError: {str(e)}",
                  "noreply@example.com",
                  [job.user.email],
                  fail_silently=True)

@shared_task
def perform_docking_task(target_id, ligand_path, config_path, job_id, email=None):
    job = Job.objects.get(job_id=job_id)
    try:
        job.status = 'RUNNING'
        job.save()
        target = Target.objects.get(id=target_id)
        receptor_pdb = f"/path/to/targets/{target.name}.pdb"
        receptor_pdbqt = f"/path/to/temp/{target.name}.pdbqt"
        convert_to_pdbqt(receptor_pdb, receptor_pdbqt)

        output_pdbqt = f"/path/to/temp/{job_id}_docked.pdbqt"
        run_vina(receptor_pdbqt, ligand_path, output_pdbqt, config_path)

        scores = parse_vina_output(output_pdbqt)
        job.results = f"Affinity scores: {scores}"
        job.status = 'COMPLETED'
        job.save()

        if email:
            send_mail("Docking Simulation Completed",
                      f"Your docking job (Job ID: {job_id}) completed.\nScores: {scores}",
                      "noreply@example.com",
                      [email],
                      fail_silently=True)
    except Exception as e:
        job.status = 'FAILED'
        job.results = str(e)
        job.save()
        if email:
            send_mail("Docking Simulation Failed",
                      f"Your docking job (Job ID: {job_id}) failed.\nError: {str(e)}",
                      "noreply@example.com",
                      [email],
                      fail_silently=True)

@shared_task
def similarity_analysis_task(smiles, threshold, job_id, email=None):
    job = Job.objects.get(job_id=job_id)
    try:
        job.status = 'RUNNING'
        job.save()
        known_drugs = KnownDrug.objects.all()
        results = calculate_tanimoto_similarity(smiles, known_drugs)
        filtered = [r for r in results if r['score'] >= threshold]
        job.results = f"Drugs above threshold: {filtered}"
        job.status = 'COMPLETED'
        job.save()

        if email:
            send_mail("Similarity Analysis Completed",
                      f"Your similarity analysis (Job ID: {job_id}) completed.\nResults: {filtered}",
                      "noreply@example.com",
                      [email],
                      fail_silently=True)
    except Exception as e:
        job.status = 'FAILED'
        job.results = str(e)
        job.save()
        if email:
            send_mail("Similarity Analysis Failed",
                      f"Your similarity analysis (Job ID: {job_id}) failed.\nError: {str(e)}",
                      "noreply@example.com",
                      [email],
                      fail_silently=True)
