from celery import shared_task
from django.core.mail import send_mail
from .models import Job, Target, KnownDrug
from .utils import (
    convert_pdb_to_pdbqt, convert_smiles_to_pdbqt, parse_vina_output, calculate_tanimoto_similarity
)
import subprocess
import uuid
import os

@shared_task
def perform_docking_task(target_id, smiles, job_id, email):
    job = Job.objects.get(job_id=job_id)
    try:
        job.status = 'RUNNING'
        job.save()

        target = Target.objects.get(id=target_id)
        receptor_pdb = f"/path/to/targets/{target.pdb_id}.pdb"
        receptor_pdbqt = f"/path/to/temp/{target.pdb_id}.pdbqt"
        convert_pdb_to_pdbqt(receptor_pdb, receptor_pdbqt)

        ligand_pdbqt = f"/path/to/temp/{uuid.uuid4()}.pdbqt"
        convert_smiles_to_pdbqt(smiles, ligand_pdbqt)

        output_pdbqt = f"/path/to/temp/{job_id}_docked.pdbqt"
        config_file = "/path/to/config.txt"

        subprocess.run([
            "vina",
            "--receptor", receptor_pdbqt,
            "--ligand", ligand_pdbqt,
            "--out", output_pdbqt,
            "--config", config_file
        ], check=True)

        affinity_scores = parse_vina_output(output_pdbqt)
        job.results = f"Affinity scores: {affinity_scores}"
        job.results_binding_pose_path = output_pdbqt
        job.status = 'COMPLETED'
        job.save()

        send_mail(
            "Docking Simulation Completed",
            f"Your docking job (Job ID: {job_id}) completed.\nAffinity Scores: {affinity_scores}",
            "noreply@example.com",
            [email],
            fail_silently=True,
        )
    except Exception as e:
        job.status = 'FAILED'
        job.results = str(e)
        job.save()
        send_mail(
            "Docking Simulation Failed",
            f"Your docking job (Job ID: {job_id}) failed.\nError: {str(e)}",
            "noreply@example.com",
            [email],
            fail_silently=True,
        )

@shared_task
def perform_similarity_analysis_task(smiles, threshold, job_id, email):
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

        send_mail(
            "Similarity Analysis Completed",
            f"Your similarity analysis (Job ID: {job_id}) completed.\nResults: {filtered}",
            "noreply@example.com",
            [email],
            fail_silently=True,
        )
    except Exception as e:
        job.status = 'FAILED'
        job.results = str(e)
        job.save()
        send_mail(
            "Similarity Analysis Failed",
            f"Your similarity analysis (Job ID: {job_id}) failed.\nError: {str(e)}",
            "noreply@example.com",
            [email],
            fail_silently=True,
        )
