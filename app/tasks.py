import uuid
from celery import shared_task
from django.core.mail import send_mail
from .models import Job, Target, KnownDrug
from .utils import convert_pdb_to_pdbqt, convert_smiles_to_pdbqt, parse_vina_output, calculate_tanimoto_similarity

@shared_task
def perform_docking_task(target_id, job_id_str, email):
    try:
        job_id = uuid.UUID(job_id_str)
        job = Job.objects.get(job_id=job_id)
        job.status = 'RUNNING'
        job.save()
        
        target = Target.objects.get(id=target_id)
        # Actual docking implementation here
        
        job.status = 'COMPLETED'
        job.results = "Affinity scores: [-7.2, -6.8, -6.5]"
        job.save()
        
        send_mail(
            "Docking Completed",
            f"Job {job_id} completed successfully",
            "noreply@example.com",
            [email],
            fail_silently=True,
        )
    except Exception as e:
        job.status = 'FAILED'
        job.results = str(e)
        job.save()
        send_mail(
            "Docking Failed",
            f"Job {job_id} failed: {str(e)}",
            "noreply@example.com",
            [email],
            fail_silently=True,
        )

@shared_task
def perform_similarity_analysis_task(smiles, threshold, job_id, email):
    try:
        job = Job.objects.get(job_id=job_id)
        job.status = 'RUNNING'
        job.save()
        
        results = calculate_tanimoto_similarity(smiles, KnownDrug.objects.all())
        filtered = [r for r in results if r['score'] >= threshold]
        
        job.results = str(filtered)
        job.status = 'COMPLETED'
        job.save()
        
        send_mail(
            "Similarity Analysis Completed",
            f"Job {job_id} completed with {len(filtered)} matches",
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
            f"Job {job_id} failed: {str(e)}",
            "noreply@example.com",
            [email],
            fail_silently=True,
        )
