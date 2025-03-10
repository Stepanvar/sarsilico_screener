from django.db import models
from django.contrib.auth.models import User
import uuid

class Job(models.Model):
    STATUS_CHOICES = [
        ('PENDING', 'Pending'),
        ('RUNNING', 'Running'),
        ('COMPLETED', 'Completed'),
        ('FAILED', 'Failed'),
    ]
    
    job_id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    user = models.ForeignKey(User, on_delete=models.CASCADE)
    status = models.CharField(max_length=10, choices=STATUS_CHOICES, default='PENDING')
    results = models.TextField(null=True, blank=True)
    results_binding_pose_path = models.CharField(max_length=500, null=True, blank=True)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    def __str__(self):
        return f"Job {self.job_id}"

class KnownDrug(models.Model):
    name = models.CharField(max_length=255)
    smiles = models.CharField(max_length=255)
    description = models.TextField()

    def __str__(self):
        return self.name

class Target(models.Model):
    name = models.CharField(max_length=255)
    pdb_id = models.CharField(max_length=10)
    description = models.TextField()

    def __str__(self):
        return f"{self.name} ({self.pdb_id})"
