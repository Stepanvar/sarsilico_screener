from django.db import models
from django.contrib.auth.models import User

class Compound(models.Model):
    id = models.BigAutoField(primary_key=True)
    STRUCTURE_CHOICES = [
        ('SMILES', 'SMILES'),
        ('PDB', 'PDB'),
    ]
    name = models.CharField(max_length=255)
    smiles = models.TextField(blank=True, null=True)
    pdb_file = models.FileField(upload_to='pdb_files/', blank=True, null=True)
    structure_type = models.CharField(max_length=50, choices=STRUCTURE_CHOICES)
    affinity_score = models.FloatField(null=True, blank=True)
    created_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return self.name

class Target(models.Model):
    name = models.CharField(max_length=255, unique=True)
    description = models.TextField(blank=True)

    def __str__(self):
        return self.name

class KnownDrug(models.Model):
    name = models.CharField(max_length=255)
    smiles = models.TextField()
    tanimoto_score = models.FloatField(null=True, blank=True)

    def __str__(self):
        return self.name

class Job(models.Model):
    STATUS_CHOICES = [
        ('PENDING', 'Pending'),
        ('RUNNING', 'Running'),
        ('COMPLETED', 'Completed'),
        ('FAILED', 'Failed'),
    ]

    user = models.ForeignKey(User, on_delete=models.CASCADE)
    job_id = models.CharField(max_length=100, unique=True)
    status = models.CharField(max_length=10, choices=STATUS_CHOICES, default='PENDING')
    results = models.TextField(blank=True)
    created_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return f"{self.job_id} - {self.status}"
