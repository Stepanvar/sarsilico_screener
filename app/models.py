# models.py
from uuid import uuid4
from django.db import models

class PredefinedProteinTarget(models.Model):
    name = models.CharField(max_length=255)
    description = models.TextField()
    pdb_file = models.FileField(upload_to='pdb_files/')

class Compound(models.Model):
    input_type = models.CharField(max_length=20)
    data = models.TextField()
    
    def save(self, *args, **kwargs):
        """Custom save logic for compound data"""
        super().save(*args, **kwargs)

class ScreeningJob(models.Model):
    id = models.UUIDField(primary_key=True, default=uuid4, editable=False)  
    target = models.ForeignKey(PredefinedProteinTarget, on_delete=models.CASCADE)
    compound = models.ForeignKey(Compound, on_delete=models.CASCADE)
    docking_info = models.JSONField()
    repurposing_info = models.JSONField()
    pdb_file = models.FileField("pdb" + str(id)[:16])
    
    def save(self, *args, **kwargs):
        """Custom save logic for screening results"""
        super().save(*args, **kwargs)
