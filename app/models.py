# app/models.py
from uuid import uuid4
from django.db import models
from django.core.validators import MaxValueValidator
from django.forms import ValidationError
from django.utils import timezone
from django.db.models import JSONField

class PredefinedProteinTarget:
    """
    Model for storing target information.
    Attributes:
        name: Target name.
        pdb_structure: Full-atomic PDB structure.
    """
    def __init__(self, name, pdb_structure):
        self.name = name
        self.pdb_structure = pdb_structure

class Target(models.Model):
    """
    Model for storing target information.
    Attributes:
        name: Target name.
        pdb_structure: Full-atomic PDB structure.
    """
    def __init__(self, name, pdb_structure):
        self.name = name
        self.pdb_structure = pdb_structure
    id = models.UUIDField(primary_key=True, default=uuid4, editable=False)  
    PROTEIN_TYPE_CHOICES = [
        ('SPIKE', 'Spike Protein'),
        ('PROTEASE', 'Protease'),
        ('POLYMERASE', 'Polymerase'),
    ]

    name = models.CharField(
        max_length=100,
        unique=True,
        help_text="Unique identifier for the target protein"
    )
    pdb_file = models.FileField(
        upload_to='targets/',
        help_text="Validated PDB structure file"
    )
    protein_type = models.CharField(
        max_length=20,
        choices=PROTEIN_TYPE_CHOICES,
        help_text="Type of viral protein"
    )
    resolution = models.FloatField(
        validators=[MaxValueValidator(3.0)],
        help_text="X-ray resolution in Ångströms (≤3.0Å required)"
    )
    created_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        ordering = ['-created_at']
        verbose_name = "Predefined Target"

    def __str__(self):
        return f"{self.name} ({self.protein_type})"

class Compound(models.Model):
    """
    Model for storing compound details.
    Attributes:
        input_type: e.g. 'SMILES', 'PDB', 'sequence', 'marvin_drawn'.
        data: Raw or converted molecular data.
    """
    
    id = models.UUIDField(primary_key=True, default=uuid4, editable=False)  
    COMPOUND_TYPE_CHOICES = [
        ('SMILES', 'SMILES (Drawn/Entered)'),
        ('PDB', 'PDB File'),
        ('PEPTIDE', 'Short Peptide (2b)'),
        ('BIOMOLECULE', 'Large Biomolecule (2c)'),
    ]
    type = models.CharField(
        max_length=20,
        choices=COMPOUND_TYPE_CHOICES,
        help_text="Input type as per section 2 requirements"
    )
    input_data = models.TextField(
        help_text="Raw input (SMILES, sequence, or PDB content)"
    )
    validation_status = models.BooleanField(
        default=False,
        help_text="Validation outcome from section 2 checks"
    )
    validation_errors = JSONField(
        blank=True,
        help_text="Validation error messages",
        default=dict  # Add default empty dict
    )

    converted_smiles = models.TextField(
        null=True,
        blank=True,
        help_text="Converted SMILES for docking (peptides/small molecules)"
    )
    converted_pdb = models.FileField(
        upload_to='converted_structures/',
        null=True,
        blank=True,
        help_text="Generated PDB from AlphaFold (biomolecules)"
    )
    alphafold_job_id = models.CharField(
        max_length=50,
        null=True,
        blank=True,
        help_text="Tracking ID for AlphaFold predictions"
    )
    created_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        ordering = ['-created_at']
        indexes = [
            models.Index(fields=['type', 'validation_status']),
        ]

    def __str__(self):
        return f"{self.type} ({self.created_at.date()})"

class ScreeningJob(models.Model):
    """
    Model for recording a screening run.
    Attributes:
        target: A PredefinedProteinTarget instance.
        compound: A Compound instance.
        docking_score: The affinity score from docking.
        repurposing_info: Additional drug repurposing info.
    """
    id = models.UUIDField(primary_key=True, default=uuid4, editable=False)  
    STATUS_CHOICES = [
        ('PENDING', 'Pending'),
        ('VALIDATING', 'Validating Input'),
        ('SIM_CHECK', 'Similarity Check (Section 3)'),
        ('DOCKING', 'Docking Execution (Section 4)'),
        ('COMPLETED', 'Completed'),
        ('FAILED', 'Failed'),
    ]
    target = models.ForeignKey(
        Target,
        on_delete=models.CASCADE,
        related_name='jobs',
        help_text="Selected target protein from section 1"
    )
    compound = models.ForeignKey(
        Compound,
        on_delete=models.CASCADE,
        related_name='jobs',
        help_text="Input compound with validated/converted data"
    )
    status = models.CharField(
        max_length=20,
        choices=STATUS_CHOICES,
        default='PENDING',
        help_text="Current workflow state per section 5 sequence diagram"
    )
    similarity_data = JSONField(
        blank=True,
        help_text="TTD/CoviDrug analysis results",
        default=dict  # Proper default for empty JSON
    )
    docking_score = models.FloatField(
        null=True,
        blank=True,
        help_text="Best AutoDock Vina affinity score (kcal/mol)"
    )
    ligand_pdbqt = models.FileField(
        upload_to='docking_files/ligands/',
        null=True,
        blank=True,
        help_text="Prepared ligand file for Vina docking"
    )
    complex_structure = models.FileField(
        upload_to='results/complexes/',
        null=True,
        blank=True,
        help_text="Final protein-ligand complex structure"
    )
    visualization_config = JSONField(
        null=True,
        blank=True,
        help_text="NGL.js visualization parameters for section 5 output"
    )
    created_at = models.DateTimeField(auto_now_add=True)
    completed_at = models.DateTimeField(
        null=True,
        blank=True,
        help_text="Timestamp of job completion/failure"
    )

    class Meta:
        ordering = ['-created_at']
        indexes = [
            models.Index(fields=['status']),
            models.Index(fields=['docking_score']),
        ]

    def __str__(self):
        return f"Job-{self.id} ({self.status})"

    def clean(self):
        """Enforce workflow rules from insilico.txt requirements"""
        # Section 1 resolution validation
        if self.target.resolution > 3.0:
            raise ValidationError("Target resolution exceeds 3.0Å limit")
        
        # Section 5 completion requirements
        if self.status == 'COMPLETED' and not self.docking_score:
            raise ValidationError("Completed jobs must have docking scores")
        
        # Section 2 validation pre-requisite
        if self.status != 'PENDING' and not self.compound.validation_status:
            raise ValidationError("Compound must be validated before processing")

    def save(self, *args, **kwargs):
        """Automatically handle completion timestamps"""
        if self.status in ['COMPLETED', 'FAILED'] and not self.completed_at:
            self.completed_at = timezone.now()
        super().save(*args, **kwargs)
        
    def __init__(self, target, compound, docking_score, repurposing_info=None):
        self.target = target
        self.compound = compound
        self.docking_score = docking_score
        self.repurposing_info = repurposing_info or {}
