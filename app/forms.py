# app/forms.py
from django import forms
from django.core.validators import MaxLengthValidator
from django.conf import settings
from rdkit import Chem
import re

class SMILESForm(forms.Form):
    """
    Validates SMILES/drawn structures (section 2a)
    Supports both direct input and MarvinJS sketches
    """
    smiles = forms.CharField(
        widget=forms.TextInput(attrs={
            'class': 'form-control',
            'placeholder': 'C1=CC=CC=C1'
        }),
        help_text="Enter valid SMILES notation or use the structure drawer"
    )
    
    marvinjs = forms.CharField(
        required=False,
        widget=forms.HiddenInput(),
        help_text="MarvinJS sketch data"
    )

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Add uniform styling
        for field in self.fields.values():
            field.widget.attrs.update({'class': 'form-control'})
            
    def clean(self):
        cleaned_data = super().clean()
        smiles = cleaned_data.get('smiles')
        marvinjs = cleaned_data.get('marvinjs')

        if not (smiles or marvinjs):
            raise forms.ValidationError(
                "Please enter a SMILES string or draw a structure",
                code='missing_input'
            )

        # Prioritize MarvinJS data if both provided
        if marvinjs:
            try:
                mol = Chem.MolFromMrvBlock(marvinjs)
                if not mol:
                    raise forms.ValidationError(
                        "Invalid chemical structure in MarvinJS sketch",
                        code='invalid_mrv'
                    )
                Chem.SanitizeMol(mol)
                cleaned_data['smiles'] = Chem.MolToSmiles(mol)
            except Exception as e:
                raise forms.ValidationError(
                    f"MarvinJS conversion failed: {str(e)}",
                    code='mrv_conversion_error'
                )
        else:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                raise forms.ValidationError(
                    "Invalid SMILES format",
                    code='invalid_smiles'
                )
            try:
                Chem.SanitizeMol(mol)
            except ValueError as e:
                raise forms.ValidationError(
                    f"Chemical validation failed: {str(e)}",
                    code='sanitization_error'
                )

        return cleaned_data

class PDBUploadForm(forms.Form):
    """
    Validates PDB files (section 2 requirements)
    Checks format, resolution, and structural validity
    """
    pdb_file = forms.FileField(
        widget=forms.ClearableFileInput(attrs={
            'accept': '.pdb',
            'class': 'pdb-upload'
        }),
        help_text="Upload PDB file (max resolution 3.0Å)"
    )

    def clean_pdb_file(self):
        pdb_file = self.cleaned_data['pdb_file']
        
        # File extension check
        if not pdb_file.name.lower().endswith('.pdb'):
            raise forms.ValidationError(
                "Invalid file format - only .pdb files accepted",
                code='invalid_format'
            )

        # Read and validate content
        try:
            content = pdb_file.read().decode()
            mol = Chem.MolFromPDBBlock(content)
            
            if not mol or mol.GetNumAtoms() == 0:
                raise forms.ValidationError(
                    "Invalid PDB structure - missing atoms",
                    code='invalid_structure'
                )

            # Extract resolution from REMARK records
            resolution = None
            for line in content.split('\n'):
                if line.startswith('REMARK   2 RESOLUTION.'):
                    try:
                        resolution = float(line.split()[-1])
                        if resolution > settings.MAX_PDB_RESOLUTION:
                            raise forms.ValidationError(
                                f"Resolution {resolution}Å exceeds maximum allowed {settings.MAX_PDB_RESOLUTION}Å",
                                code='resolution_exceeded'
                            )
                    except (IndexError, ValueError):
                        pass

            if resolution is None:
                raise forms.ValidationError(
                    "Missing resolution in PDB remarks",
                    code='missing_resolution'
                )

        except UnicodeDecodeError:
            raise forms.ValidationError(
                "Invalid file content - not a text PDB file",
                code='invalid_encoding'
            )

        return pdb_file

class SequenceForm(forms.Form):
    """
    Validates peptide/biomolecule sequences (section 2b/c)
    Handles both amino acid and nucleotide inputs
    """
    SEQUENCE_TYPE_CHOICES = [
        ('peptide', 'Short Peptide (≤20 residues)'),
        ('biomolecule', 'Large Biomolecule (≥50 residues)'),
    ]

    sequence = forms.CharField(
        widget=forms.Textarea(attrs={
            'rows': 4,
            'placeholder': 'Enter amino acid or nucleotide sequence'
        }),
        validators=[MaxLengthValidator(2000)],
        help_text="Enter sequence in single-letter code"
    )
    sequence_type = forms.ChoiceField(
        choices=SEQUENCE_TYPE_CHOICES,
        widget=forms.RadioSelect,
        initial='peptide'
    )

    def clean_sequence(self):
        sequence = self.cleaned_data['sequence'].upper().replace(' ', '')
        seq_type = self.data.get('sequence_type', 'peptide')

        # Validate characters
        if seq_type == 'peptide':
            if not re.match(r'^[ACDEFGHIKLMNPQRSTVWY]+$', sequence):
                invalid = set(sequence) - set("ACDEFGHIKLMNPQRSTVWY")
                raise forms.ValidationError(
                    f"Invalid amino acids: {', '.join(invalid)}",
                    code='invalid_amino_acids'
                )
            if len(sequence) > 20:
                raise forms.ValidationError(
                    "Peptide sequence exceeds 20 residue limit",
                    code='peptide_too_long'
                )
        else:  # Biomolecule
            if not re.match(r'^[ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy]+$', sequence) and \
               not re.match(r'^[ATCGUatcgu]+$', sequence):
                raise forms.ValidationError(
                    "Invalid sequence - must be amino acids or nucleotides",
                    code='invalid_biomolecule'
                )
            if len(sequence) < 50:
                raise forms.ValidationError(
                    "Biomolecule sequence must be ≥50 residues",
                    code='biomolecule_too_short'
                )

        return sequence

# forms.py
from django import forms

class ScreeningForm(forms.Form):
    TARGET_CHOICES = [
        # Populate with (id, display_name) tuples for PredefinedProteinTarget objects
        ('target1', 'Target 1'),
        ('target2', 'Target 2'),
    ]
    
    INPUT_TYPE_CHOICES = [
        ('marvin_drawn', 'MarvinJS Drawing'),
        ('SMILES', 'SMILES String'),
        ('PDB', 'PDB File'),
        ('sequence', 'Biomolecule Sequence'),
    ]
    
    target = forms.ChoiceField(choices=TARGET_CHOICES, label="Select Target")
    input_type = forms.ChoiceField(choices=INPUT_TYPE_CHOICES, label="Input Type")
    compound_input = forms.CharField(widget=forms.Textarea, label="Compound Input")
    
    # For file upload, you might add a FileField for PDB if necessary.
