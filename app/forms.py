from django import forms
from .models import Target
from .utils import validate_smiles, validate_pdb_file

class SMILESInputForm(forms.Form):
    smiles = forms.CharField(
        label='SMILES String',
        widget=forms.TextInput(attrs={'class': 'form-control'}),
        required=False
    )
    smiles_file = forms.FileField(
        label='SMILES File',
        widget=forms.FileInput(attrs={'class': 'form-control'}),
        required=False
    )

    def clean(self):
        cleaned_data = super().clean()
        smiles = cleaned_data.get('smiles')
        smiles_file = cleaned_data.get('smiles_file')
        
        if not smiles and not smiles_file:
            raise forms.ValidationError("Either SMILES string or file must be provided")
        
        return cleaned_data

class PDBUploadForm(forms.Form):
    pdb_file = forms.FileField(
        label='PDB File',
        widget=forms.FileInput(attrs={'class': 'form-control'})
    )

class SimilarityForm(forms.Form):
    smiles = forms.CharField(
        label='SMILES String',
        widget=forms.TextInput(attrs={'class': 'form-control'})
    )
    threshold = forms.FloatField(
        label='Similarity Threshold',
        min_value=0.0,
        max_value=1.0,
        initial=0.85,
        widget=forms.NumberInput(attrs={'class': 'form-control'})
    )

class DockingForm(forms.Form):
    target = forms.ModelChoiceField(
        queryset=Target.objects.all(),
        widget=forms.Select(attrs={'class': 'form-control'})
    )

class SequenceForm(forms.Form):
    sequence = forms.CharField(
        label='Protein Sequence',
        widget=forms.Textarea(attrs={'class': 'form-control', 'rows': 4}),
        required=False
    )
    sequence_file = forms.FileField(
        label='Sequence File',
        widget=forms.FileInput(attrs={'class': 'form-control'}),
        required=False
    )
