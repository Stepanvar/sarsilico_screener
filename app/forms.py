from django import forms
from .models import Target
from django.contrib.auth.forms import AuthenticationForm, UserCreationForm
from django.contrib.auth.models import User
from django.utils.translation import gettext_lazy as _


class SMILESInputForm(forms.Form):
    smiles = forms.CharField(
        label='SMILES String', max_length=255,
        widget=forms.TextInput(attrs={'placeholder': 'Enter SMILES string'})
    )
    smiles_file = forms.FileField(
        label='Upload SMILES File',
        required=False,
        help_text='Supported formats: .smi, .smiles, .mol2'
    )

class PDBUploadForm(forms.Form):
    pdb_file = forms.FileField(label='Upload PDB File', help_text='Supported format: .pdb')

class SimilarityForm(forms.Form):
    smiles = forms.CharField(
        label='SMILES String',
        max_length=255,
        widget=forms.TextInput(attrs={'placeholder': 'Enter SMILES string'}),
        help_text='Provide a valid SMILES string for similarity analysis.'
    )
    threshold = forms.FloatField(
        label='Similarity Threshold',
        min_value=0.0,
        max_value=1.0,
        initial=0.7,
        help_text='Set the minimum Tanimoto similarity score (e.g., 0.7).'
    )

class DockingForm(forms.Form):
    target = forms.ModelChoiceField(queryset=Target.objects.all(), label='Protein Target')
    smiles = forms.CharField(
        label='SMILES String',
        max_length=255,
        widget=forms.TextInput(attrs={'placeholder': 'Enter SMILES string'})
    )

class SequenceForm(forms.Form):
    sequence = forms.CharField(label='Sequence', max_length=500, required=False)
    sequence_file = forms.FileField(label='Sequence CIF File', required=False)

class UserRegistrationForm(UserCreationForm):
    email = forms.EmailField(required=True, help_text='Required. Enter a valid email address.')
    class Meta:
        model = User
        fields = ['username', 'email', 'password1', 'password2']

    def clean_email(self):
        email = self.cleaned_data.get('email')
        if User.objects.filter(email=email).exists():
            raise forms.ValidationError('Email address already in use.')
        return email

class TargetSelectionForm(forms.Form):
    def __init__(self, *args, **kwargs):
        targets = kwargs.pop('targets', [])
        super().__init__(*args, **kwargs)
        self.fields['target_id'] = forms.ChoiceField(
            choices=[(t.id, t.name) for t in targets],
            label='Select Target Protein',
            widget=forms.Select(attrs={'class':'form-control'})
        )