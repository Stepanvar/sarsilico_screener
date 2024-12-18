from django import forms
from django.contrib.auth.forms import AuthenticationForm, UserCreationForm
from django.utils.translation import gettext_lazy as _
from app.models import Compound

class BootstrapAuthenticationForm(AuthenticationForm):
    username = forms.CharField(max_length=254,
                               widget=forms.TextInput({
                                   'class': 'form-control',
                                   'placeholder': 'User name'}))
    password = forms.CharField(label=_("Password"),
                               widget=forms.PasswordInput({
                                   'class': 'form-control',
                                   'placeholder':'Password'}))
# class BootstrapUserCreationForm(UserCreationForm):
#     class Meta:
#         fields = ('username', 'email', 'password1', 'password2')
#     email = forms.EmailField(
# 		    max_length=254,
# 		    widget=forms.EmailInput({
# 			    'class': 'form-control',
# 			    'placeholder': 'Email address'}))

class CompoundForm(forms.ModelForm):
    class Meta:
        model = Compound
        fields = ['name', 'smiles', 'pdb_file', 'structure_type']
        widgets = {
            'structure_type': forms.RadioSelect,
        }
    def clean_smiles(self):
        from rdkit import Chem
        smiles = self.cleaned_data.get('smiles')
        if smiles and not Chem.MolFromSmiles(smiles):
            raise forms.ValidationError("Invalid SMILES string.")
        return smiles

    def clean_pdb_file(self):
        pdb_file = self.cleaned_data.get('pdb_file')
        if pdb_file and not pdb_file.name.endswith('.pdb'):
            raise forms.ValidationError("File is not a valid PDB format.")
        return pdb_file

class SMILESInputForm(forms.Form):
    smiles = forms.CharField(widget=forms.Textarea(attrs={'rows':3, 'class':'form-control'}), required=False, label='Enter SMILES')
    smiles_file = forms.FileField(required=False, label='Upload SMILES File', widget=forms.ClearableFileInput(attrs={'class':'form-control'}))

class PDBUploadForm(forms.Form):
    pdb_file = forms.FileField(required=True, label='Upload PDB File', widget=forms.ClearableFileInput(attrs={'class':'form-control'}))

class SequenceForm(forms.Form):
    sequence = forms.CharField(widget=forms.Textarea(attrs={'rows':3, 'class':'form-control'}), required=True, label='Enter Peptide/Nucleotide Sequence')

class TargetSelectionForm(forms.Form):
    def __init__(self, *args, **kwargs):
        targets = kwargs.pop('targets', [])
        super().__init__(*args, **kwargs)
        self.fields['target_id'] = forms.ChoiceField(
            choices=[(t.id, t.name) for t in targets],
            label='Select Target Protein',
            widget=forms.Select(attrs={'class':'form-control'})
        )

class SimilarityForm(forms.Form):
    smiles = forms.CharField(widget=forms.Textarea(attrs={'rows':3, 'class':'form-control'}), required=True, label='Enter SMILES')
    threshold = forms.FloatField(initial=0.7, min_value=0.5, max_value=1.0, label='Set Similarity Threshold', widget=forms.NumberInput(attrs={'class':'form-control'}))

