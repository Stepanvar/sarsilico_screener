# forms.py
from django import forms

class ScreeningForm(forms.Form):
    target = forms.ChoiceField(choices=[
        ('Nonstructural', 'Nonstructural proteins'),
        ('Structure','Structure proteins')
    ])
    input_type = forms.ChoiceField(
        choices=[
            ('jsme_drawn', 'JSME Drawing'),
            ('SMILES', 'SMILES String'),
            ('PDB', 'PDB File'),
            ('sequence', 'Biomolecule Sequence')
        ]
    )
    compound_input = forms.CharField(widget=forms.Textarea)
