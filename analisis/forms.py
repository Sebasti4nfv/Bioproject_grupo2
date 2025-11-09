from django import forms
from .models import Sequence

class SequenceUploadForm(forms.ModelForm):
    class Meta:
        model = Sequence
        fields = ["name", "fasta_file"]

    def clean_fasta_file(self):
        f = self.cleaned_data["fasta_file"]
        ok = f.name.lower().endswith((".fa",".fasta"))
        if not ok:
            raise forms.ValidationError("El archivo debe tener extensión .fa o .fasta")
        return f
