from django.db import models
from django.contrib.auth.models import User

class Sequence(models.Model):
    owner = models.ForeignKey(User, on_delete=models.CASCADE)
    name = models.CharField(max_length=200)
    fasta_file = models.FileField(upload_to="fasta/")
    length_bp = models.IntegerField(null=True, blank=True)
    created_at = models.DateTimeField(auto_now_add=True)

    def __str__(self): return self.name

class AnalysisJob(models.Model):
    MODE_CHOICES = (("DEMO","DEMO"), ("REAL","REAL"))
    sequence = models.ForeignKey(Sequence, on_delete=models.CASCADE, related_name="jobs")
    mode = models.CharField(max_length=8, choices=MODE_CHOICES, default="DEMO")
    status = models.CharField(max_length=20, default="PENDING")  # PENDING/RUNNING/DONE/ERROR
    identity_pct = models.FloatField(null=True, blank=True)
    coverage_pct = models.FloatField(null=True, blank=True)
    raw_summary = models.TextField(null=True, blank=True)
    created_at = models.DateTimeField(auto_now_add=True)
