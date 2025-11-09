from django.shortcuts import render, redirect, get_object_or_404
from django.contrib.auth.decorators import login_required
from django.contrib import messages
from django.core.paginator import Paginator
from django.db.models import Q
from django.utils.dateparse import parse_date
from django.http import HttpResponse
from .forms import SequenceUploadForm
from .models import Sequence, AnalysisJob
from .services.analyzers import analyze_demo
from .services.reports import generar_pdf
import csv

@login_required
def dashboard(request):
    seqs = Sequence.objects.filter(owner=request.user).order_by("-created_at")
    jobs = AnalysisJob.objects.filter(sequence__owner=request.user).order_by("-created_at")[:10]
    kpis = {
        "total_sequences": seqs.count(),
        "total_jobs": AnalysisJob.objects.filter(sequence__owner=request.user).count(),
        "avg_identity": round(AnalysisJob.objects.filter(sequence__owner=request.user, status="DONE").aggregate_avg("identity_pct") or 0, 2) if hasattr(AnalysisJob.objects, "aggregate_avg") else None
    }
    return render(request, "analisis/dashboard.html", {"seqs": seqs, "jobs": jobs, "kpis": kpis})

@login_required
def upload_sequence(request):
    if request.method == "POST":
        form = SequenceUploadForm(request.POST, request.FILES)
        if form.is_valid():
            seq = form.save(commit=False)
            seq.owner = request.user
            seq.save()
            messages.success(request, "✅ Secuencia cargada correctamente.")
            return redirect("sequence_detail", pk=seq.pk)  # 👈 redirección tras subir
        else:
            messages.error(request, "❌ Archivo no válido o faltan datos.")
    else:
        form = SequenceUploadForm()
    
    return render(request, "analisis/upload.html", {"form": form})
@login_required
def sequence_detail(request, pk):
    seq = get_object_or_404(Sequence, pk=pk, owner=request.user)
    jobs = seq.jobs.order_by('-created_at')

    return render(request, "analisis/detalle.html", {
        "seq": seq,
        "jobs": jobs
    })

@login_required
def run_analysis(request, pk):
    seq = get_object_or_404(Sequence, pk=pk, owner=request.user)
    job = AnalysisJob.objects.create(sequence=seq, mode="DEMO", status="RUNNING")

    try:
        length, identity, coverage, summary = analyze_demo(seq.fasta_file.path)
        seq.length_bp = length
        seq.save()

        job.identity_pct = identity
        job.coverage_pct = coverage
        job.raw_summary = summary
        job.status = "DONE"
        job.save()

        messages.success(request, f"✅ Análisis completado: {identity}% identidad, {coverage}% cobertura")
        return redirect("resultados", pk=job.pk)
    except Exception as e:
        job.status = "ERROR"
        job.raw_summary = str(e)
        job.save()
        messages.error(request, f"❌ Error en análisis: {e}")
        return redirect("sequence_detail", pk=seq.pk)
@login_required
def resultados(request, pk):
    job = get_object_or_404(AnalysisJob, pk=pk, sequence__owner=request.user)
    return render(request, "analisis/resultados.html", {"job": job})
@login_required
def exportar_pdf(request, pk):
    job = get_object_or_404(AnalysisJob, pk=pk, sequence__owner=request.user)
    context = {"job": job}
    return generar_pdf("analisis/resultados_pdf.html", context)
@login_required
def historial(request):
    """
    Lista de AnalysisJob con filtros:
    - q: busca por nombre de secuencia
    - status: PENDING/RUNNING/DONE/ERROR
    - date_from, date_to: rango de fechas (YYYY-MM-DD)
    - owner_scope: 'mine' (default) o 'all' (si is_staff)
    Paginación: 10 por página
    """
    q = request.GET.get('q', '').strip()
    status = request.GET.get('status', '').strip().upper()
    date_from = request.GET.get('date_from', '').strip()
    date_to = request.GET.get('date_to', '').strip()
    owner_scope = request.GET.get('owner_scope', 'mine')

    # Base queryset
    jobs = AnalysisJob.objects.select_related('sequence').order_by('-created_at')

    # Alcance por rol
    if request.user.is_staff and owner_scope == 'all':
        pass  # ve todo
    else:
        jobs = jobs.filter(sequence__owner=request.user)

    # Filtro por texto
    if q:
        jobs = jobs.filter(Q(sequence__name__icontains=q) | Q(raw_summary__icontains=q))

    # Filtro por estado
    if status in {'PENDING', 'RUNNING', 'DONE', 'ERROR'}:
        jobs = jobs.filter(status=status)

    # Filtro por fechas
    df = parse_date(date_from) if date_from else None
    dt = parse_date(date_to) if date_to else None
    if df:
        jobs = jobs.filter(created_at__date__gte=df)
    if dt:
        jobs = jobs.filter(created_at__date__lte=dt)

    # Paginación
    paginator = Paginator(jobs, 10)
    page_number = request.GET.get('page')
    page_obj = paginator.get_page(page_number)

    return render(request, "analisis/historial.html", {
        "page_obj": page_obj,
        "q": q,
        "status": status,
        "date_from": date_from,
        "date_to": date_to,
        "owner_scope": owner_scope,
    })


@login_required
def export_historial_csv(request):
    """
    Exporta a CSV aplicando los mismos filtros del historial.
    """
    q = request.GET.get('q', '').strip()
    status = request.GET.get('status', '').strip().upper()
    date_from = request.GET.get('date_from', '').strip()
    date_to = request.GET.get('date_to', '').strip()
    owner_scope = request.GET.get('owner_scope', 'mine')

    jobs = AnalysisJob.objects.select_related('sequence').order_by('-created_at')
    if request.user.is_staff and owner_scope == 'all':
        pass
    else:
        jobs = jobs.filter(sequence__owner=request.user)

    if q:
        jobs = jobs.filter(Q(sequence__name__icontains=q) | Q(raw_summary__icontains=q))

    if status in {'PENDING', 'RUNNING', 'DONE', 'ERROR'}:
        jobs = jobs.filter(status=status)

    df = parse_date(date_from) if date_from else None
    dt = parse_date(date_to) if date_to else None
    if df:
        jobs = jobs.filter(created_at__date__gte=df)
    if dt:
        jobs = jobs.filter(created_at__date__lte=dt)

    # Generar CSV
    response = HttpResponse(content_type='text/csv')
    response['Content-Disposition'] = 'attachment; filename="historial_analisis.csv"'
    writer = csv.writer(response)
    writer.writerow(['Fecha', 'Secuencia', 'Estado', 'Identidad (%)', 'Cobertura (%)', 'Resumen'])

    for j in jobs:
        writer.writerow([
            j.created_at.strftime("%Y-%m-%d %H:%M"),
            j.sequence.name,
            j.status,
            j.identity_pct if j.identity_pct is not None else '',
            j.coverage_pct if j.coverage_pct is not None else '',
            (j.raw_summary or '')[:200],
        ])

    return response