from django.db import models
from django.shortcuts import render, redirect, get_object_or_404
from django.contrib.auth.decorators import login_required
from django.contrib import messages
from django.core.paginator import Paginator
from django.db.models import Q
from django.utils.dateparse import parse_date
from django.http import HttpResponse
from .forms import SequenceUploadForm
from .models import Sequence, AnalysisJob, DetectedGene
from django.db.models.functions import TruncDate
from django.db.models import Count
from .services.analyzers import analyze_realistic
from .services.reports import generar_pdf
from .constants import HIGH_IDENTITY, HIGH_COVERAGE
import csv

@login_required
def dashboard(request):
    # Datos base
    seqs = Sequence.objects.filter(owner=request.user)
    jobs = AnalysisJob.objects.filter(sequence__owner=request.user)

    # KPIs
    total_sequences = seqs.count()
    total_jobs = jobs.count()
    total_genes = DetectedGene.objects.filter(job__in=jobs).count()

    avg_identity = round(
        jobs.filter(status="DONE").aggregate(models.Avg("identity_pct"))["identity_pct__avg"] or 0, 2
    )
    avg_coverage = round(
        jobs.filter(status="DONE").aggregate(models.Avg("coverage_pct"))["coverage_pct__avg"] or 0, 2
    )

    high_risk = jobs.filter(risk_level="HIGH").count()
    med_risk = jobs.filter(risk_level="MEDIUM").count()
    low_risk = jobs.filter(risk_level="LOW").count()

    # Inicializamos el diccionario antes de agregar los gráficos
    context = {
        "total_sequences": total_sequences,
        "total_jobs": total_jobs,
        "total_genes": total_genes,
        "avg_identity": avg_identity,
        "avg_coverage": avg_coverage,
        "high_risk": high_risk,
        "med_risk": med_risk,
        "low_risk": low_risk,
    }

    # 🔹 Distribución por clase antibiótica
    genes_by_class = (
        DetectedGene.objects.filter(job__in=jobs)
        .values("antibiotic_class")
        .annotate(total=models.Count("id"))
        .order_by("-total")
    )

    # 🔹 Distribución por fuente (CARD / ResFinder)
    genes_by_source = (
        DetectedGene.objects.filter(job__in=jobs)
        .values("source")
        .annotate(total=models.Count("id"))
        .order_by("-total")
    )

    # 🔹 Casos críticos por fecha
    cases_by_date = (
        jobs.filter(risk_level="HIGH")
        .annotate(date=TruncDate("created_at"))
        .values("date")
        .annotate(count=Count("id"))
        .order_by("date")
    )

    # Agregamos estos datos al contexto
    context["genes_by_class"] = list(genes_by_class)
    context["genes_by_source"] = list(genes_by_source)
    context["cases_by_date"] = list(cases_by_date)

    return render(request, "analisis/dashboard.html", context)

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
    """
    Ejecuta el análisis bioinformático principal (modo REAL).
    Basado en librerías de Biopython y patrones de genes RAM simulados.
    """
    seq = get_object_or_404(Sequence, pk=pk, owner=request.user)
    job = AnalysisJob.objects.create(sequence=seq, mode="REAL", status="RUNNING")

    try:
        from .services.analyzers import analyze_realistic
        length, identity, coverage, summary, gene_results = analyze_realistic(seq.fasta_file.path)

        seq.length_bp = length
        seq.save()

        job.identity_pct = identity
        job.coverage_pct = coverage
        job.raw_summary = summary
        job.status = "DONE"
        job.save()

        # Guardar genes y evaluar riesgo
        from .models import DetectedGene
        high_risk_found = False

        for g in gene_results:
            # g = (gene_name, matches, ident, cov, source, abx_class)
            gene_name, matches, ident, cov, source, abx_class = g

            # Clasificación por identidad (visual)
            if ident >= 90:
                level = "Alta resistencia"
            elif ident >= 60:
                level = "Resistencia moderada"
            elif ident >= 30:
                level = "Baja resistencia"
            else:
                level = "Sin resistencia"

            # Reglas de alto riesgo (anomalía)
            if ident >= HIGH_IDENTITY and cov >= HIGH_COVERAGE:
                high_risk_found = True

            DetectedGene.objects.create(
                job=job,
                gene_name=gene_name,
                source=source,
                antibiotic_class=abx_class,
                matches=matches,
                identity=ident,
                coverage=cov,
                classification=level
            )

        job.risk_level = "HIGH" if high_risk_found else ("MEDIUM" if identity>=60 else ("LOW" if identity>=30 else "NONE"))
        job.save()

        messages.success(request, f"✅ Análisis completado: {identity}% identidad, {coverage}% cobertura")
        return redirect("resultados", pk=job.pk)
    
    except Exception as e:
        # Manejo de errores
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
    
    # ⚙️ Nuevo: filtro por nivel de riesgo
    risk = request.GET.get('risk', '').upper()
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
        
    # ⚙️ Nuevo: filtrar por riesgo si fue seleccionado
    if risk in {'NONE', 'LOW', 'MEDIUM', 'HIGH'}:
        jobs = jobs.filter(risk_level=risk)
        
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
        "risk": risk, # 🔽 Nuevo
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

