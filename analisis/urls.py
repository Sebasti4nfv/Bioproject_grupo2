from django.urls import path
from . import views

urlpatterns = [
    path('', views.dashboard, name='dashboard'),
    path('upload/', views.upload_sequence, name='upload'),
    path('seq/<int:pk>/', views.sequence_detail, name='sequence_detail'),
    path('run/<int:pk>/', views.run_analysis, name='run_analysis'),
    path('resultados/<int:pk>/', views.resultados, name='resultados'),
    # 🔽 NUEVO
    path('historial/', views.historial, name='historial'),
    path('historial/export.csv', views.export_historial_csv, name='export_historial_csv'),
        # ✅ Exportar PDF
    path('exportar-pdf/<int:pk>/', views.exportar_pdf, name='exportar_pdf'),
]
