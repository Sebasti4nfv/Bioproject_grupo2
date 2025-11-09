"""
Detector de Genes de Resistencia Antimicrobiana (RAM)
Aplicación web con Streamlit - MVP
"""

from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Necesario para Streamlit Cloud
import matplotlib.pyplot as plt
from pathlib import Path
import time
import random
import io
import streamlit as st
# Configuración de la página
st.set_page_config(
    page_title="Detector RAM",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

print("✓ Módulos importados correctamente")


# ==================== BASE DE DATOS ====================
# NOTA: Secuencias dummy para demostración.
# En producción, reemplazar por CARD/ResFinder.

REF_DB = {
    'blaNDM': {
        'seq': 'ATGGAATTGCCCAATATTATGCACCCGGTCGCGAAGCTGAGCACCGCATTGGGCAATGTGCAACATACCGTGCTGCGCTTGATCGGCAACGGCGACCTGGCCATCCTGCATCCGCTGACATGGGGATGGCGTTCGTCGGTCACCAGCAGTTCGTC',
        'antibiotic_class': 'Carbapenems',
        'mechanism': 'beta-lactamase'
    },
    'blaKPC': {
        'seq': 'ATGTCACTGTATCGCCGTCTAGTTCTGCTGTCTTGTCTCTCATGGCCGCTGGCTGGCTTTTCTGCCACCGCGCTGACCAACCTCGTCGCGGAACCATTCGCTAAACTCGAACAGGACTTTGGCGGCTCCATCGGTGTGTACGCGATGGATC',
        'antibiotic_class': 'Carbapenems',
        'mechanism': 'beta-lactamase'
    },
    'mcr-1': {
        'seq': 'ATGCAGCATACTTCTGTGTGGTCTGTAGCGTCTTGGTTATTCTGGTCGCACTTGCCGCAATTATCGCGATGATTCTGGTCGTAGCTATTCTGGCGCTGGCAATGGTGATTGCCGCTATTCTGGCTTACTTCTTTCTGCTG',
        'antibiotic_class': 'Polymyxins',
        'mechanism': 'phosphoethanolamine transferase'
    },
    'qnrS': {
        'seq': 'ATGAGCAATTTCGGGAATGCACATAGTATTCCCGCAAAGCAAGTTGCTATCGATGTAGCGAAGGCTTTTCTCAAGCTCGCGAAGATGCAGGCAGCTGGCGATATCGCCAGCCTGCATAATAACGTTCTGCAGCTGCAGGCG',
        'antibiotic_class': 'Fluoroquinolones',
        'mechanism': 'target protection'
    },
    'tetO': {
        'seq': 'ATGAAATTAATAACTTTTATTGAACATCAACTTCAAGGTGTTGAAGATATTCTTCAGCGTTTGCATAGACCAATTAGTCAGGTTGTTGAGAGAAATCAGCTCAGGACTAATAGATTAATTGATTTGGATGGTTTAGAAGTTAGAGGTCAG',
        'antibiotic_class': 'Tetracyclines',
        'mechanism': 'ribosomal protection'
    },
    'ermB': {
        'seq': 'ATGAACATAAAATTAATATCACCTATACATCATTCACAGTTATTATCAATAATTCGCTTTTATGGCATTTTTGGAAAGATTTATATTATCATAGACACATTTAAACATAGAGAAACTTATTACAAACGATTTTACAAAACTTATAAAAGAC',
        'antibiotic_class': 'Macrolides',
        'mechanism': 'methyltransferase'
    },
    'gyrA_S83L': {
        'seq': 'ATGAGCGACCTTGCGAGAGAAATTACACCGGTCAACATTGAGGAAGACAGCTATGACATGATTCGCCGCCTGGGCGGCGGCACGTTCCAGTTCCGTAGCCTGACCGAGGACGGTAAACGGATGAACGTCCTGCTGCTGGAC',
        'antibiotic_class': 'Fluoroquinolones',
        'mechanism': 'target mutation'
    }
}

# ==================== FUNCIONES ====================

def read_fasta_from_string(fasta_content: str) -> str:
    """Lee contenido FASTA de un string y retorna secuencia concatenada."""
    sequence = []
    for line in fasta_content.split('\n'):
        line = line.strip()
        if not line.startswith('>') and line:
            sequence.append(line.upper())
    return ''.join(sequence)


def make_demo_fasta(ref_db: dict, genome_len: int = 4000, n_inserts: int = 3) -> tuple:
    """Genera secuencia demo con genes insertados."""
    bases = ['A', 'T', 'G', 'C']
    genome = [random.choice(bases) for _ in range(genome_len)]
    
    genes_to_insert = random.sample(list(ref_db.keys()), min(n_inserts, len(ref_db)))
    
    for gene_name in genes_to_insert:
        gene_seq = ref_db[gene_name]['seq']
        max_pos = len(genome) - len(gene_seq)
        if max_pos > 0:
            insert_pos = random.randint(0, max_pos)
            for i, base in enumerate(gene_seq):
                genome[insert_pos + i] = base
    
    genome_str = ''.join(genome)
    fasta_str = ">demo_bacterial_genome\n"
    for i in range(0, len(genome_str), 80):
        fasta_str += genome_str[i:i+80] + '\n'
    
    return fasta_str, genes_to_insert


def align_and_score(query: str, ref: str) -> dict:
    """Alinea query vs ref y calcula métricas."""
    alignments = pairwise2.align.localms(query, ref, 2, -1, -2, -1, one_alignment_only=True)
    
    if not alignments:
        return {
            'identity': 0.0,
            'coverage': 0.0,
            'start': -1,
            'end': -1,
            'alignment_length': 0
        }
    
    best = alignments[0]
    seq_a, seq_b, score, start, end = best
    
    matches = sum(1 for a, b in zip(seq_a[start:end], seq_b[start:end]) if a == b and a != '-')
    alignment_length = end - start
    
    identity = matches / alignment_length if alignment_length > 0 else 0.0
    coverage = alignment_length / len(ref) if len(ref) > 0 else 0.0
    
    return {
        'identity': identity,
        'coverage': coverage,
        'start': start,
        'end': end,
        'alignment_length': alignment_length
    }


def detect_genes(query_seq: str, ref_db: dict, id_thr: float = 0.90, cov_thr: float = 0.80) -> pd.DataFrame:
    """Detecta genes en query_seq."""
    results = []
    
    progress_bar = st.progress(0)
    status_text = st.empty()
    
    total_genes = len(ref_db)
    for idx, (gene_name, gene_data) in enumerate(ref_db.items()):
        status_text.text(f"Analizando {gene_name}...")
        progress_bar.progress((idx + 1) / total_genes)
        
        ref_seq = gene_data['seq']
        metrics = align_and_score(query_seq, ref_seq)
        
        if metrics['identity'] >= id_thr and metrics['coverage'] >= cov_thr:
            results.append({
                'gene': gene_name,
                'identity': metrics['identity'],
                'coverage': metrics['coverage'],
                'start': metrics['start'],
                'end': metrics['end'],
                'length_ref': len(ref_seq),
                'antibiotic_class': gene_data['antibiotic_class'],
                'mechanism': gene_data['mechanism']
            })
    
    progress_bar.empty()
    status_text.empty()
    
    df = pd.DataFrame(results)
    if not df.empty:
        df = df.sort_values(['identity', 'coverage'], ascending=False).reset_index(drop=True)
    
    return df


def plot_by_class(df: pd.DataFrame):
    """Genera gráfico de barras por clase de antibiótico."""
    if df.empty:
        return None
    
    class_counts = df['antibiotic_class'].value_counts()
    
    fig, ax = plt.subplots(figsize=(10, 6))
    class_counts.plot(kind='bar', ax=ax, color='steelblue')
    ax.set_xlabel('Clase de antibiótico', fontsize=12)
    ax.set_ylabel('Número de genes detectados', fontsize=12)
    ax.set_title('Genes de resistencia detectados por clase de antibiótico', fontsize=14, fontweight='bold')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    
    return fig


# ==================== INTERFAZ STREAMLIT ====================

# Header
st.title("🧬 Detector de Genes de Resistencia Antimicrobiana")
st.markdown("**MVP** - Detección rápida de genes RAM en secuencias FASTA")
st.markdown("---")

# Sidebar
with st.sidebar:
    st.header("⚙️ Configuración")
    
    st.subheader("📊 Base de Datos")
    st.info(f"**{len(REF_DB)} genes** de referencia cargados")
    
    with st.expander("Ver genes disponibles"):
        for gene, data in REF_DB.items():
            st.write(f"**{gene}**: {data['antibiotic_class']}")
    
    st.markdown("---")
    
    st.subheader("🎯 Umbrales de detección")
    identity_threshold = st.slider(
        "Identidad mínima (%)",
        min_value=50,
        max_value=100,
        value=90,
        step=5,
        help="Porcentaje de coincidencia entre secuencias"
    ) / 100
    
    coverage_threshold = st.slider(
        "Cobertura mínima (%)",
        min_value=50,
        max_value=100,
        value=80,
        step=5,
        help="Porcentaje del gen de referencia que debe estar presente"
    ) / 100
    
    st.markdown("---")
    st.caption("💡 **Nota**: Secuencias de referencia son dummy para demostración")

# Main content
tab1, tab2, tab3 = st.tabs(["📁 Análisis", "ℹ️ Información", "🧪 Demo"])

with tab1:
    st.header("Cargar secuencia FASTA")
    
    col1, col2 = st.columns([2, 1])
    
    with col1:
        uploaded_file = st.file_uploader(
            "Sube tu archivo FASTA",
            type=['fasta', 'fa', 'fna', 'txt'],
            help="Archivo con secuencia(s) en formato FASTA"
        )
    
    with col2:
        st.markdown("**O usa datos demo:**")
        use_demo = st.button("🎲 Generar secuencia demo", use_container_width=True)
    
    # Procesamiento
    if uploaded_file is not None or use_demo:
        
        if use_demo:
            st.info("🔄 Generando secuencia demo con genes insertados...")
            fasta_content, inserted_genes = make_demo_fasta(REF_DB, genome_len=4000, n_inserts=3)
            st.success(f"✅ Demo generado con: {', '.join(inserted_genes)}")
        else:
            fasta_content = uploaded_file.read().decode('utf-8')
            st.success(f"✅ Archivo cargado: {uploaded_file.name}")
        
        # Leer secuencia
        query_seq = read_fasta_from_string(fasta_content)
        
        # Mostrar info de secuencia
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Longitud", f"{len(query_seq):,} nt")
        with col2:
            gc_content = (query_seq.count('G') + query_seq.count('C')) / len(query_seq) * 100
            st.metric("GC%", f"{gc_content:.1f}%")
        with col3:
            st.metric("Umbrales", f"ID≥{identity_threshold*100:.0f}% COV≥{coverage_threshold*100:.0f}%")
        
        # Vista previa
        with st.expander("🔍 Vista previa de la secuencia"):
            st.code(query_seq[:500] + "..." if len(query_seq) > 500 else query_seq)
        
        st.markdown("---")
        
        # Botón de análisis
        if st.button("🔬 **INICIAR ANÁLISIS**", type="primary", use_container_width=True):
            
            start_time = time.time()
            
            with st.spinner("Analizando secuencia..."):
                results_df = detect_genes(query_seq, REF_DB, identity_threshold, coverage_threshold)
            
            elapsed = time.time() - start_time
            
            # Resultados
            st.markdown("---")
            st.header("📊 Resultados")
            
            if results_df.empty:
                st.warning("❌ No se detectaron genes que cumplan los umbrales especificados")
                st.info("💡 Intenta reducir los umbrales de identidad o cobertura")
            else:
                # Métricas
                col1, col2, col3, col4 = st.columns(4)
                with col1:
                    st.metric("Genes detectados", len(results_df))
                with col2:
                    st.metric("Clases antibióticas", results_df['antibiotic_class'].nunique())
                with col3:
                    avg_identity = results_df['identity'].mean() * 100
                    st.metric("Identidad promedio", f"{avg_identity:.1f}%")
                with col4:
                    st.metric("Tiempo", f"{elapsed:.2f}s")
                
                # Tabla de resultados
                st.subheader("🧬 Genes detectados")
                
                # Formatear DataFrame para mostrar
                display_df = results_df.copy()
                display_df['identity'] = display_df['identity'].apply(lambda x: f"{x*100:.1f}%")
                display_df['coverage'] = display_df['coverage'].apply(lambda x: f"{x*100:.1f}%")
                
                st.dataframe(
                    display_df,
                    use_container_width=True,
                    hide_index=True
                )
                
                # Gráfico
                st.subheader("📈 Distribución por clase de antibiótico")
                fig = plot_by_class(results_df)
                if fig:
                    st.pyplot(fig)
                
                # Descargar CSV
                csv = results_df.to_csv(index=False)
                st.download_button(
                    label="💾 Descargar resultados (CSV)",
                    data=csv,
                    file_name=f"ram_results_{time.strftime('%Y%m%d_%H%M%S')}.csv",
                    mime="text/csv",
                    use_container_width=True
                )
                
                # Genes prioritarios
                priority_genes = ['blaNDM', 'blaKPC', 'mcr-1']
                detected_priority = results_df[results_df['gene'].isin(priority_genes)]
                
                if not detected_priority.empty:
                    st.markdown("---")
                    st.error("⚠️ **ALERTA: Genes prioritarios detectados**")
                    st.dataframe(detected_priority[['gene', 'antibiotic_class']], hide_index=True)

with tab2:
    st.header("ℹ️ Información del sistema")
    
    st.markdown("""
    ### 🎯 Objetivo
    Este MVP detecta genes de resistencia antimicrobiana (RAM) en secuencias bacterianas
    mediante alineamiento local contra una base de datos de referencia.
    
    ### 🔬 Metodología
    1. **Alineamiento**: Se usa `Bio.pairwise2` (algoritmo Smith-Waterman local)
    2. **Métricas calculadas**:
       - **Identidad**: % de bases coincidentes en la región alineada
       - **Cobertura**: % del gen de referencia presente en la secuencia query
    3. **Detección**: Un gen se reporta si cumple ambos umbrales
    
    ### 📚 Base de datos actual
    - **7 genes** de referencia (secuencias dummy de 80-150 nt)
    - Clases incluidas: Carbapenems, Polymyxins, Fluoroquinolones, Tetracyclines, Macrolides
    
    ### ⚡ Limitaciones del MVP
    - ⚠️ Secuencias de referencia son **ejemplos ficticios**
    - ⚠️ No detecta mutaciones puntuales (excepto gyrA_S83L como ejemplo)
    - ⚠️ Sensibilidad limitada a variantes alélicas
    
    ### 🚀 Próxima versión
    - Integración con **CARD** y **ResFinder**
    - Uso de **BLAST+** para búsquedas más rápidas
    - **RGI** para predicción fenotípica
    - Reportes PDF con trazabilidad completa
    - Alertas automáticas para genes prioritarios
    
    ### 📖 Referencias
    - CARD: https://card.mcmaster.ca/
    - ResFinder: https://cge.food.dtu.dk/services/ResFinder/
    """)

with tab3:
    st.header("🧪 Demo interactiva")
    
    st.markdown("""
    Esta demo genera una secuencia bacteriana simulada con genes RAM insertados
    para demostrar la funcionalidad del detector.
    """)
    
    col1, col2 = st.columns(2)
    
    with col1:
        demo_length = st.number_input("Longitud del genoma (nt)", 1000, 10000, 4000, 500)
    with col2:
        demo_inserts = st.slider("Número de genes a insertar", 1, 7, 3)
    
    if st.button("🎲 Generar y analizar demo", type="primary"):
        with st.spinner("Generando secuencia demo..."):
            fasta_demo, genes_inserted = make_demo_fasta(REF_DB, demo_length, demo_inserts)
            query_demo = read_fasta_from_string(fasta_demo)
        
        st.success(f"✅ Secuencia generada: {len(query_demo)} nt")
        st.info(f"🧬 Genes insertados: **{', '.join(genes_inserted)}**")
        
        with st.spinner("Analizando..."):
            demo_results = detect_genes(query_demo, REF_DB, identity_threshold, coverage_threshold)
        
        if not demo_results.empty:
            st.success(f"✅ Detectados {len(demo_results)} genes")
            
            display_demo = demo_results.copy()
            display_demo['identity'] = display_demo['identity'].apply(lambda x: f"{x*100:.1f}%")
            display_demo['coverage'] = display_demo['coverage'].apply(lambda x: f"{x*100:.1f}%")
            
            st.dataframe(display_demo[['gene', 'identity', 'coverage', 'antibiotic_class']], hide_index=True)
            
            if len(demo_results) >= 2:
                st.balloons()
                st.success("🎉 **Validación exitosa**: Se detectaron ≥2 genes como esperado")
        else:
            st.warning("No se detectaron genes. Intenta reducir los umbrales.")

# Footer
st.markdown("---")

st.caption("🧬 Detector RAM MVP v1.0 | Desarrollado con Streamlit + Biopython")
