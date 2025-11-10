import re
from Bio import SeqIO
from analisis.models import ResistanceGene
# Genes RAM simulados del documento Avance 2
GENES_CARD = {
    "blaTEM": "ATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCTCACCCAGAAACGCTGGTGAAAGTA",
    "aac(3)-IIa": "ATGAACAAAACTATTCCTCTGGTACGTTAGCTATCCTGATGGATAAGATGGCGGGATACG",
    "tetA": "ATGTTGATAAAGCATTGGAATTACAGAGCGATCCTATCAACGAGGTTTTCTTGGAGTTGC",
    "mecA": "ATGCTAAAGTTCAAAAGAGTTCACCTTTCTTTTTAGCCAGTTTCCTTATCCTGCTGGTAA"
}

def analyze_realistic(fasta_path):
    """
    Pipeline bioinformático con datos reales de CARD / ResFinder.
    Busca coincidencias en la base de datos ResistanceGene.
    """
    record = next(SeqIO.parse(fasta_path, "fasta"))
    seq = str(record.seq).upper()
    length = len(seq)
    gc_content = round((seq.count("G") + seq.count("C")) / length * 100, 2)

    genes = ResistanceGene.objects.all()
    results = []

    for gene in genes:
        fragment = gene.sequence[:30]
        found = len(re.findall(fragment, seq))
        if found > 0:
            identity = round(min((found * len(fragment)) / length * 100, 100), 2)
            coverage = round(len(fragment) / length * 100, 2)
            results.append({
                "gene": gene.gene_name,
                "source": gene.source,
                "class": gene.antibiotic_class,
                "matches": found,
                "identity": identity,
                "coverage": coverage
            })

    if not results:
        summary = f"Secuencia sin coincidencias RAM (CARD/ResFinder). Longitud={length}bp, GC={gc_content}%"
        return length, 0.0, 0.0, summary, []

    identity_avg = round(sum(r["identity"] for r in results) / len(results), 2)
    coverage_avg = round(sum(r["coverage"] for r in results) / len(results), 2)

    summary_lines = [
        f"Secuencia analizada: {record.id}",
        f"Largo: {length} bp | GC: {gc_content}%",
        f"Genes RAM detectados:"
    ]
    for r in results:
        summary_lines.append(
            f" - {r['gene']} ({r['source']}): clase={r['class']} | identidad={r['identity']}% | cobertura={r['coverage']}%"
        )
    summary_lines.append(f"Promedio global → Identidad={identity_avg}% | Cobertura={coverage_avg}%")

    return length, identity_avg, coverage_avg, "\n".join(summary_lines), [
        (r["gene"], r["matches"], r["identity"], r["coverage"], r["source"], r["class"]) for r in results
    ]