from Bio import SeqIO, pairwise2

def analyze_demo(fasta_path):
    rec = next(SeqIO.parse(fasta_path, "fasta"))
    seq = rec.seq
    # Variante “simulada”: reverso complementario
    alignment = pairwise2.align.globalxx(seq, seq.reverse_complement(), one_alignment_only=True)[0]
    matches = alignment[2]
    identity = round(100.0 * matches / len(seq), 2)
    coverage = 100.0  # en demo asumimos cobertura completa
    summary = f"Len={len(seq)}; Matches={matches}; Identity={identity}%"
    return len(seq), identity, coverage, summary
