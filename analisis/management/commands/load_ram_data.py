from django.core.management.base import BaseCommand
from analisis.models import ResistanceGene

class Command(BaseCommand):
    help = "Carga genes RAM simulados de CARD y ResFinder"

    def handle(self, *args, **kwargs):
        genes = [
            ("CARD", "blaTEM", "ATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCTCACCCAGAAACGCTGGTGAAAGTA", "Beta-lactámicos", "Beta-lactamasa TEM-1, confiere resistencia a penicilinas"),
            ("CARD", "mecA", "ATGCTAAAGTTCAAAAGAGTTCACCTTTCTTTTTAGCCAGTTTCCTTATCCTGCTGGTAA", "Meticilina", "Gen mecA, resistencia a meticilina"),
            ("ResFinder", "tetA", "ATGTTGATAAAGCATTGGAATTACAGAGCGATCCTATCAACGAGGTTTTCTTGGAGTTGC", "Tetraciclinas", "Transportador de eflujo TetA"),
            ("ResFinder", "aac(3)-IIa", "ATGAACAAAACTATTCCTCTGGTACGTTAGCTATCCTGATGGATAAGATGGCGGGATACG", "Aminoglucósidos", "Enzima modificadora de aminoglucósidos")
        ]

        for src, name, seq, abx, desc in genes:
            ResistanceGene.objects.get_or_create(
                source=src, gene_name=name, sequence=seq,
                antibiotic_class=abx, description=desc
            )
        self.stdout.write(self.style.SUCCESS("✅ Datos CARD y ResFinder cargados correctamente."))
