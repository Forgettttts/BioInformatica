from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data.CodonTable import unambiguous_dna_by_id, unambiguous_rna_by_id
import random


def detect_sequence_type(seq):
    seq = seq.upper()
    if "U" in seq and "T" in seq:
        raise ValueError("Secuencia inválida: contiene T y U.")
    return "RNA" if "U" in seq else "DNA"


def get_table(seq_type, table_id=1):
    return unambiguous_rna_by_id[table_id] if seq_type == "RNA" else unambiguous_dna_by_id[table_id]


def altcodons(codon, table):
    if codon in table.stop_codons:
        return table.stop_codons
    try:
        aa = table.forward_table[codon]
    except KeyError:
        return []
    return [k for (k, v) in table.forward_table.items() if v == aa]


def randomize_codons(seq, table_id=1, seed=None):
    """Reemplaza cada codón por otro sinónimo elegido al azar.
       Si se pasa una seed, el resultado será reproducible."""
    if seed is not None:
        random.seed(seed)

    seq = seq.upper()
    if len(seq) % 3 != 0:
        return seq

    seq_type = detect_sequence_type(seq)
    table = get_table(seq_type, table_id)
    stop_codons = set(table.stop_codons)
    codones = [seq[i:i+3] for i in range(0, len(seq), 3)]
    nueva_seq = []

    for i, codon in enumerate(codones):
        if i == 0 and codon in ["ATG", "AUG"]:
            nueva_seq.append(codon)
            continue
        if codon in stop_codons or i == len(codones) - 1:
            nueva_seq.append(codon)
            continue

        alternativas = altcodons(codon, table)
        if alternativas:
            nueva_seq.append(random.choice(alternativas))
        else:
            nueva_seq.append(codon)

    return "".join(nueva_seq)


def gc_content(codon):
    return sum(1 for b in codon if b in ["G", "C"]) / 3


def maximize_gc_sequence(seq, table_id=1):
    """Elige siempre el codón sinónimo con mayor contenido GC."""
    seq = seq.upper()
    if len(seq) % 3 != 0:
        return seq

    seq_type = detect_sequence_type(seq)
    table = get_table(seq_type, table_id)
    stop_codons = set(table.stop_codons)
    codones = [seq[i:i+3] for i in range(0, len(seq), 3)]
    nueva_seq = []

    for i, codon in enumerate(codones):
        if i == 0 and codon in ["ATG", "AUG"]:
            nueva_seq.append(codon)
            continue
        if codon in stop_codons or i == len(codones) - 1:
            nueva_seq.append(codon)
            continue

        alternativas = altcodons(codon, table)
        if alternativas:
            mejor = max(alternativas, key=gc_content)
            nueva_seq.append(mejor)
        else:
            nueva_seq.append(codon)

    return "".join(nueva_seq)


def process_fna(input_file, output_random, output_maxgc, seed=None):
    """
    Lee todas las secuencias codificadoras de un archivo .fna
    y genera dos archivos nuevos:
      - cds_random.fna : con codones aleatorizados
      - cds_maxGC.fna  : con codones maximizados en GC
    Si se pasa una seed, los resultados aleatorios serán reproducibles.
    """
    random_records = []
    maxgc_records = []

    for record in SeqIO.parse(input_file, "fasta"):
        seq_str = str(record.seq).upper()

        random_seq = randomize_codons(seq_str, seed=seed)
        maxgc_seq = maximize_gc_sequence(seq_str)

        random_record = SeqRecord(Seq(random_seq), id=record.id, description="randomized codons")
        maxgc_record = SeqRecord(Seq(maxgc_seq), id=record.id, description="maxGC codons")

        random_records.append(random_record)
        maxgc_records.append(maxgc_record)

    SeqIO.write(random_records, output_random, "fasta")
    SeqIO.write(maxgc_records, output_maxgc, "fasta")

    print(f"Resultados guardados en:\n  - {output_random}\n  - {output_maxgc}")


if __name__ == "__main__":
    input_file = "cds.fna"
    output_random = "cds_random.fna"
    output_maxgc = "cds_maxGC.fna"
    seed = None

    process_fna(input_file, output_random, output_maxgc, seed=seed)
