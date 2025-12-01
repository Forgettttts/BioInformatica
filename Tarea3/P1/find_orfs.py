import warnings
from Bio import BiopythonWarning
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

warnings.filterwarnings("ignore", category=BiopythonWarning)


def find_orfs_with_trans(seq, trans_table=11, min_orf_len_nt=90):
    result = []
    seq_len = len(seq)

    for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
        for frame in range(3):
            trans = nuc[frame:].translate(trans_table, to_stop=False)
            trans_len = len(trans)

            aa_start = 0
            while aa_start < trans_len:
                aa_end = trans.find("*", aa_start)
                if aa_end == -1:
                    aa_end = trans_len

                orf_len_nt = (aa_end - aa_start) * 3

                if orf_len_nt >= min_orf_len_nt:
                    if strand == 1:
                        start = frame + aa_start * 3
                        end = min(seq_len, frame + aa_end * 3)
                    else:
                        start = seq_len - (frame + aa_end * 3)
                        end = seq_len - (frame + aa_start * 3)

                    if 0 <= start < end <= seq_len:
                        orf_seq = seq[start:end]
                        result.append({
                            "start": start,
                            "end": end,
                            "strand": strand,
                            "length": len(orf_seq),
                            "gc": 100 * gc_fraction(orf_seq),
                            "seq": str(orf_seq)
                        })

                aa_start = aa_end + 1

    return result


def top_and_bottom_gc_orfs(fasta_file, min_orf_len_nt=90):
    record = SeqIO.read(fasta_file, "fasta")
    orfs = find_orfs_with_trans(record.seq, min_orf_len_nt=min_orf_len_nt)

    sorted_orfs = sorted(orfs, key=lambda x: x["gc"])

    bottom3 = sorted_orfs[:3]
    top3 = sorted_orfs[-3:]
    return bottom3, top3


# ========================
# Ejecución
# ========================
bottom3, top3 = top_and_bottom_gc_orfs("plasmidoEColi.fna")

print("\n=== Top 3 GC% ===")
for o in top3:
    print(f"\nGC: {o['gc']:.2f}% | Len: {o['length']} | Pos: {o['start']}-{o['end']} | Strand: {o['strand']}")
    print("Secuencia:")
    print(o["seq"])

print("\n=== Bottom 3 GC% ===")
for o in bottom3:
    print(f"\nGC: {o['gc']:.2f}% | Len: {o['length']} | Pos: {o['start']}-{o['end']} | Strand: {o['strand']}")
    print("Secuencia:")
    print(o["seq"])
