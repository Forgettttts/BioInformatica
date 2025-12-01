import random
from collections import Counter

def count_nucleotides(seq):
    """Cuenta A, C, G, T y retorna un diccionario."""
    return Counter(seq)

def randomize_fasta(input_file, output_file, seed=None):
    if seed is not None:
        random.seed(seed)
        print(f"Seed utilizada: {seed}")

    with open(input_file, "r") as f:
        header = f.readline().rstrip("\n")
        seq = f.read().replace("\n", "").upper()

    original_count = count_nucleotides(seq)
    total_original = len(seq)

    seq_list = list(seq)
    random.shuffle(seq_list)
    randomized_seq = "".join(seq_list)

    random_count = count_nucleotides(randomized_seq)
    total_random = len(randomized_seq)

    with open(output_file, "w") as out:
        out.write(header + "\n")
        for i in range(0, len(randomized_seq), 70):
            out.write(randomized_seq[i:i+70] + "\n")

    print("\n=== Archivo original ===")
    print(f"Total nt: {total_original}")
    print(f"Conteo: {dict(original_count)}")

    print("\n=== Archivo aleatorizado ===")
    print(f"Total nt: {total_random}")
    print(f"Conteo: {dict(random_count)}")

    print("\nArchivo aleatorizado generado:", output_file)

randomize_fasta("plasmidoEColi.fna", "plasmidoEColi_random.fna", seed=1764564110)