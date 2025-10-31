import matplotlib.pyplot as plt
import subprocess
from Bio import SeqIO, AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

archivo_fasta = "sarscov2.fasta"

# Parseamos el archivo
secuencias = list(SeqIO.parse(archivo_fasta, "fasta"))

archivo_fasta = "sarscov2.fasta"
archivo_aln = "proteinas.aln"  # Nombre del archivo de salida
print(f"Cargando {archivo_fasta}...")
secuencias = list(SeqIO.parse(archivo_fasta, "fasta"))
print(f"Se cargaron {len(secuencias)} secuencias.")

comando = [
    "clustalw2",  # O "clustalw" si así lo tienes en tu PATH
    f"-INFILE={archivo_fasta}",
    f"-OUTFILE={archivo_aln}",
]

subprocess.run(comando, check=True, capture_output=True, text=True)
print(f"Alineamiento completado. Archivo guardado como '{archivo_aln}'.")

# --- Paso 3: Leer Alineamiento y Calcular Distancia ---
print("Calculando matriz de distancia...")
alineamiento = AlignIO.read(archivo_aln, "clustal")

# Usar el modelo 'blosum62' para aminoácidos
calculator = DistanceCalculator("blosum62")
dist_matrix = calculator.get_distance(alineamiento)

# --- Paso 4: Construir el Árbol (Neighbor-Joining) ---
print("Construyendo árbol...")
constructor = DistanceTreeConstructor(calculator, "nj")
tree = constructor.build_tree(alineamiento)

# Define el nombre exacto de tu outgroup (tal como aparece en el FASTA)
outgroup_name = "Wuhan Hu-1"

print(f"Enraizando el árbol con el outgroup: '{outgroup_name}'")

# 1. Buscar el nodo (hoja) del outgroup
outgroup_node = None
for leaf in tree.get_terminals():  # get_terminals() nos da todas las "hojas"
    if leaf.name == outgroup_name:
        outgroup_node = leaf
        break

# 2. Enraizar el árbol usando ese nodo
if outgroup_node:
    tree.root_with_outgroup(outgroup_node)
    print("Árbol enraizado exitosamente.")
else:
    print(f"ADVERTENCIA: No se encontró el outgroup '{outgroup_name}'.")

# Define la función de etiquetado
# n.is_terminal() devuelve True si el nodo es una hoja (hoja = secuencia)
# Si es una hoja, devuelve n.name, si no (es interno), devuelve ""
label_function = lambda n: n.name if n.is_terminal() else ""

# Opción B: Gráfico
print("Generando gráfico del árbol (limpio)...")
fig = plt.figure(figsize=(15, 10), dpi=100)
axes = fig.add_subplot(1, 1, 1)

Phylo.draw(tree, axes=axes, do_show=False, label_func=label_function)

plt.title(f"Árbol Filogenético (Enraizado con {outgroup_name})")
plt.savefig("ArbolP3-2.png")
plt.show()
