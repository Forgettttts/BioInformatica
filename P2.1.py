from Bio import SeqIO

# Nombre de tu archivo de entrada
archivo_fasta = "sarscov2.fasta"

# Parseamos el archivo
secuencias = list(SeqIO.parse(archivo_fasta, "fasta"))

print(f"Se cargaron {len(secuencias)} secuencias.")

import matplotlib.pyplot as plt
import subprocess  # <-- IMPORTANTE: Añadido
from Bio import SeqIO, AlignIO, Phylo

# Eliminamos la línea "from Bio.Align.Applications import ClustalwCommandline"
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

# --- Paso 1: Cargar FASTA ---
archivo_fasta = "sarscov2.fasta"
archivo_aln = "proteinas.aln"  # Nombre del archivo de salida
print(f"Cargando {archivo_fasta}...")
secuencias = list(SeqIO.parse(archivo_fasta, "fasta"))
print(f"Se cargaron {len(secuencias)} secuencias.")


# --- Paso 2: Correr Alineamiento (El Método Moderno) ---
print("Ejecutando ClustalW con 'subprocess'...")

# Este es el comando que escribirías en la terminal
# ClustalW usa el formato -INFILE= y -OUTFILE=
comando = [
    "clustalw2",  # O "clustalw" si así lo tienes en tu PATH
    f"-INFILE={archivo_fasta}",
    f"-OUTFILE={archivo_aln}",
]

# Ejecutamos el comando
# check=True hará que el script falle si ClustalW da un error
try:
    subprocess.run(comando, check=True, capture_output=True, text=True)
    print(f"Alineamiento completado. Archivo guardado como '{archivo_aln}'.")
except subprocess.CalledProcessError as e:
    print("¡ERROR al ejecutar ClustalW!")
    print(f"Comando: {' '.join(comando)}")
    print("Error (stderr):", e.stderr)
    exit()  # Detenemos el script si ClustalW falla
except FileNotFoundError:
    print("¡ERROR: No se encontró el ejecutable 'clustalw2'!")
    print("Asegúrate de que esté instalado y en tu PATH.")
    exit()


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
tree.root_at_midpoint()  # Enraizar el árbol


# --- Paso 5: Visualizar y Guardar (Sin etiquetas internas) ---

# Define la función de etiquetado
# n.is_terminal() devuelve True si el nodo es una hoja (hoja = secuencia)
# Si es una hoja, devuelve n.name, si no (es interno), devuelve ""
label_function = lambda n: n.name if n.is_terminal() else ""

# Opción A: ASCII
print("--- Árbol en formato ASCII ---")
Phylo.draw_ascii(tree)
print("------------------------------")

# Guardar el árbol
Phylo.write(tree, "mi_arbol_proteinas.nwk", "newick")
print("Árbol guardado en 'mi_arbol_proteinas.nwk'")

# Opción B: Gráfico
print("Generando gráfico del árbol (limpio)...")
fig = plt.figure(figsize=(15, 10), dpi=100)
axes = fig.add_subplot(1, 1, 1)

# ¡AQUÍ ESTÁ LA MODIFICACIÓN!
Phylo.draw(tree, axes=axes, do_show=False, label_func=label_function)  # <-- AÑADIDO

plt.title("Árbol Filogenético")  # Título original
plt.savefig("mi_arbol_proteinas.png")  # Nombre de archivo original
plt.show()

print("¡Proceso finalizado!")
