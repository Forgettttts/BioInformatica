import heapq


def encontrar_10_palindromos_adn(seq):
    """
    Encuentra los 10 palíndromos de ADN (pares complementarios inversos)
    más largos en una secuencia de ADN.
    """
    COMPLEMENTO = {"A": "T", "T": "A", "C": "G", "G": "C"}
    n = len(seq)
    top_10 = []
    # Iterar por cada base como un posible "centro"
    for i in range(n):

        # --- CASO 1: Palíndromo de longitud IMPAR
        l = i - 1
        r = i + 1
        print(f"Iteración: {i}, izquierda: {l}, derecha: {r}")
        while l >= 0 and r < n:
            # print(f"\tizquierda: {l}, derecha: {r}")

            # Comprobar si las bases son complementarias
            base_izquierda = seq[l]
            base_derecha = seq[r]
            print(f"\tComparando {base_izquierda} (izq) con {base_derecha} (der)")

            if (COMPLEMENTO[base_izquierda] == base_derecha) or (
                base_derecha == base_izquierda
            ):
                # Si lo son, hemos encontrado un palíndromo
                longitud = r - l + 1
                print(f"\t  Longitud palíndromo: {longitud}")
                # Usamos el heap para guardar los 10 mejores
                if len(top_10) < 10:
                    # Si el heap no está lleno, simplemente añadimos el nuevo
                    heapq.heappush(top_10, (longitud, seq[l : r + 1]))
                elif longitud > top_10[0][0]:
                    # Si está lleno, y este es MÁS LARGO que el más corto del heap...
                    # ...sacamos el más corto e insertamos este.
                    heapq.heappushpop(top_10, (longitud, seq[l : r + 1]))

                # Expandir
                l -= 1
                r += 1
            else:
                # Si no son complementarios, rompemos la expansión
                break

        # --- CASO 2: Palíndromo de longitud PAR
        l = i
        r = i + 1

        while l >= 0 and r < n:
            base_izquierda = seq[l]
            base_derecha = seq[r]

            if (COMPLEMENTO[base_izquierda] == base_derecha) or (
                base_derecha == base_izquierda
            ):
                longitud = r - l + 1

                print(f"Longitud palíndromo: {longitud}")
                if len(top_10) < 10:
                    heapq.heappush(top_10, (longitud, seq[l : r + 1]))
                elif longitud > top_10[0][0]:
                    heapq.heappushpop(top_10, (longitud, seq[l : r + 1]))

                l -= 1
                r += 1
            else:
                break
    resultados_ordenados = sorted(top_10, key=lambda x: x[0], reverse=True)
    return [palindromo for longitud, palindromo in resultados_ordenados]


def leer_secuencia_final(path="secuenciaCompleta.txt"):
    """
    Lee el contenido de secuenciaCompleta.txt y lo devuelve como una cadena. El archivo en cuestión tiene la secuencia de Genome assembly ASM18338v1 (https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000183385.1/)
    """
    contenido = open(path, "r", encoding="utf-8").read()
    return contenido


secuenciaReal = leer_secuencia_final("secuenciaCompleta.txt")
secuencia = "AGCGT"

# Recorreremos la salida de la funcion encontrar_10_palindromos_adn con la secuencia real
listita = encontrar_10_palindromos_adn(secuencia)
print("Palíndromos encontrados:")
for palindromo in listita:
    print(palindromo)
