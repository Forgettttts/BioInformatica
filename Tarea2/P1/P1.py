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
        while l >= 0 and r < n:
            # print(f"\tizquierda: {l}, derecha: {r}")

            # Comprobar si las bases son complementarias
            base_izquierda = seq[l]
            base_derecha = seq[r]

            if (COMPLEMENTO[base_izquierda] == base_derecha) or (
                base_derecha == base_izquierda
            ):
                # Si lo son, hemos encontrado un palíndromo
                longitud = r - l + 1
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

                if len(top_10) < 10:
                    heapq.heappush(top_10, (longitud, seq[l : r + 1]))
                elif longitud > top_10[0][0]:
                    heapq.heappushpop(top_10, (longitud, seq[l : r + 1]))

                l -= 1
                r += 1
            else:
                break
    resultados_ordenados = sorted(top_10, key=lambda x: x[0], reverse=True)
    top = [palindromo for longitud, palindromo in resultados_ordenados]

    print("\n TOP 10 PALÍNDROMOS MÁS LARGOS ENCONTRADOS:\n")

    for i, palindromo in enumerate(top, 1):
        print(f"{i:^8} {palindromo}")


# Recorreremos la salida de la funcion encontrar_10_palindromos_adn con la secuencia real
encontrar_10_palindromos_adn(
    open("secuenciaCompleta.txt", "r", encoding="utf-8").read()
)
