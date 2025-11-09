# 📓Tarea 2
**Fecha**: 10 de Noviembre del 2025. 
**Integrantes**:
- Alan Zapata Silva (*ROL: 201956567-2*).
- Martin Pino Cornejo (*ROL: -*).

## 📝 Requerimientos:
- Python 3.13 o superior
- Biopython
- ClustalW2 instalado y accesible desde la terminal (PATH)
- Matplotlib
- heapq

## 📐 Instrucciones de ejecución:
Cada archivo debe ejecutarse según la pregunta que corresponda:

- Pregunta 1: `P1.py`
> En esta pregunta se requiere que tanto el archivo `P1.py` como la secuencia utilizada (`secuenciaCompleta.txt`) se encuentren en el mismo directorio, por lo que en caso de mover alguno de los dos archivos, se debe mover al otro en conjunto.
- Pregunta 3:
  - 3.a) `P3.1.py`
  - 3.b) `P3.2.py`
- Pregunta 5: `process_cds_fna.py`
> En esta pregunta se requiere que tanto el archivo `process_cds_fna.py` como la secuencia utilizada (`cds.fna`) se encuentren en el mismo directorio, por lo que en caso de mover alguno de los dos archivos, se debe mover al otro en conjunto.
> 
## 📤 Salidas:
- Pregunta 1: Se imprimen en consola los 10 palíndromos más largos encontrados en la secuencia.
- Pregunta 3:
  - 3.a) Se genera un archivo de imagen `ArbolP3-1.png` con el árbol filogenético construido a partir del alineamiento de las proteínas.
  - 3.b) Se genera un archivo de imagen `ArbolP3-2.png` con el árbol filogenético enraizado utilizando la secuencia *"Wuhan Hu-1"* como outgroup.
- Pregunta 5:
  - Se genera un archivo `cds_random.fna` que contiene la secuencia codificadora con **codones aleatorizados**, donde cada codón fue reemplazado por un sinónimo escogido al azar (manteniendo el inicio y el stop).
  - Se genera un archivo `cds_maxGC.fna` que contiene la secuencia codificadora con **codones seleccionados para maximizar el contenido de G y C**.
  - Ambos archivos se guardan en formato **FASTA**.
  - En la consola se imprimen los porcentajes de contenido GC de cada versión:
    ```
    [Secuencia]
    GC% original : XX.XX%
    GC% random   : XX.XX%
    GC% maxGC    : XX.XX%
    ```
  - Los nombres de salida y la semilla aleatoria (`seed`) pueden configurarse al inicio del script `process_cds_fna.py`.
