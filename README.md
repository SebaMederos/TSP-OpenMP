# TSP-OpenMP
Trabajo Práctico Final de la asignatura Sistemas Paralelos y Distribuidos - UC 2022

El problema del agente viajero (TSP) es uno de los más estudiados en el campo de la optimización, como también en la matemática computacional. El mismo consiste en el recorrido de lugares o nodos (para entregar o recoger mercancías) con el fin de resolver problemas que impidan minimizar o maximizar algún objetivo (tiempo o costos). En el presente trabajo se presenta una implementación de la heurística de SA (Simulating Annealing) paralelizada para determinar y analizar una solución con mayor rapidez y eficiencia.

Pasos para ejecturar el proyecto
1. Abrir una ventana de terminal en la carpeta principal del
proyecto.
2. Compilar el proyecto con el comando:
    ```
    gcc -o SA SA.c -fopenmp -l
2.1 En algunos ambientes es posible que se deba ejecutar con las siguientes banderas:
    ```
    -std=c99 -D_POSIX_C_SOURCE=199309L
3. Ejecutar el proyecto con el comando:
    ```
    ./SA <archivo.tsp>

