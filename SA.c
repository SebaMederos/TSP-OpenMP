// gcc -o SA SA.c -fopenmp -lm
// macOS x86_64: /usr/local/opt/llvm/bin/clang -fopenmp -L/usr/local/opt/llvm/lib SA.c -o SA

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <omp.h>

#define NUM_THREADS 8       // numero de hilos
#define N 76 			// cuidades
#define SPEED 0.98 			// velocidad de recorrido
#define INITIAL_TEMP 10000 	// la temperatura inicial
#define MIN_TEMP 0.001	    // temperatura minima
#define L 100000

// 1: P SA TSP
// default: S SA TSP 
#define SA 2

// estructura del archivo de tsplib
typedef struct {
	char name[64];
	char type[64];
	char comment[64];
	int dimension;
	char edge_weight_type[64];
} TSPLIB;

// coordenadas 2D
typedef struct {
	int x;
	int y;
} EUC_2D;

// estructura de la solucion
typedef struct {
    double length;	// longitud de la trayectoria
    int path[N]; 	// camino
} UNIT;

// distancia entre dos puntos del plano 2d
double distance(int x1, int y1, int x2, int y2) {

	return sqrtf((double)(pow(abs(x2 - x1), 2) + pow(abs(y2 - y1), 2)));
}

// imprimir matriz de distancia
void printl(double **length_table) {

    for (unsigned long long int i = 0; i < N; i++) {
        for (unsigned long long int j = 0; j< N; j++) {
            printf("%f, ", length_table[i][j]);
        }

        printf("\n");
    }

    printf("\n\n");
}

// genera aleatoriamente una solucion inicial, 
// selecciona aleatoriamente un punto inicial para construir la ruta
void generate(UNIT *tmp, unsigned int *seed){

    int i = rand_r(seed) % N;

    for (unsigned long long int j = 0; j < N; j++) {
		// path[0] = ciudad inicial

        tmp->path[j] = i;
        i = (i + 1) % N;
    }
}

// acepta la nueva regla de solucion en dos casos
bool accept(UNIT bestone, UNIT temp, double t, unsigned int *seed){

    if (bestone.length > temp.length)
        return true;
    else {
        if ((int)(exp((bestone.length - temp.length) / t) * 100) > (rand_r(seed) % 101))
            return true;
    }

    return false;
}

// calcular la longitud del camino,
// el viaje es redondo
void calculate_length(UNIT *p, double **length_table){

    p->length = 0;

    for (unsigned long long int j = 1; j < N; j++) {
		// n = 1,2,3,4, ...
		// costo de la ciudad n a la ciudad n + 1
       p->length += length_table[p->path[j - 1]][p->path[j]];

       //printf("iteration: %lld\n", j);
       //printf("x: %d ; y %d\n", p->path[j - 1], p->path[j]);
    }

	// costo de ir a la ciudad de partida estando en la ultima ciudad
    p->length += length_table[p->path[N - 1]][p->path[0]];
}

// el metodo de intercambio genera vecinos
void get_new_solution(UNIT *p, unsigned int *seed){

    int i = rand_r(seed) % N;
    int j = rand_r(seed) % N;

    int tmp = p->path[i];

    p->path[i] = p->path[j];
    p->path[j] = tmp;
}

// imprimir matriz de soluciones
void printp(UNIT p){

    printf("la distancia del vendedor viajero es: %f \n", p.length);
    printf("las ciudades que pasan estan en orden:\n");

    for (unsigned long long int i = 0; i < N; i++) {
        printf("%06d ", p.path[i]);
        if ((i + 1) % 10 == 0)
            printf("\n");
    }

    printf("\n");
    printf("\n");
}

// distribuye la solucion de menor costo
void bestone(UNIT *p, UNIT *q) {
    UNIT tmp;
    memcpy(&tmp, &p[0], sizeof(p[0]));

    for (int i = 1; i < NUM_THREADS; i++)
        if (p[i].length < tmp.length)
            memcpy(&tmp, &p[i], sizeof(p[i]));

    for (int i = 0; i < NUM_THREADS; i++) {
        memcpy(&p[i], &tmp, sizeof(tmp));
        memcpy(&q[i], &tmp, sizeof(tmp));
    }
}

void co_operation(UNIT *p, UNIT *q) {
    for (int i = 1; i < NUM_THREADS; i++) {
        if (p[i - 1].length < p[i].length) {
            memcpy(&p[i], &p[i - i], sizeof(p[i - i]));
            memcpy(&q[i], &p[i - i], sizeof(p[i - i]));
        }
    }
}

UNIT S_SA_TST(double **length_table, unsigned int *seed) {
    UNIT initial_solution;
    generate(&initial_solution, seed);
    calculate_length(&initial_solution, length_table);

    UNIT current_solution;
    UNIT best_solution;
    memcpy(&best_solution, &initial_solution, sizeof(initial_solution));
    memcpy(&current_solution, &initial_solution, sizeof(initial_solution));

    unsigned long long int i;
    for (double t = INITIAL_TEMP; t > MIN_TEMP; t *= SPEED) {
        for (i = 0; i < L; i++) {
            get_new_solution(&current_solution, seed);
            calculate_length(&current_solution, length_table);

            if (accept(best_solution, current_solution, t, seed)) {
                memcpy(&best_solution, &current_solution, sizeof(current_solution));
            } else {
                memcpy(&current_solution, &best_solution, sizeof(best_solution));
            }
        }
    }

    return best_solution;
}

UNIT P_SA_TST(double **length_table, unsigned int *seed) {

    UNIT initial_solution;
    generate(&initial_solution, seed);
    calculate_length(&initial_solution, length_table);

    int tid;
    UNIT current_solution[NUM_THREADS];
    UNIT best_solution[NUM_THREADS];

    #pragma omp parallel num_threads(NUM_THREADS) default(none) private(tid) firstprivate(initial_solution) shared(current_solution, best_solution)
    {
        tid = omp_get_thread_num();
        memcpy(&best_solution[tid], &initial_solution, sizeof(initial_solution));
        memcpy(&current_solution[tid], &initial_solution, sizeof(initial_solution));
    }

    unsigned long long int i;
    for (double t = INITIAL_TEMP; t > MIN_TEMP; t *= SPEED) {
        // hotspot
        #pragma omp parallel num_threads(NUM_THREADS) default(none) private(i, tid) firstprivate(length_table, current_solution) shared(t, seed, best_solution)
        {
            tid = omp_get_thread_num();

            #pragma omp for
            for (i = 0; i < L; i++) {
                get_new_solution(&current_solution[tid], seed);
                calculate_length(&current_solution[tid], length_table);

                if (accept(best_solution[tid], current_solution[tid], t, seed)) {
                    memcpy(&best_solution[tid], &current_solution[tid], sizeof(current_solution[tid]));
                } else {
                    memcpy(&current_solution, &best_solution[tid], sizeof(best_solution[tid]));
                }
            }
        }

        co_operation(best_solution, current_solution);
    }

    bestone(best_solution, current_solution);
	return best_solution[0];
}

int main(int argc, char *argv[]) {

	if (argc <= 1) {
		printf("./SA <file.tsp>\n");
		return 1;
	}

	// Timer start
	struct timespec ts_start;
	clock_gettime(CLOCK_MONOTONIC, &ts_start);

    unsigned int seed = time(NULL);

	FILE *file;
	TSPLIB instance;

	file = fopen(argv[1], "r");

	fscanf(file, "NAME : %64[^\n]s", instance.name);							// nombre de la instancia TSP
	fscanf(file, "\nCOMMENT : %64[^\n]s", instance.comment);					// comentario
	fscanf(file, "\nTYPE : %64[^\n]s", instance.type);							// tipo de la instancia
	fscanf(file, "\nDIMENSION : %d", &instance.dimension);						// cantidad total de nodos
	fscanf(file, "\nEDGE_WEIGHT_TYPE : %64[^\n]s", instance.edge_weight_type);	// tipo de coordenada de los pesos
	fscanf(file, "\nNODE_COORD_SECTION");	

	if (strcmp(instance.edge_weight_type, "EUC_2D")) {	// los pesos deben ser "EUC_2D"
		return 1;
	}

	EUC_2D *coord = NULL; // coordenadas de los nodos
	if (!(coord = (EUC_2D *) malloc(instance.dimension * sizeof(EUC_2D)))) {
		return 1;
	}

	for (int i = 0; i < instance.dimension; i++) {
		fscanf(file, "\n %*[^ ] %d %d", &coord[i].x, &coord[i].y);
	}

	fclose(file);

	double **length_table = NULL; // matriz de distancia entre todas la ciudades
	if (!(length_table = (double **) malloc(instance.dimension * sizeof(double *)))) {
		return 1;
	}

	for (int i = 0; i < instance.dimension; i++) {
		length_table[i] = (double*) malloc(instance.dimension * sizeof(double));
	}

	for (int i = 0; i < instance.dimension; i++) {
		for (int j = i + 1; j < instance.dimension; j++) {
			// los costos de ida y vuelta son simetricos
			length_table[i][j] = distance(coord[i].x, coord[i].y, coord[j].x, coord[j].y);
			length_table[j][i] = length_table[i][j];
		}
	}

	free(coord);

    #if SA == 1
	UNIT bestone = P_SA_TST(length_table, &seed);
    printf("solucion P SA final:\n");
    #else
    UNIT bestone = S_SA_TST(length_table, &seed);
    printf("solucion S SA final:\n");
    #endif

	printp(bestone);

	// timer stop
	struct timespec ts_stop;
	clock_gettime(CLOCK_MONOTONIC, &ts_stop);
	double start = (double)ts_start.tv_sec + (double)ts_start.tv_nsec/1000000000.0;
	double stop = (double)ts_stop.tv_sec + (double)ts_stop.tv_nsec/1000000000.0;
	double elapsed = (stop - start);

	// display time
	printf ("Time = %fs\n", elapsed);

	for (int i = 0; i < N; i++) {
		free(length_table[i]);
	}

    free(length_table);

    return 0;
}
