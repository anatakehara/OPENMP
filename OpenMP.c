//ANA BEATRIZ GOMES TAKEHARA E LAIZA LIMA PERIA
#include <stdio.h>
#include <math.h>
#include <omp.h>

// Função a ser integrada
double func(double x, double y) 
{
    return sin(x * x + y * y);
}

// Método dos trapézios paralelo usando OpenMP
double trapezio_integral(int nx, int ny, int num_threads) {
    double x_min = 0.0, x_max = 1.5;
    double y_min = 0.0, y_max = 1.5;
    double dx = (x_max - x_min) / nx;
    double dy = (y_max - y_min) / ny;
    double sum = 0.0;

    // Paralelizar usando OpenMP
    #pragma omp parallel num_threads(num_threads)
    {
        int tid = omp_get_thread_num();
        int nthreads = omp_get_num_threads();
        
        // Sincroniza todas as threads 
        #pragma omp barrier

        #pragma omp master
        {
            // Região para imprimir o número total de threads
            printf("\nNúmero de threads = %d\n", nthreads);
        }

        // Redução paralela no cálculo da integral
        #pragma omp for reduction(+:sum)
        for (int i = 0; i < nx; i++) {
            double x = x_min + i * dx;
            for (int j = 0; j < ny; j++) {
                double y = y_min + j * dy;
                sum += func(x, y) * dx * dy;
            }
        }
    }

    return sum;
}

int main(int argc, char *argv[]) {
    // Definição de discretização e número de threads
    int nx_values[] = {1000, 10000, 100000};  // Valores de discretização em x
    int ny_values[] = {1000, 10000, 100000};  // Valores de discretização em y
    int threads_values[] = {1, 2, 4, 8};      // Quantidades de threads

    // For para testar todas as combinações de threads e discretizações
    for (int t = 0; t < 4; t++) {
        int num_threads = threads_values[t];

        // Para cada combinação de discretização de x e y
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                int nx = nx_values[i];
                int ny = ny_values[j];

                // Começa a contar o tempo
                double start_time = omp_get_wtime();

                // Chama o cálculo da integral
                double result = trapezio_integral(nx, ny, num_threads);

                // Termina a contagem do tempo
                double end_time = omp_get_wtime();
                double exec_time = end_time - start_time;

                // Imprime os resultados
                printf("nx = %d, ny = %d, Integral = %.6f, Tempo de execução = %f segundos\n", nx, ny, result, exec_time);
            }
        }
    }

    return 0;
}
