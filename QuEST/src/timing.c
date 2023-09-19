#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "timing.h"

#define WARMUP_TIME 5e-3 // 5 ms warmup time by default

// Every run must last for at least 10 ms 
// True running time about 10-100 times this
#define MIN_TIME 0.01 
#define REPLICATION 1

#define DISABLE_CPUTIME

void timeit(void (*pre)(), void (*method)(), void (*post)(), TimingResult *rst) {
    int repetition = 1, replictaion=REPLICATION;
    double tic, toc, temp;
    double wall_tarray[replictaion];
    double cpu_tarray[replictaion];
    
    pre();

    tic = omp_get_wtime();
    while ((toc = omp_get_wtime()) - tic < WARMUP_TIME)
        method();

    // Decide repetition
    while (toc - tic < MIN_TIME) {
        repetition *= 10;
        for (int i = 0; i < repetition; i++)
            method();
        
        temp = omp_get_wtime();
        tic = toc;
        toc = temp;
    }
    repetition = 1 + (double)(repetition) * 0.5 / REPLICATION / (toc - tic);

    // Measure wall time
    for (int rep = 0; rep < replictaion; rep++)
    {
        tic = omp_get_wtime();
        for (int i = 0; i < repetition; i++)
            method();
        toc = omp_get_wtime();
        wall_tarray[rep] = toc - tic;
    }

    #ifndef DISABLE_CPUTIME
    // Measure CPU time
    for (int rep = 0; rep < replictaion; rep++)
    {
        tic = (double)(clock()) / CLOCKS_PER_SEC;
        for (int i = 0; i < repetition; i++) 
            method();
        toc = (double)(clock()) / CLOCKS_PER_SEC;
        cpu_tarray[rep] = toc - tic;
    }
    #endif

    post();


    // Calculate all statistics
    _write_stat(wall_tarray, replictaion, rst->wall_stat);
    _write_stat(cpu_tarray, replictaion, rst->cpu_stat);

    rst->repetition = repetition;
    rst->replication = replictaion;     
}


int _cmp_function (const void *a, const void *b) {
   return ( *(double *)a > *(double *)b );
}


void _write_stat(double *tarray, int n, double *stat) {
    double mean=0, std=0, min, median;
    int median_counter = 0;
    for (int i = 0; i < n; i++)
        mean += tarray[i];
    
    mean /= n;

    for (int i = 0; i < n; i++)
        std += (tarray[i] - mean) * (tarray[i] - mean);
    
    std /= n;
    std = sqrt(std);

    qsort(tarray, n, sizeof(double), _cmp_function);
    min = tarray[0];
    median = tarray[n/2];

    stat[0] = mean;
    stat[1] = std;
    stat[2] = min;
    stat[3] = median;
}


void _print_time(double t) {
    char *unit;
    if (t < 1e-6) {
        unit = "ns";
        t *= 1e9;
    }
    else if (t < 1e-3) {
        unit = "μs";
        t *= 1e6;
    }
    else if (t < 1) {
        unit = "ms";
        t *= 1e3;
    }
    else {
        unit = "s";
    }

    printf("%-.4g %s", t, unit);

}


void print_timing_result(TimingResult *rst) {
    int r = rst->repetition;
    printf("Timing result with %d repetitions and %d runs\n", r, rst->replication);

    printf("[Wall time] ");
    _print_time(rst->wall_stat[0] / r);
    printf(" ± ");
    _print_time(rst->wall_stat[1] / r);
    printf("; min "); // min
    _print_time(rst->wall_stat[2] / r);
    printf("; median ");
    _print_time(rst->wall_stat[3] / r);

    printf("\n[CPU time]  ");
    _print_time(rst->cpu_stat[0] / r); // mean
    printf(" ± ");
    _print_time(rst->cpu_stat[1] / r); // std
    printf("; min "); // min
    _print_time(rst->cpu_stat[2] / r);
    printf("; median ");
    _print_time(rst->cpu_stat[3] / r);
    printf("\n");

}
