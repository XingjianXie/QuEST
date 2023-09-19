#ifndef TIMING_H
#define TIMING_H

#define NUM_STAT 4 // mean, std, min, median

typedef struct TimingResult {
    int repetition, replication;
    double wall_stat[NUM_STAT], cpu_stat[NUM_STAT];
} TimingResult;

void timeit(void (*pre)(), void (*method)(), void (*post)(), TimingResult *rst);

void print_timing_result(TimingResult *rst);

void _write_stat(double *tarray, int n, double *stat);

#endif
