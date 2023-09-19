/** @file 
 * A demo of QuEST
 *
 * @author Ania Brown
 * @author Tyson Jones
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include "QuEST.h"
#include "timing.h"

#include <omp.h>

QuESTEnv env;
Qureg qubits;
int k;
int num_qubits;

int b1, b2;

int secret;

int q;
int *s;

int numElems;
int numReps;
int solElem;

void svRandomize(Qureg qubits) {
#ifndef COMPRESS
    double sum = 0;
    ComplexArray sv = qubits.stateVec;
    for (size_t i = 0; i < qubits.numAmpsPerChunk; i++) 
    {
        sv.real[i] = (double)rand() / RAND_MAX;
        sv.imag[i] = (double)rand() / RAND_MAX;

        sum += (sv.real[i] * sv.real[i]) + (sv.imag[i] * sv.imag[i]);
    }
    
    for (size_t i = 0; i < qubits.numAmpsPerChunk; i++) 
    {
        sv.real[i] /= sum;
        sv.imag[i] /= sum;
    }
#else
    double sum = 0;
    for (size_t i = 0 ; i < qubits.compressBlockCount; i++) {
        double *content1 = malloc(qubits.compressBlockSize * sizeof(double));
        double *content2 = malloc(qubits.compressBlockSize * sizeof(double));
        for (size_t j = 0; j < qubits.compressBlockSize; j++) {
            content1[j] = (double)rand() / RAND_MAX;
            content2[j] = (double)rand() / RAND_MAX;
            sum += (content1[j] * content1[j]) + (content2[j] * content2[j]);
        }
        compressBlock(&qubits, i, content1, content2);
        free(content1);
        free(content2);
    }
    for (size_t i = 0 ; i < qubits.compressBlockCount; i++) {
        double *content1 = malloc(qubits.compressBlockSize * sizeof(double));
        double *content2 = malloc(qubits.compressBlockSize * sizeof(double));
        decompressBlock(&qubits, i, content1, content2);

        for (size_t j = 0; j < qubits.compressBlockSize; j++) {
            content1[j] /= sum;
            content2[j] /= sum;
        }
        compressBlock(&qubits, i, content1, content2);
        free(content1);
        free(content2);
    }
#endif
}

void pre_BV() {
    qubits = createQureg(k, env);
    secret = rand() % (int) pow(2, k - 1);
    initZeroState(qubits);
}

void applyOracle_BV(Qureg qureg, int numQubits, int secret) {

    int bits = secret;

    for (int q=1; q<numQubits; q++) {

        // extract the (q-1)-th bit of secret
        int bit = bits % 2;
        bits /= 2;
        
        // NOT the ancilla, controlling on the q-th qubit
        if (bit)
            controlledNot(qureg, q, 0);
    }
}

void pre_GS() {

    numElems = (int) pow(2, k);
    //numReps = ceil(M_PI/4 * sqrt(numElems));
    numReps = 1;
    solElem = rand() % numElems;
    
    // prepare |+>
    qubits = createQureg(k, env);
    initPlusState(qubits);
}

/* effect |solElem> -> -|solElem> via a 
 * multi-controlled phase flip gate 
 */
void applyOracle_GS(Qureg qureg, int numQubits, int solElem) {
    
    // apply X to transform |111> into |solElem>
    for (int q=0; q<numQubits; q++)
        if (((solElem >> q) & 1) == 0)
            pauliX(qureg, q);
        
    // effect |111> -> -|111>    
    int ctrls[numQubits];
    for (int q=0; q<numQubits; q++)
        ctrls[q] = q;
    multiControlledPhaseFlip(qureg, ctrls, numQubits);
    
    // apply X to transform |solElem> into |111>
    for (int q=0; q<numQubits; q++)
        if (((solElem >> q) & 1) == 0)
            pauliX(qureg, q);
}

/* apply 2|+><+|-I by transforming into the Hadamard basis 
 * and effecting 2|0><0|-I. We do this, by observing that 
 *   c..cZ = diag{1,..,1,-1} 
 *         = I - 2|1..1><1..1|
 * and hence 
 *   X..X c..cZ X..X = I - 2|0..0><0..0|
 * which differs from the desired 2|0><0|-I state only by 
 * the irrelevant global phase pi
 */
void applyDiffuser_GS(Qureg qureg, int numQubits) {
    
    // apply H to transform |+> into |0>
    for (int q=0; q<numQubits; q++)
        hadamard(qureg, q);

    // apply X to transform |11..1> into |00..0>
    for (int q=0; q<numQubits; q++)
        pauliX(qureg, q);
    
    // effect |11..1> -> -|11..1>
    int ctrls[numQubits];
    for (int q=0; q<numQubits; q++)
        ctrls[q] = q;
    multiControlledPhaseFlip(qureg, ctrls, numQubits);
    
    // apply X to transform |00..0> into |11..1>
    for (int q=0; q<numQubits; q++)
        pauliX(qureg, q);
    
    // apply H to transform |0> into |+>
    for (int q=0; q<numQubits; q++)
        hadamard(qureg, q);
}

void pre() {
    /*
     * PREPARE QUBIT SYSTEM
     */

    qubits = createQureg(num_qubits, env);
    svRandomize(qubits);
}

typedef struct {
    unsigned long size,resident,share,text,lib,data,dt;
} statm_t;

void read_off_memory_status(statm_t* result)
{
  unsigned long dummy;
  const char* statm_path = "/proc/self/statm";

  FILE *f = fopen(statm_path,"r");
  if(!f){
    perror(statm_path);
    abort();
  }
  if(7 != fscanf(f,"%ld %ld %ld %ld %ld %ld %ld",
    &result->size,&result->resident,&result->share,&result->text,&result->lib,&result->data,&result->dt))
  {
    perror(statm_path);
    abort();
  }
  fclose(f);
}

void post() {
    /*
     * FREE MEMORY
     */
    statm_t mem;
    //read_off_memory_status(&mem);
    //fprintf(stderr, "[Memory usage: %ld %ld %ld %ld %ld %ld %ld]\n", mem.size, mem.resident, mem.share, mem.text, mem.lib, mem.data, mem.dt);
    //system("ps ux | grep demo >> /dev/stderr");
    //fprintf(stderr, "[Done]\n\n");
    destroyQureg(qubits, env); 
}

void timing_unitary() {
    ComplexMatrix2 u = {
        .real={{0.5, 0.5}, {0.5, 0.5}},
        .imag={{0.5, -0.5}, {-0.5, 0.5}}
    };
    unitary(qubits, k, u);
}

void timing_compactUnitary() {
    compactUnitary(qubits, k, (Complex) {.real=0.5, .imag=0.5}, (Complex) {.real=0.5, .imag=-0.5});
}

void timing_pauliX() {
    pauliX(qubits, k);
}

void timing_pauliY() {
    pauliY(qubits, k);
}

void timing_pauliZ() {
    pauliZ(qubits, k);
}

void timing_phaseShift() {
    phaseShift(qubits, k, 0.5);
}

void timing_sGate() {
    sGate(qubits, k);
}

void timing_tGate() {
    tGate(qubits, k);
}

void timing_controlledUnitary() {
    ComplexMatrix2 u = {
        .real={{0.5, 0.5}, {0.5, 0.5}},
        .imag={{0.5, -0.5}, {-0.5, 0.5}}
    };
    controlledUnitary(qubits, b1, b2, u);
}

void timing_controlledNot() {
    controlledNot(qubits, b1, b2);
}

void timing_controlledPhaseShift() {
    controlledPhaseShift(qubits, b1, b2, 0.5);
}

void timing_controlledPhaseFlip() {
    controlledPhaseFlip(qubits, b1, b2);
}

void timing_multiControlledPhaseFlip() {
    multiControlledPhaseFlip(qubits, s, q);
}

void timing_multiControlledPhaseShift() {
    multiControlledPhaseShift(qubits, s, q, 0.5);
}

void timing_BV() {
    // NOT the ancilla
    pauliX(qubits, 0);

    // H all qubits, including the ancilla
    for (int q=0; q<k; q++)
        hadamard(qubits, q);

    applyOracle_BV(qubits, k, secret);

    for (int q=0; q<k; q++)
        hadamard(qubits, q);
}

void timing_GS() {
    for (int r=0; r<numReps; r++) {
        applyOracle_GS(qubits, k, solElem);
        applyDiffuser_GS(qubits, k);
    }
}

int main (int narg, char *varg[]) {

    if (narg != 3) {
        printf("usage: ./tutorial_example <num_threads> <num_qubits>\n");
        return 1;
    }

    int num_threads = atoi(varg[1]);
    omp_set_num_threads(num_threads);

    num_qubits = atoi(varg[2]);

    srand(time(NULL));

    /*
     * PREPARE QuEST environment
     * (Required only once per program)
     */

    env = createQuESTEnv();

    printf("\n=======================================================\n");

    reportQuESTEnv(env);
    char hostname[1024];
    gethostname(hostname, 1024);
    printf("Host: %s\n", hostname);

    printf("-------------------------------------------------------\n");
    printf("Running QuEST banchmark with %d qubits.\n", num_qubits);
#if (defined NEWAPI || defined COMPRESS)
    printf("\t [Prompt]\n");
#else
    printf("\t [QuEST]\n");
#endif
    printf("-------------------------------------------------------\n");

    /*
     * TIMING
     */

    TimingResult rst;

    printf("\nBV:\n");
    for (k = 1; k <= num_qubits; k++)
    {
        timeit(pre_BV, timing_BV, post, &rst);
        printf("%lf, ", rst.wall_stat[2] / rst.repetition * 1e6);
        fflush(stdout);
    }

    printf("\nGS:\n");
    for (k = 1; k <= num_qubits; k++)
    {
        timeit(pre_GS, timing_GS, post, &rst);
        printf("%lf, ", rst.wall_stat[2] / rst.repetition * 1e6);
        fflush(stdout);
    }

    // printf("\nmultiControlledPhaseFlip:\n");
    // for (b1 = 0; b1 < num_qubits; b1++)
    // {
    //     for (b2 = 0; b2 < b1; b2++)
    //         printf("N/A, ");
    //     for (b2 = b1; b2 < num_qubits; b2++)
    //     {
    //         q = b2 - b1 + 1;
    //         s = malloc(q * sizeof(int));
    //         for (int i = b1; i <= b2; i++) {
    //             s[i - b1] = i;
    //         }
    //         timeit(pre, timing_multiControlledPhaseFlip, post, &rst);
    //         printf("%f, ", rst.wall_stat[2] / rst.repetition * 1e6);
    //         free(s);
    //     }
    //     printf("\n");
    // }

    printf("\nmultiControlledPhaseShift:\n");
    for (b1 = 0; b1 < num_qubits; b1++)
    {
        for (b2 = 0; b2 < b1; b2++)
            printf("N/A, ");
        for (b2 = b1; b2 < num_qubits; b2++)
        {
            q = b2 - b1 + 1;
            s = malloc(q * sizeof(int));
            for (int i = b1; i <= b2; i++) {
                s[i - b1] = i;
            }
            timeit(pre, timing_multiControlledPhaseShift, post, &rst);
            printf("%f, ", rst.wall_stat[2] / rst.repetition * 1e6);
            fflush(stdout);
            free(s);
        }
        printf("\n");
    }

    printf("\nunitary:\n");
    for (k = 0; k < num_qubits; k++)
    {
        timeit(pre, timing_unitary, post, &rst);
        printf("%f, ", rst.wall_stat[2] / rst.repetition * 1e6);
        fflush(stdout);
    }

    // printf("\ncompactUnitary:\n");
    // for (k = 0; k < num_qubits; k++)
    // {
    //     timeit(pre, timing_compactUnitary, post, &rst);
    //     printf("%f, ", rst.wall_stat[2] / rst.repetition * 1e6);
    // }

    printf("\npauliX:\n");
    for (k = 0; k < num_qubits; k++)
    {
        timeit(pre, timing_pauliX, post, &rst);
        printf("%f, ", rst.wall_stat[2] / rst.repetition * 1e6);
        fflush(stdout);
    }

    // printf("\npauliY:\n");
    // for (k = 0; k < num_qubits; k++)
    // {
    //     timeit(pre, timing_pauliY, post, &rst);
    //     printf("%f, ", rst.wall_stat[2] / rst.repetition * 1e6);
    // }

    // printf("\npauliZ:\n");
    // for (k = 0; k < num_qubits; k++)
    // {
    //     timeit(pre, timing_pauliZ, post, &rst);
    //     printf("%f, ", rst.wall_stat[2] / rst.repetition * 1e6);
    // }

    printf("\nphaseShift:\n");
    for (k = 0; k < num_qubits; k++)
    {
        timeit(pre, timing_phaseShift, post, &rst);
        printf("%f, ", rst.wall_stat[2] / rst.repetition * 1e6);
        fflush(stdout);
    }

    // printf("\nsGate:\n");
    // for (k = 0; k < num_qubits; k++)
    // {
    //     timeit(pre, timing_sGate, post, &rst);
    //     printf("%f, ", rst.wall_stat[2] / rst.repetition * 1e6);
    // }

    // printf("\ntGate:\n");
    // for (k = 0; k < num_qubits; k++)
    // {
    //     timeit(pre, timing_tGate, post, &rst);
    //     printf("%f, ", rst.wall_stat[2] / rst.repetition * 1e6);
    // }

    printf("\ncontrolledUnitary:\n");
    for (b1 = 0; b1 < num_qubits; b1++)
    {
        for (b2 = 0; b2 < num_qubits; b2++)
        {
            if (b1 == b2) {
                printf("N/A, ");
                continue;
            }
            timeit(pre, timing_controlledUnitary, post, &rst);
            printf("%f, ", rst.wall_stat[2] / rst.repetition * 1e6);
            fflush(stdout);
        }
        printf("\n");
    }

    printf("\ncontrolledNot:\n");
    for (b1 = 0; b1 < num_qubits; b1++)
    {
        for (b2 = 0; b2 < num_qubits; b2++)
        {
            if (b1 == b2) {
                printf("N/A, ");
                continue;
            }
            timeit(pre, timing_controlledNot, post, &rst);
            printf("%f, ", rst.wall_stat[2] / rst.repetition * 1e6);
            fflush(stdout);
        }
        printf("\n");
    }

    printf("\ncontrolledPhaseShift:\n");
    for (b1 = 0; b1 < num_qubits; b1++)
    {
        for (b2 = 0; b2 < num_qubits; b2++)
        {
            if (b1 == b2) {
                printf("N/A, ");
                continue;
            }
            timeit(pre, timing_controlledPhaseShift, post, &rst);
            printf("%f, ", rst.wall_stat[2] / rst.repetition * 1e6);
            fflush(stdout);
        }
        printf("\n");
    }

    // printf("\ncontrolledPhaseFlip:\n");
    // for (b1 = 0; b1 < num_qubits; b1++)
    // {
    //     for (b2 = 0; b2 < num_qubits; b2++)
    //     {
    //         if (b1 == b2) {
    //             printf("N/A, ");
    //             continue;
    //         }
    //         timeit(pre, timing_controlledPhaseFlip, post, &rst);
    //         printf("%f, ", rst.wall_stat[2] / rst.repetition * 1e6);
    //     }
    //     printf("\n");
    // }

    finish:

    /*
     * CLOSE QUEST ENVIRONMET
     * (Required once at end of program)
     */
    destroyQuESTEnv(env);

    printf("\n\n");

    return 0;
}
