#include "QuEST.h"
#include "QuEST_internal.h"
#include "QuEST_precision.h"
#include "mt19937ar.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "zstd.h"

#define OVERRIDING
#define BLOCKSIZE (1ul << 26)
#define COMPRESS_PARALLEL 16

void compressBlockZSTD(compressed_t* compressedBlock, size_t blockSize, double* content) {
    if (compressedBlock->data != NULL) {
        free(compressedBlock->data);
        compressedBlock->data = NULL;
        compressedBlock->size = 0;
    }
    size_t compressedSize = ZSTD_compressBound(blockSize * sizeof(double));
    compressedBlock->data = malloc(compressedSize);
    size_t compressedSizeActual = ZSTD_compress(compressedBlock->data, compressedSize, content, blockSize * sizeof(double), 3);
    if (ZSTD_isError(compressedSizeActual)) {
        printf("Error compressing block: %s\n", ZSTD_getErrorName(compressedSizeActual));
        exit(1);
    }
    compressedBlock->data = realloc(compressedBlock->data, compressedSizeActual);
    compressedBlock->size = compressedSizeActual;
}

void decompressBlockZSTD(compressed_t compressedBlock, size_t blockSize, double* buffer) {
    size_t decompressedSize = blockSize * sizeof(double);
    size_t decompressedSizeActual = ZSTD_decompress(buffer, decompressedSize, compressedBlock.data, compressedBlock.size);
    if (ZSTD_isError(decompressedSizeActual)) {
        printf("Error decompressing block: %s\n", ZSTD_getErrorName(decompressedSizeActual));
        exit(1);
    }
}

void compressBlock(Qureg *qureg, size_t blockIndex, double* contentReal, double* contentImag) {
    #pragma omp parallel for shared(contentReal, contentImag)
    for (size_t j = 0; j < qureg->numSubBlocksPerBlock; j++) {
        compressBlockZSTD(&qureg->compressedRealBlocks[blockIndex * qureg->numSubBlocksPerBlock + j], qureg->subBlockSize, &contentReal[j * qureg->subBlockSize]);
        compressBlockZSTD(&qureg->compressedImagBlocks[blockIndex * qureg->numSubBlocksPerBlock + j], qureg->subBlockSize, &contentImag[j * qureg->subBlockSize]);
    }
}

void decompressBlock(Qureg *qureg, size_t blockIndex, double* bufferReal, double* bufferImag) {
    #pragma omp parallel for shared(bufferReal, bufferImag)
    for (size_t j = 0; j < qureg->numSubBlocksPerBlock; j++) {
        decompressBlockZSTD(qureg->compressedRealBlocks[blockIndex * qureg->numSubBlocksPerBlock + j], qureg->subBlockSize, &bufferReal[j * qureg->subBlockSize]);
        decompressBlockZSTD(qureg->compressedImagBlocks[blockIndex * qureg->numSubBlocksPerBlock + j], qureg->subBlockSize, &bufferImag[j * qureg->subBlockSize]);
    }
}

OVERRIDING void statevec_createQureg(Qureg *qureg, int numQubits, QuESTEnv env) {
    if (env.numRanks > 1) {
        printf("ERROR: statevec_createQureg is not implemented for distributed mode\n");
        exit(EXIT_FAILURE);
    } 
    long long int numAmps = 1LL << numQubits;
    long long int numAmpsPerRank = numAmps/env.numRanks;
    qureg->numQubitsInStateVec = numQubits;
    qureg->numAmpsTotal = numAmps;
    qureg->numAmpsPerChunk = numAmpsPerRank;
    qureg->chunkId = env.rank;
    qureg->numChunks = env.numRanks;
    qureg->isDensityMatrix = 0;

    qureg->stateVec.real = NULL;
    qureg->stateVec.imag = NULL;

    // Magic starts

    qureg->numQubits = numQubits;
    qureg->numAmps = 1ul << numQubits;
    qureg->compressBlockSize = BLOCKSIZE > qureg->numAmps ? qureg->numAmps : BLOCKSIZE;
    qureg->compressBlockCount = qureg->numAmps / qureg->compressBlockSize;
    qureg->subBlockSize = qureg->compressBlockSize > COMPRESS_PARALLEL ? qureg->compressBlockSize / COMPRESS_PARALLEL : 1;
    qureg->subBlockCount = qureg->numAmps / qureg->subBlockSize;
    qureg->numSubBlocksPerBlock = qureg->compressBlockSize / qureg->subBlockSize;
    qureg->compressedRealBlocks = (compressed_t*) malloc(qureg->subBlockCount * sizeof(compressed_t));
    qureg->compressedImagBlocks = (compressed_t*) malloc(qureg->subBlockCount * sizeof(compressed_t));
    memset(qureg->compressedRealBlocks, 0, qureg->subBlockCount * sizeof(compressed_t));
    memset(qureg->compressedImagBlocks, 0, qureg->subBlockCount * sizeof(compressed_t));
}

OVERRIDING void statevec_initZeroState(Qureg qureg)
{
    statevec_initBlankState(qureg);
    if (qureg.chunkId == 0){
        double *blankReal = (double*) malloc(qureg.compressBlockSize * sizeof(double));
        double *blankImag = (double*) malloc(qureg.compressBlockSize * sizeof(double));
        memset(blankReal, 0, qureg.compressBlockSize * sizeof(double));
        memset(blankImag, 0, qureg.compressBlockSize * sizeof(double));
        blankReal[0] = 1.0;
        compressBlock(&qureg, 0, blankReal, blankImag);
        free(blankReal);
        free(blankImag);
    }
}

OVERRIDING void statevec_initBlankState(Qureg qureg) {
    double *blank = (double*) malloc(qureg.compressBlockSize * sizeof(double));
    memset(blank, 0, qureg.compressBlockSize * sizeof(double));
    for (size_t i = 0; i < qureg.compressBlockCount; i++) {
        compressBlock(&qureg, i, blank, blank);
    }
    free(blank);
}

OVERRIDING void statevec_initPlusState(Qureg qureg) {
    double *plus = (double*) malloc(qureg.compressBlockSize * sizeof(double));
    double *blank = (double*) malloc(qureg.compressBlockSize * sizeof(double));
    memset(blank, 0, qureg.compressBlockSize * sizeof(double));
    double normFactor = 1.0/sqrt((double)qureg.numAmps);
    for (size_t i = 0; i < qureg.compressBlockSize; i++) {
        plus[i] = normFactor;
    }
    for (size_t i = 0; i < qureg.compressBlockCount; i++) {
        compressBlock(&qureg, i, plus, blank);
    }
    free(plus);
    free(blank);
}

OVERRIDING void densmatr_initPlusState(Qureg qureg) {
    double *plus = (double*) malloc(qureg.compressBlockSize * sizeof(double));
    double *blank = (double*) malloc(qureg.compressBlockSize * sizeof(double));
    memset(blank, 0, qureg.compressBlockSize * sizeof(double));
    double normFactor = 1.0/((double)(1ll << qureg.numQubitsRepresented));
    for (size_t i = 0; i < qureg.compressBlockSize; i++) {
        plus[i] = normFactor;
    }
    for (size_t i = 0; i < qureg.compressBlockCount; i++) {
        compressBlock(&qureg, i, plus, blank);
    }
    free(plus);
    free(blank);
}

OVERRIDING void statevec_initDebugState(Qureg qureg) {
    double *debug1 = (double*) malloc(qureg.compressBlockSize * sizeof(double));
    double *debug2 = (double*) malloc(qureg.compressBlockSize * sizeof(double));
    for (size_t i = 0; i < qureg.compressBlockCount; i++) {
        for (size_t j = 0; j < qureg.compressBlockSize; j++) {
            debug1[j] = ((i * qureg.compressBlockSize + j)*2.0)/10.0;
            debug2[j] = ((i * qureg.compressBlockSize + j)*2.0+1.0)/10.0;
        }
        compressBlock(&qureg, i, debug1, debug2);
    }
    free(debug1);
    free(debug2);
}

OVERRIDING void statevec_unitary(Qureg qureg, int target, ComplexMatrix2 u) {
    size_t Target = 1ul << target;
    double ar = u.real[0][0], br = u.real[0][1], cr = u.real[1][0], dr = u.real[1][1];
    double ai = u.imag[0][0], bi = u.imag[0][1], ci = u.imag[1][0], di = u.imag[1][1];
    size_t tt;
    double x_real, x_imag, y_real, y_imag;
    if (Target < qureg.compressBlockSize) {
        // in same block
        for (size_t i = 0; i < qureg.compressBlockCount; i++) {
            double* real = malloc(qureg.compressBlockSize * sizeof(double));
            double* imag = malloc(qureg.compressBlockSize * sizeof(double));

            decompressBlock(&qureg, i, real, imag);

            #pragma omp parallel for simd collapse(2) private(tt, x_real, x_imag, y_real, y_imag)
            for (size_t t = 0; t < qureg.compressBlockSize; t += 2ul * Target)
            {
                for (tt = 0; tt < Target; tt++)
                {
                    x_real = (ar * real[t+tt] - ai * imag[t+tt]) + (br * real[t+tt+Target] - bi * imag[t+tt+Target]);
                    x_imag = (ar * imag[t+tt] + ai * real[t+tt]) + (br * imag[t+tt+Target] + bi * real[t+tt+Target]);
                    y_real = (cr * real[t+tt] - ci * imag[t+tt]) + (dr * real[t+tt+Target] - di * imag[t+tt+Target]);
                    y_imag = (cr * imag[t+tt] + ci * real[t+tt]) + (dr * imag[t+tt+Target] + di * real[t+tt+Target]);
                    real[t+tt] = x_real;
                    imag[t+tt] = x_imag;
                    real[t+tt+Target] = y_real;
                    imag[t+tt+Target] = y_imag;
                }
            }
            compressBlock(&qureg, i, real, imag);
            free(real);
            free(imag);
        }
    } else {
        // in different blocks
        for (size_t i = 0; i < qureg.compressBlockCount / 2; i++) {
            size_t mask = Target / qureg.compressBlockSize - 1;
            size_t i0 = (i & mask) | ((i & ~mask) << 1);
            size_t i1 = i0 | (mask + 1);
            
            double* real0 = malloc(qureg.compressBlockSize * sizeof(double));
            double* imag0 = malloc(qureg.compressBlockSize * sizeof(double));
            double* real1 = malloc(qureg.compressBlockSize * sizeof(double));
            double* imag1 = malloc(qureg.compressBlockSize * sizeof(double));

            decompressBlock(&qureg, i0, real0, imag0);
            decompressBlock(&qureg, i1, real1, imag1);

            #pragma omp parallel for simd
            for (size_t t = 0; t < qureg.compressBlockSize; t++)
            {
                double x_real = (ar * real0[t] - ai * imag0[t]) + (br * real1[t] - bi * imag1[t]);
                double x_imag = (ar * imag0[t] + ai * real0[t]) + (br * imag1[t] + bi * real1[t]);
                double y_real = (cr * real0[t] - ci * imag0[t]) + (dr * real1[t] - di * imag1[t]);
                double y_imag = (cr * imag0[t] + ci * real0[t]) + (dr * imag1[t] + di * real1[t]);
                real0[t] = x_real;
                imag0[t] = x_imag;
                real1[t] = y_real;
                imag1[t] = y_imag;
            }
            compressBlock(&qureg, i0, real0, imag0);
            compressBlock(&qureg, i1, real1, imag1);
            free(real0);
            free(imag0);
            free(real1);
            free(imag1);
        }
    }
}

// multiControlledPhaseFlip

void statevec_pauliX(Qureg qureg, int targetQubit) {
    statevec_unitary(qureg, targetQubit, (ComplexMatrix2) {
        .real={{0,1},{1,0}}, .imag={{0,0},{0,0}}
    });
}

void statevec_hadamard(Qureg qureg, int targetQubit) {
    double recRoot2 = 1.0/sqrt(2);
    statevec_unitary(qureg, targetQubit, (ComplexMatrix2) {
        .real={{recRoot2,recRoot2},{recRoot2,-recRoot2}}, .imag={{0,0},{0,0}}
    });
}

#pragma omp declare simd
inline long long insertBits(long long num, long long mask) {
    while (mask) {
        long long lowestSetBit = mask & -mask;
        long long lowerBits = lowestSetBit - 1;
        long long higherBits = ~lowerBits;
        num = ((num & higherBits) << 1) | lowestSetBit | (num & lowerBits);
        mask &= (mask - 1);
    }
    return num;
}

void statevec_multiControlledPhaseShiftByTerm(Qureg qureg, int *controlQubits, int numControlQubits, Complex term) {
    size_t insertedIndex;
    qreal y_real, y_imag;
    qreal cosAngle = term.real, sinAngle = term.imag;
    size_t mask = getQubitBitMask(controlQubits, numControlQubits);

    size_t leftMask = (~(qureg.compressBlockSize - 1)) & mask;
    size_t rightMask = (qureg.compressBlockSize - 1) & mask;

    int numLeftMaskOnes = __builtin_popcountll(leftMask);
    int numRightMaskOnes = __builtin_popcountll(rightMask);

    for (size_t blockIndexBase = 0; blockIndexBase < (qureg.compressBlockCount >> numLeftMaskOnes); blockIndexBase++) {
        size_t blockIndex = insertBits(blockIndexBase, leftMask / qureg.compressBlockSize);
        double* real = malloc(qureg.compressBlockSize * sizeof(double));
        double* imag = malloc(qureg.compressBlockSize * sizeof(double));
        decompressBlock(&qureg, blockIndex, real, imag);
        #pragma omp parallel for simd shared (controlQubits, numControlQubits, cosAngle, sinAngle, mask) private (insertedIndex, y_real, y_imag)
        for (size_t baseIndex = 0; baseIndex < (qureg.compressBlockSize >> numRightMaskOnes); baseIndex++) {
            insertedIndex = insertBits(baseIndex, rightMask);

            y_real = cosAngle * real[insertedIndex] - sinAngle * imag[insertedIndex];
            y_imag = cosAngle * imag[insertedIndex] + sinAngle * real[insertedIndex];

            real[insertedIndex] = y_real;
            imag[insertedIndex] = y_imag;
        }
        compressBlock(&qureg, blockIndex, real, imag);
        free(real);
        free(imag);
    }
}

OVERRIDING void statevec_multiControlledPhaseFlip(Qureg qureg, int *controlQubits, int numControlQubits) {
    statevec_multiControlledPhaseShiftByTerm(qureg, controlQubits, numControlQubits, (Complex) {.real=-1, .imag=0});
}
