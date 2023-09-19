# include "QuEST.h"
# include "QuEST_internal.h"
# include "QuEST_precision.h"
# include "mt19937ar.h"
# include <stdio.h>
# include <stdlib.h>
# include <string.h>

#define OVERRIDING

// statevec functions

// void statevec_reportStateToScreen(Qureg qureg, QuESTEnv env, int reportRank) {}

// int statevec_compareStates(Qureg mq1, Qureg mq2, qreal precision) {}

// int statevec_initStateFromSingleFile(Qureg *qureg, char filename[200], QuESTEnv env) {}

// void statevec_initStateOfSingleQubit(Qureg *qureg, int qubitId, int outcome) {}

// void statevec_createQureg(Qureg *qureg, int numQubits, QuESTEnv env) {}

// void statevec_destroyQureg(Qureg qureg, QuESTEnv env) {}

// void statevec_initBlankState(Qureg qureg) {}

// void statevec_initZeroState(Qureg qureg) {}

// void statevec_initPlusState(Qureg qureg) {}

// void statevec_initDebugState(Qureg qureg) {}

// void statevec_initClassicalState(Qureg qureg, long long int stateInd) {}

// void statevec_setAmps(Qureg qureg, long long int startInd, qreal* reals, qreal* imags, long long int numAmps) {}

// void statevec_cloneQureg(Qureg targetQureg, Qureg copyQureg) {}

// qreal statevec_getRealAmp(Qureg qureg, long long int index) {}

// qreal statevec_getImagAmp(Qureg qureg, long long int index) {}

// qreal statevec_calcTotalProb(Qureg qureg) {}

// Complex statevec_calcInnerProduct(Qureg bra, Qureg ket) {}

// qreal statevec_calcProbOfOutcome(Qureg qureg, int measureQubit, int outcome) {}

// void statevec_calcProbOfAllOutcomes(qreal* retProbs, Qureg qureg, int* qubits, int numQubits) {}

// void statevec_collapseToKnownProbOutcome(Qureg qureg, int measureQubit, int outcome, qreal outcomeProb) {}

// void statevec_setWeightedQureg(Complex fac1, Qureg qureg1, Complex fac2, Qureg qureg2, Complex facOut, Qureg out) {}

// void statevec_applyDiagonalOp(Qureg qureg, DiagonalOp op) {}

// Complex statevec_calcExpecDiagonalOp(Qureg qureg, DiagonalOp op) {}

// void statevec_applyPhaseFuncOverrides(Qureg qureg, int* qubits, int numQubits, enum bitEncoding encoding, qreal* coeffs, qreal* exponents, int numTerms, long long int* overrideInds, qreal* overridePhases, int numOverrides, int conj) {}

// void statevec_applyMultiVarPhaseFuncOverrides(Qureg qureg, int* qubits, int* numQubitsPerReg, int numRegs, enum bitEncoding encoding, qreal* coeffs, qreal* exponents, int* numTermsPerReg, long long int* overrideInds, qreal* overridePhases, int numOverrides, int conj) {}

// void statevec_applyParamNamedPhaseFuncOverrides(Qureg qureg, int* qubits, int* numQubitsPerReg, int numRegs, enum bitEncoding encoding, enum phaseFunc functionNameCode, qreal* params, int numParams, long long int* overrideInds, qreal* overridePhases, int numOverrides, int conj) {}

// void statevec_copySubstateToGPU(Qureg qureg, long long int startInd, long long int numAmps) {}

// void statevec_copySubstateFromGPU(Qureg qureg, long long int startInd, long long int numAmps) {}

// void copyStateToGPU(Qureg qureg) {}

// void copyStateFromGPU(Qureg qureg) {}

// utils

// QuESTEnv createQuESTEnv(void) {}

// void destroyQuESTEnv(QuESTEnv env) {}

// void syncQuESTEnv(QuESTEnv env) {}

// int syncQuESTSuccess(int successCode) {}

// void reportQuESTEnv(QuESTEnv env) {}

// void getEnvironmentString(QuESTEnv env, char str[200]) {}

#define STRINGIFY(...) #__VA_ARGS__
#define MANY(...) __VA_ARGS__
#define REPEAT_0(X, ...)
#define REPEAT_1(X, ...) X(__VA_ARGS__)
#define REPEAT_2(X, ...) REPEAT_1(X, __VA_ARGS__) REPEAT_1(X, __VA_ARGS__)
#define REPEAT_4(X, ...) REPEAT_2(X, __VA_ARGS__) REPEAT_2(X, __VA_ARGS__)
#define REPEAT_8(X, ...) REPEAT_4(X, __VA_ARGS__) REPEAT_4(X, __VA_ARGS__)
#define REPEAT_16(X, ...) REPEAT_8(X, __VA_ARGS__) REPEAT_8(X, __VA_ARGS__)

// general unitary gates

#define COMPUTE(compute) \
    compute; \
    tt++; \

#define LOOP_CASE(n, compute, r) \
    case n: \
        _Pragma("omp for simd") \
        for (size_t t = 0; t < qureg.numAmpsPerChunk; t += (2*K)) \
        { \
            tt = 0; \
            REPEAT_##r(COMPUTE, compute) \
        } \
        break; \

#define declareGeneralUnitaryGateFunc(name, pre, shared_var, private_var, compute, ...) \
    OVERRIDING void statevec_##name(Qureg qureg, int targetQubit, ##__VA_ARGS__) { \
        size_t K = 1LL << targetQubit; \
        size_t tt; \
        pre; \
        ComplexArray amp = qureg.stateVec; \
        _Pragma(STRINGIFY(omp parallel shared (K shared_var) private (tt private_var))) \
        { \
            switch (targetQubit) { \
            LOOP_CASE(0, compute, 1) \
            LOOP_CASE(1, compute, 2) \
            LOOP_CASE(2, compute, 4) \
            LOOP_CASE(3, compute, 8) \
            LOOP_CASE(4, compute, 16) \
            default: \
                _Pragma("omp for simd collapse(2)") \
                for (size_t t = 0; t < qureg.numAmpsPerChunk; t += (2LL*K)) \
                { \
                    for (tt = 0; tt < K; tt++) \
                    { \
                        compute \
                    } \
                } \
            } \
        } \
    }

#define declareUnitaryGateFunc(name, pre, shared_var, ar, ai, br, bi, cr, ci, dr, di, ...) \
    declareGeneralUnitaryGateFunc( \
        name, \
        MANY(qreal x_real, y_real, x_imag, y_imag; pre), \
        MANY(shared_var), \
        MANY(, x_real, y_real, x_imag, y_imag), \
        { \
            x_real = (ar * amp.real[t+tt] - ai * amp.imag[t+tt]) + (br * amp.real[t+tt+K] - bi * amp.imag[t+tt+K]); \
            x_imag = (ar * amp.imag[t+tt] + ai * amp.real[t+tt]) + (br * amp.imag[t+tt+K] + bi * amp.real[t+tt+K]); \
            y_real = (cr * amp.real[t+tt] - ci * amp.imag[t+tt]) + (dr * amp.real[t+tt+K] - di * amp.imag[t+tt+K]); \
            y_imag = (cr * amp.imag[t+tt] + ci * amp.real[t+tt]) + (dr * amp.imag[t+tt+K] + di * amp.real[t+tt+K]); \
            amp.real[t+tt] = x_real; \
            amp.imag[t+tt] = x_imag; \
            amp.real[t+tt+K] = y_real; \
            amp.imag[t+tt+K] = y_imag; \
        }, \
        ##__VA_ARGS__ \
    )

#define declareUpperUnitaryGateFunc(name, pre, shared_var, dr, di, ...) \
    declareGeneralUnitaryGateFunc( \
        name, \
        MANY(qreal y_real, y_imag; pre), \
        MANY(shared_var), \
        MANY(, y_real, y_imag), \
        { \
            y_real = (dr * amp.real[t+tt+K] - di * amp.imag[t+tt+K]); \
            y_imag = (dr * amp.imag[t+tt+K] + di * amp.real[t+tt+K]); \
            amp.real[t+tt+K] = y_real; \
            amp.imag[t+tt+K] = y_imag; \
        }, \
        ##__VA_ARGS__ \
    )


// general controlled unitary gates

#define CT_LOOP_CASE(k, compute, K) \
    case k: \
        _Pragma("omp for simd collapse(2)") \
        for (size_t c = 0; c < qureg.numAmpsPerChunk; c += (2LL*KC)) \
        { \
            for (size_t t = 0; t < KC; t += (2LL*K)) \
            { \
                tt = 0; \
                REPEAT_##K(COMPUTE, compute) \
            } \
        } \
        break; \

#define TC_LOOP_CASE(kc, compute, KC) \
    case kc: \
        _Pragma("omp for simd collapse(2)") \
        for (size_t t = 0; t < qureg.numAmpsPerChunk; t += (2LL*K)) \
        { \
            for (size_t c = 0; c < K; c += (2LL*KC)) \
            { \
                tt = 0; \
                REPEAT_##KC(COMPUTE, compute) \
            } \
        } \
        break; \

#define declareGeneralControlledUnitaryGateFunc(name, pre, shared_var, private_var, compute, ...) \
    OVERRIDING void statevec_##name(Qureg qureg, int controlQubit, int targetQubit, ##__VA_ARGS__) { \
        size_t KC = 1LL << controlQubit; \
        size_t K = 1LL << targetQubit; \
        size_t tt; \
        pre; \
        ComplexArray amp = qureg.stateVec; \
        _Pragma(STRINGIFY(omp parallel shared (K, KC shared_var) private (tt private_var))) \
        { \
            if (controlQubit > targetQubit) { \
                switch(targetQubit) { \
                CT_LOOP_CASE(0, compute, 1) \
                CT_LOOP_CASE(1, compute, 2) \
                CT_LOOP_CASE(2, compute, 4) \
                default: \
                    _Pragma("omp for simd collapse(3)") \
                    for (size_t c = 0; c < qureg.numAmpsPerChunk; c += (2LL*KC)) \
                    { \
                        for (size_t t = 0; t < KC; t += (2LL*K)) \
                        { \
                            for (tt = 0; tt < K; tt++) \
                            { \
                                compute \
                            } \
                        } \
                    } \
                } \
            } else { \
                switch(controlQubit) { \
                TC_LOOP_CASE(0, compute, 1) \
                TC_LOOP_CASE(1, compute, 2) \
                TC_LOOP_CASE(2, compute, 4) \
                default: \
                    _Pragma("omp for simd collapse(3)") \
                    for (size_t t = 0; t < qureg.numAmpsPerChunk; t += (2LL*K)) \
                    { \
                        for (size_t c = 0; c < K; c += (2LL*KC)) \
                        { \
                            for (tt = 0; tt < KC; tt++) \
                            { \
                                compute \
                            } \
                        } \
                    } \
                } \
            } \
        } \
    }

#define declareControlledUnitaryGateFunc(name, pre, shared_var, ar, ai, br, bi, cr, ci, dr, di, ...) \
    declareGeneralControlledUnitaryGateFunc( \
        name, \
        MANY(qreal x_real, y_real, x_imag, y_imag, *amp_pt_r, *amp_pt_i; pre), \
        MANY(shared_var), \
        MANY(, x_real, y_real, x_imag, y_imag, amp_pt_r, amp_pt_i), \
        { \
            amp_pt_r = amp.real + c + t + tt; \
            amp_pt_i = amp.imag + c + t + tt; \
            x_real = (ar * amp_pt_r[KC] - ai * amp_pt_i[KC]) + (br * amp_pt_r[KC|K] - bi * amp_pt_i[KC|K]); \
            x_imag = (ar * amp_pt_i[KC] + ai * amp_pt_r[KC]) + (br * amp_pt_i[KC|K] + bi * amp_pt_r[KC|K]); \
            y_real = (cr * amp_pt_r[KC] - ci * amp_pt_i[KC]) + (dr * amp_pt_r[KC|K] - di * amp_pt_i[KC|K]); \
            y_imag = (cr * amp_pt_i[KC] + ci * amp_pt_r[KC]) + (dr * amp_pt_i[KC|K] + di * amp_pt_r[KC|K]); \
            amp_pt_r[KC] = x_real; \
            amp_pt_i[KC] = x_imag; \
            amp_pt_r[KC|K] = y_real; \
            amp_pt_i[KC|K] = y_imag; \
        }, \
        ##__VA_ARGS__ \
    )

#define declareControlledUpperUnitaryGateFunc(name, pre, shared_var, dr, di, ...) \
    declareGeneralControlledUnitaryGateFunc( \
        name, \
        MANY(qreal y_real, y_imag, *amp_pt_r, *amp_pt_i; pre), \
        MANY(shared_var), \
        MANY(, y_real, y_imag, amp_pt_r, amp_pt_i), \
        { \
            amp_pt_r = amp.real + c + t + tt; \
            amp_pt_i = amp.imag + c + t + tt; \
            y_real = (dr * amp_pt_r[KC|K] - di * amp_pt_i[KC|K]); \
            y_imag = (dr * amp_pt_i[KC|K] + di * amp_pt_r[KC|K]); \
            amp_pt_r[KC|K] = y_real; \
            amp_pt_i[KC|K] = y_imag; \
        }, \
        ##__VA_ARGS__ \
    )

// controlled and uncontrolled unitary gates

#define declareGateFuncs(type, name1, name2, ...) \
    declare##type##GateFunc(name1, ##__VA_ARGS__) \
    declareControlled##type##GateFunc(name2, ##__VA_ARGS__)

// overridings

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

void statevec_controlledPhaseShiftByTerm(Qureg, int, int, Complex);

void statevec_multiControlledPhaseShiftByTerm(Qureg qureg, int *controlQubits, int numControlQubits, Complex term) {

    if (numControlQubits == 1) {
        statevec_phaseShiftByTerm(qureg, controlQubits[0], term);
        return;
    } else if (numControlQubits == 2 && controlQubits[0] >= 3 && controlQubits[1] >= 3) {
        statevec_controlledPhaseShiftByTerm(qureg, controlQubits[0], controlQubits[1], term);
        return;
    }

    ComplexArray amp = qureg.stateVec;
    size_t insertedIndex;
    qreal y_real, y_imag;
    qreal cosAngle = term.real, sinAngle = term.imag;
    long long mask = getQubitBitMask(controlQubits, numControlQubits);
    
    _Pragma("omp parallel shared (controlQubits, numControlQubits, cosAngle, sinAngle, mask) private (insertedIndex, y_real, y_imag)")
    {
        _Pragma("omp for simd")
        for (size_t baseIndex = 0; baseIndex < (qureg.numAmpsPerChunk >> numControlQubits); baseIndex++) {
            insertedIndex = insertBits(baseIndex, mask);

            y_real = cosAngle * amp.real[insertedIndex] - sinAngle * amp.imag[insertedIndex];
            y_imag = cosAngle * amp.imag[insertedIndex] + sinAngle * amp.real[insertedIndex];

            amp.real[insertedIndex] = y_real;
            amp.imag[insertedIndex] = y_imag;
        }
    }
}

OVERRIDING void statevec_multiControlledPhaseFlip(Qureg qureg, int *controlQubits, int numControlQubits) {
    statevec_multiControlledPhaseShiftByTerm(qureg, controlQubits, numControlQubits, (Complex) {.real=-1, .imag=0});
}

declareControlledUpperUnitaryGateFunc(controlledPhaseFlip,
    MANY(),
    MANY(),
    -1, 0
);

declareGateFuncs(UpperUnitary, phaseShiftByTerm, controlledPhaseShiftByTerm,
    MANY(qreal cosAngle = term.real, sinAngle = term.imag),
    MANY(, cosAngle, sinAngle),
    cosAngle, sinAngle,
    Complex term
);

OVERRIDING void statevec_controlledPhaseShift(Qureg qureg, int idQubit1, int idQubit2, qreal angle) {
    Complex term; 
    term.real = cos(angle); 
    term.imag = sin(angle);
    statevec_controlledPhaseShiftByTerm(qureg, idQubit1, idQubit2, term);
};

OVERRIDING void statevec_multiControlledPhaseShift(Qureg qureg, int *controlQubits, int numControlQubits, qreal angle) {
    statevec_multiControlledPhaseShiftByTerm(qureg, controlQubits, numControlQubits, (Complex) {.real=cos(angle), .imag=sin(angle)});
};

declareGateFuncs(Unitary, pauliX, controlledNot,
    MANY(),
    MANY(),
    0, 0, 1, 0, 1, 0, 0, 0
);

declareGateFuncs(Unitary, pauliY, controlledPauliY,
    MANY(),
    MANY(),
    0, 0, 0, -1, 0, 1, 0, 0
);

declareGateFuncs(Unitary, pauliYConj, controlledPauliYConj,
    MANY(),
    MANY(),
    0, 0, 0, 1, 0, -1, 0, 0
);
 
declareGateFuncs(Unitary, compactUnitary, controlledCompactUnitary,
    MANY(qreal ar = alpha.real, br = beta.real, ai = alpha.imag, bi = beta.imag),
    MANY(, ar, ai, br, bi),
    ar, ai, -br, bi, br, bi, ar, -ai,
    Complex alpha, Complex beta
);

declareGateFuncs(Unitary, unitary, controlledUnitary,
    MANY(qreal ar = u.real[0][0], br = u.real[0][1], cr = u.real[1][0], dr = u.real[1][1],
               ai = u.imag[0][0], bi = u.imag[0][1], ci = u.imag[1][0], di = u.imag[1][1]),
    MANY(, ar, ai, br, bi, cr, ci, dr, di),
    ar, ai, br, bi, cr, ci, dr, di,
    ComplexMatrix2 u
);

// OVERRIDING void statevec_multiControlledTwoQubitUnitary(Qureg qureg, long long int ctrlMask, int targetQubit1, int targetQubit2, ComplexMatrix4 u) {}

// long long removeBits(long long original, unsigned long long mask) {
//     long long result = 0;
//     long long bitSelector = 1;
//     mask = ~mask;
//     for (long long i = 1; mask; i <<= 1, mask >>= 1) {
//         if (mask & 1) {
//             if (original & i) {
//                 result |= bitSelector;
//             }
//             bitSelector <<= 1;
//         }
//     }
//     return result;
// }

// #pragma omp declare simd
// inline long long insertBits2(long long num, long long mask, long long content) {
//     while (mask) {
//         long long lowestSetBit = mask & -mask;
//         long long lowerBits = lowestSetBit - 1;
//         long long higherBits = ~lowerBits;
//         num = ((num & higherBits) << 1) | (lowestSetBit & content) | (num & lowerBits);
//         mask &= (mask - 1);
//     }
//     return num;
// }

// #pragma omp declare simd
// inline long long convertMaskContent(size_t i, int* targs) {
//     long long content = 0;
//     size_t cur = 0;
//     while (i) {
//         if (i & 1)
//             content |= 1LL << targs[cur];
//         i >>= 1;
//         cur++;
//     }
//     return content;
// }

// OVERRIDING void statevec_multiControlledMultiQubitUnitary(Qureg qureg, long long int ctrlMask, int* targs, int numTargs, ComplexMatrixN u) {
//     long long targetMask = getQubitBitMask(targs, numTargs);
//     long long mask = ctrlMask | targetMask;
//     int numMaskedQubits = __builtin_popcountll(mask);
//     long long ctrlMask2 = removeBits(ctrlMask, targetMask);

//     ComplexArray amp = qureg.stateVec;
//     size_t baseIndex2;
//     qreal *result_real;
//     qreal *result_imag;
//     size_t *insertedIndeces;

//     _Pragma("omp parallel shared (ctrlMask, ctrlMask2, targetMask, mask, numMaskedQubits, u) private (baseIndex2, result_real, result_imag, insertedIndeces)")
//     {
//         result_real = malloc(sizeof(*result_real) * (1LL << numTargs));
//         result_imag = malloc(sizeof(*result_imag) * (1LL << numTargs));
//         insertedIndeces = malloc(sizeof(*insertedIndeces) * (1LL << numTargs));

//         _Pragma("omp for simd")
//         for (size_t baseIndex = 0; baseIndex < (qureg.numAmpsPerChunk >> numMaskedQubits); baseIndex++) {
//             baseIndex2 = insertBits(baseIndex, ctrlMask2);
//             for (size_t i = 0; i < (1LL << numTargs); i++) {
//                 result_real[i] = result_imag[i] = 0;
//                 insertedIndeces[i] = insertBits2(baseIndex2, targetMask, convertMaskContent(i, targs));
//             }
//             for (size_t i = 0; i < (1LL << numTargs); i++) {
//                 for (size_t j = 0; j < (1LL << numTargs); j++) {
//                     result_real[i] += u.real[i][j] * amp.real[insertedIndeces[j]] - u.imag[i][j] * amp.imag[insertedIndeces[j]];
//                     result_imag[i] += u.real[i][j] * amp.imag[insertedIndeces[j]] + u.imag[i][j] * amp.real[insertedIndeces[j]];
//                 }
//             }
//             for (size_t i = 0; i < (1LL << numTargs); i++) {
//                 amp.real[insertedIndeces[i]] = result_real[i];
//                 amp.imag[insertedIndeces[i]] = result_imag[i];
//             }
//         }

//         free(result_real);
//         free(result_imag);
//         free(insertedIndeces);
//     }
// }

// OVERRIDING void statevec_multiControlledUnitary(Qureg qureg, long long int ctrlQubitsMask, long long int ctrlFlipMask, int targetQubit, ComplexMatrix2 u) {}

declareUnitaryGateFunc(hadamard,
    MANY(qreal recRoot2 = 1.0/sqrt(2)),
    MANY(, recRoot2),
    recRoot2, 0, recRoot2, 0, recRoot2, 0, -recRoot2, 0
);

// OVERRIDING void statevec_multiControlledMultiQubitNot(Qureg qureg, int ctrlMask, int targMask) {}

// OVERRIDING void statevec_swapQubitAmps(Qureg qureg, int qb1, int qb2) {}

// OVERRIDING void statevec_multiRotateZ(Qureg qureg, long long int mask, qreal angle) {}

// OVERRIDING void statevec_multiControlledMultiRotateZ(Qureg qureg, long long int ctrlMask, long long int targMask, qreal angle) {}
