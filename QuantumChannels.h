#ifndef QuantumChannels_H
#define QuantumChannels_H

#include <stdio.h>
#include "QuEST.h"

typedef struct OneQubitKrausOperator
{
	qreal real[2][2];
	qreal imag[2][2];
} OneQubitKrausOperator;

typedef struct OneQubitSuperOperator
{
	qreal real[4][4];
	qreal imag[4][4];
	int isComplex;
} OneQubitSuperOperator;

typedef struct TwoQubitSuperOperator
{
	qreal real[16][16];
	qreal imag[16][16];
	int isComplex;
} TwoQubitSuperOperator;


// This calculates the superoperator from the Kraus operators A and B and adds it to the superoperator C
void KrausOperator2SuperOperator(OneQubitKrausOperator *A, OneQubitKrausOperator *B, OneQubitSuperOperator *C);

void ApplyOneQubitChannel_local(Qureg qureg, const int targetQubit, OneQubitSuperOperator supop);

void ApplyOneQubitKrausMap(Qureg qureg, const int targetQubit, OneQubitKrausOperator *operators, int numberOfOperators);

void ApplyOneQubitPauliChannel(Qureg qureg, const int targetQubit, qreal probX, qreal probY, qreal probZ);

void ApplyTwoQubitChannel_local(Qureg qureg, const int qubit1, const int qubit2, TwoQubitSuperOperator supop);

#endif