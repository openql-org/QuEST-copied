#ifndef QuantumChannels_H
#define QuantumChannels_H

#include <stdio.h>
#include "QuEST.h"

typedef struct OneQubitKraussOperator
{
	qreal real[2][2];
	qreal imag[2][2];
} OneQubitKraussOperator;

typedef struct OneQubitSuperOperator
{
	qreal real[4][4];
	qreal imag[4][4];
	int isComplex;
} OneQubitSuperOperator;


// This calculates the superoperator from the Krauss operators A and B and adds it to the superoperator C
void KraussOperator2SuperOperator(OneQubitKraussOperator *A, OneQubitKraussOperator *B, OneQubitSuperOperator *C);

void ApplyOneQubitChannel_local(Qureg qureg, const int targetQubit, OneQubitSuperOperator supop);

void ApplyArbitraryKraussMap(Qureg qureg, const int targetQubit, OneQubitKraussOperator *operators, int numberOfOperators);

#endif