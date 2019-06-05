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
} OneQubitSuperOperator;

typedef struct OneQubitSparseSuperOperator
{
	int length;
	int indexes[16][2];
	qreal valuesRe[16];
	qreal valuesIm[16];
} OneQubitSparseSuperOperator;


typedef struct OneQubitChannel
{	
	OneQubitSuperOperator SupOp;
	OneQubitSparseSuperOperator SparseSupOp;
} OneQubitChannel;

// This calculates the superoperator from the Krauss operators A and B and adds it to the superoperator C
void KraussOperator2SuperOperator(OneQubitKraussOperator *A, OneQubitKraussOperator *B, OneQubitSuperOperator *C);

// This calculates the sparse represenation of the superoperator
void CalculateOneQubitSparseSuperOperator(OneQubitChannel *thisChannel);

void ApplyOneQubitChannel_local(Qureg qureg, const int targetQubit, OneQubitSparseSuperOperator supop);

void ApplyArbitraryKraussMap(Qureg qureg, const int targetQubit, OneQubitKraussOperator operators[4]);

#endif