#ifndef ExampleChannels_H
#define ExampleChannels_H

#include "../../QuantumChannels.h"
#include "QuEST.h"
#include <stdio.h>

void ApplyOneQubitDephaseChannel(Qureg qureg, const int targetQubit, qreal prob);

void ApplyOneQubitDepolariseChannel(Qureg qureg, const int targetQubit, qreal prob);

void ApplyOneQubitDampingChannel(Qureg qureg, const int targetQubit, qreal prob);

void ApplyTwoQubitDephaseChannel(Qureg qureg, const int qubit1, const int qubit2, qreal prob);

void rotXtest(Qureg qureg, const int targetQubit);

#endif