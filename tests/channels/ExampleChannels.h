#ifndef ExampleChannels_H
#define ExampleChannels_H

#include "../../QuantumChannels.h"
#include "QuEST.h"
#include <stdio.h>

void ApplyOneQubitDephaseChannel(Qureg qureg, const int targetQubit, qreal prob);

void ApplyOneQubitDepolariseChannel(Qureg qureg, const int targetQubit, qreal prob);

void ApplyOneQubitDampingChannel(Qureg qureg, const int targetQubit, qreal prob);

void ApplyOneQubitUnitalChannel(Qureg qureg, const int targetQubit, qreal prefactors[4]);

void ApplyArbitraryKraussMap(Qureg qureg, const int targetQubit, OneQubitKraussOperator operators[4]);

#endif