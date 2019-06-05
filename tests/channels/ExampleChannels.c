#include "ExampleChannels.h"


void ApplyOneQubitDephaseChannel(Qureg qureg, const int targetQubit, qreal prob)
{	
	//DO THE CHECKS HERE -- e.g., validate p

	//Define the channel via its superoperator
	OneQubitChannel thisChannel = {
		.SupOp.real = {
			{1, 0, 0, 0}, {0, 1 - 2*prob, 0, 0},
			{0, 0, 1 - 2*prob, 0}, {0, 0, 0, 1}
			}
	};
	
	//Calculate the sparse representation of the superoperator
	CalculateOneQubitSparseSuperOperator(&thisChannel);
	
	//Apply channel to the qubit
	ApplyOneQubitChannel_local(qureg, targetQubit, thisChannel.SparseSupOp);	
}	


void ApplyOneQubitDepolariseChannel(Qureg qureg, const int targetQubit, qreal prob)
{	
	//DO THE CHECKS HERE -- e.g., validate p

	//Define the channel via its superoperator
	OneQubitChannel thisChannel = {
		.SupOp.real = {
			{1 - (2*prob)/3, 0, 0, (2*prob)/3},
			{0, 1 - (4*prob)/3, 0, 0},
			{0, 0, 1 - (4*prob)/3, 0}, 
			{(2*prob)/3, 0, 0, 1 - (2*prob)/3}
			}
	};
	
	//Calculate the sparse representation of the superoperator
	CalculateOneQubitSparseSuperOperator(&thisChannel);
	
	//Apply channel to the qubit
	ApplyOneQubitChannel_local(qureg, targetQubit, thisChannel.SparseSupOp);	
}	


void ApplyOneQubitDampingChannel(Qureg qureg, const int targetQubit, qreal prob)
{	
	//DO THE CHECKS HERE -- e.g., validate p

	//Define the channel via its superoperator
	OneQubitChannel thisChannel = {
		.SupOp.real = {
			{1, 0, 0, prob}, {0, sqrt(1 - prob), 0, 0}, 
			{0, 0, sqrt(1 - prob), 0}, {0, 0, 0, 1 - prob}
			}
	};
	
	//Calculate the sparse representation of the superoperator
	CalculateOneQubitSparseSuperOperator(&thisChannel);
	
	//Apply channel to the qubit
	ApplyOneQubitChannel_local(qureg, targetQubit, thisChannel.SparseSupOp);	
}



void ApplyOneQubitUnitalChannel(Qureg qureg, const int targetQubit, qreal prefactors[4])
{
	//DO the checks on the prefactros to verify that the channel is unital and completely positive
	OneQubitKraussOperator Pauli0 = {.real = {{prefactors[0] * 1, 0},{0, prefactors[0] * 1}}, .imag = {0}};
	OneQubitKraussOperator Pauli1 = {.real = {{0, prefactors[1] * 1},{prefactors[1] * 1, 0}}, .imag = {0}};
	OneQubitKraussOperator Pauli2 = {.imag = {{0, prefactors[2] * -1},{prefactors[2] * 1, 0}}, .real = {0}};
	OneQubitKraussOperator Pauli3 = {.real = {{prefactors[2] * 1, 0},{0, prefactors[2] * -1}}, .imag = {0}};
	
	
	OneQubitKraussOperator operators[4] = {Pauli0, Pauli1, Pauli2, Pauli3};
	
	ApplyArbitraryKraussMap(qureg, targetQubit, operators);
	
}