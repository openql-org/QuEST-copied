#include "ExampleChannels.h"


void ApplyOneQubitDephaseChannel(Qureg qureg, const int targetQubit, qreal prob)
{	
	//DO THE CHECKS HERE -- e.g., validate p

	//Define the channel via its superoperator
	OneQubitSuperOperator supop = {
		.real = {
			{1, 0, 0, 0}, {0, 1 - 2*prob, 0, 0},
			{0, 0, 1 - 2*prob, 0}, {0, 0, 0, 1}
			},
		.isComplex = 0
	};
	
	//Calculate the sparse representation of the superoperator
	//CalculateOneQubitSparseSuperOperator(&thisChannel);
	
	//Apply channel to the qubit
	ApplyOneQubitChannel_local(qureg, targetQubit, supop);	
}	


void ApplyOneQubitDepolariseChannel(Qureg qureg, const int targetQubit, qreal prob)
{	
	//DO THE CHECKS HERE -- e.g., validate p

	//Define the channel via its superoperator
	OneQubitSuperOperator supop = {
		.real = {
			{1 - (2*prob)/3, 0, 0, (2*prob)/3},
			{0, 1 - (4*prob)/3, 0, 0},
			{0, 0, 1 - (4*prob)/3, 0}, 
			{(2*prob)/3, 0, 0, 1 - (2*prob)/3}
			},
		.isComplex = 0
	};
	
	//ApplyOneQubitChannel_local(qureg, targetQubit, thisChannel.SparseSupOp);	
	ApplyOneQubitChannel_local(qureg, targetQubit, supop);	
}	


void ApplyOneQubitDampingChannel(Qureg qureg, const int targetQubit, qreal prob)
{	
	//DO THE CHECKS HERE -- e.g., validate p

	//Define the channel via its superoperator
	OneQubitSuperOperator supop = {
		.real = {
			{1, 0, 0, prob}, {0, sqrt(1 - prob), 0, 0}, 
			{0, 0, sqrt(1 - prob), 0}, {0, 0, 0, 1 - prob}
			},
		.isComplex = 0
	};
	
	//Apply channel to the qubit
	ApplyOneQubitChannel_local(qureg, targetQubit, supop);	
}
