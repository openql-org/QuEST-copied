#include "QuantumChannels.h"

void ApplyTwoQubitChannel_local(Qureg qureg, const int qubit1, const int qubit2, TwoQubitSuperOperator supop) {

    const long long int numTasks = qureg.numAmpsPerChunk;
    const long long int innerMask1 = 1LL << qubit1;
    const long long int outerMask1 = 1LL << (qubit1 + (qureg.numQubitsRepresented));
	const long long int innerMask2 = 1LL << qubit2;
    const long long int outerMask2 = 1LL << (qubit2 + (qureg.numQubitsRepresented));
    const long long int totMask = outerMask1 | outerMask2 | innerMask1 | innerMask2;
	
	const long long int mask1 = innerMask1;
	const long long int mask2 = innerMask2;
	const long long int mask3 = innerMask1 | innerMask2;
	const long long int mask4 = outerMask1;
	const long long int mask5 = outerMask1 | innerMask1;
	const long long int mask6 = outerMask1 | innerMask2;
	const long long int mask7 = outerMask1 | innerMask1 | innerMask2;
	const long long int mask8 = outerMask2;
	const long long int mask9 = outerMask2 | innerMask1;
	const long long int mask10 = outerMask2 | innerMask2;
	const long long int mask11 = outerMask2 | innerMask1 | innerMask2;
	const long long int mask12 = outerMask1 | outerMask2;
	const long long int mask13 = outerMask1 | outerMask2 | innerMask1;
	const long long int mask14 = outerMask1 | outerMask2 | innerMask2;
	const long long int mask15 = totMask;
	
	long long int thisTask;
	long long int indexes[16];
	
	qreal rhoRe[16], rhoIm[16], newrhoRe[16], newrhoIm[16];
	
	const TwoQubitSuperOperator thisSupOp = supop;
	
    for (thisTask=0; thisTask<numTasks; thisTask++){
		if ((thisTask&totMask)==0){ //this element relates to both qubits in state 0 -- upper diagonal
			indexes[0] = thisTask;
			indexes[1] = thisTask | mask1;
			indexes[2] = thisTask | mask2;
			indexes[3] = thisTask | mask3;
			indexes[4] = thisTask | mask4;
			indexes[5] = thisTask | mask5;
			indexes[6] = thisTask | mask6;
			indexes[7] = thisTask | mask7;
			indexes[8] = thisTask | mask8;
			indexes[9] = thisTask | mask9;
			indexes[10] = thisTask | mask10;
			indexes[11] = thisTask | mask11;
			indexes[12] = thisTask | mask12;
			indexes[13] = thisTask | mask13;
			indexes[14] = thisTask | mask14;
			indexes[15] = thisTask | mask15;
			
			
			
			//store current values of the density matrix
			//and set new values to 0
			for(int i = 0; i < 16; i++){
				rhoRe[i] = qureg.stateVec.real[indexes[i]];
				newrhoRe[i] = 0.;
				rhoIm[i] = qureg.stateVec.imag[indexes[i]];
				newrhoIm[i] = 0.;
			}
			
			for(int i = 0; i < 16; i++){
				for(int j = 0; j < 16; j++){
					newrhoRe[i] += thisSupOp.real[i][j]*rhoRe[j];
					newrhoIm[i] += thisSupOp.real[i][j]*rhoIm[j];
					if (supop.isComplex){
						newrhoRe[i] -= thisSupOp.imag[i][j]*rhoIm[j];
						newrhoIm[i] += thisSupOp.imag[i][j]*rhoRe[j];
					}
				}
			}

			//Update the density matrix
			for(int i = 0; i < 16; i++){
				qureg.stateVec.real[indexes[i]] = newrhoRe[i];
				qureg.stateVec.imag[indexes[i]] = newrhoIm[i];
			}
			
		}
	}
}



void ApplyOneQubitChannel_local(Qureg qureg, const int targetQubit, OneQubitSuperOperator supop) {

    const long long int numTasks = qureg.numAmpsPerChunk;
    const long long int innerMask = 1LL << targetQubit;
    const long long int outerMask = 1LL << (targetQubit + (qureg.numQubitsRepresented));
    const long long int totMask = innerMask|outerMask;

    long long int thisTask;
	
	
	qreal rhoRe0, rhoRe1, rhoRe2, rhoRe3;
	qreal rhoIm0, rhoIm1, rhoIm2, rhoIm3;
	long long int ind0, ind1, ind2, ind3;
	
	//Unpack primitives of the super operator
	const int isComplex = supop.isComplex;
	const qreal
	r0c0=supop.real[0][0], r0c1=supop.real[0][1], r0c2=supop.real[0][2], r0c3=supop.real[0][3], 
	r1c0=supop.real[1][0], r1c1=supop.real[1][1], r1c2=supop.real[1][2], r1c3=supop.real[1][3], 
	r2c0=supop.real[2][0], r2c1=supop.real[2][1], r2c2=supop.real[2][2], r2c3=supop.real[2][3], 
	r3c0=supop.real[3][0], r3c1=supop.real[3][1], r3c2=supop.real[3][2], r3c3=supop.real[3][3]; 

	const qreal
	Imr0c0=supop.imag[0][0], Imr0c1=supop.imag[0][1], Imr0c2=supop.imag[0][2], Imr0c3=supop.imag[0][3], 
	Imr1c0=supop.imag[1][0], Imr1c1=supop.imag[1][1], Imr1c2=supop.imag[1][2], Imr1c3=supop.imag[1][3], 
	Imr2c0=supop.imag[2][0], Imr2c1=supop.imag[2][1], Imr2c2=supop.imag[2][2], Imr2c3=supop.imag[2][3], 
	Imr3c0=supop.imag[3][0], Imr3c1=supop.imag[3][1], Imr3c2=supop.imag[3][2], Imr3c3=supop.imag[3][3]; 
	
    for (thisTask=0; thisTask<numTasks; thisTask++){
		if ((thisTask&totMask)==0){ //this element relates to targetQubit in state 0 -- upper diagonal
			//store indexes
			ind0 = thisTask;					// element A -- upper left
			ind1 = thisTask | innerMask;		// element B -- lower left
			ind2 = thisTask | outerMask;		// element C -- upper right
			ind3 = thisTask | totMask;			// element D -- lower right	
			
			//store current values
			//i.e., copy elements of the density matrix into a vector |rho>
			rhoRe0 = qureg.stateVec.real[ind0];
			rhoRe1 = qureg.stateVec.real[ind1];
			rhoRe2 = qureg.stateVec.real[ind2];
			rhoRe3 = qureg.stateVec.real[ind3];
			
			rhoIm0 = qureg.stateVec.imag[ind0];
			rhoIm1 = qureg.stateVec.imag[ind1];
			rhoIm2 = qureg.stateVec.imag[ind2];
			rhoIm3 = qureg.stateVec.imag[ind3];
		
			// apply the superoperator matrix to the vector |rho>
			qureg.stateVec.real[ind0] = r0c0*rhoRe0 + r0c1*rhoRe1 + r0c2*rhoRe2 + r0c3*rhoRe3; 
			qureg.stateVec.real[ind1] = r1c0*rhoRe0 + r1c1*rhoRe1 + r1c2*rhoRe2 + r1c3*rhoRe3; 
			qureg.stateVec.real[ind2] = r2c0*rhoRe0 + r2c1*rhoRe1 + r2c2*rhoRe2 + r2c3*rhoRe3; 
			qureg.stateVec.real[ind3] = r3c0*rhoRe0 + r3c1*rhoRe1 + r3c2*rhoRe2 + r3c3*rhoRe3; 

			qureg.stateVec.imag[ind0] = r0c0*rhoIm0 + r0c1*rhoIm1 + r0c2*rhoIm2 + r0c3*rhoIm3; 
			qureg.stateVec.imag[ind1] = r1c0*rhoIm0 + r1c1*rhoIm1 + r1c2*rhoIm2 + r1c3*rhoIm3; 
			qureg.stateVec.imag[ind2] = r2c0*rhoIm0 + r2c1*rhoIm1 + r2c2*rhoIm2 + r2c3*rhoIm3; 
			qureg.stateVec.imag[ind3] = r3c0*rhoIm0 + r3c1*rhoIm1 + r3c2*rhoIm2 + r3c3*rhoIm3; 

			if (isComplex == 1){
			qureg.stateVec.real[ind0] -= Imr0c0*rhoIm0 + Imr0c1*rhoIm1 + Imr0c2*rhoIm2 + Imr0c3*rhoIm3; 
			qureg.stateVec.real[ind1] -= Imr1c0*rhoIm0 + Imr1c1*rhoIm1 + Imr1c2*rhoIm2 + Imr1c3*rhoIm3; 
			qureg.stateVec.real[ind2] -= Imr2c0*rhoIm0 + Imr2c1*rhoIm1 + Imr2c2*rhoIm2 + Imr2c3*rhoIm3; 
			qureg.stateVec.real[ind3] -= Imr3c0*rhoIm0 + Imr3c1*rhoIm1 + Imr3c2*rhoIm2 + Imr3c3*rhoIm3; 

			qureg.stateVec.imag[ind0] += Imr0c0*rhoRe0 + Imr0c1*rhoRe1 + Imr0c2*rhoRe2 + Imr0c3*rhoRe3; 
			qureg.stateVec.imag[ind1] += Imr1c0*rhoRe0 + Imr1c1*rhoRe1 + Imr1c2*rhoRe2 + Imr1c3*rhoRe3; 
			qureg.stateVec.imag[ind2] += Imr2c0*rhoRe0 + Imr2c1*rhoRe1 + Imr2c2*rhoRe2 + Imr2c3*rhoRe3; 
			qureg.stateVec.imag[ind3] += Imr3c0*rhoRe0 + Imr3c1*rhoRe1 + Imr3c2*rhoRe2 + Imr3c3*rhoRe3; 
			}
			
		}
    }  
}



void KrausOperatorMultiply(OneQubitKrausOperator *A, OneQubitKrausOperator *B, OneQubitKrausOperator *C) 
{ 	// This calculates the matrix product C += ConjugateTranspose(A) x B and adds the result to C
	const int N = 2;
	qreal tempAr, tempAi, tempBr, tempBi; 
    for (int i = 0; i < N; i++) { 
        for (int j = 0; j < N; j++) { 
            for (int k = 0; k < N; k++) {
				tempAr = A->real[k][i];  //because it is the conjugate transpose of A
				tempAi = -A->imag[k][i]; //because it is the conjugate transpose of A
				tempBr = B->real[k][j];
				tempBi = B->imag[k][j];
				C->real[i][j] += tempAr*tempBr - tempAi*tempBi;
				C->imag[i][j] += tempAr*tempBi + tempAi*tempBr;
			}	 
		}
	}
} 

int ValidateOneQubitKrausMap(OneQubitKrausOperator *operators, int numberOfOperators)
{
		const int N = 2;
		OneQubitKrausOperator result = {.real = {0}, .imag = {0}};
		OneQubitKrausOperator id2 = {.real = {{1.,0.},{0.,1.}}, .imag = {0}};
		
		for (int i = 0; i < numberOfOperators; i++) {
			KrausOperatorMultiply(&operators[i], &operators[i], &result);
		}
		
		qreal distance = 0., re, im;
		// calculate the Hilbert-Schmidt distance between the
		// result and the identity matrix
		for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
			re = result.real[i][j] - id2.real[i][j];
			im = result.imag[i][j] - id2.imag[i][j];
			distance += fabs(re)*fabs(re) + fabs(im)*fabs(im);
		}
		}
		
		return(distance != 0.);
}


void KrausOperator2SuperOperator(OneQubitKrausOperator *A, OneQubitKrausOperator *B, OneQubitSuperOperator *C)
{ // This calculates the tensor product      C += conjugate(A) (x) B   and adds the result to the superopertor C
    qreal tempAr, tempAi, tempBr, tempBi;
	const int N = 2;
	
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) { 
			for (int k = 0; k < N; k++) { 
                for (int l = 0; l < N; l++) { 
					tempAr = A->real[i][j];
					tempAi = -A->imag[i][j]; // minus -- conjugate of A
					tempBr = B->real[k][l];
					tempBi = B->imag[k][l];
                    C->real[i*2 + k][j*2 + l] +=  tempAr * tempBr - tempAi * tempBi;
					C->imag[i*2 + k][j*2 + l] +=  tempAi * tempBr + tempAr * tempBi;
					if ( (C->imag[i*2 + k][j*2 + l] != 0.) && C->isComplex == 0)
					{
						C->isComplex = 1;
					}
                } 
            } 
        } 
    } 
}



void ApplyOneQubitKrausMap(Qureg qureg, const int targetQubit, OneQubitKrausOperator *operators, int numberOfOperators)
{   
	// Validate the Kraus map
	if ( ValidateOneQubitKrausMap(operators, numberOfOperators) )
	{ printf("The specified Kraus map is not a completely positive, trace preserving map\n");}
    // This if the test is not passed, it could still be allowed
	// to do the calculation --> the resulting rho is not physical
	
	//Initialize the channel with 0 superoperator
	OneQubitSuperOperator supop = {.real = {0}, .imag = {0}, .isComplex = 0 };
	
	//turn the Kraus operators into a superoperator
	for (int i = 0; i < numberOfOperators; i++) {
		KrausOperator2SuperOperator(&operators[i], &operators[i], &supop);
	}

	//Apply the superoperator to the qubit
	ApplyOneQubitChannel_local(qureg, targetQubit, supop);	
	
}


void ApplyOneQubitUnitalChannel(Qureg qureg, const int targetQubit, qreal probX, qreal probY, qreal probZ)
{
	// validate the probabilities here
	// not necessary to valideate -- the resulting Kraus operators will be validated anyway
	
	// Turn the probabilities into Kraus operators
	qreal prefactors[4] = {
		sqrt(1-(probX + probY + probZ)),
		sqrt(probX),
		sqrt(probY),
		sqrt(probZ)
	};
	OneQubitKrausOperator Pauli0 = {.real = {{prefactors[0] * 1, 0},{0, prefactors[0] * 1}}, .imag = {0}};
	OneQubitKrausOperator Pauli1 = {.real = {{0, prefactors[1] * 1},{prefactors[1] * 1, 0}}, .imag = {0}};
	OneQubitKrausOperator Pauli2 = {.imag = {{0, prefactors[2] * -1},{prefactors[2] * 1, 0}}, .real = {0}};
	OneQubitKrausOperator Pauli3 = {.real = {{prefactors[2] * 1, 0},{0, prefactors[2] * -1}}, .imag = {0}};
	
	OneQubitKrausOperator operators[4] = {Pauli0, Pauli1, Pauli2, Pauli3};
	ApplyOneQubitKrausMap(qureg, targetQubit, operators, 4);
	
}
