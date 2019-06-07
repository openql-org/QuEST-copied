#include "QuantumChannels.h"



void ApplyOneQubitChannel_local(Qureg qureg, const int targetQubit, OneQubitSuperOperator supop) {

    const long long int numTasks = qureg.numAmpsPerChunk;
    const long long int innerMask = 1LL << targetQubit;
    const long long int outerMask = 1LL << (targetQubit + (qureg.numQubitsRepresented));
    const long long int totMask = innerMask|outerMask;

    long long int thisTask;
	
	
	qreal rhoReA, rhoReB, rhoReC, rhoReD;
	qreal rhoImA, rhoImB, rhoImC, rhoImD;
	long long int indA, indB, indC, indD;
	
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
			indA = thisTask;					// element A -- upper left
			indB = thisTask | innerMask;		// element B -- lower left
			indC = thisTask | outerMask;		// element C -- upper right
			indD = thisTask | totMask;			// element D -- lower right	
			
			//store current values
			//i.e., copy elements of the density matrix into a vector |rho>
			rhoReA = qureg.stateVec.real[indA];
			rhoReB = qureg.stateVec.real[indB];
			rhoReC = qureg.stateVec.real[indC];
			rhoReD = qureg.stateVec.real[indD];
			
			rhoImA = qureg.stateVec.imag[indA];
			rhoImB = qureg.stateVec.imag[indB];
			rhoImC = qureg.stateVec.imag[indC];
			rhoImD = qureg.stateVec.imag[indD];
		
			// apply the superoperator matrix to the vector |rho>
			qureg.stateVec.real[indA] = r0c0*rhoReA + r0c1*rhoReB + r0c2*rhoReC + r0c3*rhoReD;
			qureg.stateVec.imag[indA] = r0c0*rhoImA + r0c1*rhoImB + r0c2*rhoImC + r0c3*rhoImD;
			
			qureg.stateVec.real[indB] = r1c0*rhoReA + r1c1*rhoReB + r1c2*rhoReC + r1c3*rhoReD;
			qureg.stateVec.imag[indB] = r1c0*rhoImA + r1c1*rhoImB + r1c2*rhoImC + r1c3*rhoImD;
			
			qureg.stateVec.real[indC] = r2c0*rhoReA + r2c1*rhoReB + r2c2*rhoReC + r2c3*rhoReD;
			qureg.stateVec.imag[indC] = r2c0*rhoImA + r2c1*rhoImB + r2c2*rhoImC + r2c3*rhoImD;
					
			qureg.stateVec.real[indD] = r3c0*rhoReA + r3c1*rhoReB + r3c2*rhoReC + r3c3*rhoReD;
			qureg.stateVec.imag[indD] = r3c0*rhoImA + r3c1*rhoImB + r3c2*rhoImC + r3c3*rhoImD;
			
			if (isComplex == 1){
				qureg.stateVec.real[indA] -= Imr0c0*rhoImA + Imr0c1*rhoImB + Imr0c2*rhoImC + Imr0c3*rhoImD;
				qureg.stateVec.imag[indA] += Imr0c0*rhoReA + Imr0c1*rhoReB + Imr0c2*rhoReC + Imr0c3*rhoReD;
				
				qureg.stateVec.real[indB] -= Imr1c0*rhoImA + Imr1c1*rhoImB + Imr1c2*rhoImC + Imr1c3*rhoImD;
				qureg.stateVec.imag[indB] += Imr1c0*rhoReA + Imr1c1*rhoReB + Imr1c2*rhoReC + Imr1c3*rhoReD;
				
				qureg.stateVec.real[indC] -= Imr2c0*rhoImA + Imr2c1*rhoImB + Imr2c2*rhoImC + Imr2c3*rhoImD;
				qureg.stateVec.imag[indC] += Imr2c0*rhoReA + Imr2c1*rhoReB + Imr2c2*rhoReC + Imr2c3*rhoReD;
				
				qureg.stateVec.real[indD] -= Imr3c0*rhoImA + Imr3c1*rhoImB + Imr3c2*rhoImC + Imr3c3*rhoImD;
				qureg.stateVec.imag[indD] += Imr3c0*rhoReA + Imr3c1*rhoReB + Imr3c2*rhoReC + Imr3c3*rhoReD;
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
