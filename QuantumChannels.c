#include "QuantumChannels.h"



void CalculateOneQubitSparseSuperOperator(OneQubitChannel *thisChannel)
{
	int length = 0;
	for (int m = 0; m < 4; m++){
		for (int n = 0; n < 4; n++){
			if (thisChannel->SupOp.real[m][n] != 0.){
				thisChannel->SparseSupOp.indexes[length][0] = m;
				thisChannel->SparseSupOp.indexes[length][1] = n;
				//Currently implemented only for real superoperators -- straightforward to extend
				thisChannel->SparseSupOp.valuesRe[length] = thisChannel->SupOp.real[m][n];
				/*printf("%d,    %d, %d,      %f \n",length,
					thisChannel->SparseSupOp.indexes[length][0],
					thisChannel->SparseSupOp.indexes[length][1],
					thisChannel->SparseSupOp.valuesRe[length]
				);*/
				length += 1;
			}
		}
	}
	
	thisChannel->SparseSupOp.length = length;
}

void ApplyOneQubitChannel_local(Qureg qureg, const int targetQubit, OneQubitSparseSuperOperator supop) {

    const long long int numTasks = qureg.numAmpsPerChunk;
    long long int innerMask = 1LL << targetQubit;
    long long int outerMask = 1LL << (targetQubit + (qureg.numQubitsRepresented));
    long long int totMask = innerMask|outerMask;

    long long int thisTask;
	

// Elements of the density matrix are denoted as ((A,C),(B,D))	
	long long int thisIndexes[4];
	qreal rhoRe[4], rhoIm[4], newrhoRe[4], newrhoIm[4];

	int currentRow, currentCol;

    for (thisTask=0; thisTask<numTasks; thisTask++){
		if ((thisTask&totMask)==0){ //this element relates to targetQubit in state 0 -- upper diagonal
			//Copy elements of the density matrix into a vector |rho>
			thisIndexes[0] = thisTask;					// element A -- upper left
			thisIndexes[1] = thisTask | innerMask;		// element B -- lower left
			thisIndexes[2] = thisTask | outerMask;		// element C -- upper right
			thisIndexes[3] = thisTask | totMask;		// element D -- lower right	
			
			//store current values of the density matrix
			//and set new values to 0
			for(int i = 0; i < 4; i++){
				rhoRe[i] = qureg.stateVec.real[thisIndexes[i]];
				newrhoRe[i] = 0.;
				rhoIm[i] = qureg.stateVec.imag[thisIndexes[i]];
				newrhoIm[i] = 0.;
			}
					
			for(int i = 0; i < supop.length; i++){
				currentCol = supop.indexes[i][0];
				currentRow = supop.indexes[i][1];
				//IMPLEMENTED ONLY FOR REAL SUPEROPERATORS -- straightforward to extend
				newrhoRe[currentCol] += supop.valuesRe[i]*rhoRe[currentRow];
				newrhoIm[currentCol] += supop.valuesRe[i]*rhoIm[currentRow];
			}
					
			//Update the density matrix
			for(int i = 0; i < 4; i++){
				qureg.stateVec.real[thisIndexes[i]] = newrhoRe[i];
				qureg.stateVec.imag[thisIndexes[i]] = newrhoIm[i];
			}
		}
    }  
}




void ApplyOneQubitChannel_nonsparse_local(Qureg qureg, const int targetQubit, OneQubitSuperOperator supop) {

    const long long int numTasks = qureg.numAmpsPerChunk;
    long long int innerMask = 1LL << targetQubit;
    long long int outerMask = 1LL << (targetQubit + (qureg.numQubitsRepresented));
    long long int totMask = innerMask|outerMask;

    long long int thisTask;
	

// Elements of the density matrix are denoted as ((A,C),(B,D))	
	long long int thisIndexes[4];
	qreal rhoRe[4], rhoIm[4], newrhoRe[4], newrhoIm[4];

    for (thisTask=0; thisTask<numTasks; thisTask++){
		if ((thisTask&totMask)==0){ //this element relates to targetQubit in state 0 -- upper diagonal
			//Copy elements of the density matrix into a vector |rho>
			thisIndexes[0] = thisTask;					// element A -- upper left
			thisIndexes[1] = thisTask | innerMask;		// element B -- lower left
			thisIndexes[2] = thisTask | outerMask;		// element C -- upper right
			thisIndexes[3] = thisTask | totMask;		// element D -- lower right	
			
			//store current values of the density matrix
			//and set new values to 0
			for(int i = 0; i < 4; i++){
				rhoRe[i] = qureg.stateVec.real[thisIndexes[i]];
				newrhoRe[i] = 0.;
				rhoIm[i] = qureg.stateVec.imag[thisIndexes[i]];
				newrhoIm[i] = 0.;
			}
					
			for(int i = 0; i < 4; i++){
				for(int j = 0; j < 4; j++){
					//IMPLEMENTED ONLY FOR REAL SUPEROPERATORS -- straightforward to extend
					newrhoRe[i] += supop.real[i][j]*rhoRe[j];
					newrhoIm[i] += supop.real[i][j]*rhoIm[j];
				}
			}
					
			//Update the density matrix
			for(int i = 0; i < 4; i++){
				qureg.stateVec.real[thisIndexes[i]] = newrhoRe[i];
				qureg.stateVec.imag[thisIndexes[i]] = newrhoIm[i];
			}
		}
    }  
}


void KraussOperator2SuperOperator(OneQubitKraussOperator *A, OneQubitKraussOperator *B, OneQubitSuperOperator *C)
{ // This calculates the tensor product      C += conjugate(A) (x) B   and adds the result to the superopertor C
    qreal tempAr, tempAi, tempBr, tempBi;
  
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) { 
			for (int k = 0; k < 2; k++) { 
                for (int l = 0; l < 2; l++) { 
                    // Each element of matrix A is 
                    // multiplied by the matrix B 
                    // and stored as Matrix C 
					tempAr = A->real[i][j];
					tempAi = A->imag[i][j];
					// this calculates the conjugate of A
					tempAi = -tempAi;
					tempBr = B->real[k][l];
					tempBi = B->imag[k][l];
                    C->real[i*2 + k][j*2 + l] +=  tempAr * tempBr - tempAi * tempBi;
					C->imag[i*2 + k][j*2 + l] +=  tempAi * tempBr + tempAr * tempBi;
					if (C->imag[i*2 + k][j*2 + l] != 0)
					{
						printf("ERROR:  This channel has a complex super operator -- NOT YET IMPLEMENTED!!!!");
					}
                } 
            } 
        } 
    } 
}


void ApplyArbitraryKraussMap(Qureg qureg, const int targetQubit, OneQubitKraussOperator operators[4])
{
	//DO the checks on the Krauss operators
	
	//Initialize the channel with 0 superoperator
	OneQubitChannel thisChannel = {.SupOp.real = {0}, .SupOp.imag = {0}};
	
    KraussOperator2SuperOperator(&operators[0], &operators[0], &thisChannel.SupOp);
	KraussOperator2SuperOperator(&operators[1], &operators[1], &thisChannel.SupOp);
	KraussOperator2SuperOperator(&operators[2], &operators[2], &thisChannel.SupOp);
	KraussOperator2SuperOperator(&operators[3], &operators[3], &thisChannel.SupOp);

	/*printf("The super operator of the channel is:\n");
	for (int i = 0; i < 4; i++) {
		for (int k = 0; k < 4; k++) { 
			printf("%f %f    ", thisChannel.SupOp.real[i][k], thisChannel.SupOp.imag[i][k]);
		}
		printf("\n");
	}*/
	
	//Calculate the sparse representation of the superoperator
	CalculateOneQubitSparseSuperOperator(&thisChannel);
	
	//Apply channel to the qubit
	ApplyOneQubitChannel_local(qureg, targetQubit, thisChannel.SparseSupOp);	
	
}
