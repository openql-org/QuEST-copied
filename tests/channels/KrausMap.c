#include "../../QuantumChannels.h"
#include "ExampleChannels.h"
#include <time.h> 



void printrho(int N, Qureg qubits) {
	qreal re,im;
	
	for(int i=0; i<pow(2,N); i++) {
	for(int j=0; j<pow(2,N); j++) {
		re = getDensityAmp(qubits, i, j).real;
		im = getDensityAmp(qubits, i, j).imag; 
		printf("(%.3f %.3f)  ", re, im);
	}
	printf("\n");
	}
	printf("\n");
}

void HilbertSchmidtDistance(int N, Qureg qubitsA, Qureg qubitsB) {
	qreal re,im;
	qreal result = 0.;
	for(int i=0; i<pow(2,N); i++) {
	for(int j=0; j<pow(2,N); j++) {
		re = getDensityAmp(qubitsA, i, j).real-getDensityAmp(qubitsB, i, j).real;
		im = getDensityAmp(qubitsA, i, j).imag-getDensityAmp(qubitsB, i, j).imag; 
		result += fabs(re)*fabs(re) + fabs(im)*fabs(im);
	}
	}
	printf("\nHilbert-Schmidt distance of the two density matrices is: %f\n", result);
}







int main (int narg, char *varg[]) {
	
	QuESTEnv env = createQuESTEnv();

	const int N = 11;
    Qureg qubitsA = createDensityQureg(N, env);
	Qureg qubitsB = createDensityQureg(N, env);
		
	for(int i = 0; i < N; i++)
	{
		rotateX (qubitsA, i, 0.12134234*i);
		rotateX (qubitsB, i, 0.12134234*i);
	}
	
	qreal prob = 0.3;
	
	clock_t t; 
	double time_taken;
	
t = clock();
	for(int i = 0; i < N; i++)
	{
		//ApplyOneQubitDepolariseChannel(qubitsA, i, prob);
		//ApplyOneQubitDephaseChannel(qubitsA, i, prob);
		//ApplyOneQubitDampingChannel(qubitsA, i, prob);
		ApplyOneQubitPauliChannel(qubitsA, i, prob/3, prob/3, prob/3);
		//rotXtest(qubitsA, i);
		for(int j = 0; j < N; j++) {
			//if (i!=j) ApplyTwoQubitDephaseChannel(qubitsA, i, j,0.1);
		}
	} 
t = clock() - t; 
    time_taken = ((double)t)/CLOCKS_PER_SEC;
    printf("Time with superopertators: %f \n", time_taken); 

t = clock(); 
	for(int i = 0; i < N; i++)
	{
		applyOneQubitDepolariseError(qubitsB, i, prob);
		//applyOneQubitDephaseError(qubitsB, i, prob);
		//applyOneQubitDampingError(qubitsB, i, prob);
		//applyOneQubitDepolariseError(qubitsB, i, prob);
		//rotateX (qubitsB, i, 1.);
		for(int j = 0; j < N; j++) {
			//if (i!=j) applyTwoQubitDephaseError(qubitsB, i, j,0.1);
		}
	} 
t = clock() - t; 
	time_taken = ((double)t)/CLOCKS_PER_SEC;
	printf("Time with the previous implementation: %f \n", time_taken); 
	
	//printrho(N, qubitsA);
	//printrho(N, qubitsB);
	
	HilbertSchmidtDistance(N, qubitsA, qubitsB);
    
    destroyQureg(qubitsA, env); 
	destroyQureg(qubitsB, env); 
    destroyQuESTEnv(env);
    return 0;
}



















