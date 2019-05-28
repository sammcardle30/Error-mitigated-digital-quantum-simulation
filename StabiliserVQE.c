//------------------------------------------------------------------------------------------------------------------------
//----------------------------------------- Error-mitigated digital quantum simulation  ----------------------------------
//------------------------------------------------------------------------------------------------------------------------


/* 
 * This program carries out the numerical simulations detailed in https://arxiv.org/abs/1807.02467.
 * It computes the ground state energy of the di-hydrogen molecule, using the Hamiltonian averaging 
 * approach, along with a unitary coupled cluster ansatz. The circuit is subject to both stochastic
 * depolarising noise, and temporally correlated over/under-rotations. These effects are mitigated
 * using the stabiliser VQE technique, and linear extrapolation.
 */






#include <stdio.h>
# include <stdlib.h>
# include <time.h>
#include <sys/time.h>
# include <math.h>
# include <unistd.h>
# include <string.h>
# include <omp.h>

#include "QuEST.h"



const long double PI = 3.14159265358979323846264338327950288419716939937510;



int nParams = 3;


const int NUMTERMS = 15;
const int TERMLENGTH = 6;
const int CORRLEN = 6;










//// DATA FUNCTIONS ////////////////////////////////////////


void readData(double hamiltonianData[NUMTERMS][TERMLENGTH])
	{
		// Read in the Hamiltonian data.
	FILE *fp;

	fp = fopen("H2at075.txt", "r");

	for(int i=0; i<NUMTERMS; i++)
		{
		for(int j=0; j<TERMLENGTH; j++)
			{
			fscanf(fp, "%lf ", &hamiltonianData[i][j]);
			}
		fscanf(fp, "\n");
		}

	fclose(fp);	
	}




//// CIRCUIT FUNCTIONS /////////////////////////////////////////
	
 
void oneQubitNoise(MultiQubit qReg, double gateNoise, int qubit, double errorArray[2])
	{
		// Adds depolarising noise after a single qubit gate.
	double singleErrorProb = gateNoise/10;
	double probNoError = 1 - singleErrorProb;
	int noiseChanceInt = rand();      
	double noiseChance = (double)noiseChanceInt/(double)RAND_MAX;

	if(noiseChance < probNoError)
		{
		// If the random number is below the error threshold, do nothing to the qubit.
		}
	else if(probNoError < noiseChance && noiseChance < (probNoError + singleErrorProb/3) )
		{
		sigmaX(qReg, qubit);
		errorArray[0] = errorArray[0] + 1;
		}
	else if( (probNoError + singleErrorProb/3) < noiseChance && noiseChance < (probNoError + 2*singleErrorProb/3) )
		{
		sigmaY(qReg, qubit);
		errorArray[0] = errorArray[0] + 1;
		}
	else
		{
		sigmaZ(qReg, qubit);
		errorArray[0] = errorArray[0] + 1;
		}

	return;
	}






void twoQubitNoise(MultiQubit qReg, double gateNoise, int qubit1, int qubit2, double errorArray[2])
	{
		// Adds depolarising noise after a two qubit gate.
	double probNoError = 1 - gateNoise;
	int noiseChanceInt = rand();      
	double noiseChance = (double)noiseChanceInt/(double)RAND_MAX;

	if(noiseChance < probNoError)
		{
		// If the random number is below the error threshold, do nothing to the qubits.
		}

	else if(probNoError < noiseChance && noiseChance < (probNoError + 1*gateNoise/15) )
		{
		sigmaX(qReg, qubit1);
		errorArray[0] = errorArray[0] + 1;
		}

	else if((probNoError + 1*gateNoise/15) < noiseChance && noiseChance < (probNoError + 2*gateNoise/15) )
		{
		sigmaY(qReg, qubit1);
		errorArray[0] = errorArray[0] + 1;
		}

	else if((probNoError + 2*gateNoise/15) < noiseChance && noiseChance < (probNoError + 3*gateNoise/15) )
		{
		sigmaZ(qReg, qubit1);
		errorArray[0] = errorArray[0] + 1;
		}

	else if((probNoError + 3*gateNoise/15) < noiseChance && noiseChance < (probNoError + 4*gateNoise/15) )
		{
		sigmaX(qReg, qubit2);
		errorArray[0] = errorArray[0] + 1;
		}

	else if((probNoError + 4*gateNoise/15) < noiseChance && noiseChance < (probNoError + 5*gateNoise/15) )
		{
		sigmaX(qReg, qubit1);
		sigmaX(qReg, qubit2);
		errorArray[0] = errorArray[0] + 1;
		}

	else if((probNoError + 5*gateNoise/15) < noiseChance && noiseChance < (probNoError + 6*gateNoise/15) )
		{
		sigmaY(qReg, qubit1);
		sigmaX(qReg, qubit2);
		errorArray[0] = errorArray[0] + 1;
		}

	else if((probNoError + 6*gateNoise/15) < noiseChance && noiseChance < (probNoError + 7*gateNoise/15) )
		{
		sigmaZ(qReg, qubit1);
		sigmaX(qReg, qubit2);
		errorArray[0] = errorArray[0] + 1;
		}

	else if((probNoError + 7*gateNoise/15) < noiseChance && noiseChance < (probNoError + 8*gateNoise/15) )
		{
		sigmaY(qReg, qubit2);
		errorArray[0] = errorArray[0] + 1;
		}

	else if((probNoError + 8*gateNoise/15) < noiseChance && noiseChance < (probNoError + 9*gateNoise/15) )
		{
		sigmaX(qReg, qubit1);
		sigmaY(qReg, qubit2);
		errorArray[0] = errorArray[0] + 1;
		}

	else if((probNoError + 9*gateNoise/15) < noiseChance && noiseChance < (probNoError + 10*gateNoise/15) )
		{
		sigmaY(qReg, qubit1);
		sigmaY(qReg, qubit2);
		errorArray[0] = errorArray[0] + 1;
		}

	else if((probNoError + 10*gateNoise/15) < noiseChance && noiseChance < (probNoError + 11*gateNoise/15) )
		{
		sigmaZ(qReg, qubit1);
		sigmaY(qReg, qubit2);
		errorArray[0] = errorArray[0] + 1;
		}

	else if((probNoError + 11*gateNoise/15) < noiseChance && noiseChance < (probNoError + 12*gateNoise/15) )
		{
		sigmaZ(qReg, qubit2);
		errorArray[0] = errorArray[0] + 1;
		}

	else if((probNoError + 12*gateNoise/15) < noiseChance && noiseChance < (probNoError + 13*gateNoise/15) )
		{
		sigmaX(qReg, qubit1);
		sigmaZ(qReg, qubit2);
		errorArray[0] = errorArray[0] + 1;
		}

	else if((probNoError + 13*gateNoise/15) < noiseChance && noiseChance < (probNoError + 14*gateNoise/15) )
		{
		sigmaY(qReg, qubit1);
		sigmaZ(qReg, qubit2);
		errorArray[0] = errorArray[0] + 1;
		}

	else 
		{
		sigmaZ(qReg, qubit1);
		sigmaZ(qReg, qubit2);
		errorArray[0] = errorArray[0] + 1;
		}
	
	return;
	}





void fillCorrArray(double corrArray[CORRLEN], double corrBound)
	{
		// Fills the correlation array with a random seed between +corrBound and -corrBound.

		for (int i=0; i<CORRLEN; i++)
		{
			double randRot = corrBound * ((double) rand() / (RAND_MAX));
			double sign = ((double) rand() / (RAND_MAX));
			if (sign < 0.5)
			{
				corrArray[i] = 1.0 - randRot;
			}
			else
			{
				corrArray[i] = 1.0 + randRot;
			}			
		}

		return;
	}	



void parametrisedCNOT(MultiQubit qReg, int controlQubit, int targetQubit, double correlationArray[CORRLEN])
	{
		// Performs a parametrised CNOT gate given by Control-(iRx(theta)).
	
		double rotationAngle = correlationArray[targetQubit] * PI;
		
		Complex r0c0 = {0, cos(rotationAngle/2)};
		Complex r0c1 = {sin(rotationAngle/2), 0};
		Complex r1c0 = {sin(rotationAngle/2), 0};
		Complex r1c1 = {0, cos(rotationAngle/2)};

		ComplexMatrix2 rotateMatrix = {r0c0, r0c1, r1c0, r1c1};
	
		controlledUnitary(qReg, controlQubit, targetQubit, rotateMatrix);
	
		return;
	}



void parametrisedHadamard(MultiQubit qReg, int targetQubit, double correlationArray[CORRLEN])
	{
		// Performs a parametrised Hadamard gate.
	
		double rotationAngle = (1.0 + ((correlationArray[targetQubit] - 1.0)/10.0)) * (PI/2);
	
		Complex r0c0 = {cos(rotationAngle/2), 0};
		Complex r0c1 = {sin(rotationAngle/2), 0};
		Complex r1c0 = {sin(rotationAngle/2), 0};
		Complex r1c1 = {-cos(rotationAngle/2), 0};
	
		ComplexMatrix2 rotateMatrix = {r0c0, r0c1, r1c0, r1c1};
	
		unitary(qReg, targetQubit, rotateMatrix);
	
		return;		
	}






//// ANSATZ FUNCTIONS ////////////////////////////////////////////
	

void hydrogenUCCterm0a(double parameter, double gateNoise, MultiQubit qReg, double corrArray[CORRLEN], double errorArray[2])
	{
		// Applies the UCC term exp(i/2 X1 Y0) to the qubit register.
		rotateX(qReg, 0, (PI/2) * (1.0 + ((corrArray[0] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 0, errorArray);
		parametrisedHadamard(qReg, 1, corrArray);
		oneQubitNoise(qReg, gateNoise, 1, errorArray);
		parametrisedCNOT(qReg, 0, 1, corrArray);
		twoQubitNoise(qReg, gateNoise, 0, 1, errorArray);

		// Note - sign as exponent has + sign.
		rotateZ(qReg, 1, -2*parameter * (1.0 + ((corrArray[1] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 1, errorArray);

		parametrisedCNOT(qReg, 0, 1, corrArray);
		twoQubitNoise(qReg, gateNoise, 0, 1, errorArray);
		parametrisedHadamard(qReg, 1, corrArray);
		oneQubitNoise(qReg, gateNoise, 1, errorArray);
		rotateX(qReg, 0, -(PI/2) * (1.0 + ((corrArray[0] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 0, errorArray);

		return;
	}


void hydrogenUCCterm0b(double parameter, double gateNoise, MultiQubit qReg, double corrArray[CORRLEN], double errorArray[2])
	{
		// Applies the UCC term exp(-i/2 Y1 X0) to the qubit register.
		rotateX(qReg, 1, (PI/2) * (1.0 + ((corrArray[1] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 1, errorArray);
		parametrisedHadamard(qReg, 0, corrArray);
		oneQubitNoise(qReg, gateNoise, 0, errorArray);
		parametrisedCNOT(qReg, 0, 1, corrArray);
		twoQubitNoise(qReg, gateNoise, 0, 1, errorArray);

		// Note + sign as exponent has - sign.
		rotateZ(qReg, 1, 2*parameter * (1.0 + ((corrArray[1] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 1, errorArray);

		parametrisedCNOT(qReg, 0, 1, corrArray);
		twoQubitNoise(qReg, gateNoise, 0, 1, errorArray);
		parametrisedHadamard(qReg, 0, corrArray);
		oneQubitNoise(qReg, gateNoise, 0, errorArray);
		rotateX(qReg, 1, -(PI/2) * (1.0 + ((corrArray[1] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 1, errorArray);

		return;
	}


void hydrogenUCCterm1a(double parameter, double gateNoise, MultiQubit qReg, double corrArray[CORRLEN], double errorArray[2])
	{
		// Applies the UCC term exp(i/2 X3 Y2) to the qubit register.
		rotateX(qReg, 2, (PI/2) * (1.0 + ((corrArray[2] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 2, errorArray);
		parametrisedHadamard(qReg, 3, corrArray);
		oneQubitNoise(qReg, gateNoise, 3, errorArray);
		parametrisedCNOT(qReg, 2, 3, corrArray);
		twoQubitNoise(qReg, gateNoise, 2, 3, errorArray);

		// Note - sign as exponent has + sign.
		rotateZ(qReg, 3, -2*parameter * (1.0 + ((corrArray[3] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 3, errorArray);

		parametrisedCNOT(qReg, 2, 3, corrArray);
		twoQubitNoise(qReg, gateNoise, 2, 3, errorArray);
		parametrisedHadamard(qReg, 3, corrArray);
		oneQubitNoise(qReg, gateNoise, 3, errorArray);
		rotateX(qReg, 2, -(PI/2) * (1.0 + ((corrArray[2] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 2, errorArray);

		return;
	}


void hydrogenUCCterm1b(double parameter, double gateNoise, MultiQubit qReg, double corrArray[CORRLEN], double errorArray[2])
	{
		// Applies the UCC term exp(-i/2 Y3 X2) to the qubit register.
		rotateX(qReg, 3, (PI/2) * (1.0 + ((corrArray[3] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 3, errorArray);
		parametrisedHadamard(qReg, 2, corrArray);
		oneQubitNoise(qReg, gateNoise, 2, errorArray);
		parametrisedCNOT(qReg, 2, 3, corrArray);
		twoQubitNoise(qReg, gateNoise, 2, 3, errorArray);

		// Note + sign as exponent has - sign.
		rotateZ(qReg, 3, 2*parameter * (1.0 + ((corrArray[3] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 3, errorArray);

		parametrisedCNOT(qReg, 2, 3, corrArray);
		twoQubitNoise(qReg, gateNoise, 2, 3, errorArray);
		parametrisedHadamard(qReg, 2, corrArray);
		oneQubitNoise(qReg, gateNoise, 2, errorArray);
		rotateX(qReg, 3, -(PI/2) * (1.0 + ((corrArray[3] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 3, errorArray);

		return;
	}

	
	
void hydrogenUCCterm2a(double parameter, double gateNoise, MultiQubit qReg, double corrArray[CORRLEN], double errorArray[2])
	{
		// Applies the UCC term exp(i $\theta$ X3 X2 Y1 X0) to the qubit register.
		parametrisedHadamard(qReg, 0, corrArray);
		oneQubitNoise(qReg, gateNoise, 0, errorArray);
		rotateX(qReg, 1, (PI/2) * (1.0 + ((corrArray[1] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 1, errorArray);
		parametrisedHadamard(qReg, 2, corrArray);
		oneQubitNoise(qReg, gateNoise, 2, errorArray);
		parametrisedHadamard(qReg, 3, corrArray);
		oneQubitNoise(qReg, gateNoise, 3, errorArray);

		parametrisedCNOT(qReg, 0, 1, corrArray);
		twoQubitNoise(qReg, gateNoise, 0, 1, errorArray);
		parametrisedCNOT(qReg, 1, 2, corrArray);
		twoQubitNoise(qReg, gateNoise, 1, 2, errorArray);
		parametrisedCNOT(qReg, 2, 3, corrArray);
		twoQubitNoise(qReg, gateNoise, 2, 3, errorArray);

		// Note - sign as exponent has + sign.
		rotateZ(qReg, 3, -2*parameter * (1.0 + ((corrArray[3] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 3, errorArray);

		parametrisedCNOT(qReg, 2, 3, corrArray);
		twoQubitNoise(qReg, gateNoise, 2, 3, errorArray);
		parametrisedCNOT(qReg, 1, 2, corrArray);
		twoQubitNoise(qReg, gateNoise, 1, 2, errorArray);
		parametrisedCNOT(qReg, 0, 1, corrArray);
		twoQubitNoise(qReg, gateNoise, 0, 1, errorArray);

		parametrisedHadamard(qReg, 3, corrArray);
		oneQubitNoise(qReg, gateNoise, 3, errorArray);
		parametrisedHadamard(qReg, 2, corrArray);
		oneQubitNoise(qReg, gateNoise, 2, errorArray);
		rotateX(qReg, 1, -(PI/2) * (1.0 + ((corrArray[1] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 1, errorArray);
		parametrisedHadamard(qReg, 0, corrArray);
		oneQubitNoise(qReg, gateNoise, 0, errorArray);

		return;
	}
	

void hydrogenUCCterm2b(double parameter, double gateNoise, MultiQubit qReg, double corrArray[CORRLEN], double errorArray[2])
	{
		// Applies the UCC term exp(i $\theta$ Y3 X2 X1 X0) to the qubit register.
		parametrisedHadamard(qReg, 0, corrArray);
		oneQubitNoise(qReg, gateNoise, 0, errorArray);
		parametrisedHadamard(qReg, 1, corrArray);
		oneQubitNoise(qReg, gateNoise, 1, errorArray);
		parametrisedHadamard(qReg, 2, corrArray);
		oneQubitNoise(qReg, gateNoise, 2, errorArray);
		rotateX(qReg, 3, (PI/2) * (1.0 + ((corrArray[3] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 3, errorArray);

		parametrisedCNOT(qReg, 0, 1, corrArray);
		twoQubitNoise(qReg, gateNoise, 0, 1, errorArray);
		parametrisedCNOT(qReg, 1, 2, corrArray);
		twoQubitNoise(qReg, gateNoise, 1, 2, errorArray);
		parametrisedCNOT(qReg, 2, 3, corrArray);
		twoQubitNoise(qReg, gateNoise, 2, 3, errorArray);

		// Note - sign as exponent has + sign.
		rotateZ(qReg, 3, -2*parameter * (1.0 + ((corrArray[3] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 3, errorArray);

		parametrisedCNOT(qReg, 2, 3, corrArray);
		twoQubitNoise(qReg, gateNoise, 2, 3, errorArray);
		parametrisedCNOT(qReg, 1, 2, corrArray);
		twoQubitNoise(qReg, gateNoise, 1, 2, errorArray);
		parametrisedCNOT(qReg, 0, 1, corrArray);
		twoQubitNoise(qReg, gateNoise, 0, 1, errorArray);

		rotateX(qReg, 3, -(PI/2) * (1.0 + ((corrArray[3] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 3, errorArray);
		parametrisedHadamard(qReg, 2, corrArray);
		oneQubitNoise(qReg, gateNoise, 2, errorArray);
		parametrisedHadamard(qReg, 1, corrArray);
		oneQubitNoise(qReg, gateNoise, 1, errorArray);
		parametrisedHadamard(qReg, 0, corrArray);
		oneQubitNoise(qReg, gateNoise, 0, errorArray);

		return;
	}
	
void hydrogenUCCterm2c(double parameter, double gateNoise, MultiQubit qReg, double corrArray[CORRLEN], double errorArray[2])
	{
		// Applies the UCC term exp(i $\theta$ Y3 Y2 Y1 X0) to the qubit register.
		parametrisedHadamard(qReg, 0, corrArray);
		oneQubitNoise(qReg, gateNoise, 0, errorArray);
		rotateX(qReg, 1, (PI/2) * (1.0 + ((corrArray[1] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 1, errorArray);
		rotateX(qReg, 2, (PI/2) * (1.0 + ((corrArray[2] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 2, errorArray);
		rotateX(qReg, 3, (PI/2) * (1.0 + ((corrArray[3] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 3, errorArray);

		parametrisedCNOT(qReg, 0, 1, corrArray);
		twoQubitNoise(qReg, gateNoise, 0, 1, errorArray);
		parametrisedCNOT(qReg, 1, 2, corrArray);
		twoQubitNoise(qReg, gateNoise, 1, 2, errorArray);
		parametrisedCNOT(qReg, 2, 3, corrArray);
		twoQubitNoise(qReg, gateNoise, 2, 3, errorArray);

		// Note - sign as exponent has + sign.
		rotateZ(qReg, 3, -2*parameter * (1.0 + ((corrArray[3] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 3, errorArray);

		parametrisedCNOT(qReg, 2, 3, corrArray);
		twoQubitNoise(qReg, gateNoise, 2, 3, errorArray);
		parametrisedCNOT(qReg, 1, 2, corrArray);
		twoQubitNoise(qReg, gateNoise, 1, 2, errorArray);
		parametrisedCNOT(qReg, 0, 1, corrArray);
		twoQubitNoise(qReg, gateNoise, 0, 1, errorArray);

		rotateX(qReg, 3, -(PI/2) * (1.0 + ((corrArray[3] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 3, errorArray);
		rotateX(qReg, 2, -(PI/2) * (1.0 + ((corrArray[2] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 2, errorArray);
		rotateX(qReg, 1, -(PI/2) * (1.0 + ((corrArray[1] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 1, errorArray);
		parametrisedHadamard(qReg, 0, corrArray);
		oneQubitNoise(qReg, gateNoise, 0, errorArray);

		return;
	}


void hydrogenUCCterm2d(double parameter, double gateNoise, MultiQubit qReg, double corrArray[CORRLEN], double errorArray[2])
	{
		// Applies the UCC term exp(i $\theta$ Y3 X2 Y1 Y0) to the qubit register.
		rotateX(qReg, 0, (PI/2) * (1.0 + ((corrArray[0] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 0, errorArray);
		rotateX(qReg, 1, (PI/2) * (1.0 + ((corrArray[1] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 1, errorArray);
		parametrisedHadamard(qReg, 2, corrArray);
		oneQubitNoise(qReg, gateNoise, 2, errorArray);
		rotateX(qReg, 3, (PI/2) * (1.0 + ((corrArray[3] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 3, errorArray);

		parametrisedCNOT(qReg, 0, 1, corrArray);
		twoQubitNoise(qReg, gateNoise, 0, 1, errorArray);
		parametrisedCNOT(qReg, 1, 2, corrArray);
		twoQubitNoise(qReg, gateNoise, 1, 2, errorArray);
		parametrisedCNOT(qReg, 2, 3, corrArray);
		twoQubitNoise(qReg, gateNoise, 2, 3, errorArray);

		// Note - sign as exponent has + sign.
		rotateZ(qReg, 3, -2*parameter * (1.0 + ((corrArray[3] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 3, errorArray);

		parametrisedCNOT(qReg, 2, 3, corrArray);
		twoQubitNoise(qReg, gateNoise, 2, 3, errorArray);
		parametrisedCNOT(qReg, 1, 2, corrArray);
		twoQubitNoise(qReg, gateNoise, 1, 2, errorArray);
		parametrisedCNOT(qReg, 0, 1, corrArray);
		twoQubitNoise(qReg, gateNoise, 0, 1, errorArray);

		rotateX(qReg, 3, -(PI/2) * (1.0 + ((corrArray[3] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 3, errorArray);
		parametrisedHadamard(qReg, 2, corrArray);
		oneQubitNoise(qReg, gateNoise, 2, errorArray);
		rotateX(qReg, 1, -(PI/2) * (1.0 + ((corrArray[1] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 1, errorArray);
		rotateX(qReg, 0, -(PI/2) * (1.0 + ((corrArray[0] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 0, errorArray);

		return;
	}
	
void hydrogenUCCterm2e(double parameter, double gateNoise, MultiQubit qReg, double corrArray[CORRLEN], double errorArray[2])
	{
		// Applies the UCC term exp(-i $\theta$ X3 Y2 Y1 Y0) to the qubit register.
		rotateX(qReg, 0, (PI/2) * (1.0 + ((corrArray[0] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 0, errorArray);
		rotateX(qReg, 1, (PI/2) * (1.0 + ((corrArray[1] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 1, errorArray);
		rotateX(qReg, 2, (PI/2) * (1.0 + ((corrArray[2] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 2, errorArray);
		parametrisedHadamard(qReg, 3, corrArray);
		oneQubitNoise(qReg, gateNoise, 3, errorArray);

		parametrisedCNOT(qReg, 0, 1, corrArray);
		twoQubitNoise(qReg, gateNoise, 0, 1, errorArray);
		parametrisedCNOT(qReg, 1, 2, corrArray);
		twoQubitNoise(qReg, gateNoise, 1, 2, errorArray);
		parametrisedCNOT(qReg, 2, 3, corrArray);
		twoQubitNoise(qReg, gateNoise, 2, 3, errorArray);

		// Note + sign as exponent has - sign.
		rotateZ(qReg, 3, 2*parameter * (1.0 + ((corrArray[3] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 3, errorArray);

		parametrisedCNOT(qReg, 2, 3, corrArray);
		twoQubitNoise(qReg, gateNoise, 2, 3, errorArray);
		parametrisedCNOT(qReg, 1, 2, corrArray);
		twoQubitNoise(qReg, gateNoise, 1, 2, errorArray);
		parametrisedCNOT(qReg, 0, 1, corrArray);
		twoQubitNoise(qReg, gateNoise, 0, 1, errorArray);

		parametrisedHadamard(qReg, 3, corrArray);
		oneQubitNoise(qReg, gateNoise, 3, errorArray);
		rotateX(qReg, 2, -(PI/2) * (1.0 + ((corrArray[2] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 2, errorArray);
		rotateX(qReg, 1, -(PI/2) * (1.0 + ((corrArray[1] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 1, errorArray);
		rotateX(qReg, 0, -(PI/2) * (1.0 + ((corrArray[0] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 0, errorArray);

		return;
	}


void hydrogenUCCterm2f(double parameter, double gateNoise, MultiQubit qReg, double corrArray[CORRLEN], double errorArray[2])
	{
		// Applies the UCC term exp(-i $\theta$ Y3 Y2 X1 Y0) to the qubit register.
		rotateX(qReg, 0, (PI/2) * (1.0 + ((corrArray[0] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 0, errorArray);
		parametrisedHadamard(qReg, 1, corrArray);
		oneQubitNoise(qReg, gateNoise, 1, errorArray);
		rotateX(qReg, 2, (PI/2) * (1.0 + ((corrArray[2] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 2, errorArray);
		rotateX(qReg, 3, (PI/2) * (1.0 + ((corrArray[3] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 3, errorArray);

		parametrisedCNOT(qReg, 0, 1, corrArray);
		twoQubitNoise(qReg, gateNoise, 0, 1, errorArray);
		parametrisedCNOT(qReg, 1, 2, corrArray);
		twoQubitNoise(qReg, gateNoise, 1, 2, errorArray);
		parametrisedCNOT(qReg, 2, 3, corrArray);
		twoQubitNoise(qReg, gateNoise, 2, 3, errorArray);

		// Note + sign as exponent has - sign.
		rotateZ(qReg, 3, 2*parameter * (1.0 + ((corrArray[3] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 3, errorArray);

		parametrisedCNOT(qReg, 2, 3, corrArray);
		twoQubitNoise(qReg, gateNoise, 2, 3, errorArray);
		parametrisedCNOT(qReg, 1, 2, corrArray);
		twoQubitNoise(qReg, gateNoise, 1, 2, errorArray);
		parametrisedCNOT(qReg, 0, 1, corrArray);
		twoQubitNoise(qReg, gateNoise, 0, 1, errorArray);

		rotateX(qReg, 3, -(PI/2) * (1.0 + ((corrArray[3] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 3, errorArray);
		rotateX(qReg, 2, -(PI/2) * (1.0 + ((corrArray[2] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 2, errorArray);
		parametrisedHadamard(qReg, 1, corrArray);
		oneQubitNoise(qReg, gateNoise, 1, errorArray);
		rotateX(qReg, 0, -(PI/2) * (1.0 + ((corrArray[0] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 0, errorArray);

		return;
	}

void hydrogenUCCterm2g(double parameter, double gateNoise, MultiQubit qReg, double corrArray[CORRLEN], double errorArray[2])
	{
		// Applies the UCC term exp(-i $\theta$ X3 Y2 X1 X0) to the qubit register.
		parametrisedHadamard(qReg, 0, corrArray);
		oneQubitNoise(qReg, gateNoise, 0, errorArray);
		parametrisedHadamard(qReg, 1, corrArray);
		oneQubitNoise(qReg, gateNoise, 1, errorArray);
		rotateX(qReg, 2, (PI/2) * (1.0 + ((corrArray[2] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 2, errorArray);
		parametrisedHadamard(qReg, 3, corrArray);
		oneQubitNoise(qReg, gateNoise, 3, errorArray);

		parametrisedCNOT(qReg, 0, 1, corrArray);
		twoQubitNoise(qReg, gateNoise, 0, 1, errorArray);
		parametrisedCNOT(qReg, 1, 2, corrArray);
		twoQubitNoise(qReg, gateNoise, 1, 2, errorArray);
		parametrisedCNOT(qReg, 2, 3, corrArray);
		twoQubitNoise(qReg, gateNoise, 2, 3, errorArray);

		// Note + sign as exponent has - sign.
		rotateZ(qReg, 3, 2*parameter * (1.0 + ((corrArray[3] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 3, errorArray);

		parametrisedCNOT(qReg, 2, 3, corrArray);
		twoQubitNoise(qReg, gateNoise, 2, 3, errorArray);
		parametrisedCNOT(qReg, 1, 2, corrArray);
		twoQubitNoise(qReg, gateNoise, 1, 2, errorArray);
		parametrisedCNOT(qReg, 0, 1, corrArray);
		twoQubitNoise(qReg, gateNoise, 0, 1, errorArray);

		parametrisedHadamard(qReg, 3, corrArray);
		oneQubitNoise(qReg, gateNoise, 3, errorArray);
		rotateX(qReg, 2, -(PI/2) * (1.0 + ((corrArray[2] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 2, errorArray);
		parametrisedHadamard(qReg, 1, corrArray);
		oneQubitNoise(qReg, gateNoise, 1, errorArray);
		parametrisedHadamard(qReg, 0, corrArray);
		oneQubitNoise(qReg, gateNoise, 0, errorArray);

		return;
	}


void hydrogenUCCterm2h(double parameter, double gateNoise, MultiQubit qReg, double corrArray[CORRLEN], double errorArray[2])
	{
		// Applies the UCC term exp(-i $\theta$ X3 X2 X1 Y0) to the qubit register.
		rotateX(qReg, 0, (PI/2) * (1.0 + ((corrArray[0] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 0, errorArray);
		parametrisedHadamard(qReg, 1, corrArray);
		oneQubitNoise(qReg, gateNoise, 1, errorArray);
		parametrisedHadamard(qReg, 2, corrArray);
		oneQubitNoise(qReg, gateNoise, 2, errorArray);
		parametrisedHadamard(qReg, 3, corrArray);
		oneQubitNoise(qReg, gateNoise, 3, errorArray);

		parametrisedCNOT(qReg, 0, 1, corrArray);
		twoQubitNoise(qReg, gateNoise, 0, 1, errorArray);
		parametrisedCNOT(qReg, 1, 2, corrArray);
		twoQubitNoise(qReg, gateNoise, 1, 2, errorArray);
		parametrisedCNOT(qReg, 2, 3, corrArray);
		twoQubitNoise(qReg, gateNoise, 2, 3, errorArray);

		// Note + sign as exponent has - sign.
		rotateZ(qReg, 3, 2*parameter * (1.0 + ((corrArray[3] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 3, errorArray);

		parametrisedCNOT(qReg, 2, 3, corrArray);
		twoQubitNoise(qReg, gateNoise, 2, 3, errorArray);
		parametrisedCNOT(qReg, 1, 2, corrArray);
		twoQubitNoise(qReg, gateNoise, 1, 2, errorArray);
		parametrisedCNOT(qReg, 0, 1, corrArray);
		twoQubitNoise(qReg, gateNoise, 0, 1, errorArray);

		parametrisedHadamard(qReg, 3, corrArray);
		oneQubitNoise(qReg, gateNoise, 3, errorArray);
		parametrisedHadamard(qReg, 2, corrArray);
		oneQubitNoise(qReg, gateNoise, 2, errorArray);
		parametrisedHadamard(qReg, 1, corrArray);
		oneQubitNoise(qReg, gateNoise, 1, errorArray);
		rotateX(qReg, 0, -(PI/2) * (1.0 + ((corrArray[0] - 1.0)/10.0)));
		oneQubitNoise(qReg, gateNoise, 0, errorArray);

		return;
	}


			
void makeAnsatz(double parameters[nParams], double gateNoise, MultiQubit qReg, double corrArray[CORRLEN], double errorArray[2])
{
	// Creates the UCC ansatz applied to H2, starting from a HF state.

	initStateZero(&qReg);

	// Create Hartree-Fock state.
	sigmaX(qReg, 0);
	sigmaX(qReg, 2);

	// Apply UCC circuit.
	hydrogenUCCterm2a(parameters[2], gateNoise, qReg, corrArray, errorArray);
	hydrogenUCCterm2b(parameters[2], gateNoise, qReg, corrArray, errorArray);
	hydrogenUCCterm2c(parameters[2], gateNoise, qReg, corrArray, errorArray);
	hydrogenUCCterm2d(parameters[2], gateNoise, qReg, corrArray, errorArray);
	hydrogenUCCterm2e(parameters[2], gateNoise, qReg, corrArray, errorArray);
	hydrogenUCCterm2f(parameters[2], gateNoise, qReg, corrArray, errorArray);
	hydrogenUCCterm2g(parameters[2], gateNoise, qReg, corrArray, errorArray);
	hydrogenUCCterm2h(parameters[2], gateNoise, qReg, corrArray, errorArray);


	hydrogenUCCterm1a(parameters[1], gateNoise, qReg, corrArray, errorArray);
	hydrogenUCCterm1b(parameters[1], gateNoise, qReg, corrArray, errorArray);

	hydrogenUCCterm0a(parameters[0], gateNoise, qReg, corrArray, errorArray);
	hydrogenUCCterm0b(parameters[0], gateNoise, qReg, corrArray, errorArray);

	return;
}				
			
				





//// MEASUREMENT FUNCTIONS ///////////////////////////////////////////////

		
void measureStateOnce(int numMeasured, double measuredReg[numMeasured], MultiQubit qReg)
	{
		// Measures the register such that it collapses to a single state in the computational basis.
		
		for (int i=0; i<numMeasured; i++)
		{
			measuredReg[i] = measure(qReg, i);			
		}
	
		return;
	}	


double getSubtermMeasurement(int numQubits, double parameters[nParams], double gateNoise,
	MultiQubit qReg, double hamilTerm[TERMLENGTH], double errorArray[2], int numUp, int numParticles, double corrArray[CORRLEN])
	{
		// Performs a single measurement of one of the Pauli strings in the Hamiltonian.
	
		int nRegister = numQubits - 2;
		makeAnsatz(parameters, gateNoise, qReg, corrArray, errorArray);
	
		// Parity check.
		// First check spin up parity.
		parametrisedCNOT(qReg, 1, 0, corrArray);
		twoQubitNoise(qReg, gateNoise, 1, 0, errorArray);
	
		parametrisedCNOT(qReg, 0, 4, corrArray);
		twoQubitNoise(qReg, gateNoise, 0, 4, errorArray);
	
		int parityUpMeasure = measure(qReg, 4);
		double parityUpValue = pow(-1, parityUpMeasure);
	
		parametrisedCNOT(qReg, 0, 4, corrArray);
		twoQubitNoise(qReg, gateNoise, 0, 4, errorArray);
	
		parametrisedCNOT(qReg, 1, 0, corrArray);
		twoQubitNoise(qReg, gateNoise, 1, 0, errorArray);
	

		// Check spin down parity.
		parametrisedCNOT(qReg, 2, 3, corrArray);
		twoQubitNoise(qReg, gateNoise, 2, 3, errorArray);
	
		parametrisedCNOT(qReg, 3, 5, corrArray);
		twoQubitNoise(qReg, gateNoise, 3, 5, errorArray);
	
		int parityDownMeasure = measure(qReg, 5);
		double parityDownValue = pow(-1, parityDownMeasure);
	
		parametrisedCNOT(qReg, 3, 5, corrArray);
		twoQubitNoise(qReg, gateNoise, 3, 5, errorArray);
	
		parametrisedCNOT(qReg, 2, 3, corrArray);
		twoQubitNoise(qReg, gateNoise, 2, 3, errorArray);
	

		// Check if parity has been violated. Return zero if so.
		// If not, carry out the measurement of h_i.
		int numDown = numParticles - numUp;
		double parityUp = pow(-1, numUp);
		double parityDown = pow(-1, numDown);
	
		// Avoid comparing doubles in c.
		double upDiff = fabs(parityUpValue - parityUp);
		double downDiff = fabs(parityDownValue - parityDown);
	
		if (upDiff > 0.1 || downDiff > 0.1)
		// If either parity is not equal to true value, discard measurement.
		{
			return 0.0;
		}

		else
		// Carry out the measurement of h_i.
		{
			// Change basis such that measuring Pauli term becomes equivalent to measuring in the computational basis.
			for (int k=0; k<nRegister; k++)
			{
				double gate = hamilTerm[k+2];
				if (gate > 0.98 && gate < 1.02)
				{
					hadamard(qReg, k);
				}
				else if (gate > 1.98 && gate < 2.02)
				{
					rotateX(qReg, k, PI/2);
				}
				else
				{
					// Do nothing for I or Z.
				}
			}
			// Measure the main register in the Z basis.
			double reg[nRegister];
			measureStateOnce(nRegister, reg, qReg);
		
			// Calculate the parity of relevant qubits in the term.
			double paritySum = 0.0;
			for (int l=0; l<nRegister; l++)
			{
				double gate = hamilTerm[l+2];
				if (gate > 0.02)
				{
					paritySum = paritySum + reg[l];
				}
			}
			double parityMeasurement = pow(-1, paritySum);
		
			return parityMeasurement;
		}
	}


double getSubtermExpectation(int numQubits, double parameters[nParams], double gateNoise, 
	MultiQubit qReg, double hamilTerm[TERMLENGTH], int nRepeats, double errorArray[2], int numUp, int numParticles, double corrBound)
	{
		// Performs multiple state initialisations and measurements to measure the expectation value of a single Pauli term in the Hamiltonian.
		
		int errorCounter = 0;
		double expectationSum = 0.0;
		double corrArray[CORRLEN];

		for (int i=0; i<nRepeats; i++)
		{
			fillCorrArray(corrArray, corrBound);
			double measuredValue = getSubtermMeasurement(numQubits, parameters, gateNoise, qReg, hamilTerm,
				errorArray, numUp, numParticles, corrArray);
			if (measuredValue < 0.01 && measuredValue > -0.01)
			// Returned values are +/- 1 for no error detected, 0 for error detected, so a value between +/- 0.01 implies an error detection. 
			{
				errorCounter = errorCounter + 1;
			}
			expectationSum = expectationSum + measuredValue;
		}
		
		double nTrials = nRepeats - errorCounter;
		double expectation = expectationSum/nTrials;
		errorArray[1] = errorArray[1] + errorCounter;
		
		return expectation;
	}



double getExtrapolatedExpectation(int numQubits, double parameters[nParams], double gateNoise, MultiQubit qReg, 
	double hamilTerm[TERMLENGTH], int nRepeats, double errorArray[2], int numUp, int numParticles, double corrBound)
	{
		// Performs linear error extrapolation of a single Hamiltonian term expectation value.
		double stretchFactor = 1.0 + 0.0005/gateNoise;
		double stretchOverRot = pow(stretchFactor, 0.5);
		
		double expectationNormalError = getSubtermExpectation(numQubits, parameters, gateNoise, qReg, 
			hamilTerm, nRepeats/2, errorArray, numUp, numParticles, corrBound);
		double expectationStretchedError = getSubtermExpectation(numQubits, parameters, stretchFactor*gateNoise, qReg, 
			hamilTerm, nRepeats/2, errorArray, numUp, numParticles, stretchOverRot*corrBound);
		
		double gradient = (expectationStretchedError - expectationNormalError)/((stretchFactor - 1)*(gateNoise));
	    double extrapolatedValue = expectationNormalError - gradient*gateNoise;
		
		return extrapolatedValue;
	}



double getRepititionRatio(double desiredPrecision, double hamiltonian[NUMTERMS][TERMLENGTH])
	{
		// Calculates the repitition ratio required to achieve a specified precision in the final energy value.

		// Find the max value in the Hamiltonian, and get the sum of the terms in the Hamiltonian.
		double termSum = 0.0;
		double maxValue = 0.0;
		for (int i=1; i<NUMTERMS; i++)
		{
			if (fabs(hamiltonian[i][1]) > maxValue)
			{
				maxValue = fabs(hamiltonian[i][1]);
			}
			termSum = termSum + fabs(hamiltonian[i][1]);
		}

		// Get the repitition ratio, i.e. the maximum number of measurements any term gets, divided by the maximum Hamiltonian coefficient. 
		double maxReps = (maxValue*termSum)/(desiredPrecision*desiredPrecision);
		double repititionRatio = maxReps/maxValue;

		return repititionRatio;
	} 


	
double measureEnergy(int numQubits, double parameters[nParams], double gateNoise, MultiQubit qReg, double hamiltonian[NUMTERMS][TERMLENGTH], double precision,
	 double errorArray[2], int numUp, int numParticles, 
	 int extrapolation, double corrBound)
	{
		// Measures the energy of the state defined by the parameters by performing many state initialisations and measurements.
		
		double energy = 0.0;
		
		double nuclearRepulsion = hamiltonian[0][1];
		energy = energy + nuclearRepulsion;

		double repRatio = getRepititionRatio(precision, hamiltonian);
		
		for (int i=1; i<NUMTERMS; i++)
		{
			double termCoefficient = hamiltonian[i][1];
			double termExpectation;
			int termReps = floor(repRatio * fabs(termCoefficient));
			
			if (extrapolation == 1)
			{
				termExpectation = getExtrapolatedExpectation(numQubits, parameters, gateNoise, qReg, 
					hamiltonian[i], termReps, errorArray, numUp, numParticles, corrBound);
				// Stop over-extrapolation at points close to the boundaries.
				if (termExpectation < -1.0)
				{
					termExpectation = -1.0;
				}
				else if (termExpectation > 1.0)
				{
					termExpectation = 1.0;
				}
			}
			else
			{
				termExpectation = getSubtermExpectation(numQubits, parameters, gateNoise, qReg, hamiltonian[i], 
				termReps, errorArray, numUp, numParticles, corrBound);
			}
			double termEnergy = termExpectation * termCoefficient;
			printf("Term expectation:%.14f\n", termExpectation);
			energy = energy + termEnergy;
		}

		return energy;
	}


/* === Hash function for seeding Rand(). === */
unsigned long hashParams(int numQubits, float distance, long unsigned jobNumber) 
	{
	float dist = 1000*distance;
	int dist2 = (int) dist;
	
	unsigned long hash = 5381;
	hash = ((hash << 6) + hash) + dist2;
	hash = ((hash << 3) + hash) + numQubits;
	hash = ((hash << 4) + hash) + jobNumber;
	hash = ((hash << 2) + hash) + dist2;
	
	return hash;
	}

















//--------------------------------------------------------------
//---------------------- START OF main()  ----------------------
//--------------------------------------------------------------
int main (int narg, char** varg) {

	// ==== INITIALISE ENVIRONMENT
	QuESTEnv env;
	initQuESTEnv(&env);

	
	// Get command line arguments.	
	float distance = atof(varg[1]);	
	char **restrict extraParam;
	unsigned long jobID = strtoul(varg[2], extraParam, 10);

	int numQubits = 6;
	MultiQubit qubits; 
	createMultiQubit(&qubits, numQubits, env);
    initStateZero(&qubits);
	
	// Turn randomness on.	
	unsigned long int seedArray[2] = {jobID, numQubits};
	QuESTSeedRandom	(seedArray, 2);
	
	unsigned long myHash = hashParams(numQubits, distance, jobID);
	srand(myHash);
	printf("Seed: %li\n", myHash);
	printf("\n");
	
	
	// END INITIALISATION






	// ==== MAIN PROGRAM
	
	// Start timing.
	struct timeval timeInstance;
	gettimeofday(&timeInstance, NULL);
	long double startTime = (timeInstance.tv_sec + (long double) timeInstance.tv_usec/pow(10,6));
	
	

	
	
	double hamiltonian[NUMTERMS][TERMLENGTH];
	readData(hamiltonian);
	
	double parameters[3] = {0.0, 0.0, 0.114833226384/8};
	double gateNoise = 0.001;
	double overRot = 0.01;
	double chosenPrecision = 0.0002;
	double errorArray[2] = {0.0, 0.0};
	int numElectrons = 2;
	int numSpinUp = 1;
	int extrap = 1;
	
	double energy = measureEnergy(numQubits, parameters, gateNoise, qubits, hamiltonian, chosenPrecision, errorArray,
		 numSpinUp, numElectrons, extrap, overRot);

	printf("Energy:%.14f\n", energy);
	printf("Number of errors:%.2f\n", errorArray[0]);
	printf("Number of detections:%.2f\n", errorArray[1]);

	
	

	// Finish timing.
	gettimeofday(&timeInstance, NULL);
	long double endTime = (timeInstance.tv_sec + (long double) timeInstance.tv_usec/pow(10,6));
	long double duration = endTime - startTime;
	if (env.rank==0) printf("Time taken (s) : %.14Lf\n", duration);
	
	
	// END MAIN PROGRAM






	

	
	// ==== MEMORY AND ENVIRONMENT CLEANUP
	
	destroyMultiQubit(qubits, env);
	closeQuESTEnv(env);



	return EXIT_SUCCESS;
}
