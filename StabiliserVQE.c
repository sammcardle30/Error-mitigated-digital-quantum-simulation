
//------------------------------------------------------------------------------------------------------------------------
//----------------------------------------- Error mitigated digital quantum simulation -----------------------------------
//------------------------------------------------------------------------------------------------------------------------


/* 
 * This program implements the stabilier-VQE technique introduced in arXiv:1807.02467.  
/*



/*
 * The QuEST library can be found at https://quest.qtechtheory.org/
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



//! 1: print end qubit state to file, 0: don't print.
# define REPORT_STATE 1

const long double PI = 3.14159265358979323846264338327950288419716939937510;



int nParams = 3;


const int NUMTERMS = 15;
const int TERMLENGTH = 6;










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


void controlledY(MultiQubit qReg, int controlQubit, int targetQubit)
{
	// Performs a controlled-Y gate.
	
	controlledRotateX(qReg, controlQubit, targetQubit, PI);
	hadamard(qReg, targetQubit);
	controlledNot(qReg, controlQubit, targetQubit);
	hadamard(qReg, targetQubit);
	
	return;	
}


void controlledZ(MultiQubit qReg, int controlQubit, int targetQubit)
{
	// Performs a controlled-Z gate.
	
	hadamard(qReg, targetQubit);
	controlledNot(qReg, controlQubit, targetQubit);
	hadamard(qReg, targetQubit);
	
	return;	
}


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


void hydrogenUCCterm0a(double parameter, double gateNoise, MultiQubit qReg, double errorArray[2])
		{
			// Applies the UCC term exp(i/2 $\theta$ X3 Z2 Y1) to the qubit register.
			rotateX(qReg, 1, PI/2);
			oneQubitNoise(qReg, gateNoise, 1, errorArray);
			hadamard(qReg, 3);
			oneQubitNoise(qReg, gateNoise, 3, errorArray);
			controlledNot(qReg, 1, 2);
			twoQubitNoise(qReg, gateNoise, 1, 2, errorArray);
			controlledNot(qReg, 2, 3);
			twoQubitNoise(qReg, gateNoise, 2, 3, errorArray);
			
			// Note - sign as exponent has + sign.
			rotateZ(qReg, 3, -2*parameter);
			oneQubitNoise(qReg, gateNoise, 3, errorArray);
			
			controlledNot(qReg, 2, 3);
			twoQubitNoise(qReg, gateNoise, 2, 3, errorArray);
			controlledNot(qReg, 1, 2);
			twoQubitNoise(qReg, gateNoise, 1, 2, errorArray);
			hadamard(qReg, 3);
			oneQubitNoise(qReg, gateNoise, 3, errorArray);
			rotateX(qReg, 1, -PI/2);
			oneQubitNoise(qReg, gateNoise, 1, errorArray);
			
			return;
		}


void hydrogenUCCterm0b(double parameter, double gateNoise, MultiQubit qReg, double errorArray[2])
		{
			// Applies the UCC term exp(-i/2 $\theta$ Y3 Z2 X1) to the qubit register.
			rotateX(qReg, 3, PI/2);
			oneQubitNoise(qReg, gateNoise, 3, errorArray);
			hadamard(qReg, 1);
			oneQubitNoise(qReg, gateNoise, 1, errorArray);
			controlledNot(qReg, 1, 2);
			twoQubitNoise(qReg, gateNoise, 1, 2, errorArray);
			controlledNot(qReg, 2, 3);
			twoQubitNoise(qReg, gateNoise, 2, 3, errorArray);
	
			// Note + sign as exponent has - sign.
			rotateZ(qReg, 3, 2*parameter);
			oneQubitNoise(qReg, gateNoise, 3, errorArray);
	
			controlledNot(qReg, 2, 3);
			twoQubitNoise(qReg, gateNoise, 2, 3, errorArray);
			controlledNot(qReg, 1, 2);
			twoQubitNoise(qReg, gateNoise, 1, 2, errorArray);
			hadamard(qReg, 1);
			oneQubitNoise(qReg, gateNoise, 1, errorArray);
			rotateX(qReg, 3, -PI/2);
			oneQubitNoise(qReg, gateNoise, 3, errorArray);
	
			return;
		}


void hydrogenUCCterm1a(double parameter, double gateNoise, MultiQubit qReg, double errorArray[2])
		{
			// Applies the UCC term exp(i/2 $\theta$ X2 Z1 Y0) to the qubit register.
			rotateX(qReg, 0, PI/2);
			oneQubitNoise(qReg, gateNoise, 0, errorArray);
			hadamard(qReg, 2);
			oneQubitNoise(qReg, gateNoise, 2, errorArray);
			controlledNot(qReg, 0, 1);
			twoQubitNoise(qReg, gateNoise, 0, 1, errorArray);
			controlledNot(qReg, 1, 2);
			twoQubitNoise(qReg, gateNoise, 1, 2, errorArray);
	
			// Note - sign as exponent has + sign.
			rotateZ(qReg, 2, -2*parameter);
			oneQubitNoise(qReg, gateNoise, 2, errorArray);
	
			controlledNot(qReg, 1, 2);
			twoQubitNoise(qReg, gateNoise, 1, 2, errorArray);
			controlledNot(qReg, 0, 1);
			twoQubitNoise(qReg, gateNoise, 0, 1, errorArray);
			hadamard(qReg, 2);
			oneQubitNoise(qReg, gateNoise, 2, errorArray);
			rotateX(qReg, 0, -PI/2);
			oneQubitNoise(qReg, gateNoise, 0, errorArray);
	
			return;
		}


void hydrogenUCCterm1b(double parameter, double gateNoise, MultiQubit qReg, double errorArray[2])
		{
			// Applies the UCC term exp(-i/2 $\theta$ Y2 Z1 X0) to the qubit register.
			rotateX(qReg, 2, PI/2);
			oneQubitNoise(qReg, gateNoise, 2, errorArray);
			hadamard(qReg, 0);
			oneQubitNoise(qReg, gateNoise, 0, errorArray);
			controlledNot(qReg, 0, 1);
			twoQubitNoise(qReg, gateNoise, 0, 1, errorArray);
			controlledNot(qReg, 1, 2);
			twoQubitNoise(qReg, gateNoise, 1, 2, errorArray);

			// Note + sign as exponent has - sign.
			rotateZ(qReg, 2, 2*parameter);
			oneQubitNoise(qReg, gateNoise, 2, errorArray);

			controlledNot(qReg, 1, 2);
			twoQubitNoise(qReg, gateNoise, 1, 2, errorArray);
			controlledNot(qReg, 0, 1);
			twoQubitNoise(qReg, gateNoise, 0, 1, errorArray);
			hadamard(qReg, 0);
			oneQubitNoise(qReg, gateNoise, 0, errorArray);
			rotateX(qReg, 2, -PI/2);
			oneQubitNoise(qReg, gateNoise, 2, errorArray);

			return;
		}
				
				
void hydrogenUCCterm2(double parameter, double gateNoise, MultiQubit qReg, double errorArray[2])
		{
			// Applies the UCC term exp(-i/2 $\theta$ X3 X2 X1 Y0) to the qubit register.
			rotateX(qReg, 0, PI/2);
			oneQubitNoise(qReg, gateNoise, 0, errorArray);
			hadamard(qReg, 1);
			oneQubitNoise(qReg, gateNoise, 1, errorArray);
			hadamard(qReg, 2);
			oneQubitNoise(qReg, gateNoise, 2, errorArray);
			hadamard(qReg, 3);
			oneQubitNoise(qReg, gateNoise, 3, errorArray);
			
			controlledNot(qReg, 0, 1);
			twoQubitNoise(qReg, gateNoise, 0, 1, errorArray);
			controlledNot(qReg, 1, 2);
			twoQubitNoise(qReg, gateNoise, 1, 2, errorArray);
			controlledNot(qReg, 2, 3);
			twoQubitNoise(qReg, gateNoise, 2, 3, errorArray);

			// Note + sign as exponent has - sign.
			rotateZ(qReg, 3, 2*parameter);
			oneQubitNoise(qReg, gateNoise, 3, errorArray);

			controlledNot(qReg, 2, 3);
			twoQubitNoise(qReg, gateNoise, 2, 3, errorArray);
			controlledNot(qReg, 1, 2);
			twoQubitNoise(qReg, gateNoise, 1, 2, errorArray);
			controlledNot(qReg, 0, 1);
			twoQubitNoise(qReg, gateNoise, 0, 1, errorArray);
			
			hadamard(qReg, 3);
			oneQubitNoise(qReg, gateNoise, 3, errorArray);
			hadamard(qReg, 2);
			oneQubitNoise(qReg, gateNoise, 2, errorArray);
			hadamard(qReg, 1);
			oneQubitNoise(qReg, gateNoise, 1, errorArray);
			rotateX(qReg, 0, -PI/2);
			oneQubitNoise(qReg, gateNoise, 0, errorArray);

			return;
		}
				


				
				
void makeAnsatz(int numQubits, double parameters[nParams], double gateNoise, MultiQubit qReg, double errorArray[2])
	{
		// Creates the UCC ansatz applied to H2, starting from a HF state.
		initStateZero(&qReg);
		
		// Create Hartree-Fock state.
		sigmaX(qReg, 0);
		sigmaX(qReg, 1);
		
		// Apply UCC circuit.
		hydrogenUCCterm2(parameters[2], gateNoise, qReg, errorArray);
		hydrogenUCCterm1a(parameters[1], gateNoise, qReg, errorArray);
		hydrogenUCCterm1b(parameters[1], gateNoise, qReg, errorArray);
		hydrogenUCCterm0a(parameters[0], gateNoise, qReg, errorArray);
		hydrogenUCCterm0b(parameters[0], gateNoise, qReg, errorArray);
		
		return;
	}				
				
				
				
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
	MultiQubit qReg, double hamilTerm[TERMLENGTH], double errorArray[2], int numUp, int numParticles)
	{
		// Performs a single measurement of one of the Pauli strings in the Hamiltonian.
		int nRegister = numQubits - 1;
		makeAnsatz(numQubits, parameters, gateNoise, qReg, errorArray);
		
		// Parity check.
		// First check spin up parity.
		for (int i=0; i<nRegister; i++)
		{
			if (i%2 == 0)
			{
				controlledNot(qReg, i, numQubits-1);
				twoQubitNoise(qReg, gateNoise, i, numQubits-1, errorArray);
			}
		}
		int parityUpMeasure = measure(qReg, numQubits-1);
		double parityUpValue = pow(-1, parityUpMeasure);

		// Reinitialise ancilla in |0>.
		if (parityUpMeasure == 1)
		{
			sigmaX(qReg, numQubits-1);
		}

		// Check spin down parity.
		for (int j=0; j<nRegister; j++)
		{
			if (j%2 == 1)
			{
				controlledNot(qReg, j, numQubits-1);
				twoQubitNoise(qReg, gateNoise, j, numQubits-1, errorArray);
			}
		}
		int parityDownMeasure = measure(qReg, numQubits-1);
		double parityDownValue = pow(-1, parityDownMeasure);
		
		// Check if parity has been violated. Return zero if so.
		// If not, carry out the measurement of h_i.
		int numDown = numParticles - numUp;
		double parityUp = pow(-1, numUp);
		double parityDown = pow(-1, numDown);
		
		// Avoid comparing doubles in c.
		double upDiff = fabs(parityUpValue - parityUp);
		double downDiff = fabs(parityDownValue - parityDown);
		
		if (upDiff > 0.1 || downDiff > 0.1)
		{
			return 0.0;
		}
		// Carry out the measurement of h_i.
		else
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
	MultiQubit qReg, double hamilTerm[TERMLENGTH], int nRepeats, double errorArray[2], int numUp, int numParticles)
	{
		// Performs multiple state initialisations and measurements to measure the expectation value of a single Pauli term in the Hamiltonian.
		int errorCounter = 0;
		double expectationSum = 0.0;
		for (int i=0; i<nRepeats; i++)
		{
			double measuredValue = getSubtermMeasurement(numQubits, parameters, gateNoise, qReg, hamilTerm,
				errorArray, numUp, numParticles);
			if (measuredValue < 0.01 && measuredValue > -0.01)
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


double getExtrapolatedExpectation(int numQubits, double parameters[nParams], double gateNoise, MultiQubit qReg, double hamilTerm[TERMLENGTH], int nRepeats,
	double errorArray[2], int numUp, int numParticles)
	{
		// Performs linear error extrapolation of a single Hamiltonian term expectation value.
		double stretchFactor = 1.0 + 0.005/gateNoise;

		double expectationNormalError = getSubtermExpectation(numQubits, parameters, gateNoise, qReg, hamilTerm, nRepeats/2,
			errorArray, numUp, numParticles);
		double expectationStretchedError =  getSubtermExpectation(numQubits, parameters, stretchFactor*gateNoise, qReg, hamilTerm, nRepeats/2,
			errorArray, numUp, numParticles);

		double gradient = (expectationStretchedError - expectationNormalError)/((stretchFactor - 1)*(gateNoise));
	    double extrapolatedValue = expectationNormalError - gradient*gateNoise;

		return extrapolatedValue;
	}


double measureEnergy(int numQubits, double parameters[nParams], double gateNoise, MultiQubit qReg, 
	double hamiltonian[NUMTERMS][TERMLENGTH], int nRepeats, double errorArray[2], int numUp, int numParticles, 
	int extrapolation)
	{
		// Measures the energy of the state defined by the parameters by performing many state initialisations and measurements.
		double energy = 0.0;
		
		for (int i=0; i<NUMTERMS; i++)
		{
			double termCoefficient = hamiltonian[i][1];
			double termExpectation;
			if (extrapolation == 1)
			{
				termExpectation = getExtrapolatedExpectation(numQubits, parameters, gateNoise, qReg, 
					hamiltonian[i], nRepeats, errorArray, numUp, numParticles);
			}
			else
			{
				termExpectation = getSubtermExpectation(numQubits, parameters, gateNoise, qReg, hamiltonian[i], 
				nRepeats, errorArray, numUp, numParticles);
			}
			double termEnergy = termExpectation * termCoefficient;
			
			energy = energy + termEnergy;
		}
		return energy;
	}


unsigned long hashParams(int numQubits, float distance, long unsigned jobNumber) 
	{
	// Hash function for seeding Rand().

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

	int numQubits = 5;
	MultiQubit qubits; 
	createMultiQubit(&qubits, numQubits, env);
    initStateZero(&qubits);
	
	//unsigned long int seedArray[2] = {jobID, numQubits};
	
	//QuESTSeedRandom	(seedArray, 2);
	QuESTSeedRandomDefault();
	// END INITIALISATION






	// ==== MAIN PROGRAM

	// Turn randomness on.	
	unsigned long myHash = hashParams(numQubits, distance, jobID);
	srand(myHash);
	printf("Seed: %li\n", myHash);
	printf("\n");
	
	// Start timing.
	struct timeval timeInstance;
	gettimeofday(&timeInstance, NULL);
	long double startTime = (timeInstance.tv_sec + (long double) timeInstance.tv_usec/pow(10,6));
	
	




	
	
	double hamiltonian[NUMTERMS][TERMLENGTH];
	readData(hamiltonian);
	
	double parameters[3] = {0.05, -0.07, 0.1};
	double gateNoise = 0.005;
	int nRepeats = 1000000;
	double errorArray[2] = {0.0, 0.0};
	int numElectrons = 2;
	int numSpinUp = 1;
	int extrap = 0; // Set to 1 if extrapolation is desired.
	
	double energy = measureEnergy(numQubits, parameters, gateNoise, qubits, hamiltonian, nRepeats, errorArray,
		 numSpinUp, numElectrons, extrap);

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

//--------------------------------------------------------------
//----------------------- END OF main()  -----------------------
//--------------------------------------------------------------
