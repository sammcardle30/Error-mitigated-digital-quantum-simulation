# Error mitigated digital quantum simulation
Simulation code generated in preparation of the work 'Error mitigated digital quantum simulation' https://arxiv.org/abs/1807.02467.



# Setup

First download the QuEST folder, and ensure that the C file and makefile are in the same directory as this folder (not within
the folder itself). QuEST has been updated since this research was carried out, so the current version of QuEST available from 
https://www.quest.qtechtheory.org/ will not work with this program. 

In order to run the file, open the command prompt and navigate to the directory containing the C file. Type:

'make clean'

'make'

'./VQE Distance JobID'

to run the file. (Distance is the interatomic separation, and JobID is an identifier that is different for every simulation.
I used the Slurm Job ID for my simulations. These variables ensure random seeding for each run.)


# Hamiltonian data

The Hamiltonian data files give the Jordan-Wigner encoded Hamiltonian for the Hydrogen molecule at a range of interatomic
separations. Each row denotes a term in the Hamiltonian. The first column gives the interatomic separation, which is also 
given in the filename. The second column gives the coefficient of that term. The other columns represent Pauli operators 
acting on the corresponding qubit (3rd column= zeroth qubit, sixth column = last qubit). The notation is as follows:

0 = I

1 = X

2 = Y

3 = Z

For example, the first two terms of the Hamiltonian at 0.75 A are: -0.109731 IIII + 0.169885 ZIII + ...


# Running StabiliserVQE.c

The simulation program StabiliserVQE.c takes several variables, which are updated within the C code itself. These are:

double parameters[3] : The three parameters for the UCCSD ansatz.

double gateNoise : The 2 qubit gate noise (A number between 0 and 1).

double overRot : The percentage over/under-rotation in the temporally correlated noise model.

double chosenPrecision : An upper bound on the desired shot-noise in the result, when error mitigation is not applied.
	
int numElectrons = The number of electrons in the H2 molecule ground state.

int numSpinUp = The number of spin-up electrons in the H2 molecule ground state.
  
int extrap :  Whether error extrapolation is on or off; set as 0 for off, set as 1 for on.

