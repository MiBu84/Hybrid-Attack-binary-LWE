I. Source Code Structure
The source code should contain the following folders:
	src			: contains implementations (.cpp files)
	include			: contains headers (.h files)
	testcases		: contains input bases for HybridAttack
	reduced_basis		: contains reduced bases
	GenerateBinaryError	: code for generating input basis
	GenerateRandomError	: code for generating input basis
	MPI-output		: an output folder is generated for each MPI run.
					Each folder contains as many text files as the number of MPI processes.
Notes:
If the folder "MPI-output"is not present, it must be created otherwise the output of MPI is lost.
The program HybridAttack has two main routines called "precomputing" and "hybridAttack" which are implemented in "src/Attack.cpp".

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
II. Libraries for Compiling HybridAttack

In order to compile and run the program, the following libraries / compiler are needed:
GCC (g++)
OpenMPI (mpic++)
NTL (thread-safe compiled)
GMP
Intel TBB

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
III. Building the Application
1. Running Modes of HybridAttack
To run HybridAttack in different modes, choose

- Mode 1:
	(1) Basis is reduced without permutation and randomization.
	(2) The reduced basis is fed to a test for necessity condition. The test checks whether the hybrid attack could succeed on this given basis.
		If the test is passed, return this basis and go to (3). Else repeat step (1).
	(3) Running attack based on the reduced basis.

- Mode 2: 
	(1) Basis is permutated and multiplied with an unimodular matrix before being reduced.
	(2) The reduced basis is fed to a test for necessity condition. The test checks whether the hybrid attack could succeed on this given basis.
		If the test is passed, go to (3). Else repeat step (1).
	(3) Running attack on the basis returned from (2).

Notes: 
This mode is used for an input basis with a random binary error, which is generated by GenerateInputHybridAttackRandomError.
For example, try the input basis named test17.txt found in folder "testcases".

- Mode 3:
	Each process reads a reduced basis and goes directly to the attack step.

- Mode 4:
	MPI processes with even IDs reduce input basis before executing attack.
	MPI processes with odd IDs ignore the BKZ reduction, simply read a reduced basis and begin with attack.

Notes:
For mode 3 and 4: the naming convention for a reduced basis is: [name_of_testcase(without .txt)]_r_[value_of_r]_beta_[value_of_beta].txt
For example, if the arguments of HybridAttack are -input testcases/test13.txt -m 160 -beta 20 -r 20 then the name of the reduced basis must be 
test13_r_20_beta_20.txt. This reduced basis must be found in folder "reduced_basis".

- Mode 5:
	(1) Basis is multiplied with an unimodular matrix before being reduced.
		There is no test for necessity condition at the end of the precomputing step.
	(2) Each MPI process begins the attack step with the reduced basis from step (1).
		If timeout (given in seconds, see below) is exceeded, process returns to step (1).

Notes:
For mode 5, an argument -timeout must be given in order to determine the timeout to stop the attack process.
For example, HybridAttack is started with
./HybridAttack -input testcases/test13.txt -m 160 -beta 20 -r 20 -c 5 -numthread 8 -timeout 600

2. Output Reduced Bases to an Output File
In order to output the reduced basis to an output and choosing mode 1 remove comments in the following lines
	output_PrecomputedBasis_ToFile(B, B_1_reduced_transposed, B_reduced_gs,euclid_norm);
	return;
These lines are found at the end of the procedure "precomputing" in "src/Attack.cpp".
Recompiling the program with the chosen mode and running the program with 1 MPI process.
A reduced basis following the naming convention mentioned above is created in folder "reduced_basis".

Notes:
If there is no need to output reduced bases to text files, these two lines mentioned above should be commented out.

3. Building the Program
The program HybridAttack is compiled by running a script called build.sh, which is found in the source code folder,
with the following arguments:
	m			: as described in paper
	r			: as described in paper
	mode			: as described above	

Additionally, running build.sh also builds two executables named GenerateInputHybridAttackBinaryError and GenerateInputHybridAttackRandomError
for generating input basis. Code for these two executables GenerateInputHybridAttackBinaryError and GenerateInputHybridAttackRandomError
reside in two folders GenerateBinaryError and GenerateRandomError, correspondingly.
In summary, building the program with build.sh generates three executables, namely:
- HybridAttack
- GenerateInputHybridAttackBinaryError
- GenerateInputHybridAttackRandomError
For the usage of GenerateInputHybridAttackBinaryError and GenerateInputHybridAttackRandomError, see Section IV.

For example, running 
build.sh -m 160 -r 20 -mode 1 
will compile the program for m = 160, r = 20 in mode 1.

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
IV. Generate Input Basis for HybridAttack

1. Alternating Binary Error
To generate an input basis with a binary error vector in which 0 and 1 alternate, i.e. [0 1 (0 1)], run ./GenerateInputHybridAttackBinaryError with the following arguments:
	m			: as described in paper
	n			: as described in paper
	q			: as described in paper
	output			: name (inclusive path) for output
For example, running with
./GenerateInputHybridAttackBinaryError -m 160 -n 80 -q 521 -output test13.txt
will generate an output named test13.txt containing a matrix A, a vector b, a vector error in form [0 1 (0 1)].

2. Random Binary Error 
To generate an input basis with a random binary error, run ./GenerateInputHybridAttackRandomError with arguments equal 
to those mentioned above for the GenerateInputHybridAttackBinaryError.

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
V. Recompiling Code for Different Values of m, r
For a different value of r, code needs to be recompiled. 
Run build.sh with the desired values of m,r and mode.

For example, run build.sh -m 160 -r 20 -mode 1

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
VI. Running HybridAttack
HybridAttack needs the following arguments:
	m			: as described in paper
	r			: as described in paper
	c			: as described in paper
	input			: path to basis generated by GenerateInputHybridAttack
	output			: path to output folder (for MPI run)
	numthread		: the number of OMP threads
	beta			: blocksize for BKZ
For example, let's consider this script:
	
	currentDate=$(date '+%d_%m_%Y_%T.%6N')
	mpiOutput="MPI-output/$currentDate"
	mkdir $mpiOutput
	instance="srun ./HybridAttack -input testcases/test13.txt -m 160 -beta 20 -r 20 -c 5 -numthread 8 -output $currentDate"
	$instance

Before executing HybridAttack with MPI, for example, by calling srun in slurm, this script creates a folder whose name is set as the current system time for separately writing output stream of each MPI process.
In this example, HybridAttack gets an input basis named test13.txt and runs for m=160, beta = 20, r = 20, and c = 5.
The attack step is executed by 8 OpenMP threads within each MPI process.

Notes:
Output of MPI processes are found in folder "MPI-output". As mentioned before, the subfolder containing the MPI's output is the system time at
time of execution. The naming convention for the output of a MPI process with process ID "procID" is "output_proc_procID.txt".

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
VII. Demo
- For demo of mode 1,3, and 4, compile and run HybridAttack with the arguments: -input testcases/test13.txt -m 160 -beta 20 -r 20 -c 5
(for MPI and OMP, arguments "numthread" and "output" must also be set) 
- For demo of mode 2, compile and run HybridAttack with the arguments: -input testcases/test17.txt -m 100 -beta 3 -r 4 -c 1
(for MPI and OMP, arguments "numthread" and "output" must also be set)
